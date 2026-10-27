# -*- coding: utf-8 -*-
"""
Tests for the incoming-wave split in nest2_sfincs_in_sfincs (filter_incoming=True).

A progressive long wave is prescribed at the overall observation points with
q = c * eta along a chosen direction, so the Guza split is exact: if that direction
points into the detail model the whole signal is incoming (bzi = eta, bzs = mean);
if it points out of the model nothing is incoming (bzi = 0).
"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from conftest import EPSG, make_mock_overall_sfincs
from shapely.geometry import Point

from cht_nesting.nest2.nest2_sfincs_in_sfincs import nest2_sfincs_in_sfincs

pytest.importorskip("hydromt_sfincs")

G = 9.81
DEPTH = 20.0
AMP = 0.5
PERIOD = 600.0
PREFIX = "nest"


def _times(n=241, dt=30.0):
    return pd.date_range("2020-01-01", periods=n, freq=f"{int(dt)}s")


def _wave(times):
    t = (times - times[0]).total_seconds().values
    # cosine over whole periods: zero mean and zero linear trend, so the split is exact
    return AMP * np.cos(2.0 * np.pi * t / PERIOD)


def _write_his(path: Path, times, xy, direction, zb_1d=True, mean_level=0.3):
    """his-file with a wave travelling in `direction` (unit vector, x/y) at every point."""
    eta = _wave(times)
    hh = DEPTH + eta
    c = np.sqrt(G * hh)
    speed = c * eta / hh  # q = u * hh = c * eta -> exact progressive wave
    n_t, n_p = len(times), len(xy)
    zs = np.repeat((mean_level + eta)[:, None], n_p, axis=1)
    u = np.repeat((speed * direction[0])[:, None], n_p, axis=1)
    v = np.repeat((speed * direction[1])[:, None], n_p, axis=1)
    zb = np.full(n_p, -DEPTH + mean_level)
    names = np.array(
        [f"{PREFIX}_{i + 1:04d}".encode() for i in range(n_p)], dtype="S40"
    )
    ds = xr.Dataset(
        {
            "point_zs": (["time", "stations"], zs),
            "point_u": (["time", "stations"], u),
            "point_v": (["time", "stations"], v),
            "point_zb": (["stations"], zb)
            if zb_1d
            else (["time", "stations"], np.repeat(zb[None, :], n_t, axis=0)),
            "station_name": (["stations"], names),
            "station_x": (["stations"], xy[:, 0]),
            "station_y": (["stations"], xy[:, 1]),
            "crs": ([], EPSG),
        },
        coords={"time": times},
    )
    ds.to_netcdf(path / "sfincs_his.nc")


def _detail(tmp_path, xy, name="detail"):
    from hydromt_sfincs import SfincsModel

    root = tmp_path / name
    root.mkdir()
    mod = SfincsModel(root=str(root), mode="w+")
    mod.grid.create(
        x0=0.0, y0=0.0, dx=1000.0, dy=1000.0, nmax=5, mmax=5, rotation=0.0, epsg=EPSG
    )
    gdf = gpd.GeoDataFrame(
        {"name": [f"b{i}" for i in range(len(xy))]},
        geometry=[Point(x, y) for x, y in xy],
        crs=f"EPSG:{EPSG}",
    )
    mod.water_level.set_locations(gdf)
    return mod


# West boundary (x = 0) of a 5 x 5 km grid: inward = +x.  South boundary (y = 0): inward = +y.
WEST = np.array([[0.0, 1000.0], [0.0, 2000.0], [0.0, 3000.0], [0.0, 4000.0]])
SOUTH = np.array([[1000.0, 0.0], [2000.0, 0.0], [3000.0, 0.0], [4000.0, 0.0]])


def _split(mod):
    data = mod.water_level.data
    return data["bzs"].transpose("time", ...).values, data["bzi"].transpose(
        "time", ...
    ).values


@pytest.mark.parametrize(
    "points, direction",
    [(WEST, (1.0, 0.0)), (SOUTH, (0.0, 1.0))],
    ids=["west_boundary_u", "south_boundary_v"],
)
def test_incoming_wave_goes_to_bzi(tmp_path, points, direction):
    """A wave travelling into the model ends up in bzi; bzs keeps the slow level."""
    times = _times()
    _write_his(tmp_path, times, points, np.array(direction))
    mod = _detail(tmp_path, points)
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
        filter_incoming=True,
    )
    bzs, bzi = _split(mod)
    eta = _wave(times)
    assert np.allclose(bzi, eta[:, None], atol=0.03 * AMP)
    assert np.allclose(bzs, 0.3, atol=0.03 * AMP)


@pytest.mark.parametrize(
    "points, direction",
    [(WEST, (-1.0, 0.0)), (SOUTH, (0.0, -1.0))],
    ids=["west_boundary_outgoing", "south_boundary_outgoing"],
)
def test_outgoing_wave_not_in_bzi(tmp_path, points, direction):
    """A wave travelling out of the model gives (almost) no incoming signal."""
    times = _times()
    _write_his(tmp_path, times, points, np.array(direction))
    mod = _detail(tmp_path, points)
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
        filter_incoming=True,
    )
    _, bzi = _split(mod)
    assert np.abs(bzi).max() < 0.03 * AMP


def test_alongshore_velocity_ignored(tmp_path):
    """Velocity parallel to the boundary does not create an incoming signal."""
    times = _times()
    _write_his(tmp_path, times, WEST, np.array([0.0, 1.0]))
    mod = _detail(tmp_path, WEST)
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
        filter_incoming=True,
    )
    _, bzi = _split(mod)
    # no normal flux: incoming = half the elevation (standing-wave split)
    assert np.allclose(bzi, 0.5 * _wave(times)[:, None], atol=0.03 * AMP)


@pytest.mark.parametrize("zb_1d", [True, False], ids=["point_zb_1d", "point_zb_2d"])
def test_point_zb_shapes(tmp_path, zb_1d):
    """point_zb may be (stations,) or (time, stations)."""
    times = _times()
    _write_his(tmp_path, times, WEST, np.array([1.0, 0.0]), zb_1d=zb_1d)
    mod = _detail(tmp_path, WEST)
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
        filter_incoming=True,
    )
    assert "bzi" in mod.water_level.data


def test_unfiltered_has_no_bzi(tmp_path):
    """Without filter_incoming the total level goes to bzs and there is no bzi."""
    times = _times()
    _write_his(tmp_path, times, WEST, np.array([1.0, 0.0]))
    mod = _detail(tmp_path, WEST)
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
    )
    data = mod.water_level.data
    assert "bzi" not in data
    assert np.allclose(
        data["bzs"].transpose("time", ...).values[:, 0], 0.3 + _wave(times)
    )


def test_write_ascii_bzi(tmp_path):
    """ASCII route: bzs and bzi files are written and bzifile is set in the config."""
    times = _times()
    _write_his(tmp_path, times, WEST, np.array([1.0, 0.0]))
    mod = _detail(tmp_path, WEST)
    mod.config.set("tref", times[0].to_pydatetime())
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
        filter_incoming=True,
        bc_path=str(mod.root.path),
    )
    root = Path(mod.root.path)
    assert mod.config.get("bzifile")
    bzi = np.loadtxt(root / Path(mod.config.get("bzifile")).name)
    assert bzi.shape == (len(times), 1 + len(WEST))
    assert np.allclose(bzi[:, 1], _wave(times), atol=0.03 * AMP + 1e-3)


def test_write_netcdf_zi(tmp_path):
    """NetCDF route: zs (slow level) and zi (incoming) end up in netbndbzsbzifile."""
    times = _times()
    _write_his(tmp_path, times, WEST, np.array([1.0, 0.0]))
    mod = _detail(tmp_path, WEST)
    mod.config.set("tref", times[0].to_pydatetime())
    mod.config.set("netbndbzsbzifile", "sfincs_netbndbzsbzifile.nc")
    nest2_sfincs_in_sfincs(
        make_mock_overall_sfincs(tmp_path),
        mod,
        obs_point_prefix=PREFIX,
        output_path=str(tmp_path),
        filter_incoming=True,
        bc_path=str(mod.root.path),
    )
    root = Path(mod.root.path)
    ds = xr.open_dataset(root / mod.config.get("netbndbzsbzifile"), decode_times=False)
    for var in ("zs", "zi", "x", "y", "time"):
        assert var in ds, f"{var} missing in netbndbzsbzi file: {list(ds.variables)}"
    zi = ds["zi"].transpose("time", ...).values
    assert np.allclose(zi[:, 0], _wave(times), atol=0.03 * AMP + 1e-3)
    assert not (root / "sfincs.bzi").exists()
