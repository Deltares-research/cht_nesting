"""
Shared pytest fixtures for cht_nesting tests.

Synthetic model data is generated in-memory so no external files are required.
The only external dependency is hydromt_sfincs (for SfincsModel detail fixtures).
"""

from pathlib import Path
from unittest.mock import MagicMock

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import xarray as xr
from shapely.geometry import Point

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

OBS_PREFIX = "nest"
N_POINTS = 2
N_TIMES = 12
EPSG = 32633  # UTM 33N — matches hydromt_sfincs test convention

# Boundary point names on the detail model (stored in GeoDataFrame "name" column)
DETAIL_POINT_NAMES = ["bnd_0001", "bnd_0002"]

# Corresponding station names in the overall model output files.
# The implementation builds names as <prefix>_<zfill4(gdf_index + 1)>,
# so for a 0-based integer GDF index the first two points become
# "nest_0001" and "nest_0002".
OVERALL_STATION_NAMES = [
    f"{OBS_PREFIX}_{str(i + 1).zfill(4)}" for i in range(N_POINTS)
]


# ---------------------------------------------------------------------------
# Time and water level helpers
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def sim_times() -> pd.DatetimeIndex:
    """Hourly datetime index used across all synthetic output files."""
    return pd.date_range("2020-01-01", periods=N_TIMES, freq="h")


@pytest.fixture(scope="session")
def wl_values(sim_times) -> np.ndarray:
    """Reproducible synthetic water level array (time × stations)."""
    rng = np.random.default_rng(42)
    return rng.uniform(-0.5, 2.0, size=(len(sim_times), N_POINTS))


# ---------------------------------------------------------------------------
# Synthetic sfincs_his.nc (used by nest2_sfincsmodel_in_sfincsmodel)
# ---------------------------------------------------------------------------


@pytest.fixture
def sfincs_his_path(tmp_path, sim_times, wl_values) -> Path:
    """
    Write a minimal ``sfincs_his.nc`` to *tmp_path* and return the path.

    Station names are stored as bytes (``|S`` dtype) to match what the
    SFINCS Fortran kernel produces.  All velocity variables are zero;
    bed levels are -5 m so that the incoming-wave filter leaves water
    levels unchanged in the basic test.
    """
    n_t = len(sim_times)
    # Encode station names as fixed-length bytes (Fortran convention)
    station_names = np.array(
        [s.encode("utf-8") for s in OVERALL_STATION_NAMES], dtype="S40"
    )

    ds = xr.Dataset(
        {
            "point_zs": (["time", "stations"], wl_values),
            "point_u": (["time", "stations"], np.zeros((n_t, N_POINTS))),
            "point_v": (["time", "stations"], np.zeros((n_t, N_POINTS))),
            "point_zb": (["time", "stations"], np.full((n_t, N_POINTS), -5.0)),
            "station_name": (["stations"], station_names),
            "station_x": (["stations"], np.array([0.0, 0.0])),
            "station_y": (["stations"], np.array([1000.0, 3000.0])),
            "crs": ([], EPSG),
        },
        coords={"time": sim_times},
    )
    path = tmp_path / "sfincs_his.nc"
    ds.to_netcdf(path)
    return path


# ---------------------------------------------------------------------------
# Synthetic hurrywave_his.nc (used by nest2_hurrywave_in_sfincsmodel)
# ---------------------------------------------------------------------------


@pytest.fixture
def hurrywave_his_path(tmp_path, sim_times) -> Path:
    """
    Write a minimal ``hurrywave_his.nc`` with wave parameters to *tmp_path*.

    Values are deterministic:
    - Hm0:  1.5 m  (constant)
    - Tp:   8.0 s  (constant)
    - Dir: 270.0 ° (constant)
    - Spr:  20.0 ° (constant, within the valid 4–55 ° clipping range)
    """
    n_t = len(sim_times)
    n_s = N_POINTS

    station_names = np.array(
        [s.encode("utf-8") for s in OVERALL_STATION_NAMES], dtype="S40"
    )

    ds = xr.Dataset(
        {
            "point_hm0": (["time", "stations"], np.full((n_t, n_s), 1.5)),
            "point_tp": (["time", "stations"], np.full((n_t, n_s), 8.0)),
            "point_wavdir": (["time", "stations"], np.full((n_t, n_s), 270.0)),
            "point_dirspr": (["time", "stations"], np.full((n_t, n_s), 20.0)),
            "station_name": (["stations"], station_names),
        },
        coords={"time": sim_times},
    )
    path = tmp_path / "hurrywave_his.nc"
    ds.to_netcdf(path)
    return path


# ---------------------------------------------------------------------------
# Synthetic beware_his.nc (used by nest2_beware_in_sfincsmodel)
# ---------------------------------------------------------------------------


@pytest.fixture
def beware_his_path(tmp_path, sim_times) -> Path:
    """
    Write a minimal ``beware_his.nc`` with BEWARE offshore locations and
    water level / setup to *tmp_path*.

    Offshore points are placed at UTM 33N coordinates that fall inside the
    bounding box of the detail model fixture (0–5000 m × 0–5000 m).
    """
    n_t = len(sim_times)
    n_s = N_POINTS

    # Offshore locations inside the detail model bbox
    x_off = np.array([500.0, 4500.0])
    y_off = np.array([1000.0, 3000.0])

    tref = pd.Timestamp("1970-01-01", tz="UTC")
    tsec = ((sim_times - tref) / pd.Timedelta("1s")).values.astype(np.int64)

    ds = xr.Dataset(
        {
            "x_off": (["stations"], x_off),
            "y_off": (["stations"], y_off),
            "WL": (["stations", "time"], np.full((n_s, n_t), 0.5)),
            "R2_setup": (["stations", "time"], np.full((n_s, n_t), 0.1)),
            "time": (["time"], tsec),
        },
    )
    path = tmp_path / "beware_his.nc"
    ds.to_netcdf(path)
    return path


# ---------------------------------------------------------------------------
# Detail SfincsModel fixture (water level BCs)
# ---------------------------------------------------------------------------


@pytest.fixture
def detail_model_wl(tmp_path):
    """
    A minimal hydromt_sfincs ``SfincsModel`` with two water level boundary
    points pre-set.  Points are placed along the western edge of a 5 × 5 km
    UTM 33N grid.

    Yields the model object; the caller is responsible for calling
    ``water_level.set_timeseries`` and checking the result.
    """
    from hydromt_sfincs import SfincsModel

    root = tmp_path / "detail_wl"
    root.mkdir()
    mod = SfincsModel(root=str(root), mode="w+")

    # Minimal regular grid (sets CRS to EPSG 32633)
    mod.grid.create(
        mmax=5,
        nmax=5,
        dx=1000.0,
        dy=1000.0,
        x0=0.0,
        y0=0.0,
        rotation=0.0,
        epsg=EPSG,
    )

    # Add two water level boundary points
    gdf = gpd.GeoDataFrame(
        {"name": DETAIL_POINT_NAMES},
        geometry=[Point(0.0, 1000.0), Point(0.0, 3000.0)],
        crs=f"EPSG:{EPSG}",
    )
    mod.water_level.set_locations(gdf)
    return mod


# ---------------------------------------------------------------------------
# Detail SfincsModel fixture (SnapWave BCs)
# ---------------------------------------------------------------------------


@pytest.fixture
def detail_model_sw(tmp_path):
    """
    A minimal hydromt_sfincs ``SfincsModel`` with two SnapWave boundary
    points pre-set.
    """
    from hydromt_sfincs import SfincsModel

    root = tmp_path / "detail_sw"
    root.mkdir()
    mod = SfincsModel(root=str(root), mode="w+")

    mod.grid.create(
        mmax=5,
        nmax=5,
        dx=1000.0,
        dy=1000.0,
        x0=0.0,
        y0=0.0,
        rotation=0.0,
        epsg=EPSG,
    )

    gdf = gpd.GeoDataFrame(
        {"name": DETAIL_POINT_NAMES},
        geometry=[Point(0.0, 1000.0), Point(0.0, 3000.0)],
        crs=f"EPSG:{EPSG}",
    )
    mod.snapwave_boundary_conditions.set_locations(gdf)
    return mod


# ---------------------------------------------------------------------------
# Mock overall model helpers
# ---------------------------------------------------------------------------


def make_mock_overall_sfincs(root_path: Path) -> MagicMock:
    """Return a mock SfincsModel whose ``root.path`` points to *root_path*."""
    mock = MagicMock(name="overall_sfincsmodel")
    mock.root.path = root_path
    mock.name = "sfincs"
    return mock


def make_mock_delft3dfm(output_path: Path, bzs_df: pd.DataFrame) -> MagicMock:
    """
    Return a mock Delft3DFM model whose ``read_timeseries_output`` returns
    *bzs_df* regardless of the arguments passed.
    """
    mock = MagicMock(name="overall_delft3dfm")
    mock.read_timeseries_output.return_value = bzs_df
    return mock


def make_mock_hurrywave(output_path: Path) -> MagicMock:
    """Return a mock HurryWave model whose ``path`` points to *output_path*."""
    mock = MagicMock(name="overall_hurrywave")
    mock.path = str(output_path)
    return mock


def make_mock_beware(output_path: Path) -> MagicMock:
    """Return a mock BEWARE model whose ``path`` points to *output_path*."""
    from pyproj import CRS

    mock = MagicMock(name="overall_beware")
    mock.path = str(output_path)
    mock.crs = CRS.from_epsg(EPSG)
    return mock
