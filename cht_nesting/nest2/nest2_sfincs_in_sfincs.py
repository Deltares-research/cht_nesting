"""
Nest 2 script for nesting a hydromt_sfincs SfincsModel within another SfincsModel.

Reads water level time series from the overall SfincsModel's ``sfincs_his.nc``
output and sets them as water level boundary conditions on the detail SfincsModel.
Station names in the his-file are stored as fixed-length bytes by the SFINCS
Fortran kernel and are decoded automatically.
"""

import os
from typing import Any, List, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr
from cht_tide.tide_predict import predict
from cht_utils.fileio.deltares_ini import IniStruct
from cht_utils.physics.waves import split_waves_guza
from pyproj import CRS, Transformer


def nest2_sfincs_in_sfincs(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    boundary_water_level_correction: float = 0.0,
    return_maximum: bool = False,
    filter_incoming: bool = False,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> Union[Any, pd.Series]:
    """Nest a hydromt_sfincs SfincsModel within another hydromt_sfincs SfincsModel.

    Reads water level time series (``point_zs``) from the overall model's
    ``sfincs_his.nc`` output and applies them as water level boundary conditions
    on the detail SfincsModel via ``detail.water_level.set_timeseries()``.

    An optional incoming-wave filter (Guza method) can be applied using the
    velocity and bed-level variables also present in the his-file.

    Parameters
    ----------
    overall : hydromt_sfincs.SfincsModel
        The coarse SfincsModel whose ``sfincs_his.nc`` is read.
    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel that receives the water level boundary conditions.
    obs_point_prefix : str, optional
        Prefix prepended to each detail boundary point name when looking up
        stations in the overall his-file.  Defaults to ``detail.name``.
    output_path : str, optional
        Directory containing the overall model's output files.  Defaults to
        ``str(overall.root.path)``.
    output_file : str, optional
        Name of the SFINCS his output file.  Defaults to ``"sfincs_his.nc"``.
    boundary_water_level_correction : float, optional
        Uniform offset (metres) added to every boundary time series.  Useful
        for datum corrections.  Defaults to ``0.0``.
    return_maximum : bool, optional
        When ``True``, return a :class:`~pandas.Series` of peak (maximum)
        water levels across all boundary points instead of the BC component.
        Defaults to ``False``.
    filter_incoming : bool, optional
        When ``True``, split the overall water level into a slowly varying level
        (mean + linear trend, written as ``bzs``) and the incoming wave (``bzi``)
        with the Guza method, using the velocity component normal to the detail
        boundary (positive into the detail model). SFINCS then applies a
        weakly-reflective (Riemann) boundary that lets outgoing waves leave.
        Requires ``point_u``, ``point_v`` and ``point_zb`` in the his-file; the
        boundary normals come from the detail boundary points and its active
        cells. Defaults to ``False``.
    bc_path : str, optional
        If provided, write the boundary conditions to disk by calling
        ``detail.water_level.write()`` after setting.  The value of *bc_path*
        is not used directly — writing goes to the model's configured root.
        With ``filter_incoming``: if the detail model uses ``netbndbzsbzifile``,
        ``bzi`` is written as variable ``zi`` in that netcdf file; otherwise an
        ASCII bzi file is written and ``bzifile`` is set in the detail config
        (write the config yourself).
    **kwargs : Any
        Extra keyword arguments are accepted and silently ignored so that the
        function can be called uniformly through the :func:`nest2` dispatcher.

    Returns
    -------
    hydromt_sfincs water_level component or pandas.Series
        The ``detail.water_level`` component (default) or, when
        *return_maximum* is ``True``, a :class:`~pandas.Series` with the
        peak water level for each time step across all boundary points.

    Raises
    ------
    ValueError
        If a required observation point name is not found among the stations
        in the overall his-file.
    """
    if not output_path:
        output_path = str(overall.root.path)
    if not obs_point_prefix:
        obs_point_prefix = detail.name
    if not output_file:
        output_file = "sfincs_his.nc"

    # The name is the index with 4 leading zeros
    point_names: List[str] = [
        f"{obs_point_prefix}_{str(ind + 1).zfill(4)}"
        for ind, _row in detail.water_level.gdf.iterrows()
    ]

    # Open the overall his-file and decode station names (bytes or str).
    fn_his: str = os.path.join(output_path, output_file)
    ds = xr.open_dataset(fn_his)
    all_stations: List[str] = [
        s.decode("utf-8").strip()
        if isinstance(s, (bytes, np.bytes_))
        else str(s).strip()
        for s in ds["station_name"].values
    ]

    # Map each required point name to a station index (case-insensitive).
    ireq: List[int] = []
    for pname in point_names:
        matched = False
        for ist, st in enumerate(all_stations):
            if pname.lower() == st.lower():
                ireq.append(ist)
                matched = True
                break
        if not matched:
            raise ValueError(
                f"Observation point '{pname}' not found in '{fn_his}'. "
                f"Available stations: {all_stations}"
            )

    times = pd.DatetimeIndex(ds["time"].values)
    bzs = pd.DataFrame(ds["point_zs"].values[:, ireq], index=times)
    bzi: Optional[np.ndarray] = None

    # Optional Guza incoming-wave filter. The total water level is split into a
    # slowly varying level (mean + linear trend, -> bzs) and the incoming wave
    # (-> bzi), using the velocity component normal to the detail boundary,
    # positive into the detail model.
    if filter_incoming:
        zb_all = np.asarray(ds["point_zb"].values)
        # point_zb is (stations,) in current SFINCS his files, (time, stations) in older ones
        bzb = (zb_all[0] if zb_all.ndim == 2 else zb_all)[ireq]
        u = ds["point_u"].values[:, ireq]
        v = ds["point_v"].values[:, ireq]
        gdf = detail.water_level.gdf
        xy = np.column_stack([gdf.geometry.x.values, gdf.geometry.y.values])
        normals = _inward_normals(detail, xy)
        normals = _normals_in_overall_frame(xy, normals, detail.crs, _his_crs(ds))

        bzi = np.zeros((times.size, len(ireq)))
        for icol in range(len(ireq)):
            un = u[:, icol] * normals[icol, 0] + v[:, icol] * normals[icol, 1]
            z = bzs.iloc[:, icol].values
            df_col = pd.DataFrame({"time": times, "z": z, "u": un})
            df_col = split_waves_guza(df_col, float(bzb[icol]))
            # zin = incoming + slow, zout = outgoing + slow, z = incoming + outgoing + slow
            slow = df_col["zin"].values + df_col["zout"].values - z
            bzs.iloc[:, icol] = slow
            bzi[:, icol] = df_col["zin"].values - slow

    ds.close()

    # Astronomic correction from an optional .cor file.
    vcor: Union[float, pd.Series] = 0.0
    corfile: Optional[str] = detail.config.get("corfile")
    if corfile:
        vcor = get_vcor(os.path.join(str(detail.root.path), corfile), times)

    # Build the BC DataFrame: integer GDF-index columns matched to time index.
    df_bcs = pd.DataFrame(index=times)
    for icol, (ind, _row) in enumerate(detail.water_level.gdf.iterrows()):
        df_bcs[ind] = bzs.iloc[:, icol].values + boundary_water_level_correction + vcor

    detail.water_level.set_timeseries(df_bcs)

    data = detail.water_level.data
    if bzi is not None:
        data["bzi"] = (("time", "index"), bzi)
        data["bzi"] = data["bzi"].transpose(*data["bzs"].dims)
    elif "bzi" in data:
        # unfiltered nesting: no incoming-wave signal
        detail.water_level._data = data.drop_vars("bzi")

    if bc_path is not None:
        detail.water_level.write()
        if bzi is not None and not detail.config.get("netbndbzsbzifile"):
            # ASCII route: hydromt_sfincs writes bnd/bzs only, the bzi file is written here.
            # (With netbndbzsbzifile, bzi is written as variable zi in the netcdf file.)
            _write_bzi_ascii(detail)

    if return_maximum:
        total = bzs if bzi is None else bzs + bzi
        return (total + boundary_water_level_correction + vcor).max(axis=1)
    return detail.water_level


def _write_bzi_ascii(detail: Any) -> None:
    """Write the incoming-wave signal ``bzi`` of ``detail.water_level`` to the ASCII bzi file.

    Uses the same writer and time convention (seconds since tref) as the bzs file and sets
    ``bzifile`` in the detail config (default ``sfincs.bzi``).
    """
    from hydromt_sfincs import writers

    path = detail.config.get_set_file_variable("bzifile", default="sfincs.bzi")
    df = detail.water_level.data["bzi"].transpose("time", ...).to_pandas()
    writers.write_timeseries(path, df, detail.config.get("tref"))


def _his_crs(ds: xr.Dataset) -> Optional[CRS]:
    """CRS of the overall model from the ``crs`` variable of its his file (None if absent)."""
    if "crs" not in ds:
        return None
    var = ds["crs"]
    for key in ("epsg_code", "epsg"):
        if key in var.attrs:
            return CRS.from_user_input(var.attrs[key])
    try:
        return CRS.from_epsg(int(var.values))
    except Exception:  # noqa: BLE001
        return None


def _active_centre(detail: Any, xy: np.ndarray) -> np.ndarray:
    """A point inside the detail model: centroid of its active cells (fallbacks: bounds, xy)."""
    try:
        if detail.config.get("qtrfile") is not None:
            data = detail.quadtree_grid.data
            cells = np.asarray(data.grid.face_coordinates)
            mask = np.asarray(data["mask"].values) > 0
        else:
            data = detail.grid.data
            xx, yy = np.meshgrid(data["x"].values, data["y"].values)
            cells = np.column_stack([xx.ravel(), yy.ravel()])
            mask = np.asarray(data["msk"].values).ravel() > 0
        if mask.any():
            return cells[mask].mean(axis=0)
    except Exception:  # noqa: BLE001
        pass
    try:
        x0, y0, x1, y1 = detail.bounds
        return np.array([0.5 * (x0 + x1), 0.5 * (y0 + y1)])
    except Exception:  # noqa: BLE001
        return xy.mean(axis=0)


def _inward_normals(detail: Any, xy: np.ndarray) -> np.ndarray:
    """Unit normals of the boundary polyline at each point, pointing into the detail model.

    The points are split into separate open-boundary segments where the spacing jumps
    (> 3 x the median spacing). Per segment the normal is perpendicular to the local
    tangent; its sign is chosen so that, summed over the segment, it points towards the
    centre of the active detail domain. A single isolated point gets the direction to
    that centre.
    """
    centre = _active_centre(detail, xy)
    normals = np.zeros_like(xy, dtype=float)
    if len(xy) == 1:
        d = centre - xy[0]
        normals[0] = d / max(np.linalg.norm(d), 1e-12)
        return normals
    spacing = np.hypot(*np.diff(xy, axis=0).T)
    gap = spacing > 3.0 * np.median(spacing)
    segment = np.concatenate([[0], np.cumsum(gap)])
    for s in np.unique(segment):
        idx = np.where(segment == s)[0]
        if idx.size == 1:
            d = centre - xy[idx[0]]
            normals[idx[0]] = d / max(np.linalg.norm(d), 1e-12)
            continue
        t = np.gradient(xy[idx], axis=0)
        t /= np.linalg.norm(t, axis=1)[:, None]
        n = np.column_stack([-t[:, 1], t[:, 0]])
        if np.sum((n * (centre - xy[idx])).sum(axis=1)) < 0.0:
            n = -n
        normals[idx] = n
    return normals


def _normals_in_overall_frame(
    xy: np.ndarray, normals: np.ndarray, detail_crs: Any, overall_crs: Optional[CRS]
) -> np.ndarray:
    """Express normals (detail CRS) in the overall model's velocity frame.

    SFINCS his velocities are x/y components of the overall grid: east/north for a
    geographic model. The detail normal is mapped by transforming a short step along it,
    which accounts for grid convergence between the two CRSs.
    """
    if overall_crs is None or detail_crs is None or CRS(detail_crs) == overall_crs:
        return normals
    step = 10.0 if not CRS(detail_crs).is_geographic else 1.0e-4
    tr = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)
    x0, y0 = tr.transform(xy[:, 0], xy[:, 1])
    x1, y1 = tr.transform(
        xy[:, 0] + step * normals[:, 0], xy[:, 1] + step * normals[:, 1]
    )
    dx, dy = np.asarray(x1) - x0, np.asarray(y1) - y0
    if overall_crs.is_geographic:
        dx = dx * np.cos(np.radians(y0))
    out = np.column_stack([dx, dy])
    return out / np.linalg.norm(out, axis=1)[:, None]


def get_vcor(corfile: str, times: pd.DatetimeIndex) -> pd.Series:
    """Compute an astronomic water level correction from a ``.cor`` file.

    Parameters
    ----------
    corfile : str
        Path to the ``.cor`` file (Deltares INI format) containing tidal
        component amplitudes and phases.
    times : pandas.DatetimeIndex
        Time axis for which to evaluate the tidal prediction.

    Returns
    -------
    pandas.Series
        Predicted astronomic water levels at each time in *times*.
    """
    d = IniStruct(filename=corfile)
    astro = d.section[0].data
    names: List[str] = []
    amp: List[float] = []
    phi: List[float] = []
    for icmp, cmp in enumerate(astro.index):
        names.append(cmp)
        amp.append(astro[1][icmp])
        phi.append(astro[2][icmp])
    df = pd.DataFrame(
        {"amplitude": pd.Series(amp), "phase": pd.Series(phi)},
        index=pd.Series(names, name="component"),
    )
    return predict(df, times)
