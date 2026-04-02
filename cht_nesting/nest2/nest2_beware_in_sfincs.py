# -*- coding: utf-8 -*-
"""
Nest 2 script for nesting a hydromt_sfincs SfincsModel within a BEWARE model.

Reads water level (WL + R2_setup) time series from BEWARE output and sets them
as water level boundary conditions on the detail SfincsModel.

Only the ``"flow"`` option is implemented; BEWARE wave (infragravity) output is
not applicable to SfincsModel SnapWave boundary conditions.
"""

import datetime
import os
from typing import Any, List, Optional, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer
from shapely.geometry import Point


def nest2_beware_in_sfincs(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    boundary_water_level_correction: float = 0.0,
    option: Optional[str] = "flow",
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a hydromt_sfincs SfincsModel within a BEWARE model (flow BCs only).

    Finds BEWARE offshore output locations inside (or near) the detail model's
    bounding box, extracts WL + R2_setup time series, and sets them as water
    level boundary conditions on the detail SfincsModel via
    ``detail.water_level.set_timeseries()``.

    The bounding box of the detail model is expanded by 10 % in each direction
    before filtering BEWARE stations, so that points lying on or just outside
    the exact boundary are included.

    Parameters
    ----------
    overall : BEWARE
        The coarse BEWARE model.  Must expose:

        * ``path`` – directory containing the BEWARE output file.
        * ``crs`` – :class:`~pyproj.CRS` of the BEWARE coordinate system.

    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel that receives the water level boundary conditions.
        Must expose:

        * ``bounds`` – ``(xmin, ymin, xmax, ymax)`` tuple in the model's CRS.
        * ``crs`` – :class:`~pyproj.CRS` of the detail model.
        * ``water_level`` – boundary condition component with ``clear()``,
          ``set_locations()``, ``set_timeseries()``, and ``write()`` methods.

    obs_point_prefix : str, optional
        Unused for BEWARE nesting (boundary points are named by index).
        Accepted for interface consistency.
    output_path : str, optional
        Directory containing the BEWARE output files.  Defaults to
        ``overall.path``.
    output_file : str, optional
        Name of the BEWARE output file.  Defaults to ``"beware_his.nc"``.
    boundary_water_level_correction : float, optional
        Uniform offset (metres) added to every boundary time series.  Useful
        for datum corrections.  Defaults to ``0.0``.
    option : str, optional
        Nesting option; only ``"flow"`` is supported.  Defaults to ``"flow"``.
    bc_path : str, optional
        If provided, write the boundary conditions to disk by calling
        ``detail.water_level.write()`` after setting.  The value of *bc_path*
        is not used directly — writing goes to the model's configured root.
    **kwargs : Any
        Extra keyword arguments are accepted and silently ignored so that the
        function can be called uniformly through the :func:`nest2` dispatcher.

    Returns
    -------
    None

    Raises
    ------
    NotImplementedError
        If *option* is anything other than ``"flow"``.
    ValueError
        If no BEWARE output locations are found within the (expanded) bounding
        box of the detail model.
    """
    if not output_path:
        output_path = overall.path
    if not output_file:
        output_file = "beware_his.nc"

    if option != "flow":
        raise NotImplementedError(
            "Only option='flow' is implemented for nest2_beware_in_sfincsmodel. "
            "BEWARE infragravity wave output is not applicable to SfincsModel SnapWave BCs."
        )

    # Transform detail model bounding box to BEWARE CRS.
    xmin, ymin, xmax, ymax = detail.bounds
    transformer_to_overall = Transformer.from_crs(
        detail.crs, overall.crs, always_xy=True
    )
    corners_x: Tuple[float, ...] = transformer_to_overall.transform(
        [xmin, xmax, xmin, xmax], [ymin, ymin, ymax, ymax]
    )[0]
    corners_y: Tuple[float, ...] = transformer_to_overall.transform(
        [xmin, xmax, xmin, xmax], [ymin, ymin, ymax, ymax]
    )[1]
    x_range: List[float] = [min(corners_x), max(corners_x)]
    y_range: List[float] = [min(corners_y), max(corners_y)]

    # Expand bounding box by 10 % to include stations on or near the boundary.
    dx = (x_range[1] - x_range[0]) / 10.0
    dy = (y_range[1] - y_range[0]) / 10.0
    x_range[0] -= dx
    x_range[1] += dx
    y_range[0] -= dy
    y_range[1] += dy

    # Open BEWARE output and find stations within the expanded bounding box.
    file_name: str = os.path.join(output_path, output_file)
    ddd = xr.open_dataset(file_name)

    xb: np.ndarray = ddd.x_off.values
    yb: np.ndarray = ddd.y_off.values
    inear: np.ndarray = np.where(
        (xb > x_range[0]) & (xb < x_range[1]) & (yb > y_range[0]) & (yb < y_range[1])
    )[0]
    xb = xb[inear]
    yb = yb[inear]
    nb: int = xb.size

    if nb == 0:
        ddd.close()
        raise ValueError(
            "No BEWARE output locations found within the detail model bounding box."
        )

    # Transform selected BEWARE locations to the detail model CRS.
    transformer_to_detail = Transformer.from_crs(
        overall.crs, detail.crs, always_xy=True
    )
    xs: np.ndarray
    ys: np.ndarray
    xs, ys = transformer_to_detail.transform(xb, yb)

    # Replace existing water level boundary locations on the detail model.
    names: List[str] = [str(ip + 1).zfill(4) for ip in range(nb)]
    gdf_new = gpd.GeoDataFrame(
        {"name": names},
        geometry=[Point(x, y) for x, y in zip(xs, ys)],
        crs=detail.crs,
    )
    detail.water_level.clear()
    detail.water_level.set_locations(gdf_new)

    # Parse BEWARE time axis (seconds since 1970-01-01 UTC).
    tref = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)
    tsec: np.ndarray = ddd.time.values  # int64 array of epoch seconds
    times = pd.DatetimeIndex([tref + datetime.timedelta(seconds=int(t)) for t in tsec])

    # Build the BC DataFrame: integer GDF-index columns matched to time index.
    df_bcs = pd.DataFrame(index=times)
    for ip in range(nb):
        wl: np.ndarray = ddd.WL.values[inear[ip], :]
        setup: np.ndarray = ddd.R2_setup.values[inear[ip], :]
        df_bcs[ip] = wl + setup + boundary_water_level_correction

    ddd.close()

    detail.water_level.set_timeseries(df_bcs)

    if bc_path is not None:
        detail.water_level.write()
