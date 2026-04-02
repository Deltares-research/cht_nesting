"""Nest 2 script for nesting an XBeach model within a hydromt_sfincs SfincsModel.

Reads water level time series from the overall SfincsModel's ``sfincs_his.nc``
output and sets them as flow boundary conditions on the detail XBeach model.
"""

import os
from typing import Any, List, Optional

import numpy as np
import pandas as pd
import xarray as xr


def nest2_xbeach_in_sfincs(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    boundary_water_level_correction: float = 0.0,
    return_maximum: bool = False,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    """Nest an XBeach model within a hydromt_sfincs SfincsModel.

    Parameters
    ----------
    overall : hydromt_sfincs.SfincsModel
        The coarse SfincsModel whose sfincs_his.nc is read.
    detail : XBeach
        The fine XBeach model that receives the water level boundary conditions.
    obs_point_prefix : str, optional
        Prefix for observation point names.
    output_path : str, optional
        Directory containing the overall model's output files.
    output_file : str, optional
        Name of the SFINCS his output file. Defaults to "sfincs_his.nc".
    boundary_water_level_correction : float, optional
        Uniform offset added to every boundary time series.
    return_maximum : bool, optional
        When True, return peak water levels instead of boundary points.
    bc_path : str, optional
        If provided, write boundary conditions to disk.
    """
    if not output_path:
        output_path = str(overall.root.path)
    if not obs_point_prefix:
        obs_point_prefix = detail.name
    if not output_file:
        output_file = "sfincs_his.nc"

    point_names: List[str] = [
        f"{obs_point_prefix}_{point.name}" for point in detail.flow_boundary_point
    ]
    fn_his = os.path.join(output_path, output_file)

    # Open the overall his-file and decode station names.
    ds = xr.open_dataset(fn_his)
    all_stations: List[str] = [
        s.decode("utf-8").strip()
        if isinstance(s, (bytes, np.bytes_))
        else str(s).strip()
        for s in ds["station_name"].values
    ]

    ireq: List[int] = []
    for pname in point_names:
        for ist, st in enumerate(all_stations):
            if pname.lower() == st.lower():
                ireq.append(ist)
                break

    times = pd.DatetimeIndex(ds["time"].values)
    bzs = pd.DataFrame(ds["point_zs"].values[:, ireq], index=times)
    ds.close()

    # Interpolate on desired format for XBeach forcing
    bzs_resampled = bzs.resample("10min").mean()
    bzs_interpolated = bzs_resampled.interpolate(method="linear")
    bzs_filtered = bzs_interpolated[detail.tref : detail.tstop]

    ts = bzs_filtered.index
    for icol, point in enumerate(detail.flow_boundary_point):
        point.data = (
            pd.Series(bzs_filtered.iloc[:, icol].values, index=ts)
            + boundary_water_level_correction
        )

    # Write boundary conditions
    if bc_path is not None:
        detail.write_flow_boundary_conditions()

    if return_maximum:
        zmax = -999.0
        if len(detail.flow_boundary_point) <= 2:
            for icol, point in enumerate(detail.flow_boundary_point):
                if icol == 1:
                    break
                zx = point.data.max()
                if zx > zmax:
                    zs = point.data
                    zmax = zx
        elif len(detail.flow_boundary_point) == 4:
            for icol, point in enumerate(detail.flow_boundary_point):
                if icol == 2:
                    break
                zx = point.data.max()
                if zx > zmax:
                    zs = point.data
                    zmax = zx
        return zs
    else:
        return detail.flow_boundary_point
