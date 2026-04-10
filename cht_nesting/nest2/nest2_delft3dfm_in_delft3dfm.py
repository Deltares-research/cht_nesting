"""Nest 2 script for nesting Delft3D-FM within another Delft3D-FM model.

Reads water level time series from the overall Delft3D-FM ``flow_his.nc`` output
and sets them as boundary conditions on the detail Delft3D-FM model.
"""

from typing import Any, Optional

import pandas as pd


def nest2_delft3dfm_in_delft3dfm(
    overall: Any,
    detail: Any,
    output_path: Optional[str] = None,
    output_file: str = "flow_his.nc",
    boundary_water_level_correction: float = 0.0,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a Delft3DFM model within another Delft3DFM model.

    Reads water level time series from the overall model's ``flow_his.nc``
    output via ``overall.read_timeseries_output()`` and applies them as
    boundary conditions on each boundary of the detail model.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3D-FM model.  Must expose
        ``read_timeseries_output(name_list, path, file_name)``.
    detail : Delft3DFM
        The fine Delft3D-FM model that receives the boundary conditions.
    output_path : str, optional
        Directory containing the overall model's output files.
    output_file : str, optional
        Name of the overall his output file.  Defaults to ``"flow_his.nc"``.
    boundary_water_level_correction : float, optional
        Uniform offset (metres) added to every boundary time series.
        Defaults to ``0.0``.
    bc_path : str, optional
        If provided, write boundary conditions to disk.
    **kwargs : Any
        Extra keyword arguments are silently ignored.

    Returns
    -------
    None
    """
    for ind, bnd in enumerate(detail.boundary):
        point_names = []
        for ip, point in enumerate(bnd.point):
            point_names.append(f"{detail.name}_{point.name}")

        # Return DataFrame bzs
        bzs = overall.read_timeseries_output(
            name_list=point_names, path=output_path, file_name=output_file
        )

        ts = bzs.index
        for ip, point in enumerate(bnd.point):
            point.data = (
                pd.Series(bzs.iloc[:, ip].values, index=ts)
                + boundary_water_level_correction
            )

    if bc_path is not None:
        detail.write_flow_boundary_conditions(path=bc_path)
