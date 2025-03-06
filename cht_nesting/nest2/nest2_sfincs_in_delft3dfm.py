# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_sfincs_in_delft3dfm function for handling nesting of SFINCS models within Delft3DFM models.
"""

import os
import pandas as pd
from typing import Optional, Any

def nest2_sfincs_in_delft3dfm(overall: Any,
                              detail: Any,
                              output_path: Optional[str] = None,
                              output_file: Optional[str] = None,
                              boundary_water_level_correction: float = 0,
                              bc_path: Optional[str] = None,
                              **kwargs) -> None:
    """
    Nest a SFINCS model within a Delft3DFM model.

    Parameters:
    overall (Any): The overall Delft3DFM model.
    detail (Any): The detailed SFINCS model.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    boundary_water_level_correction (float): The correction to apply to the boundary water levels. Default is 0.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.
    """
    if not output_file:
        output_file = "flow_his.nc"

    point_names = [detail.name + "_" + point.name for point in detail.flow_boundary_point]
    output_file = os.path.join(output_path, output_file)

    # Return DataFrame bzs
    bzs = overall.read_timeseries_output(name_list=point_names,
                                         path=output_path,
                                         file_name=output_file)

    ts = bzs.index
    for icol, point in enumerate(detail.flow_boundary_point):
        point.data = pd.Series(bzs.iloc[:, icol].values, index=ts) + boundary_water_level_correction

    if bc_path is not None:
        detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))