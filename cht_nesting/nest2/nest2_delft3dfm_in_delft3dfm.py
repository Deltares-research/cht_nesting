# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_delft3dfm_in_delft3dfm function for handling nesting in Delft3DFM models.
"""

import os
import pandas as pd
from typing import Optional, Any

def nest2_delft3dfm_in_delft3dfm(overall: Any,
                                 detail: Any,
                                 output_path: Optional[str] = None,
                                 output_file: str = "flow_his.nc",
                                 boundary_water_level_correction: float = 0.0,
                                 bc_path: Optional[str] = None,
                                 **kwargs) -> None:
    """
    Nest a Delft3DFM model within another Delft3DFM model.

    Parameters:
    overall (Any): The overall Delft3DFM model.
    detail (Any): The detailed Delft3DFM model.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (str): The name of the output file. Default is "flow_his.nc".
    boundary_water_level_correction (float): The correction to apply to the boundary water levels. Default is 0.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.
    """
    for ind, bnd in enumerate(detail.boundary):
        point_names = []
        for ip, point in enumerate(bnd.point):
            point_names.append(detail.name + "_" + point.name)

        # Return DataFrame bzs
        bzs = overall.read_timeseries_output(name_list=point_names,
                                             path=output_path,
                                             file_name=output_file)

        ts = bzs.index
        for ip, point in enumerate(bnd.point):
            point.data = pd.Series(bzs.iloc[:, ip].values, index=ts) + boundary_water_level_correction

    if bc_path is not None:
        detail.write_flow_boundary_conditions(path=bc_path)