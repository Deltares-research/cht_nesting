# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_xbeach_in_sfincs function for handling nesting of XBeach models within SFINCS models.
"""

import os
import pandas as pd
from typing import Optional, Any

def nest2_xbeach_in_sfincs(overall: Any,
                           detail: Any,
                           output_path: Optional[str] = None,
                           output_file: Optional[str] = None,
                           boundary_water_level_correction: float = 0,
                           return_maximum: bool = False,
                           bc_path: Optional[str] = None,
                           **kwargs) -> Any:
    """
    Nest an XBeach model within a SFINCS model.

    Parameters:
    overall (Any): The overall SFINCS model.
    detail (Any): The detailed XBeach model.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    boundary_water_level_correction (float): The correction to apply to the boundary water levels. Default is 0.
    return_maximum (bool): Whether to return the maximum values. Default is False.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.

    Returns:
    Any: The boundary points or maximum values.
    """
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path
        
    if overall.input.outputformat[0:3] == "bin":
        # ascii output        
        if not output_file:
            output_file = "zst.txt"
    else:
        # netcdf        
        if not output_file:
            output_file = "sfincs_his.nc"
    
    point_names = [detail.name + "_" + point.name for point in detail.flow_boundary_point]
    zstfile = os.path.join(output_path, output_file)

    # Return DataFrame bzs
    bzs = overall.output.read_his_file(station=point_names,
                                       parameter="point_zs",
                                       file_name=zstfile)

    # Interpolate on desired format for XBeach forcing
    bzs_resampled = bzs.resample('10min').mean()
    bzs_interpolated = bzs_resampled.interpolate(method='linear')
    bzs_filtered = bzs_interpolated[detail.tref:detail.tstop]
    
    ts = bzs_filtered.index
    for icol, point in enumerate(detail.flow_boundary_point):
        point.data = pd.Series(bzs_filtered.iloc[:, icol].values, index=ts) + boundary_water_level_correction

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