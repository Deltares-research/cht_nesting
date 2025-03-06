# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_sfincs_in_beware function for handling nesting of SFINCS models within BEWARE models.
"""

import os
import xarray as xr
import pandas as pd
import numpy as np
import datetime
from pyproj import Transformer
from typing import Optional, Any
from cht_sfincs.sfincs import FlowBoundaryPoint, WaveMakerForcingPoint

def nest2_sfincs_in_beware(overall: Any,
                           detail: Any,
                           output_path: Optional[str] = None,
                           output_file: Optional[str] = None,
                           boundary_water_level_correction: float = 0,
                           option: Optional[str] = None,
                           bc_path: Optional[str] = None,
                           **kwargs) -> None:
    """
    Nest a SFINCS model within a BEWARE model.

    Parameters:
    overall (Any): The overall BEWARE model.
    detail (Any): The detailed SFINCS model.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    boundary_water_level_correction (float): The correction to apply to the boundary water levels. Default is 0.
    option (Optional[str]): The option for nesting ("flow" or "wave"). Default is None.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.
    """
    if not output_file:
        output_file = "beware_his.nc"
    # Get bounding box for sfincs model
    # Convert bbox to beware crs

    x_range, y_range = detail.bounding_box(crs=overall.crs)
    dx = (x_range[1] - x_range[0]) / 10
    dy = (y_range[1] - y_range[0]) / 10
    x_range[0] = x_range[0] - dx
    x_range[1] = x_range[1] + dx
    y_range[0] = y_range[0] - dy
    y_range[1] = y_range[1] + dy
    
    # Read BEWARE offshore locations
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path
        
    file_name = os.path.join(output_path, output_file)

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    
    if option == "flow":
    
        xb = ddd.x_off.values
        yb = ddd.y_off.values

        # Find beware locations in bounding box
        inear = np.where((xb > x_range[0]) & (xb < x_range[1]) & (yb > y_range[0]) & (yb < y_range[1]))    
        xb = xb[inear]
        yb = yb[inear]
        nb = xb.size
        
        # Clear existing flow boundary points
        detail.flow_boundary_point = []
    
        # Convert to coordinate system of detail model
        transformer = Transformer.from_crs(overall.crs,
                                           detail.crs,
                                           always_xy=True)
        
        for ip in range(nb):
            name = str(ip + 1).zfill(4)        
            x, y = transformer.transform(xb[ip], yb[ip])
            point = FlowBoundaryPoint(x,
                                      y,
                                      name=name)
            detail.flow_boundary_point.append(point)
        
        # Extract data and set water level boundary conditions
        tref = datetime.datetime(1970, 1, 1)
        tsec = ddd.time.values  # array of int64
        times = tref + tsec * datetime.timedelta(seconds=1)

        for ip, point in enumerate(detail.flow_boundary_point):
            point.data = pd.Series(ddd.WL.values[inear[0][ip], :], index=times) + pd.Series(ddd.R2_setup.values[inear[0][ip], :], index=times) + boundary_water_level_correction

        if bc_path is not None:    
            detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))

    elif option == "wave":

        xb = ddd.x_coast.values
        yb = ddd.y_coast.values

        # Find beware locations in bounding box
        inear = np.where((xb > x_range[0]) & (xb < x_range[1]) & (yb > y_range[0]) & (yb < y_range[1]))    
        xb = xb[inear]
        yb = yb[inear]
        nb = xb.size

        # Clear existing flow boundary points
        detail.wavemaker_forcing_point = []
    
        # Find wave boundary forcing at intersection wfp file and beware transects
        # Load sfincs.wfp, find intersect with wfp and beware transect.

        # Convert to coordinate system of detail model
        transformer = Transformer.from_crs(overall.crs,
                                           detail.crs,
                                           always_xy=True)
        
        for ip in range(nb):
            name = str(ip + 1).zfill(4)        
            x, y = transformer.transform(xb[ip], yb[ip])
            point = WaveMakerForcingPoint(x,
                                          y,
                                          name=name)
            detail.wavemaker_forcing_point.append(point)
        
        # Extract data and set water level boundary conditions
        tref = datetime.datetime(1970, 1, 1)
        tsec = ddd.time.values  # array of int64
        times = tref + tsec * datetime.timedelta(seconds=1)
        
        for ip, point in enumerate(detail.wavemaker_forcing_point):

            df = pd.DataFrame()
            df["hm0_ig"] = ddd.hm0_ig.values[inear[0][ip], :]
            df["tp_ig"] = ddd.tp_ig.values[inear[0][ip], :]
            df["setup"] = ddd.setup.values[inear[0][ip], :]
            df["time"] = times
            df = df.set_index("time")
            
            df["hm0_ig"] = df["hm0_ig"].replace(np.nan, 0.1)
            df["tp_ig"] = df["tp_ig"].replace(np.nan, 60.0)
            df["setup"] = df["setup"].replace(np.nan, 0.0)
            
            point.data = df

    ddd.close()

    if bc_path is not None:
        detail.write_whi_file(file_name=os.path.join(bc_path, detail.input.whifile))
        detail.write_wti_file(file_name=os.path.join(bc_path, detail.input.wtifile))