# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_beware_in_sfincs function for handling nesting of BEWARE models within SFINCS models.
"""

import os
import pandas as pd
import xarray as xr
import numpy as np
import datetime
from pyproj import Transformer
from typing import Optional, Any

def nest2_beware_in_sfincs(overall: Any,
                           detail: Any,
                           output_path: Optional[str] = None,
                           output_file: Optional[str] = None,
                           boundary_water_level_correction: float = 0,
                           option: Optional[str] = None,
                           bc_path: Optional[str] = None,
                           **kwargs) -> None:
    """
    Nest a BEWARE model within a SFINCS model.

    Parameters:
    overall (Any): The overall SFINCS model.
    detail (Any): The detailed BEWARE model.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    boundary_water_level_correction (float): The correction to apply to the boundary water levels. Default is 0.
    option (Optional[str]): The option for nesting ("flow" or "wave"). Default is None.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.
    """
    if option == 'flow':
        if not output_file:
            output_file = "sfincs_his.nc"
        
        point_names = [detail.name + "_" + point.name for point in detail.flow_boundary_point]

        bzs = overall.output.read_his_file(station=point_names,
                                           parameter="point_zs",
                                           file_name=os.path.join(output_path, output_file))

        # Replace -999.0 with zeros. This should not happen, but easy fix for now.
        bzs = bzs.replace(-999.0, 0.0)
        for icol, point in enumerate(detail.flow_boundary_point):
            point.data = pd.Series(bzs.iloc[:, icol].values, index=bzs.index) + boundary_water_level_correction    
            
        if bc_path is not None:
            detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))

    elif option == 'wave':
        if not output_path:
            output_path = overall.path
            
        if not output_file:
            file_name = os.path.join(output_path, "wavh-wave-nest.nc")
        else:
            file_name = os.path.join(output_path, output_file)
        
        file_name_dfm = os.path.join(output_path, "flow_his.nc")
        ddd = xr.open_dataset(file_name_dfm)
        stations = ddd.station_name.values
        all_stations = [str(st.strip())[2:-1] for st in stations]

        point_names = [detail.name + "_" + point.name for point in detail.wave_boundary_point]

        ddd = xr.load_dataset(file_name)
        times = ddd.Hsig.time.values

        for ip, point in enumerate(detail.wave_boundary_point):
            for ist, st in enumerate(all_stations):
                ireq = -1
                if point_names[ip] == st.lower():
                    ireq = ist
                    break

            if ireq > -1:
                hm0 = ddd.Hsig.values[:, ireq]
                tp = ddd.RTpeak.values[:, ireq]
                wavdir = ddd.Dir.values[:, ireq]
                dirspr = ddd.Dspr.values[:, ireq]

                df = pd.DataFrame(index=times)
                df.insert(0, "hm0", hm0)
                df.insert(1, "tp", tp)
                df.insert(2, "wavdir", wavdir)
                df.insert(3, "dirspr", dirspr)

                point.data = df

        if bc_path is not None:
            detail.write_wave_boundary_conditions(path=bc_path)