# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_beware_in_hurrywave function for handling nesting of BEWARE models within HurryWave models.
"""

import os
import xarray as xr
import pandas as pd
from typing import Optional, Any

def nest2_beware_in_hurrywave(overall: Any,
                              detail: Any,
                              output_path: Optional[str] = None,
                              output_file: Optional[str] = None,
                              bc_path: Optional[str] = None,
                              **kwargs) -> None:
    """
    Nest a BEWARE model within a HurryWave model.

    Parameters:
    overall (Any): The overall HurryWave model.
    detail (Any): The detailed BEWARE model.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.
    """
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path
    if not output_file:
        output_file = "hurrywave_his.nc"
        
    file_name = os.path.join(output_path, output_file)

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    stations = ddd.station_name.values
    all_stations = [str(st.strip())[2:-1] for st in stations]

    point_names = []
    if detail.wave_boundary_point:
        # Find required boundary points        
        for point in detail.wave_boundary_point:
            point_names.append(detail.name + "_" + point.name)                    
    else:
        point_names = all_stations.copy()
        
    times = ddd.point_hm0.coords["time"].values

    ireq = []
    for ip, point in enumerate(point_names):
        for ist, st in enumerate(all_stations):
            if point.lower() == st.lower():
                ireq.append(ist)
                break

    for ip, point in enumerate(detail.wave_boundary_point):
        hm0 = ddd.point_hm0.values[:, ireq[ip]]
        tp = ddd.point_tp.values[:, ireq[ip]]
        wavdir = ddd.point_wavdir.values[:, ireq[ip]]
        dirspr = ddd.point_dirspr.values[:, ireq[ip]]

        df = pd.DataFrame(index=times)
        df.insert(0, "hm0", hm0)
        df.insert(1, "tp", tp)
        df.insert(2, "wavdir", wavdir)
        df.insert(3, "dirspr", dirspr)

        point.data = df

    if bc_path is not None:
        detail.write_wave_boundary_conditions(path=bc_path)