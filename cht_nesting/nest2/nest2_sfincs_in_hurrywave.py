# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_sfincs_in_hurrywave function for handling nesting of SFINCS models within HurryWave models.
"""

import os
import xarray as xr
import pandas as pd
from typing import Optional, Any

def nest2_sfincs_in_hurrywave(overall: Any,
                              detail: Any,
                              obs_point_prefix: Optional[str] = None,
                              output_path: Optional[str] = None,
                              output_file: Optional[str] = None,
                              bc_path: Optional[str] = None,
                              **kwargs) -> None:
    """
    Nest a SFINCS model within a HurryWave model.

    Parameters:
    overall (Any): The overall HurryWave model.
    detail (Any): The detailed SFINCS model.
    obs_point_prefix (Optional[str]): The prefix for observation points. Default is None.
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
    print("Nesting in " + file_name)

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    stations = ddd.station_name.values
    all_stations = [str(st.strip())[2:-1] for st in stations]

    point_names = []
    if len(detail.snapwave.boundary_conditions.gdf) > 0:
        # Find required boundary points        
        for ind, point in detail.snapwave.boundary_conditions.gdf.iterrows():
            point_names.append(obs_point_prefix + "_" + point["name"])        
    else:
        point_names = all_stations.copy()
        
    times = ddd.point_hm0.coords["time"].values

    ireq = []
    for ip, point in enumerate(point_names):
        for ist, st in enumerate(all_stations):
            if point.lower() == st.lower():
                ireq.append(ist)
                break

    for ip, point in detail.snapwave.boundary_conditions.gdf.iterrows():
        hm0 = ddd.point_hm0.values[:, ireq[ip]]
        tp = ddd.point_tp.values[:, ireq[ip]]
        wavdir = ddd.point_wavdir.values[:, ireq[ip]]
        dirspr = ddd.point_dirspr.values[:, ireq[ip]]
        # Limit direction spreading to values between 4.0 and 55.0 degrees
        dirspr = np.clip(dirspr, 4.0, 55.0)
        # Limit period to values between 1.0 and 25.0 seconds
        tp = np.clip(tp, 1.0, 25.0)

        df = pd.DataFrame(index=times)
        df.insert(0, "hs", hm0)
        df.insert(1, "tp", tp)
        df.insert(2, "wd", wavdir)
        df.insert(3, "ds", dirspr)

        detail.snapwave.boundary_conditions.gdf.at[ip, "timeseries"] = df

    if bc_path is not None:
        detail.snapwave.boundary_conditions.write_boundary_conditions_timeseries()