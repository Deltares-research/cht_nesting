# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_hurrywave_in_hurrywave function for handling nesting of HurryWave models within HurryWave models.
"""

import os
import xarray as xr
import pandas as pd
from typing import Optional, Any

def nest2_hurrywave_in_hurrywave(overall: Any,
                                 detail: Any,
                                 obs_point_prefix: Optional[str] = None,
                                 output_path: Optional[str] = None,
                                 output_file: Optional[str] = None,
                                 bc_path: Optional[str] = None,
                                 **kwargs) -> None:
    """
    Nest a HurryWave model within another HurryWave model.

    Parameters:
    overall (Any): The overall HurryWave model.
    detail (Any): The detailed HurryWave model.
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
        output_file = "hurrywave_sp2.nc"

    file_name = os.path.join(output_path, output_file)
    
    detail.boundary_conditions.forcing = "spectra"

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    stations = ddd.station_name.values
    all_stations = [str(st.strip())[2:-1] for st in stations]

    if len(detail.boundary_conditions.gdf) == 0:
        detail.input.variables.bndfile = "hurrywave.bnd"
        detail.boundary_conditions.read_boundary_points()

    point_names = []
    if len(detail.boundary_conditions.gdf) > 0:
        for ind, row in detail.boundary_conditions.gdf.iterrows():
            # Find required boundary points        
            point_names.append(obs_point_prefix + "_" + row["name"])                    
    else:
        point_names = all_stations.copy()
        
    times = ddd.point_spectrum2d.coords["time"].values
    sigma = ddd.point_spectrum2d.coords["sigma"].values
    theta = ddd.point_spectrum2d.coords["theta"].values

    ireq = []    
    for ip, point in enumerate(point_names):
        for ist, st in enumerate(all_stations):
            if point.lower() == st.lower():
                ireq.append(ist)            
                break

    for ind, row in detail.boundary_conditions.gdf.iterrows():
        sp2 = ddd.point_spectrum2d.values[:, ireq[ind], :, :]

        ds = xr.Dataset(
            data_vars=dict(point_spectrum2d=(["time", "theta", "sigma"], sp2)),
            coords=dict(time=times, theta=theta, sigma=sigma)
        )
        detail.boundary_conditions.gdf.loc[ind, "spectra"] = ds.to_array()

    if bc_path is not None:
        detail.boundary_conditions.write_boundary_conditions_spectra(file_name=os.path.join(bc_path, detail.input.variables.bspfile))