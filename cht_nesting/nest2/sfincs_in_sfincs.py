# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

@author: ormondt
"""

import os
# from pyproj import CRS
# from pyproj import Transformer
import pandas as pd
# import xarray as xr
# import numpy as np
# import glob
# import datetime 
from cht_utils.deltares_ini import IniStruct
from cht_tide.tide_predict import predict
from cht_physics.waves import split_waves_guza

def sfincs_in_sfincs(overall,
                     detail,
                     obs_point_prefix=None,
                     output_path=None,
                     output_file=None,
                     boundary_water_level_correction=0.0,
                     return_maximum=False,
                     filter_incoming=False,
                     bctype="waterlevel",
                     bc_path=None):

    if not output_path:
        # Path of the overall output time series
        output_path = overall.path

    if not obs_point_prefix:
        # Name of obs points
        obs_point_prefix = detail.name    
        
    if overall.input.variables.outputformat[0:3] == "bin":
        # ascii output        
        if not output_file:
            output_file = "zst.txt"
    else:
        # netcdf        
        if not output_file:
            output_file = "sfincs_his.nc"
    
    point_names = []
    # Loop through gdf of boundary conditions
    for ind, point in detail.boundary_conditions.gdf.iterrows():
    # for point in detail.flow_boundary_point:
        point_names.append(obs_point_prefix + "_" + point["name"])                    
    zstfile = os.path.join(output_path, output_file)

    # Return DataFrame bzs
    bzs = overall.output.read_his_file(station=point_names,
                                parameter="zs",
                                file_name=zstfile)
    ts  = bzs.index

    # Check if filtering is needed
    if filter_incoming:

        # Read velocities
        bzu = overall.output.read_his_file(station=point_names,
                                    parameter="point_u",
                                    file_name=zstfile)
        bzv = overall.output.read_his_file(station=point_names,
                                    parameter="point_v",
                                    file_name=zstfile)

        # Read bed levels
        bzb = overall.output.read_his_file(station=point_names,
                                    parameter="zb",
                                    file_name=zstfile)

        # Loop over all points
        for icol, point in detail.boundary_conditions.gdf.iterrows():
            # Make new DataFrame with water level and velocity time series for this point
            df = pd.DataFrame()
            df["time"] = ts
            df["z"]    = bzs.iloc[:,icol].values
            df["u"]    = bzu.iloc[:,icol].values
            df["v"]    = bzv.iloc[:,icol].values
            # Get the bed level at this point
            zb = bzb.iloc[0, icol]
            # Filter out the outgoing waves
            df = split_waves_guza(df, zb)
            # Update water levels in bzs
            bzs.iloc[:,icol] = df["zin"].values
            # Do the same for the velocity


    # Astro correction
    vcor = 0.0
    if detail.input.variables.corfile:        
        vcor = get_vcor(os.path.join(detail.path, detail.input.corfile), ts)

    for icol, point in detail.boundary_conditions.gdf.iterrows():
        df = pd.DataFrame()     
        df["time"] = ts
        df["wl"] = bzs.iloc[:,icol].values + boundary_water_level_correction + vcor
        df = df.set_index("time")
        detail.boundary_conditions.gdf.at[icol, "timeseries"] = df

    # Write bzs file
    if bc_path is not None:
        detail.boundary_conditions.write_boundary_conditions_timeseries()

    # zs_maximum for clustering
    if return_maximum:
        zmax = -999.0
        for icol, point in detail.boundary_conditions.gdf.iterrows():
            zx = point.data.max()
            if zx>zmax:
                zs = point.data
                zmax = zx
        return zs                            
    else:    
        return detail.boundary_conditions

def get_vcor(corfile, times):
    # Add astronomic correction to time series        
    # Read cor file
#        corfile = os.path.join(detail.path, detail.input.corfile)
    d = IniStruct(filename=corfile)
    astro = d.section[0].data
#    times = self.domain.flow_boundary_point[0].data.index
    names = []
    amp   = []
    phi   = []        
    for icmp, cmp in enumerate(astro.index):                
        names.append(cmp)
        amp.append(astro[1][icmp])
        phi.append(astro[2][icmp])        
    df = pd.DataFrame()
    df["component"] = pd.Series(names) 
    df["amplitude"] = pd.Series(amp) 
    df["phase"]     = pd.Series(phi) 
    df = df.set_index("component")
    return predict(df, times)
    