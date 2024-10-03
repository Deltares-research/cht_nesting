# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

@author: ormondt
"""

import os
from pyproj import CRS
from pyproj import Transformer
import pandas as pd
import xarray as xr
import numpy as np
import glob
import datetime 

def nest1(overall, detail, option=None):
    
    # Returns a list with observation point objects
    
    if overall.type.lower() == "delft3dfm":
        if detail.type.lower() == "delft3dfm":
            nest1_delft3dfm_in_delft3dfm(overall, detail)
        elif detail.type.lower() == "sfincs":
            nest1_sfincs_in_delft3dfm(overall, detail)
        elif detail.type.lower() == "beware":
            nest1_beware_in_delft3dfm(overall, detail)
            
    elif overall.type.lower() == "sfincs":
        if detail.type.lower() == "sfincs":
            nest1_sfincs_in_sfincs(overall, detail)
        elif detail.type.lower() == "xbeach":
            nest1_xbeach_in_sfincs(overall, detail)
        elif detail.type.lower() == "beware":
            nest1_beware_in_sfincs(overall, detail)

    elif overall.type.lower() == "hurrywave":
        if detail.type.lower() == "hurrywave":
            nest1_hurrywave_in_hurrywave(overall, detail)
        elif detail.type.lower() == "xbeach":    
            nest1_xbeach_in_hurrywave(overall, detail)
        elif detail.type.lower() == "sfincs":    
            nest1_sfincs_in_hurrywave(overall, detail)
        elif detail.type.lower() == "beware":
            nest1_beware_in_hurrywave(overall, detail)

    elif overall.type.lower() == "beware":
        if detail.type.lower() == "sfincs":
            # No need to do anything here. BEWARE output points are fixed
            pass

#        elif detail.type == "delft3dfm":
#            obs = nest1_delft3dfm_in_sfincs(overall, detail)

#    return obs    


def nest1_delft3dfm_in_delft3dfm(overall, detail):
    
#    from delft3dfm import ObservationPoint as obspoint
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, bnd in enumerate(detail.boundary):        
        for ip, point in enumerate(bnd.point):
            x, y = transformer.transform(point.geometry.x,
                                         point.geometry.y)
            overall.add_observation_point(x, y, detail.name + "_" + point.name)
    
def nest1_sfincs_in_delft3dfm(overall, detail):
    
#    from delft3dfm import ObservationPoint as obspoint
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.flow_boundary_point):

        name = detail.name + "_" + str(ind + 1).zfill(4)
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
#        obs_list.append(obspoint(x, y, name, crs=overall.crs))
        overall.add_observation_point(x, y, name)

def nest1_beware_in_delft3dfm(overall, detail):
        
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.flow_boundary_point):

        name = detail.name + "_" + point.name
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
#        obs_list.append(obspoint(x, y, name, crs=overall.crs))
        overall.add_observation_point(x, y, name)
    
def nest1_sfincs_in_sfincs(overall, detail):
    
#    from sfincs import ObservationPoint as obspoint

    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    # Loop over all points in the gdf of the boundary conditions
    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = detail.name + "_" + str(ind + 1).zfill(4)
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points.add_point(x, y, name)
    
#     for ind, point in enumerate(detail.flow_boundary_point):

#         name = detail.name + "_" + str(ind + 1).zfill(4)
#         x, y = transformer.transform(point.geometry.x,
#                                      point.geometry.y)
# #        obs_list.append(obspoint(x, y, name, crs=overall.crs))
#         overall.add_observation_point(x, y, name)

def nest1_xbeach_in_sfincs(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.flow_boundary_point):

        name = detail.name + "_" + str(ind + 1).zfill(4)
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
#        obs_list.append(obspoint(x, y, name, crs=overall.crs))
        overall.add_observation_point(x, y, name)

def nest1_beware_in_sfincs(overall, detail):
        
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.flow_boundary_point):

        name = detail.name + "_" + point.name
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
#        obs_list.append(obspoint(x, y, name, crs=overall.crs))
        overall.add_observation_point(x, y, name)

def nest1_hurrywave_in_hurrywave(overall, detail):

    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)

    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = detail.name + "_" + str(ind + 1).zfill(4)
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points_sp2.add_point(x, y, name)    

def nest1_xbeach_in_hurrywave(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.wave_boundary_point):

        name = detail.name + "_" + str(ind + 1).zfill(4)
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)

        overall.observation_points_sp2.add_point(x, y, name)
        overall.observation_points_regular.add_point(x, y, name)

def nest1_sfincs_in_hurrywave(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)

    for ind, point in enumerate(detail.wave_boundary_point):

        name = detail.name + "_" + str(ind + 1).zfill(4)
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
        overall.observation_points_regular.add_point(x, y, name)

def nest1_beware_in_hurrywave(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)

    for ind, point in enumerate(detail.wave_boundary_point):

        
        name = detail.name + "_" + point.name
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
        overall.observation_points_regular.add_point(x, y, name)
