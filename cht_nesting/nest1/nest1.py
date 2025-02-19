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

def nest1(overall, detail, option=None, obs_point_prefix=None):
    
    # Returns a list with observation point objects

    # Check if detail model has attribute name
    if obs_point_prefix is None:
        if not hasattr(detail, "name"):
            detail.name = "nest"
    else:
        detail.name = obs_point_prefix       

    # Get the types of overall and detail classes
    overall_type = overall.__class__.__name__.lower()
    detail_type = detail.__class__.__name__.lower()

    if overall_type == "delft3dfm":
        if detail_type == "delft3dfm":
            nest1_delft3dfm_in_delft3dfm(overall, detail)
        elif detail_type == "sfincs":
            nest1_sfincs_in_delft3dfm(overall, detail)
        elif detail_type == "beware":
            nest1_beware_in_delft3dfm(overall, detail)
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False
            
    elif overall_type == "sfincs":
        if detail_type == "sfincs":
            nest1_sfincs_in_sfincs(overall, detail)
        elif detail_type == "xbeach":
            nest1_xbeach_in_sfincs(overall, detail)
        elif detail_type == "beware":
            nest1_beware_in_sfincs(overall, detail)
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    elif overall_type == "hurrywave":
        if detail_type == "hurrywave":
            nest1_hurrywave_in_hurrywave(overall, detail)
        elif detail_type == "xbeach":    
            nest1_xbeach_in_hurrywave(overall, detail)
        elif detail_type == "sfincs":    
            nest1_sfincs_in_hurrywave(overall, detail)
        elif detail_type == "beware":
            nest1_beware_in_hurrywave(overall, detail)
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    elif overall_type == "beware":
        if detail_type == "sfincs":
            # No need to do anything here. BEWARE output points are fixed
            pass
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    return True    

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
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)

    # Check if overall model is in geographic coordinates. If so, make sure that longitude is correct. Either 0-360 or -180-180.
    overall_degrees_west = False
    if overall.crs.is_geographic:
        # Get max longitude of the overall model
        if overall.grid.data is None:
            overall.grid.read()
        x_max = overall.grid.data.grid.face_coordinates[:, 0].max()
        if x_max > 180.0:
            overall_degrees_west = True

    # Get list of names of the observation points
    overall_names = overall.observation_points.list_names()
    # Loop over all points in the gdf of the boundary conditions
    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = detail.name + "_" + str(ind + 1).zfill(4)
        if name in overall_names:
            print(f"Observation point {name} already exists in the model")
            continue
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        if x < 0 and overall_degrees_west:
            x += 360
        overall.observation_points.add_point(x, y, name)

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
        # overall.add_observation_point(x, y, name)
        overall.observation_points.add_point(x, y, name)

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

    for ind, point in detail.snapwave.boundary_conditions.gdf.iterrows():

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
