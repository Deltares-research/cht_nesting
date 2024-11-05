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
from cht_utils.deltares_ini import IniStruct
from cht_tide.tide_predict import predict
from .sfincs_in_sfincs import sfincs_in_sfincs

def nest2(overall,
          detail,
          obs_point_prefix=None,
          output_path=None,
          output_file=None,
          bc_path=None,
          bc_file=None,
          overall_crs=None,
          detail_crs=None,
	      option=None,
          boundary_water_level_correction=0.0,
          filter_incoming=False,
          bctype="waterlevel",          
          return_maximum=False):

    if type(overall) == str:
        # Overall is a string, so we need to instantiate the class
        if overall == "sfincs":
            from cht_sfincs.sfincs import SFINCS
            overall = SFINCS()
        elif overall == "hurrywave":
            from cht_hurrywave.hurrywave import HurryWave
            overall = HurryWave()
        elif overall == "xbeach":
            from cht_xbeach.xbeach import XBeach
            overall = XBeach()
        elif overall == "beware":
            from cht_beware.beware import BEWARE
            overall = BEWARE()
        elif overall == "delft3dfm":
            from cht_delft3dfm.delft3dfm import Delft3DFM
            overall = Delft3DFM()           

    # Get the types of overall and detail classes
    overall_type = overall.__class__.__name__.lower()
    detail_type = detail.__class__.__name__.lower()

    # Check if detail model has attribute name
    if obs_point_prefix is None:
        if hasattr(detail, "name"):
            obs_point_prefix = detail.name
        else:
            obs_point_prefix = "nest"

    if overall_crs is not None:
        overall.crs = CRS(overall_crs)
    if detail_crs is not None:
        detail.crs = CRS(detail_crs)

    if overall_type == "delft3dfm":

        if detail_type == "delft3dfm":
            nest2_delft3dfm_in_delft3dfm(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            boundary_water_level_correction=boundary_water_level_correction)

        elif detail_type == "sfincs":
            nest2_sfincs_in_delft3dfm(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            boundary_water_level_correction=boundary_water_level_correction)

        elif detail_type == "beware":
            nest2_beware_in_delft3dfm(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            option=option,
                                            boundary_water_level_correction=boundary_water_level_correction)
            
    elif overall_type == "sfincs":

        if detail_type == "sfincs":
            zs = sfincs_in_sfincs(overall,    
                                  detail,
                                  obs_point_prefix=obs_point_prefix,
                                  output_path=output_path,
                                  output_file=output_file,
                                  bc_path=bc_path,
                                  boundary_water_level_correction=boundary_water_level_correction,
                                  return_maximum=return_maximum,
                                  filter_incoming=filter_incoming,
                                  bctype=bctype          
                                 ) 
            return zs
        elif detail_type == "xbeach":
            bc = nest2_xbeach_in_sfincs(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            boundary_water_level_correction=boundary_water_level_correction,
                                            return_maximum=return_maximum)
            return bc
        elif detail_type == "beware":
            nest2_beware_in_sfincs(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            option=option,
                                            boundary_water_level_correction=boundary_water_level_correction)

    elif overall_type == "hurrywave":

        if detail_type == "hurrywave":
            nest2_hurrywave_in_hurrywave(overall,
                                            detail,
                                            obs_point_prefix=obs_point_prefix,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path)
            
        elif detail_type == "xbeach":
            bc = nest2_xbeach_in_hurrywave(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            option=option,
                                            return_maximum=return_maximum)
            return bc
        elif detail_type == "sfincs":
            nest2_sfincs_in_hurrywave(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path)

        elif detail_type == "beware":
            nest2_beware_in_hurrywave(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path)

    elif overall_type == "beware":

        if detail_type == "sfincs":
            nest2_sfincs_in_beware(overall,
                                            detail,
                                            output_path=output_path,
                                            output_file=output_file,
                                            bc_path=bc_path,
                                            option=option,
                                            boundary_water_level_correction=boundary_water_level_correction)

    return True

def nest2_delft3dfm_in_delft3dfm(overall,
                                 detail,
                                 output_path=None,
                                 output_file= "flow_his.nc",
                                 boundary_water_level_correction=0,
                                 bc_path= None):

    for ind, bnd in enumerate(detail.boundary):        
        point_names = []
        for ip, point in enumerate(bnd.point):
            point_names.append(detail.name + "_" + point.name)

        # Return DataFrame bzs
        bzs = overall.read_timeseries_output(name_list=point_names,
                                             path=output_path,
                                             file_name=output_file)
        
        ts  = bzs.index
        for ip, point in enumerate(bnd.point):
            point.data = pd.Series(bzs.iloc[:,ip].values, index=ts) + boundary_water_level_correction    
    
    if bc_path is not None:
        detail.write_flow_boundary_conditions(path=bc_path) # check if this works

def nest2_sfincs_in_delft3dfm(overall,
                              detail,
                              output_path=None,
                              output_file=None,
                              boundary_water_level_correction=0,
                              bc_path=None):
    
    if not output_file:
        output_file = "flow_his.nc"

    point_names = []
    for point in detail.flow_boundary_point:
        point_names.append(detail.name + "_" + point.name)                    
    output_file = os.path.join(output_path, output_file)

    # Return DataFrame bzs
    bzs = overall.read_timeseries_output(name_list=point_names,
                                         path=output_path,
                                         file_name=output_file)

    ts  = bzs.index
    for icol, point in enumerate(detail.flow_boundary_point):
        point.data = pd.Series(bzs.iloc[:,icol].values, index=ts) + boundary_water_level_correction  
    
    # Write bzs file
    if bc_path is not None:
        detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))

def nest2_beware_in_delft3dfm(overall,
                              detail,
                              output_path=None,
                              output_file=None,
                              boundary_water_level_correction = 0,
                              option=None,
                              bc_path=None):

    if option == 'flow':

        if not output_file:
            # Path of the overall output time series
            output_file = "flow_his.nc"
        
        point_names = []
        for point in detail.flow_boundary_point:
            point_names.append(detail.name + "_" + point.name) #reopen                   
            # point_names.append(point.name)                    
        # output_file = os.path.join(output_path, output_file)

        # Return DataFrame bzs
        bzs = overall.read_timeseries_output(name_list=point_names,
                                             path=output_path,
                                             file_name=output_file)
        ts  = bzs.index
        for icol, point in enumerate(detail.flow_boundary_point):
            point.data = pd.Series(bzs.iloc[:,icol].values, index=ts) + boundary_water_level_correction  
        
        if bc_path is not None:
            detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))
            
    if option == 'wave':

        if not output_path:
            # Path of the overall output time series
            output_path = overall.path
            
        if not output_file:
            file_name = os.path.join(output_path, "wavh-wave-nest.nc")
        else:
            file_name = os.path.join(output_path, output_file)
        
        # Open netcdf file
        file_name_dfm = os.path.join(output_path, "flow_his.nc")
        ddd = xr.open_dataset(file_name_dfm)
        stations=ddd.station_name.values
        all_stations = []
        for ist, st in enumerate(stations):
            st=str(st.strip())[2:-1] #reopen
            # st=str(st.strip())[2:6]
            all_stations.append(st)

        point_names = []
        for point in detail.wave_boundary_point:
            point_names.append(detail.name + "_" + point.name) #reopen                   
            # point_names.append(point.name)  
        # ireq = []    
        # for ip, point in enumerate(point_names):
        #     for ist,st in enumerate(all_stations):
        #         if point.lower() == st.lower():
        #             ireq.append(ist)            
        #             break

        ddd = xr.open_dataset(file_name)
        times   = ddd.Hsig.time.values

        for ip, point in enumerate(detail.wave_boundary_point):
            for ist,st in enumerate(all_stations):
                ireq = -1
                if point_names[ip] == st.lower():
                    ireq=ist
                    break

            if ireq>-1:
                hm0     = ddd.Hsig.values[:,ireq]
                tp      = ddd.RTpeak.values[:,ireq]
                wavdir  = ddd.Dir.values[:,ireq]
                dirspr  = ddd.Dspr.values[:,ireq]
    
                df = pd.DataFrame(index=times)
                df.insert(0,"hm0",hm0)
                df.insert(1,"tp",tp)
                df.insert(2,"wavdir",wavdir)
                df.insert(3,"dirspr",dirspr)
    
                point.data = df

        if bc_path is not None:
            detail.write_wave_boundary_conditions(path=bc_path)




def nest2_xbeach_in_sfincs(overall,
                           detail,
                           output_path=None,
                           output_file=None,
                           boundary_water_level_correction=0,
                           return_maximum=False,
                           bc_path=None):

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
    
    point_names = []
    for point in detail.flow_boundary_point:
        point_names.append(detail.name + "_" + point.name)                    
    zstfile = os.path.join(output_path, output_file)

    # Return DataFrame bzs
    bzs = overall.output.read_his_file(station=point_names,
                                       parameter="point_zs",
                                       file_name=zstfile)

    # Interpolate on desired format for XBeach forcing
    bzs_resampled = bzs.resample('10min').mean()
    bzs_interpolated = bzs_resampled.interpolate(method='linear')
    bzs_filtered = bzs_interpolated[detail.tref:detail.tstop]
    
    ts  = bzs_filtered.index
    for icol, point in enumerate(detail.flow_boundary_point):
        point.data = pd.Series(bzs_filtered.iloc[:,icol].values, index=ts) + boundary_water_level_correction

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
                if zx>zmax:
                    zs = point.data
                    zmax = zx
                    
        elif len(detail.flow_boundary_point) == 4:
            for icol, point in enumerate(detail.flow_boundary_point):
                if icol == 2:
                    break
                zx = point.data.max()
                if zx>zmax:
                    zs = point.data
                    zmax = zx

        return zs                          
    else:    
        return detail.flow_boundary_point
    
    # Write boundary conditions
    if bc_path is not None:
        detail.write_flow_boundary_conditions() # check Roel

def nest2_beware_in_sfincs(overall,
                           detail,
                           output_path=None,
                           output_file=None,
                           boundary_water_level_correction=0,
                           option=None,
                           bc_path=None):

    if option == 'flow':
        if not output_file:
            # Path of the overall output time series
            output_file = "sfincs_his.nc"
        
        point_names = []
        for point in detail.flow_boundary_point:
            point_names.append(detail.name + "_" + point.name) #reopen                   
            # point_names.append(point.name)                    
        # output_file = os.path.join(output_path, output_file)

        # Return DataFrame bzs
        bzs = overall.output.read_his_file(station=point_names,
                                        parameter="point_zs",
                                        file_name=os.path.join(output_path, output_file))

        # Replace -999.0 with zeros. This should not happen, but easy fix for now.
        bzs = bzs.replace(-999.0,0.0)
        for icol, point in enumerate(detail.flow_boundary_point):
            point.data = pd.Series(bzs.iloc[:,icol].values, index=bzs.index) + boundary_water_level_correction    
            
        if bc_path is not None:
            detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))

def nest2_hurrywave_in_hurrywave(overall,
                                 detail,
                                 obs_point_prefix=None,
                                 output_path=None,
                                 output_file=None,
                                 bc_path=None):
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path
    if not output_file:
        output_file = "hurrywave_sp2.nc"

    file_name = os.path.join(output_path, output_file)
    
    detail.boundary_conditions.forcing = "spectra"

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    stations=ddd.station_name.values
    all_stations = []
    for ist, st in enumerate(stations):
        st=str(st.strip())[2:-1]
        all_stations.append(st)

    if len(detail.boundary_conditions.gdf)==0:
        detail.input.variables.bndfile = "hurrywave.bnd"
        detail.boundary_conditions.read_boundary_points()

    point_names = []
    if len(detail.boundary_conditions.gdf)>0:
        for ind, row in detail.boundary_conditions.gdf.iterrows():
        # Find required boundary points        
            point_names.append(obs_point_prefix + "_" + row["name"])                    
        
    else:
        point_names = all_stations.copy()
        
    times   = ddd.point_spectrum2d.coords["time"].values
    sigma   = ddd.point_spectrum2d.coords["sigma"].values
    theta   = ddd.point_spectrum2d.coords["theta"].values

    ireq = []    
    for ip, point in enumerate(point_names):
        for ist,st in enumerate(all_stations):
            if point.lower() == st.lower():
                ireq.append(ist)            
                break

    for ind, row in detail.boundary_conditions.gdf.iterrows():

        sp2 = ddd.point_spectrum2d.values[:,ireq[ind],:,:]

        ds = xr.Dataset(
                data_vars = dict(point_spectrum2d=(["time", "theta", "sigma"], sp2)),
                coords    = dict(time=times,
                                 theta=theta,
                                 sigma=sigma)
                )
        detail.boundary_conditions.gdf.loc[ind, "spectra"] = ds.to_array()

    if bc_path is not None:
        detail.boundary_conditions.write_boundary_conditions_spectra(file_name=os.path.join(bc_path, detail.input.variables.bspfile))

     
def nest2_xbeach_in_hurrywave(overall,
                              detail,
                              output_path=None,
                              output_file=None,
                              option=None,
                              return_maximum=False,
                              bc_path=None):
    
    from cht_physics.deshoal import deshoal
    from scipy import interpolate
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path

    if option == "sp2":    
        if not output_file:
            output_file = "hurrywave_sp2.nc"
    
        file_name = os.path.join(output_path, output_file)


    
        # Open netcdf file
        ddd = xr.open_dataset(file_name)
        stations=ddd.station_name.values
        all_stations = []
        for ist, st in enumerate(stations):
            st=str(st.strip())[2:-1]
            all_stations.append(st)
    
        point_names = []    
        if detail.wave_boundary_point:
            # Find required boundary points        
            for point in detail.wave_boundary_point:
                point_names.append(detail.name + "_" + point.name)                    
            
        else:
            point_names = all_stations.copy()
            
        times   = ddd.point_spectrum2d.coords["time"].values
        sigma   = ddd.point_spectrum2d.coords["sigma"].values
        theta   = ddd.point_spectrum2d.coords["theta"].values
    
        ireq = []    
        for ip, point in enumerate(point_names):
            for ist,st in enumerate(all_stations):
                if point.lower() == st.lower():
                    ireq.append(ist)            
                    break
    
        for ip, point in enumerate(detail.wave_boundary_point):
    
            sp2 = ddd.point_spectrum2d.values[:,ireq[ip],:,:]
    
            ds = xr.Dataset(
                    data_vars = dict(point_spectrum2d=(["time", "theta", "sigma"], sp2)),
                    coords    = dict(time=times,
                                     theta=theta,
                                     sigma=sigma)
                    )
            
            point.data = ds

    elif option == "timeseries":

        if not output_file:
            output_file = "hurrywave_his.nc"

        file_name = os.path.join(output_path, output_file)
    
        # Open netcdf file
        ddd = xr.open_dataset(file_name)
        stations=ddd.station_name.values
        all_stations = []
        for ist, st in enumerate(stations):
            st=str(st.strip())[2:-1]
            all_stations.append(st)
    
        point_names = []    
        if detail.wave_boundary_point:
            # Find required boundary points        
            for point in detail.wave_boundary_point:
                point_names.append(detail.name + "_" + point.name)                    
            
        else:
            point_names = all_stations.copy()
            
        times   = ddd.point_hm0.coords["time"].values

        ireq = []    
        for ip, point in enumerate(point_names):
            for ist,st in enumerate(all_stations):
                if point.lower() == st.lower():
                    ireq.append(ist)            
                    break

        for ip, point in enumerate(detail.wave_boundary_point):
            hm0     = ddd.point_hm0.values[:,ireq[ip]]
            tp      = ddd.point_tp.values[:,ireq[ip]]
            zb_point = ddd.station_z.values[ip]

            #deshoal waveheights to offshore boundary
            hm0_deshoal = []
            if detail.zb_deshoal:
                try:
                    zs = detail.flow_boundary_point[0].data
                    # Interpolate to wave timeseries
                    wave_secs = times.astype(float)
                    flow_secs = zs.index.values.astype(float)
                    f = interpolate.interp1d(flow_secs, zs,fill_value=0,bounds_error=False)
                    zs = f(wave_secs)
                except:
                    zs = 0*hm0
                    
                for ih, h in enumerate(hm0):
                    hm0_deshoal.append(deshoal(h, tp[ih], abs(zb_point)+zs[ih],abs(detail.zb_deshoal))[0]+zs[ih])
                    
                hm0 = hm0_deshoal
            
            #set wavedir such that waves are forced perpendicular to coast instead of real direction
            wavdir =  np.mean([detail.params["thetamin"], detail.params["thetamax"]])
            #wavdir  = ddd.point_wavdir.values[:,ireq[ip]]
            
            #convert directional spread in degrees to xbeach spreading parameter
            dirspr  = ddd.point_dirspr.values[:,ireq[ip]]
            s = 2/(dirspr*np.pi/180)**2 - 1
            
    
            df = pd.DataFrame(index=times)
            df.insert(0,"hm0",hm0)
            df.insert(1,"tp",tp)
            df.insert(2,"wavdir",wavdir)
            df.insert(3,'gammajsp', 3.3)
            df.insert(4,"s",s)
    
            #resample to half-hourly data
            df_resampled = df.resample('30min').max()
            df_interpolated = df_resampled.interpolate(method='linear')
            mask = (df_interpolated.index >= detail.tref) & (df_interpolated.index <= detail.tstop)
            df_filtered = df_interpolated[mask]
    
            df_filtered.insert(5,'duration', 1800)
            df_filtered.insert(6,"dtbc",1)
            
            point.data = df_filtered

        if bc_path is not None:
            detail.write_wave_boundary_conditions(option=option) # check roel

    if return_maximum:
        hmax = -999.0
        for icol, point in enumerate(detail.wave_boundary_point):
            hx = point.data["hm0"].max()
            if hx>hmax:
                hmx = point.data["hm0"]
                hmax = hx
        return hmx        
                    
    else:    
        return detail.wave_boundary_point


def nest2_sfincs_in_hurrywave(overall,
                              detail,
                              output_path=None,
                              output_file=None,
                              bc_path=None):
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path
    if not output_file:
        output_file = "hurrywave_his.nc"
        
    file_name = os.path.join(output_path, output_file)
    print("Nesting in " + file_name)

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    stations=ddd.station_name.values
    all_stations = []
    for ist, st in enumerate(stations):
        st=str(st.strip())[2:-1]
        all_stations.append(st)

    point_names = []    
    if detail.wave_boundary_point:
        # Find required boundary points        
        for point in detail.wave_boundary_point:
            point_names.append(detail.name + "_" + point.name)                    
        
    else:
        point_names = all_stations.copy()
        
    times   = ddd.point_hm0.coords["time"].values

    ireq = []    
    for ip, point in enumerate(point_names):
        for ist,st in enumerate(all_stations):
            if point.lower() == st.lower():
                ireq.append(ist)            
                break

    for ip, point in enumerate(detail.wave_boundary_point):

        hm0     = ddd.point_hm0.values[:,ireq[ip]]
        tp      = ddd.point_tp.values[:,ireq[ip]]
        wavdir  = ddd.point_wavdir.values[:,ireq[ip]]
        dirspr  = ddd.point_dirspr.values[:,ireq[ip]]

        df = pd.DataFrame(index=times)
        df.insert(0,"hm0",hm0)
        df.insert(1,"tp",tp)
        df.insert(2,"wavdir",wavdir)
        df.insert(3,"dirspr",dirspr)

        point.data = df

    if bc_path is not None:
        detail.write_wave_boundary_conditions(path=bc_path)

def nest2_beware_in_hurrywave(overall,
                              detail,
                              output_path=None,
                              output_file=None,
                              bc_path=None):

    if not output_path:
        # Path of the overall output time series
        output_path = overall.path
    if not output_file:
        output_file = "hurrywave_his.nc"
        
    file_name = os.path.join(output_path, output_file)

    # Open netcdf file
    ddd = xr.open_dataset(file_name)
    stations=ddd.station_name.values
    all_stations = []
    for ist, st in enumerate(stations):
        st=str(st.strip())[2:-1]
        all_stations.append(st)

    point_names = []    
    if detail.wave_boundary_point:
        # Find required boundary points        
        for point in detail.wave_boundary_point:
            point_names.append(detail.name + "_" + point.name)                    
        
    else:
        point_names = all_stations.copy()
        
    times   = ddd.point_hm0.coords["time"].values

    ireq = []    
    for ip, point in enumerate(point_names):
        for ist,st in enumerate(all_stations):
            if point.lower() == st.lower():
                ireq.append(ist)            
                break

    for ip, point in enumerate(detail.wave_boundary_point):

        hm0     = ddd.point_hm0.values[:,ireq[ip]]
        tp      = ddd.point_tp.values[:,ireq[ip]]
        wavdir  = ddd.point_wavdir.values[:,ireq[ip]]
        dirspr  = ddd.point_dirspr.values[:,ireq[ip]]

        df = pd.DataFrame(index=times)
        df.insert(0,"hm0",hm0)
        df.insert(1,"tp",tp)
        df.insert(2,"wavdir",wavdir)
        df.insert(3,"dirspr",dirspr)

        point.data = df

    if bc_path is not None:
        detail.write_wave_boundary_conditions(path=bc_path)


def nest2_sfincs_in_beware(overall,
                           detail,
                           output_path=None,
                           output_file=None,
                           boundary_water_level_correction=0,
                           option=None,
                           bc_path=None):

    from cht_sfincs.sfincs import FlowBoundaryPoint
    from cht_sfincs.sfincs import WaveMakerForcingPoint

    if not output_file:
        output_file = "beware_his.nc"
    # Get bounding box for sfincs model
    # Convert bbox to beware crs

    x_range, y_range = detail.bounding_box(crs=overall.crs)
    dx = (x_range[1] - x_range[0])/10
    dy = (y_range[1] - y_range[0])/10
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
        inear = np.where((xb>x_range[0]) & (xb<x_range[1]) & (yb>y_range[0]) & (yb<y_range[1]))    
        xb=xb[inear]
        yb=yb[inear]
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
        tref = datetime.datetime(1970,1,1)
        tsec = ddd.time.values # array of int64
        times = tref + tsec*datetime.timedelta(seconds=1)

        for ip, point in enumerate(detail.flow_boundary_point):
            point.data = pd.Series(ddd.WL.values[inear[0][ip],:], index=times) + pd.Series(ddd.R2_setup.values[inear[0][ip],:], index=times) + boundary_water_level_correction

        if bc_path is not None:    
            detail.write_flow_boundary_conditions(file_name=os.path.join(bc_path, detail.input.bzsfile))

    elif option == "wave":

        xb = ddd.x_coast.values
        yb = ddd.y_coast.values

        # Find beware locations in bounding box
        inear = np.where((xb>x_range[0]) & (xb<x_range[1]) & (yb>y_range[0]) & (yb<y_range[1]))    
        xb=xb[inear]
        yb=yb[inear]
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
        tref = datetime.datetime(1970,1,1)
        tsec = ddd.time.values # array of int64
        times = tref + tsec*datetime.timedelta(seconds=1)
        
        for ip, point in enumerate(detail.wavemaker_forcing_point):

            df = pd.DataFrame()
            df["hm0_ig"] = ddd.hm0_ig.values[inear[0][ip],:]
            df["tp_ig"]  = ddd.tp_ig.values[inear[0][ip],:]
            df["setup"]  = ddd.setup.values[inear[0][ip],:]
            df["time"]   = times
            df = df.set_index("time")
            
            df["hm0_ig"]=df["hm0_ig"].replace(np.nan, 0.1)
            df["tp_ig"]=df["tp_ig"].replace(np.nan, 60.0)
            df["setup"]=df["setup"].replace(np.nan, 0.0)
            
            point.data = df

    ddd.close()

    if bc_path is not None:
        detail.write_whi_file(file_name=os.path.join(bc_path, detail.input.whifile))
        detail.write_wti_file(file_name=os.path.join(bc_path, detail.input.wtifile))


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
    