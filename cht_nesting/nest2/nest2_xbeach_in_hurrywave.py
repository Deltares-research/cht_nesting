# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2_xbeach_in_hurrywave function for handling nesting of XBeach models within HurryWave models.
"""

import os
import xarray as xr
import pandas as pd
import numpy as np
from scipy import interpolate
from typing import Optional, Any

def nest2_xbeach_in_hurrywave(overall: Any,
                              detail: Any,
                              obs_point_prefix: Optional[str] = None,
                              output_path: Optional[str] = None,
                              output_file: Optional[str] = None,
                              option: Optional[str] = None,
                              return_maximum: bool = False,
                              bc_path: Optional[str] = None,
                              **kwargs) -> Any:
    """
    Nest an XBeach model within a HurryWave model.

    Parameters:
    overall (Any): The overall HurryWave model.
    detail (Any): The detailed XBeach model.
    obs_point_prefix (Optional[str]): The prefix for observation points. Default is None.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    option (Optional[str]): The option for nesting ("sp2" or "timeseries"). Default is None.
    return_maximum (bool): Whether to return the maximum values. Default is False.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    **kwargs: Additional keyword arguments.

    Returns:
    Any: The boundary points or maximum values.
    """
    if not output_path:
        # Path of the overall output time series
        output_path = overall.path

    if option == "sp2":
        if not output_file:
            output_file = "hurrywave_sp2.nc"

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

        times = ddd.point_spectrum2d.coords["time"].values
        sigma = ddd.point_spectrum2d.coords["sigma"].values
        theta = ddd.point_spectrum2d.coords["theta"].values

        ireq = []
        for ip, point in enumerate(point_names):
            for ist, st in enumerate(all_stations):
                if point.lower() == st.lower():
                    ireq.append(ist)
                    break

        for ip, point in enumerate(detail.wave_boundary_point):
            sp2 = ddd.point_spectrum2d.values[:, ireq[ip], :, :]

            ds = xr.Dataset(
                data_vars=dict(point_spectrum2d=(["time", "theta", "sigma"], sp2)),
                coords=dict(time=times, theta=theta, sigma=sigma)
            )

            point.data = ds

    elif option == "timeseries":
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
            zb_point = ddd.station_z.values[ip]

            # Deshoal wave heights to offshore boundary
            hm0_deshoal = []
            if detail.zb_deshoal:
                try:
                    zs = detail.flow_boundary_point[0].data
                    # Interpolate to wave timeseries
                    wave_secs = times.astype(float)
                    flow_secs = zs.index.values.astype(float)
                    f = interpolate.interp1d(flow_secs, zs, fill_value=0, bounds_error=False)
                    zs = f(wave_secs)
                except:
                    zs = 0 * hm0

                for ih, h in enumerate(hm0):
                    hm0_deshoal.append(deshoal(h, tp[ih], abs(zb_point) + zs[ih], abs(detail.zb_deshoal))[0] + zs[ih])

                hm0 = hm0_deshoal

            # Set wave direction such that waves are forced perpendicular to coast instead of real direction
            wavdir = np.mean([detail.params["thetamin"], detail.params["thetamax"]])
            # wavdir = ddd.point_wavdir.values[:, ireq[ip]]

            # Convert directional spread in degrees to XBeach spreading parameter
            dirspr = ddd.point_dirspr.values[:, ireq[ip]]
            s = 2 / (dirspr * np.pi / 180) ** 2 - 1

            df = pd.DataFrame(index=times)
            df.insert(0, "hm0", hm0)
            df.insert(1, "tp", tp)
            df.insert(2, "wavdir", wavdir)
            df.insert(3, 'gammajsp', 3.3)
            df.insert(4, "s", s)

            # Resample to half-hourly data
            df_resampled = df.resample('30min').max()
            df_interpolated = df_resampled.interpolate(method='linear')
            mask = (df_interpolated.index >= detail.tref) & (df_interpolated.index <= detail.tstop)
            df_filtered = df_interpolated[mask]

            df_filtered.insert(5, 'duration', 1800)
            df_filtered.insert(6, "dtbc", 1)

            point.data = df_filtered

        if bc_path is not None:
            detail.write_wave_boundary_conditions(option=option)

    if return_maximum:
        hmax = -999.0
        for icol, point in enumerate(detail.wave_boundary_point):
            hx = point.data["hm0"].max()
            if hx > hmax:
                hmx = point.data["hm0"]
                hmax = hx
        return hmx
    else:
        return detail.wave_boundary_point