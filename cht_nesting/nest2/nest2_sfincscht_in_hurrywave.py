"""Nest 2 script for nesting cht_sfincs SFINCS within hydromt_hurrywave HurrywaveModel.

Reads wave parameter time series from the HurrywaveModel's hurrywave_his.nc output
and sets them as SnapWave boundary conditions on the detail SFINCS model (cht_sfincs API).
"""

import os
from typing import Any, Optional

import numpy as np
import pandas as pd
import xarray as xr


def nest2_sfincscht_in_hurrywave(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs,
) -> None:
    """Nest a cht_sfincs SFINCS model within a HurrywaveModel.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel whose his output is read.
    detail : cht_sfincs.SFINCS
        The fine SFINCS model that receives the SnapWave boundary conditions.
    obs_point_prefix : str, optional
        Prefix for observation point names.
    output_path : str, optional
        Directory containing the HurrywaveModel output files.
    output_file : str, optional
        Name of the his output file. Defaults to "hurrywave_his.nc".
    bc_path : str, optional
        If provided, write boundary conditions to disk.
    """
    if not output_path:
        output_path = str(overall.root.path)
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
        dirspr = np.clip(dirspr, 4.0, 55.0)
        tp = np.clip(tp, 1.0, 25.0)

        df = pd.DataFrame(index=times)
        df.insert(0, "hs", hm0)
        df.insert(1, "tp", tp)
        df.insert(2, "wd", wavdir)
        df.insert(3, "ds", dirspr)

        detail.snapwave.boundary_conditions.gdf.loc[ip, "timeseries"] = df

    if bc_path is not None:
        detail.snapwave.boundary_conditions.write_boundary_conditions_timeseries()
