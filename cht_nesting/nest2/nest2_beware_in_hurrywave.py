"""Nest 2 script for nesting BEWARE within hydromt_hurrywave HurrywaveModel.

Reads wave parameter time series from the HurrywaveModel's hurrywave_his.nc output
and sets them as wave boundary conditions on the detail BEWARE model.
"""

import os
from typing import Any, Optional

import pandas as pd
import xarray as xr


def nest2_beware_in_hurrywave(
    overall: Any,
    detail: Any,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a BEWARE model within a HurrywaveModel.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel whose his output is read.
    detail : BEWARE
        The fine BEWARE model that receives wave boundary conditions.
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

    ddd = xr.open_dataset(file_name)
    stations = ddd.station_name.values
    all_stations = [str(st.strip())[2:-1] for st in stations]

    point_names = []
    if detail.wave_boundary_point:
        for point in detail.wave_boundary_point:
            point_names.append(f"{detail.name}_{point.name}")
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
