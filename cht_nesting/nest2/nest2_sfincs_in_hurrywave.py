"""Nest 2 script for nesting hydromt_sfincs SfincsModel within hydromt_hurrywave HurrywaveModel.

Reads wave parameter time series from the HurrywaveModel's hurrywave_his.nc output
and sets them as SnapWave boundary conditions on the detail SfincsModel.
"""

import os
from typing import Any, List, Optional

import numpy as np
import pandas as pd
import xarray as xr


def nest2_sfincs_in_hurrywave(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a SfincsModel within a HurrywaveModel (SnapWave BCs).

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel whose his output is read.
    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel that receives the SnapWave boundary conditions.
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

    fn_his: str = os.path.join(output_path, output_file)
    ddd = xr.open_dataset(fn_his)

    # Decode station names
    all_stations: List[str] = [
        st.decode("utf-8").strip()
        if isinstance(st, (bytes, np.bytes_))
        else str(st).strip()
        for st in ddd.station_name.values
    ]

    point_names: List[str] = [
        f"{obs_point_prefix}_{str(ind + 1).zfill(4)}"
        for ind, _row in detail.snapwave_boundary_conditions.gdf.iterrows()
    ]

    times = pd.DatetimeIndex(ddd.point_hm0.coords["time"].values)

    ireq: List[int] = []
    for pname in point_names:
        matched = False
        for ist, st in enumerate(all_stations):
            if pname.lower() == st.lower():
                ireq.append(ist)
                matched = True
                break
        if not matched:
            raise ValueError(
                f"Observation point '{pname}' not found in '{fn_his}'. "
                f"Available stations: {all_stations}"
            )

    gdf_indices: List[Any] = [
        ind for ind, _row in detail.snapwave_boundary_conditions.gdf.iterrows()
    ]

    # Extract and clip wave parameters to valid physical ranges.
    hm0_all: np.ndarray = ddd.point_hm0.values[:, ireq]
    tp_all: np.ndarray = np.clip(ddd.point_tp.values[:, ireq], 1.0, 25.0)
    wavdir_all: np.ndarray = ddd.point_wavdir.values[:, ireq]
    dirspr_all: np.ndarray = np.clip(ddd.point_dirspr.values[:, ireq], 4.0, 55.0)

    ddd.close()

    ds_bcs = xr.Dataset(
        data_vars={
            "hs": (["time", "index"], hm0_all),
            "tp": (["time", "index"], tp_all),
            "wd": (["time", "index"], wavdir_all),
            "ds": (["time", "index"], dirspr_all),
        },
        coords={
            "time": times,
            "index": gdf_indices,
        },
    )

    detail.snapwave_boundary_conditions.set_timeseries(ds_bcs)

    if bc_path is not None:
        detail.snapwave_boundary_conditions.write()
