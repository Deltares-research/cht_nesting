# -*- coding: utf-8 -*-
"""
Nest 2 script for nesting a hydromt_sfincs SfincsModel within a HurryWave model.

Reads wave parameter time series from the HurryWave ``hurrywave_his.nc`` output
and sets them as SnapWave boundary conditions on the detail SfincsModel.
Station names in the his-file are stored as fixed-length bytes by the HurryWave
Fortran kernel and are decoded automatically.
"""

import os
from typing import Any, List, Optional

import numpy as np
import pandas as pd
import xarray as xr


def nest2_hurrywave_in_sfincsmodel(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a hydromt_sfincs SfincsModel within a HurryWave model (SnapWave BCs).

    Reads significant wave height (Hm0), peak period (Tp), wave direction and
    directional spread from the HurryWave his output and applies them as SnapWave
    boundary conditions on the detail SfincsModel via
    ``detail.snapwave_boundary_conditions.set_timeseries()``.

    Wave parameters are clipped to physically realistic ranges before being set:

    * Tp is clipped to [1, 25] seconds.
    * Directional spread is clipped to [4, 55] degrees.

    Parameters
    ----------
    overall : HurryWave
        The coarse HurryWave model.  Must expose a ``path`` attribute pointing
        to the directory that contains the his output file.
    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel that receives the SnapWave boundary conditions.
    obs_point_prefix : str, optional
        Prefix prepended to each detail SnapWave boundary point name when
        looking up stations in the overall his-file.
    output_path : str, optional
        Directory containing the HurryWave output files.  Defaults to
        ``overall.path``.
    output_file : str, optional
        Name of the HurryWave his output file.  Defaults to
        ``"hurrywave_his.nc"``.
    bc_path : str, optional
        If provided, write the boundary conditions to disk by calling
        ``detail.snapwave_boundary_conditions.write()`` after setting.  The
        value of *bc_path* is not used directly — writing goes to the model's
        configured root.
    **kwargs : Any
        Extra keyword arguments are accepted and silently ignored so that the
        function can be called uniformly through the :func:`nest2` dispatcher.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If a required observation point name is not found among the stations
        in the overall his-file.
    """
    if not output_path:
        output_path = overall.path
    if not output_file:
        output_file = "hurrywave_his.nc"

    fn_his: str = os.path.join(output_path, output_file)
    ddd = xr.open_dataset(fn_his)

    # Decode station names — HurryWave Fortran kernel stores them as bytes.
    all_stations: List[str] = [
        st.decode("utf-8").strip()
        if isinstance(st, (bytes, np.bytes_))
        else str(st).strip()
        for st in ddd.station_name.values
    ]

    # Build required SnapWave boundary point names from the detail model's GDF.
    point_names: List[str] = [
        f"{obs_point_prefix}_{row['name']}"
        for _ind, row in detail.snapwave_boundary_conditions.gdf.iterrows()
    ]

    times = pd.DatetimeIndex(ddd.point_hm0.coords["time"].values)

    # Map each required point name to a station index (case-insensitive).
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

    npts: int = len(ireq)
    gdf_indices: List[Any] = [
        ind for ind, _row in detail.snapwave_boundary_conditions.gdf.iterrows()
    ]

    # Extract and clip wave parameters to valid physical ranges.
    hm0_all: np.ndarray = ddd.point_hm0.values[:, ireq]  # (nt, npts)
    tp_all: np.ndarray = np.clip(ddd.point_tp.values[:, ireq], 1.0, 25.0)  # seconds
    wavdir_all: np.ndarray = ddd.point_wavdir.values[:, ireq]  # degrees
    dirspr_all: np.ndarray = np.clip(
        ddd.point_dirspr.values[:, ireq], 4.0, 55.0
    )  # degrees

    ddd.close()

    # Build an xarray Dataset with all wave variables for set_timeseries.
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
