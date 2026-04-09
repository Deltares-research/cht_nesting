"""Nest 2 script for nesting hydromt_hurrywave HurrywaveModel within HurrywaveModel.

Reads spectral output from the overall HurrywaveModel and sets spectral
boundary conditions on the detail HurrywaveModel.
"""

import os
from typing import Any, List, Optional

import numpy as np
import xarray as xr


def nest2_hurrywave_in_hurrywave(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs,
) -> None:
    """Nest a HurrywaveModel within another HurrywaveModel.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel whose spectral output is read.
    detail : hydromt_hurrywave.HurrywaveModel
        The fine HurrywaveModel that receives the spectral boundary conditions.
    obs_point_prefix : str, optional
        Prefix for observation point names.
    output_path : str, optional
        Directory containing the overall model's output files.
    output_file : str, optional
        Name of the spectral output file. Defaults to "hurrywave_sp2.nc".
    bc_path : str, optional
        If provided, write boundary conditions to disk.
    """
    if not output_path:
        output_path = str(overall.root.path)
    if not output_file:
        output_file = "hurrywave_sp2.nc"

    file_name = os.path.join(output_path, output_file)

    # Read boundary points from detail model
    bndfile = detail.config.get("bndfile")
    if bndfile is None:
        detail.config.set("bndfile", "hurrywave.bnd")
    detail.boundary_conditions._read_boundary_points()
    bnd_gdf = detail.boundary_conditions.gdf

    if len(bnd_gdf) == 0:
        print("Warning: no boundary points found in detail model")
        return

    # Open the overall spectral output
    ddd = xr.open_dataset(file_name)
    stations = ddd.station_name.values
    all_stations: List[str] = [
        st.decode("utf-8").strip()
        if isinstance(st, (bytes, np.bytes_))
        else str(st).strip()
        for st in stations
    ]

    # Build required point names
    point_names: List[str] = [
        f"{obs_point_prefix}_{row['name']}" for _, row in bnd_gdf.iterrows()
    ]

    # Map point names to station indices
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
                f"Observation point '{pname}' not found in '{file_name}'. "
                f"Available stations: {all_stations}"
            )

    times = ddd.point_spectrum2d.coords["time"].values
    sigma = ddd.point_spectrum2d.coords["sigma"].values
    theta = ddd.point_spectrum2d.coords["theta"].values

    # Extract spectra for required stations
    sp2 = ddd.point_spectrum2d.values[:, ireq, :, :]  # (time, stations, theta, sigma)

    # Get boundary point coordinates
    xs = np.array([row.geometry.x for _, row in bnd_gdf.iterrows()])
    ys = np.array([row.geometry.y for _, row in bnd_gdf.iterrows()])
    station_names = [row["name"] for _, row in bnd_gdf.iterrows()]

    ddd.close()

    # Build a proper spectra dataset for set_spectra()
    ds_spectra = xr.Dataset(
        data_vars={
            "point_spectrum2d": (["time", "stations", "theta", "sigma"], sp2),
            "station_x": ("stations", xs),
            "station_y": ("stations", ys),
        },
        coords={
            "time": times,
            "stations": station_names,
            "theta": theta,
            "sigma": sigma,
        },
    )

    detail.boundary_conditions.set_spectra(ds_spectra)

    if bc_path is not None:
        detail.boundary_conditions._write_boundary_spectra()
