"""Nest 2 script for nesting a SFINCS model within a BEWARE model.

Reads water level (flow) or infragravity wave forcing (wave) from BEWARE output
and sets them as boundary conditions on the detail SFINCS model.
"""

import datetime
import os
from typing import Any, Optional

import numpy as np
import pandas as pd
import xarray as xr
from cht_sfincs.sfincs import FlowBoundaryPoint, WaveMakerForcingPoint
from pyproj import Transformer


def nest2_sfincs_in_beware(
    overall: Any,
    detail: Any,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    boundary_water_level_correction: float = 0,
    option: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a SFINCS model within a BEWARE model.

    Finds BEWARE offshore (flow) or coastal (wave) output locations inside
    the detail model's bounding box, extracts the relevant time series, and
    sets them as boundary conditions on the SFINCS model.

    Parameters
    ----------
    overall : BEWARE
        The coarse BEWARE model.  Must expose ``path`` and ``crs``.
    detail : SFINCS
        The fine SFINCS model that receives boundary conditions.
    output_path : str, optional
        Directory containing the BEWARE output files.
        Defaults to ``overall.path``.
    output_file : str, optional
        Name of the BEWARE output file.  Defaults to ``"beware_his.nc"``.
    boundary_water_level_correction : float, optional
        Uniform offset (metres) added to flow boundary water levels.
        Defaults to ``0``.
    option : str, optional
        Nesting option: ``"flow"`` or ``"wave"``.
    bc_path : str, optional
        If provided, write boundary conditions to disk.
    **kwargs : Any
        Extra keyword arguments are silently ignored.

    Returns
    -------
    None
    """
    if not output_file:
        output_file = "beware_his.nc"
    # Get bounding box for sfincs model
    # Convert bbox to beware crs

    x_range, y_range = detail.bounding_box(crs=overall.crs)
    dx = (x_range[1] - x_range[0]) / 10
    dy = (y_range[1] - y_range[0]) / 10
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
        inear = np.where(
            (xb > x_range[0])
            & (xb < x_range[1])
            & (yb > y_range[0])
            & (yb < y_range[1])
        )
        xb = xb[inear]
        yb = yb[inear]
        nb = xb.size

        # Clear existing flow boundary points
        detail.flow_boundary_point = []

        # Convert to coordinate system of detail model
        transformer = Transformer.from_crs(overall.crs, detail.crs, always_xy=True)

        for ip in range(nb):
            name = str(ip + 1).zfill(4)
            x, y = transformer.transform(xb[ip], yb[ip])
            point = FlowBoundaryPoint(x, y, name=name)
            detail.flow_boundary_point.append(point)

        # Extract data and set water level boundary conditions
        tref = datetime.datetime(1970, 1, 1)
        tsec = ddd.time.values  # array of int64
        times = tref + tsec * datetime.timedelta(seconds=1)

        for ip, point in enumerate(detail.flow_boundary_point):
            point.data = (
                pd.Series(ddd.WL.values[inear[0][ip], :], index=times)
                + pd.Series(ddd.R2_setup.values[inear[0][ip], :], index=times)
                + boundary_water_level_correction
            )

        if bc_path is not None:
            detail.write_flow_boundary_conditions(
                file_name=os.path.join(bc_path, detail.input.bzsfile)
            )

    elif option == "wave":
        xb = ddd.x_coast.values
        yb = ddd.y_coast.values

        # Find beware locations in bounding box
        inear = np.where(
            (xb > x_range[0])
            & (xb < x_range[1])
            & (yb > y_range[0])
            & (yb < y_range[1])
        )
        xb = xb[inear]
        yb = yb[inear]
        nb = xb.size

        # Clear existing flow boundary points
        detail.wavemaker_forcing_point = []

        # Find wave boundary forcing at intersection wfp file and beware transects
        # Load sfincs.wfp, find intersect with wfp and beware transect.

        # Convert to coordinate system of detail model
        transformer = Transformer.from_crs(overall.crs, detail.crs, always_xy=True)

        for ip in range(nb):
            name = str(ip + 1).zfill(4)
            x, y = transformer.transform(xb[ip], yb[ip])
            point = WaveMakerForcingPoint(x, y, name=name)
            detail.wavemaker_forcing_point.append(point)

        # Extract data and set water level boundary conditions
        tref = datetime.datetime(1970, 1, 1)
        tsec = ddd.time.values  # array of int64
        times = tref + tsec * datetime.timedelta(seconds=1)

        for ip, point in enumerate(detail.wavemaker_forcing_point):
            df = pd.DataFrame()
            df["hm0_ig"] = ddd.hm0_ig.values[inear[0][ip], :]
            df["tp_ig"] = ddd.tp_ig.values[inear[0][ip], :]
            df["setup"] = ddd.setup.values[inear[0][ip], :]
            df["time"] = times
            df = df.set_index("time")

            df["hm0_ig"] = df["hm0_ig"].replace(np.nan, 0.1)
            df["tp_ig"] = df["tp_ig"].replace(np.nan, 60.0)
            df["setup"] = df["setup"].replace(np.nan, 0.0)

            point.data = df

    ddd.close()

    if bc_path is not None:
        detail.write_whi_file(file_name=os.path.join(bc_path, detail.input.whifile))
        detail.write_wti_file(file_name=os.path.join(bc_path, detail.input.wtifile))
