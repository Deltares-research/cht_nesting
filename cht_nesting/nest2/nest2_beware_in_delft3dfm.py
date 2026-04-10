"""Nest 2 script for nesting BEWARE within Delft3D-FM.

Reads water level and wave time series from the Delft3D-FM model output
and sets them as boundary conditions on the detail BEWARE model.
"""

import os
from typing import Any, Optional

import pandas as pd
import xarray as xr


def nest2_beware_in_delft3dfm(
    overall: Any,
    detail: Any,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    boundary_water_level_correction: float = 0,
    option: Optional[str] = None,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a BEWARE model within a Delft3DFM model.

    Reads water level (flow) or wave parameter (wave) time series from the
    Delft3D-FM output and sets them as boundary conditions on the detail
    BEWARE model.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3D-FM model.
    detail : BEWARE
        The fine BEWARE model that receives boundary conditions.
    output_path : str, optional
        Directory containing the Delft3D-FM output files.
    output_file : str, optional
        Name of the output file. Defaults depend on *option*.
    boundary_water_level_correction : float, optional
        Uniform offset (metres) added to flow boundary water levels.
        Defaults to 0.
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
    if option == "flow":
        if not output_file:
            output_file = "flow_his.nc"

        point_names = [
            f"{detail.name}_{point.name}" for point in detail.flow_boundary_point
        ]

        bzs = overall.read_timeseries_output(
            name_list=point_names, path=output_path, file_name=output_file
        )
        ts = bzs.index
        for icol, point in enumerate(detail.flow_boundary_point):
            point.data = (
                pd.Series(bzs.iloc[:, icol].values, index=ts)
                + boundary_water_level_correction
            )

        if bc_path is not None:
            detail.write_flow_boundary_conditions(
                file_name=os.path.join(bc_path, detail.input.bzsfile)
            )

    if option == "wave":
        if not output_path:
            output_path = overall.path

        if not output_file:
            file_name = os.path.join(output_path, "wavh-wave-nest.nc")
        else:
            file_name = os.path.join(output_path, output_file)

        file_name_dfm = os.path.join(output_path, "flow_his.nc")
        ddd = xr.open_dataset(file_name_dfm)
        stations = ddd.station_name.values
        all_stations = [str(st.strip())[2:-1] for st in stations]

        point_names = [
            f"{detail.name}_{point.name}" for point in detail.wave_boundary_point
        ]

        ddd = xr.load_dataset(file_name)
        times = ddd.Hsig.time.values

        for ip, point in enumerate(detail.wave_boundary_point):
            for ist, st in enumerate(all_stations):
                ireq = -1
                if point_names[ip] == st.lower():
                    ireq = ist
                    break

            if ireq > -1:
                hm0 = ddd.Hsig.values[:, ireq]
                tp = ddd.RTpeak.values[:, ireq]
                wavdir = ddd.Dir.values[:, ireq]
                dirspr = ddd.Dspr.values[:, ireq]

                df = pd.DataFrame(index=times)
                df.insert(0, "hm0", hm0)
                df.insert(1, "tp", tp)
                df.insert(2, "wavdir", wavdir)
                df.insert(3, "dirspr", dirspr)

                point.data = df

        if bc_path is not None:
            detail.write_wave_boundary_conditions(path=bc_path)
