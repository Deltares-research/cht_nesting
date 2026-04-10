"""
Nest 2 script for nesting a hydromt_sfincs SfincsModel within a Delft3DFM model.

Reads water level time series from the Delft3DFM ``flow_his.nc`` output via
``overall.read_timeseries_output()`` and sets them as water level boundary
conditions on the detail SfincsModel.
"""

import os
from typing import Any, List, Optional

import pandas as pd


def nest2_sfincs_in_delft3dfm(
    overall: Any,
    detail: Any,
    obs_point_prefix: Optional[str] = None,
    output_path: Optional[str] = None,
    output_file: Optional[str] = None,
    boundary_water_level_correction: float = 0.0,
    bc_path: Optional[str] = None,
    **kwargs: Any,
) -> None:
    """Nest a hydromt_sfincs SfincsModel within a Delft3DFM model.

    Reads water level time series from the Delft3DFM his output via
    ``overall.read_timeseries_output()`` and applies them as water level
    boundary conditions on the detail SfincsModel via
    ``detail.water_level.set_timeseries()``.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3DFM model.  Must expose a
        ``read_timeseries_output(name_list, path, file_name)`` method that
        returns a :class:`~pandas.DataFrame` with a
        :class:`~pandas.DatetimeIndex` and one column per station.
    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel that receives the water level boundary conditions.
    obs_point_prefix : str, optional
        Prefix prepended to each detail boundary point name when requesting
        stations from the overall model output.
    output_path : str, optional
        Directory containing the Delft3DFM output files.
    output_file : str, optional
        Name of the Delft3DFM his output file.  Defaults to ``"flow_his.nc"``.
    boundary_water_level_correction : float, optional
        Uniform offset (metres) added to every boundary time series.  Useful
        for datum corrections.  Defaults to ``0.0``.
    bc_path : str, optional
        If provided, write the boundary conditions to disk by calling
        ``detail.water_level.write()`` after setting.  The value of *bc_path*
        is not used directly — writing goes to the model's configured root.
    **kwargs : Any
        Extra keyword arguments are accepted and silently ignored so that the
        function can be called uniformly through the :func:`nest2` dispatcher.

    Returns
    -------
    None

    Raises
    ------
    KeyError
        If a required point name is not present in the Delft3DFM output.
    """
    if not output_file:
        output_file = "flow_his.nc"
    if not obs_point_prefix:
        obs_point_prefix = detail.name

    # The name is the index with 4 leading zeros
    point_names: List[str] = [
        f"{obs_point_prefix}_{str(ind + 1).zfill(4)}"
        for ind, _row in detail.water_level.gdf.iterrows()
    ]

    output_filepath: str = os.path.join(output_path, output_file)
    bzs: pd.DataFrame = overall.read_timeseries_output(
        name_list=point_names,
        path=output_path,
        file_name=output_filepath,
    )
    ts = bzs.index

    # Build the BC DataFrame: integer GDF-index columns matched to time index.
    df_bcs = pd.DataFrame(index=ts)
    for icol, (ind, _row) in enumerate(detail.water_level.gdf.iterrows()):
        df_bcs[ind] = bzs.iloc[:, icol].values + boundary_water_level_correction

    detail.water_level.set_timeseries(df_bcs)

    if bc_path is not None:
        detail.water_level.write()
