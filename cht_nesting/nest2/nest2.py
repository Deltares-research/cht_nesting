# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:56 2021

This module defines the nest2 function for handling nesting in various models.
"""

from pyproj import CRS
from typing import Optional, Any

from .nest2_delft3dfm_in_delft3dfm import nest2_delft3dfm_in_delft3dfm
from .nest2_sfincs_in_delft3dfm import nest2_sfincs_in_delft3dfm
from .nest2_beware_in_delft3dfm import nest2_beware_in_delft3dfm
from .nest2_sfincs_in_sfincs import nest2_sfincs_in_sfincs
from .nest2_xbeach_in_sfincs import nest2_xbeach_in_sfincs
from .nest2_beware_in_sfincs import nest2_beware_in_sfincs
from .nest2_hurrywave_in_hurrywave import nest2_hurrywave_in_hurrywave
from .nest2_xbeach_in_hurrywave import nest2_xbeach_in_hurrywave
from .nest2_sfincs_in_hurrywave import nest2_sfincs_in_hurrywave
from .nest2_beware_in_hurrywave import nest2_beware_in_hurrywave
# from .nest2_sfincs_in_beware import nest2_sfincs_in_beware # This function is no longer compatible with the current version of the cht_sfincs


def nest2(overall: Any,
          detail: Any,
          obs_point_prefix: Optional[str] = None,
          output_path: Optional[str] = None,
          output_file: Optional[str] = None,
          bc_path: Optional[str] = None,
          bc_file: Optional[str] = None,
          overall_crs: Optional[str] = None,
          detail_crs: Optional[str] = None,
          option: Optional[str] = None,
          boundary_water_level_correction: float = 0.0,
          filter_incoming: bool = False,
          bctype: str = "waterlevel",
          return_maximum: bool = False) -> Any:
    """
    Nest a detailed model within an overall model.

    Parameters:
    overall (Any): The overall model.
    detail (Any): The detailed model.
    obs_point_prefix (Optional[str]): The prefix for observation points. Default is None.
    output_path (Optional[str]): The path to the output files. Default is None.
    output_file (Optional[str]): The name of the output file. Default is None.
    bc_path (Optional[str]): The path to the boundary conditions files. Default is None.
    bc_file (Optional[str]): The name of the boundary conditions file. Default is None.
    overall_crs (Optional[str]): The coordinate reference system of the overall model. Default is None.
    detail_crs (Optional[str]): The coordinate reference system of the detailed model. Default is None.
    option (Optional[str]): The option for nesting. Default is None.
    boundary_water_level_correction (float): The correction to apply to the boundary water levels. Default is 0.0.
    filter_incoming (bool): Whether to filter incoming data. Default is False.
    bctype (str): The type of boundary condition. Default is "waterlevel".
    return_maximum (bool): Whether to return the maximum values. Default is False.

    Returns:
    Any: The result of the nesting function.
    """
    # Overall can be a string, because we may not have the overall model as an object
    # We may only know the location of the output files
    # In that case, we need to instantiate the class
    if isinstance(overall, str):
        # Overall is a string, so we need to instantiate the class
        if overall == "sfincs":
            from cht_sfincs import SFINCS
            overall = SFINCS()
        elif overall == "hydromt_sfincs":
            from hydromt_sfincs import SfincsModel
            overall = SfincsModel()
        elif overall == "hurrywave":
            from cht_hurrywave import HurryWave
            overall = HurryWave()
        elif overall == "xbeach":
            from cht_xbeach import XBeach
            overall = XBeach()
        elif overall == "beware":
            from cht_beware import BEWARE
            overall = BEWARE()
        elif overall == "delft3dfm":
            from cht_delft3dfm import Delft3DFM
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

    nest2_fcn = None

    kwargs = {
        "output_path": output_path,
        "output_file": output_file,
        "bc_path": bc_path,
        "boundary_water_level_correction": boundary_water_level_correction,
        "option": option,
        "return_maximum": return_maximum,
        "filter_incoming": filter_incoming,
        "bctype": bctype,
        "obs_point_prefix": obs_point_prefix
    }

    if overall_type == "delft3dfm":
        if detail_type == "delft3dfm":
            nest2_fcn = nest2_delft3dfm_in_delft3dfm
        elif detail_type == "sfincs":
            nest2_fcn = nest2_sfincs_in_delft3dfm
        elif detail_type == "beware":
            nest2_fcn = nest2_beware_in_delft3dfm
    elif overall_type == "sfincs":
        if detail_type == "sfincs":
            nest2_fcn = nest2_sfincs_in_sfincs
        elif detail_type == "xbeach":
            nest2_fcn = nest2_xbeach_in_sfincs
        elif detail_type == "beware":
            nest2_fcn = nest2_beware_in_sfincs
    elif overall_type == "hurrywave":
        if detail_type == "hurrywave":
            nest2_fcn = nest2_hurrywave_in_hurrywave
        elif detail_type == "xbeach":
            nest2_fcn = nest2_xbeach_in_hurrywave
        elif detail_type == "sfincs":
            nest2_fcn = nest2_sfincs_in_hurrywave
        elif detail_type == "beware":
            nest2_fcn = nest2_beware_in_hurrywave
    elif overall_type == "beware":
        if detail_type == "sfincs":
            nest2_fcn = nest2_sfincs_in_beware

    output = nest2_fcn(overall, detail, **kwargs)

    if output is None:        
        output = True

    return output
