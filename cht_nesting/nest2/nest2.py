"""Nesting step 2 dispatcher.

Reads output from the overall (coarse) model at observation point locations
and sets boundary conditions on the detail (fine) model.
"""

from typing import Any, Optional

from pyproj import CRS


def nest2(
    overall: Any,
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
    return_maximum: bool = False,
) -> Any:
    """Nest a detailed model within an overall model.

    The correct sub-function is selected automatically based on the class names
    of *overall* and *detail*.

    The class name ``SfincsModel`` (from hydromt_sfincs) is mapped to type
    ``"sfincs"``, while ``SFINCS`` (from cht_sfincs) is mapped to ``"sfincscht"``.

    Parameters
    ----------
    overall : Any
        The overall model, or a string type name (e.g. "sfincs", "hurrywave").
    detail : Any
        The detailed model.
    obs_point_prefix : str, optional
        The prefix for observation points.
    output_path : str, optional
        The path to the output files.
    output_file : str, optional
        The name of the output file.
    bc_path : str, optional
        The path to the boundary conditions files.
    bc_file : str, optional
        The name of the boundary conditions file.
    overall_crs : str, optional
        The coordinate reference system of the overall model.
    detail_crs : str, optional
        The coordinate reference system of the detailed model.
    option : str, optional
        The option for nesting.
    boundary_water_level_correction : float
        The correction to apply to the boundary water levels.
    filter_incoming : bool
        Whether to filter incoming data.
    bctype : str
        The type of boundary condition.
    return_maximum : bool
        Whether to return the maximum values.

    Returns
    -------
    Any
        The result of the nesting function.
    """
    # Overall can be a string, because we may not have the overall model as an object

    if isinstance(overall, str):
        if overall == "sfincs":
            # When passed as string "sfincs", assume hydromt-sfincs
            from hydromt_sfincs import SfincsModel
            overall = SfincsModel()
        elif overall == "cht_sfincs":
            from cht_sfincs import SFINCS
            overall = SFINCS()
        elif overall == "hydromt_sfincs":
            from hydromt_sfincs import SfincsModel
            overall = SfincsModel()
        elif overall == "hurrywave" or overall == "hydromt_hurrywave":            
            from hydromt_hurrywave import HurrywaveModel
            overall = HurrywaveModel()
        elif overall == "cht_hurrywave":
            # When passed as string "hurrywave", assume cht_hurrywave for backwards compat
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

    # Resolve types
    overall_type = _resolve_type(overall)
    detail_type = _resolve_type(detail)

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
        "obs_point_prefix": obs_point_prefix,
    }

    if overall_type == "delft3dfm":
        if detail_type == "delft3dfm":
            from .nest2_delft3dfm_in_delft3dfm import nest2_delft3dfm_in_delft3dfm
            nest2_fcn = nest2_delft3dfm_in_delft3dfm
        elif detail_type == "sfincscht":
            from .nest2_sfincscht_in_delft3dfm import nest2_sfincscht_in_delft3dfm
            nest2_fcn = nest2_sfincscht_in_delft3dfm
        elif detail_type == "sfincs":
            from .nest2_sfincs_in_delft3dfm import nest2_sfincs_in_delft3dfm
            nest2_fcn = nest2_sfincs_in_delft3dfm
        elif detail_type == "beware":
            from .nest2_beware_in_delft3dfm import nest2_beware_in_delft3dfm
            nest2_fcn = nest2_beware_in_delft3dfm

    elif overall_type == "sfincs":
        if detail_type == "sfincs":
            from .nest2_sfincs_in_sfincs import nest2_sfincs_in_sfincs
            nest2_fcn = nest2_sfincs_in_sfincs
        elif detail_type == "xbeach":
            from .nest2_xbeach_in_sfincs import nest2_xbeach_in_sfincs
            nest2_fcn = nest2_xbeach_in_sfincs
        elif detail_type == "beware":
            from .nest2_beware_in_sfincs import nest2_beware_in_sfincs
            nest2_fcn = nest2_beware_in_sfincs

    elif overall_type == "sfincscht":
        if detail_type == "sfincscht":
            from .nest2_sfincscht_in_sfincscht import nest2_sfincscht_in_sfincscht
            nest2_fcn = nest2_sfincscht_in_sfincscht
        elif detail_type == "sfincs":
            # SfincsModel nested in cht_sfincs — use hydromt variant
            from .nest2_sfincs_in_sfincs import nest2_sfincs_in_sfincs
            nest2_fcn = nest2_sfincs_in_sfincs
        elif detail_type == "xbeach":
            from .nest2_xbeach_in_sfincscht import nest2_xbeach_in_sfincscht
            nest2_fcn = nest2_xbeach_in_sfincscht
        elif detail_type == "beware":
            from .nest2_beware_in_sfincscht import nest2_beware_in_sfincscht
            nest2_fcn = nest2_beware_in_sfincscht

    elif overall_type == "hurrywave":
        if detail_type == "hurrywave":
            from .nest2_hurrywave_in_hurrywave import nest2_hurrywave_in_hurrywave
            nest2_fcn = nest2_hurrywave_in_hurrywave
        elif detail_type == "xbeach":
            from .nest2_xbeach_in_hurrywave import nest2_xbeach_in_hurrywave
            nest2_fcn = nest2_xbeach_in_hurrywave
        elif detail_type == "sfincscht":
            from .nest2_sfincscht_in_hurrywave import nest2_sfincscht_in_hurrywave
            nest2_fcn = nest2_sfincscht_in_hurrywave
        elif detail_type == "sfincs":
            from .nest2_sfincs_in_hurrywave import nest2_sfincs_in_hurrywave
            nest2_fcn = nest2_sfincs_in_hurrywave
        elif detail_type == "beware":
            from .nest2_beware_in_hurrywave import nest2_beware_in_hurrywave
            nest2_fcn = nest2_beware_in_hurrywave

    elif overall_type == "hurrywavecht":
        if detail_type == "hurrywavecht":
            from .nest2_hurrywavecht_in_hurrywavecht import nest2_hurrywavecht_in_hurrywavecht
            nest2_fcn = nest2_hurrywavecht_in_hurrywavecht
        elif detail_type == "xbeach":
            from .nest2_xbeach_in_hurrywavecht import nest2_xbeach_in_hurrywavecht
            nest2_fcn = nest2_xbeach_in_hurrywavecht
        elif detail_type == "sfincscht":
            from .nest2_sfincscht_in_hurrywavecht import nest2_sfincscht_in_hurrywavecht
            nest2_fcn = nest2_sfincscht_in_hurrywavecht
        elif detail_type == "sfincs":
            from .nest2_sfincs_in_hurrywavecht import nest2_sfincs_in_hurrywavecht
            nest2_fcn = nest2_sfincs_in_hurrywavecht
        elif detail_type == "beware":
            from .nest2_beware_in_hurrywavecht import nest2_beware_in_hurrywavecht
            nest2_fcn = nest2_beware_in_hurrywavecht

    elif overall_type == "beware":
        if detail_type == "sfincscht":
            from .nest2_sfincscht_in_beware import nest2_sfincscht_in_beware
            nest2_fcn = nest2_sfincscht_in_beware
        elif detail_type == "sfincs":
            from .nest2_beware_in_sfincs import nest2_beware_in_sfincs
            nest2_fcn = nest2_beware_in_sfincs

    output = nest2_fcn(overall, detail, **kwargs)

    if output is None:
        output = True

    return output


def _resolve_type(model: Any) -> str:
    """Map a model class name to a nesting type string."""
    class_name = model.__class__.__name__
    if class_name == "SfincsModel":
        return "sfincs"
    elif class_name == "SFINCS":
        return "sfincscht"
    elif class_name == "HurrywaveModel":
        return "hurrywave"
    elif class_name == "HurryWave":
        return "hurrywavecht"
    else:
        return class_name.lower()
