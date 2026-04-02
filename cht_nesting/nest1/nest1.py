"""Nesting step 1 dispatcher.

Adds observation points to the overall model at the boundary point locations
of the detail model, so that the overall model writes output at those positions
for use in nesting step 2.
"""

from typing import Any, Optional


def nest1(
    overall: Any,
    detail: Any,
    option: Optional[str] = None,
    obs_point_prefix: Optional[str] = None,
) -> bool:
    """Add observation points to the overall model at detail model boundary locations.

    The correct sub-function is selected automatically based on the class names of
    *overall* and *detail*.

    The class name ``SfincsModel`` (from hydromt_sfincs) is mapped to type
    ``"sfincs"``, while ``SFINCS`` (from cht_sfincs) is mapped to ``"sfincscht"``.

    Parameters
    ----------
    overall : Any
        The coarse model that receives the new observation points.
    detail : Any
        The fine model whose boundary points are used as observation point
        locations.
    option : str, optional
        Reserved for future use; currently unused.
    obs_point_prefix : str, optional
        Prefix for observation point names.

    Returns
    -------
    bool
        ``True`` if observation points were added successfully, ``False`` if the
        overall/detail model combination is not supported.
    """
    # Resolve observation point prefix.
    if obs_point_prefix is None:
        if not hasattr(detail, "name"):
            detail.name = "nest"
    else:
        detail.name = obs_point_prefix

    overall_type: str = _resolve_type(overall)
    detail_type: str = _resolve_type(detail)

    nest1_fcn = None

    if overall_type == "delft3dfm":
        if detail_type == "delft3dfm":
            from .nest1_delft3dfm_in_delft3dfm import nest1_delft3dfm_in_delft3dfm
            nest1_fcn = nest1_delft3dfm_in_delft3dfm
        elif detail_type == "sfincs":
            from .nest1_sfincs_in_delft3dfm import nest1_sfincs_in_delft3dfm
            nest1_fcn = nest1_sfincs_in_delft3dfm
        elif detail_type == "sfincscht":
            from .nest1_sfincscht_in_delft3dfm import nest1_sfincscht_in_delft3dfm
            nest1_fcn = nest1_sfincscht_in_delft3dfm
        elif detail_type == "beware":
            from .nest1_beware_in_delft3dfm import nest1_beware_in_delft3dfm
            nest1_fcn = nest1_beware_in_delft3dfm
        elif detail_type == "xbeach":
            from .nest1_xbeach_in_delft3dfm import nest1_xbeach_in_delft3dfm
            nest1_fcn = nest1_xbeach_in_delft3dfm

    elif overall_type == "sfincs":
        if detail_type == "sfincs":
            from .nest1_sfincs_in_sfincs import nest1_sfincs_in_sfincs
            nest1_fcn = nest1_sfincs_in_sfincs
        elif detail_type == "xbeach":
            from .nest1_xbeach_in_sfincs import nest1_xbeach_in_sfincs
            nest1_fcn = nest1_xbeach_in_sfincs
        elif detail_type == "beware":
            from .nest1_beware_in_sfincs import nest1_beware_in_sfincs
            nest1_fcn = nest1_beware_in_sfincs

    elif overall_type == "sfincscht":
        if detail_type == "sfincscht":
            from .nest1_sfincscht_in_sfincscht import nest1_sfincscht_in_sfincscht
            nest1_fcn = nest1_sfincscht_in_sfincscht
        elif detail_type == "xbeach":
            from .nest1_xbeach_in_sfincscht import nest1_xbeach_in_sfincscht
            nest1_fcn = nest1_xbeach_in_sfincscht
        elif detail_type == "beware":
            from .nest1_beware_in_sfincscht import nest1_beware_in_sfincscht
            nest1_fcn = nest1_beware_in_sfincscht

    elif overall_type == "hurrywave":
        if detail_type == "hurrywave":
            from .nest1_hurrywave_in_hurrywave import nest1_hurrywave_in_hurrywave
            nest1_fcn = nest1_hurrywave_in_hurrywave
        elif detail_type == "xbeach":
            from .nest1_xbeach_in_hurrywave import nest1_xbeach_in_hurrywave
            nest1_fcn = nest1_xbeach_in_hurrywave
        elif detail_type == "sfincs":
            from .nest1_sfincs_in_hurrywave import nest1_sfincs_in_hurrywave
            nest1_fcn = nest1_sfincs_in_hurrywave
        elif detail_type == "sfincscht":
            from .nest1_sfincscht_in_hurrywave import nest1_sfincscht_in_hurrywave
            nest1_fcn = nest1_sfincscht_in_hurrywave
        elif detail_type == "beware":
            from .nest1_beware_in_hurrywave import nest1_beware_in_hurrywave
            nest1_fcn = nest1_beware_in_hurrywave

    elif overall_type == "hurrywavecht":
        if detail_type == "hurrywavecht":
            from .nest1_hurrywavecht_in_hurrywavecht import nest1_hurrywavecht_in_hurrywavecht
            nest1_fcn = nest1_hurrywavecht_in_hurrywavecht
        elif detail_type == "xbeach":
            from .nest1_xbeach_in_hurrywavecht import nest1_xbeach_in_hurrywavecht
            nest1_fcn = nest1_xbeach_in_hurrywavecht
        elif detail_type == "sfincs":
            from .nest1_sfincs_in_hurrywavecht import nest1_sfincs_in_hurrywavecht
            nest1_fcn = nest1_sfincs_in_hurrywavecht
        elif detail_type == "sfincscht":
            from .nest1_sfincscht_in_hurrywavecht import nest1_sfincscht_in_hurrywavecht
            nest1_fcn = nest1_sfincscht_in_hurrywavecht
        elif detail_type == "beware":
            from .nest1_beware_in_hurrywavecht import nest1_beware_in_hurrywavecht
            nest1_fcn = nest1_beware_in_hurrywavecht

    elif overall_type == "beware":
        if detail_type in ("sfincs", "sfincscht"):
            # BEWARE output points are fixed; nothing to add.
            return True

    if nest1_fcn is None:
        print("Nesting step 1 not implemented for this combination of models")
        return False

    nest1_fcn(overall, detail)
    return True


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
