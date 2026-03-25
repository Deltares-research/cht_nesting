# -*- coding: utf-8 -*-
"""
Nesting step 1 dispatcher.

Adds observation points to the overall model at the boundary point locations
of the detail model, so that the overall model writes output at those positions
for use in nesting step 2.
"""

from typing import Any, Optional

from .nest1_beware_in_delft3dfm import nest1_beware_in_delft3dfm
from .nest1_beware_in_hurrywave import nest1_beware_in_hurrywave
from .nest1_beware_in_sfincs import nest1_beware_in_sfincs
from .nest1_delft3dfm_in_delft3dfm import nest1_delft3dfm_in_delft3dfm
from .nest1_hurrywave_in_hurrywave import nest1_hurrywave_in_hurrywave
from .nest1_sfincs_in_delft3dfm import nest1_sfincs_in_delft3dfm
from .nest1_sfincs_in_hurrywave import nest1_sfincs_in_hurrywave
from .nest1_sfincs_in_sfincs import nest1_sfincs_in_sfincs
from .nest1_xbeach_in_delft3dfm import nest1_xbeach_in_delft3dfm
from .nest1_xbeach_in_hurrywave import nest1_xbeach_in_hurrywave
from .nest1_xbeach_in_sfincs import nest1_xbeach_in_sfincs


def nest1(
    overall: Any,
    detail: Any,
    option: Optional[str] = None,
    obs_point_prefix: Optional[str] = None,
) -> bool:
    """Add observation points to the overall model at detail model boundary locations.

    This is nesting step 1.  For each boundary point of the detail (fine) model,
    a corresponding observation point is added to the overall (coarse) model using
    the name ``<obs_point_prefix>_<index>``.  The overall model will then write
    output at those positions, which is consumed by nesting step 2 (:func:`nest2`).

    The correct sub-function is selected automatically based on the class names of
    *overall* and *detail*.

    Parameters
    ----------
    overall : Any
        The coarse model that receives the new observation points.  Supported
        types: ``SFINCS`` / ``SfincsModel``, ``HurryWave``, ``Delft3DFM``,
        ``BEWARE``.
    detail : Any
        The fine model whose boundary points are used as observation point
        locations.  Supported types: ``SFINCS`` / ``SfincsModel``, ``XBeach``,
        ``HurryWave``, ``Delft3DFM``, ``BEWARE``.
    option : str, optional
        Reserved for future use; currently unused.
    obs_point_prefix : str, optional
        Prefix for observation point names.  If provided, sets ``detail.name``
        to this value.  If omitted, ``detail.name`` is used as-is; if *detail*
        has no ``name`` attribute it is set to ``"nest"``.

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

    overall_type: str = overall.__class__.__name__.lower()
    detail_type: str = detail.__class__.__name__.lower()

    if overall_type == "delft3dfm":
        if detail_type == "delft3dfm":
            nest1_delft3dfm_in_delft3dfm(overall, detail)
        elif detail_type == "sfincs":
            nest1_sfincs_in_delft3dfm(overall, detail)
        elif detail_type == "beware":
            nest1_beware_in_delft3dfm(overall, detail)
        elif detail_type == "xbeach":
            nest1_xbeach_in_delft3dfm(overall, detail)
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    elif overall_type in ("sfincs", "sfincsmodel"):
        if detail_type in ("sfincs", "sfincsmodel"):
            nest1_sfincs_in_sfincs(overall, detail)
        elif detail_type == "xbeach":
            nest1_xbeach_in_sfincs(overall, detail)
        elif detail_type == "beware":
            nest1_beware_in_sfincs(overall, detail)
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    elif overall_type == "hurrywave":
        if detail_type == "hurrywave":
            nest1_hurrywave_in_hurrywave(overall, detail)
        elif detail_type == "xbeach":
            nest1_xbeach_in_hurrywave(overall, detail)
        elif detail_type == "sfincs":
            nest1_sfincs_in_hurrywave(overall, detail)
        elif detail_type == "beware":
            nest1_beware_in_hurrywave(overall, detail)
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    elif overall_type == "beware":
        if detail_type == "sfincs":
            # BEWARE output points are fixed; nothing to add.
            pass
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    return True
