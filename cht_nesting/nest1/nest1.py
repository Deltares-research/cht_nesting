# -*- coding: utf-8 -*-
"""
Nesting step 1 script

Adds observation points in overall model at boundary points of detail model

"""

from .nest1_beware_in_delft3dfm import nest1_beware_in_delft3dfm
from .nest1_beware_in_hurrywave import nest1_beware_in_hurrywave
from .nest1_beware_in_sfincs import nest1_beware_in_sfincs
from .nest1_delft3dfm_in_delft3dfm import nest1_delft3dfm_in_delft3dfm
from .nest1_hurrywave_in_hurrywave import nest1_hurrywave_in_hurrywave
from .nest1_sfincs_in_delft3dfm import nest1_sfincs_in_delft3dfm
from .nest1_sfincs_in_hurrywave import nest1_sfincs_in_hurrywave
from .nest1_sfincs_in_sfincs import nest1_sfincs_in_sfincs
from .nest1_xbeach_in_hurrywave import nest1_xbeach_in_hurrywave
from .nest1_xbeach_in_sfincs import nest1_xbeach_in_sfincs
from .nest1_xbeach_in_delft3dfm import nest1_xbeach_in_delft3dfm


def nest1(overall, detail, option=None, obs_point_prefix=None):
    
    # Returns a list with observation point objects

    # Check if detail model has attribute name
    if obs_point_prefix is None:
        if not hasattr(detail, "name"):
            detail.name = "nest"
    else:
        detail.name = obs_point_prefix       

    # Get the types of overall and detail classes
    overall_type = overall.__class__.__name__.lower()
    detail_type = detail.__class__.__name__.lower()

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
            
    elif overall_type == "sfincs":
        if detail_type == "sfincs":
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
            # No need to do anything here. BEWARE output points are fixed
            pass
        else:
            print("Nesting step 1 not implemented for this combination of models")
            return False

    return True    
