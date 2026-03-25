# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting XBeach within SFINCS.

Adds observation points to the overall SFINCS model at the flow boundary point
locations of the detail XBeach model.
"""

from typing import Any

from pyproj import Transformer


def nest1_xbeach_in_sfincs(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SFINCS model at XBeach flow boundary locations.

    Iterates over the flow boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as observation points.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded sequential index.

    Parameters
    ----------
    overall : SFINCS
        The coarse SFINCS model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``observation_points.add_point(x, y, name)`` — adds a point.

    detail : XBeach
        The fine XBeach model whose flow boundary points are used as
        observation locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``flow_boundary_point`` — iterable of point objects with a
          ``geometry`` attribute (Shapely ``Point``).

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.flow_boundary_point):
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points.add_point(x, y, name)
