# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting BEWARE within HurryWave.

Adds observation points to the overall HurryWave model at the wave boundary
point locations of the detail BEWARE model.
"""

from typing import Any

from pyproj import Transformer


def nest1_beware_in_hurrywavecht(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurryWave model at BEWARE wave boundary locations.

    Iterates over the wave boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as regular observation
    points to the HurryWave model.

    Each point is named ``<detail.name>_<point.name>``.

    Parameters
    ----------
    overall : HurryWave
        The coarse HurryWave model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``observation_points_regular.add_point(x, y, name)`` — adds a point.

    detail : BEWARE
        The fine BEWARE model whose wave boundary points are used as
        observation locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``wave_boundary_point`` — iterable of point objects, each with a
          ``name`` attribute and a ``geometry`` attribute (Shapely ``Point``).

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.wave_boundary_point):
        name = f"{detail.name}_{point.name}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points_regular.add_point(x, y, name)
