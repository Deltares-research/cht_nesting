# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting XBeach within HurryWave.

Adds observation points to the overall HurryWave model at the wave boundary
point locations of the detail XBeach model.  Points are added to both the
spectral (sp2) and regular observation point collections so that the overall
model writes spectral and bulk wave output at those positions.
"""

from typing import Any

from pyproj import Transformer


def nest1_xbeach_in_hurrywavecht(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurryWave model at XBeach wave boundary locations.

    Iterates over the wave boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them to both
    ``observation_points_sp2`` (spectral output) and
    ``observation_points_regular`` (bulk output) of the overall HurryWave model.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded sequential index.

    Parameters
    ----------
    overall : HurryWave
        The coarse HurryWave model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``observation_points_sp2.add_point(x, y, name)`` — adds a spectral point.
        * ``observation_points_regular.add_point(x, y, name)`` — adds a bulk point.

    detail : XBeach
        The fine XBeach model whose wave boundary points are used as
        observation locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``wave_boundary_point`` — iterable of point objects with a
          ``geometry`` attribute (Shapely ``Point``).

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.wave_boundary_point):
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points_sp2.add_point(x, y, name)
        overall.observation_points_regular.add_point(x, y, name)
