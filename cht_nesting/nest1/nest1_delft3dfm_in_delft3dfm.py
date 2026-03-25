# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting Delft3D-FM within Delft3D-FM.

Adds observation points to the overall Delft3D-FM model at the boundary point
locations of the detail Delft3D-FM model.
"""

from typing import Any

from pyproj import Transformer


def nest1_delft3dfm_in_delft3dfm(overall: Any, detail: Any) -> None:
    """Add observation points to the overall Delft3D-FM model at detail boundary locations.

    Iterates over all boundaries and their points in *detail*, transforms the
    coordinates to the CRS of *overall*, and adds them as observation points.

    Each point is named ``<detail.name>_<point.name>``.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3D-FM model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``add_observation_point_gdf(x, y, name)`` — adds a point.

    detail : Delft3DFM
        The fine Delft3D-FM model whose boundary points are used as observation
        locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``boundary`` — iterable of boundary objects, each with a ``point``
          attribute that is an iterable of point objects with a ``name``
          attribute and a ``geometry`` attribute (Shapely ``Point``).

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, bnd in enumerate(detail.boundary):
        for ip, point in enumerate(bnd.point):
            name = f"{detail.name}_{point.name}"
            x, y = transformer.transform(point.geometry.x, point.geometry.y)
            overall.add_observation_point_gdf(x, y, name)
