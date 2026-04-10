"""Nest 1 script for nesting XBeach within Delft3D-FM.

Adds observation points to the overall Delft3D-FM model at the flow boundary
point locations of the detail XBeach model.
"""

from typing import Any

from pyproj import Transformer


def nest1_xbeach_in_delft3dfm(overall: Any, detail: Any) -> None:
    """Add observation points to the overall Delft3D-FM model at XBeach flow boundary locations.

    Iterates over the flow boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as observation points
    using ``overall.add_observation_point_gdf()``.

    Each point is named ``<detail.name>_<point.name>``.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3D-FM model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``add_observation_point_gdf(x, y, name)`` — adds a point.

    detail : XBeach
        The fine XBeach model whose flow boundary points are used as
        observation locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``flow_boundary_point`` — iterable of point objects, each with a
          ``name`` attribute and a ``geometry`` attribute (Shapely ``Point``).

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.flow_boundary_point):
        name = f"{detail.name}_{point.name}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.add_observation_point_gdf(x, y, name)
