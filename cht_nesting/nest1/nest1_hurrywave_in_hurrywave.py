# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting HurryWave within HurryWave.

Adds spectral (sp2) observation points to the overall HurryWave model at the
boundary point locations of the detail HurryWave model.
"""

from typing import Any

from pyproj import Transformer


def nest1_hurrywave_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurryWave model at detail boundary locations.

    Iterates over the boundary condition points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as spectral (sp2)
    observation points.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded GDF row index.

    Parameters
    ----------
    overall : HurryWave
        The coarse HurryWave model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``observation_points_sp2.add_point(x, y, name)`` — adds a spectral point.

    detail : HurryWave
        The fine HurryWave model whose boundary condition points are used as
        observation locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``boundary_conditions.gdf`` — :class:`~geopandas.GeoDataFrame` of
          boundary points with ``geometry`` column.

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points_sp2.add_point(x, y, name)
