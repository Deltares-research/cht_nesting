# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting SFINCS within HurryWave.

Adds observation points to the overall HurryWave model at the SnapWave boundary
point locations of the detail SFINCS model.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincs_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurryWave model at SFINCS SnapWave boundary locations.

    Iterates over the SnapWave boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as regular observation
    points to the HurryWave model.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded GDF row index.

    Parameters
    ----------
    overall : HurryWave
        The coarse HurryWave model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``observation_points_regular.add_point(x, y, name)`` — adds a point.

    detail : SFINCS
        The fine SFINCS model whose SnapWave boundary points are used as
        observation locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``snapwave.boundary_conditions.gdf`` —
          :class:`~geopandas.GeoDataFrame` of SnapWave boundary points with
          ``geometry`` column.

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in detail.snapwave.boundary_conditions.gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x = point["geometry"].coords[0][0]
        y = point["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points_regular.add_point(x, y, name)
