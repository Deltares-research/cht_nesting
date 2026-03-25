# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting SFINCS within Delft3D-FM.

Adds observation points to the overall Delft3D-FM model at the water level
boundary point locations of the detail SFINCS model.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincs_in_delft3dfm(overall: Any, detail: Any) -> None:
    """Add observation points to the overall Delft3D-FM model at SFINCS boundary locations.

    Iterates over the water level boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as observation points.
    Points that already exist in *overall* are silently skipped.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded GDF row index.

    When *overall* uses a geographic CRS and its grid extends east of 180 °
    (i.e. uses 0–360 ° longitudes), any transformed longitude < 0 is shifted
    by +360 ° before the point is added.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3D-FM model that receives the new observation points.
        Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``list_observation_names()`` — returns existing observation point names.
        * ``add_observation_point_gdf(x, y, name)`` — adds a point.
        * ``grid.data`` / ``grid.read()`` — used to detect 0–360 ° grids.

    detail : SFINCS
        The fine SFINCS model whose boundary points are used as observation
        locations.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``name`` — prefix string for observation point names.
        * ``boundary_conditions.gdf`` — :class:`~geopandas.GeoDataFrame` of
          boundary points with ``geometry`` column.

    Returns
    -------
    None
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    # Detect whether the overall grid uses 0–360 ° longitudes.
    overall_degrees_west = False
    if overall.crs.is_geographic:
        if overall.grid.data is None:
            overall.grid.read()
        x_max = overall.grid.data.grid.face_coordinates[:, 0].max()
        if x_max > 180.0:
            overall_degrees_west = True

    overall_names = overall.list_observation_names()

    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        if name in overall_names:
            print(f"Observation point {name} already exists in the model")
            continue
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        if x < 0 and overall_degrees_west:
            x += 360.0
        overall.add_observation_point_gdf(x, y, name)
