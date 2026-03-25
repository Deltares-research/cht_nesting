# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting SFINCS (or SfincsModel) within SFINCS (or SfincsModel).

Adds observation points to the overall SFINCS model at the water level boundary
point locations of the detail SFINCS model.  Both ``cht_sfincs.SFINCS`` and
``hydromt_sfincs.SfincsModel`` are supported for either role.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincs_in_sfincs(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SFINCS model at detail boundary locations.

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
    overall : SFINCS or SfincsModel
        The coarse model that receives the new observation points.  Must expose:

        * ``crs`` — :class:`~pyproj.CRS` of the model.
        * ``observation_points.list_names()`` — returns existing point names.
        * ``observation_points.add_point(x, y, name)`` — adds a point.
        * ``grid.data`` / ``grid.read()`` — used to detect 0–360 ° grids.

    detail : SFINCS or SfincsModel
        The fine model whose boundary points are used as observation locations.
        Must expose:

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
        try:
            grid_data = overall.grid.data
            if grid_data is None:
                overall.grid.read()
                grid_data = overall.grid.data
            x_max = grid_data.grid.face_coordinates[:, 0].max()
            if x_max > 180.0:
                overall_degrees_west = True
        except Exception as e:
            print(f"Could not determine if overall model is in degrees west: {e}")

    overall_names = overall.observation_points.list_names()

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
        overall.observation_points.add_point(x, y, name)
