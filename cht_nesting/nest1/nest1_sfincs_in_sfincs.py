"""Nest 1 script for nesting hydromt_sfincs SfincsModel within SfincsModel.

Adds observation points to the overall SfincsModel at the water level boundary
point locations of the detail SfincsModel.
"""

from typing import Any

import geopandas as gpd
import numpy as np
from pyproj import Transformer


def nest1_sfincs_in_sfincs(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SfincsModel at detail boundary locations.

    Iterates over the water level boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as observation points.
    Points that already exist in *overall* are silently skipped.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded GDF row index.

    Parameters
    ----------
    overall : hydromt_sfincs.SfincsModel
        The coarse model that receives the new observation points.
    detail : hydromt_sfincs.SfincsModel
        The fine model whose boundary points are used as observation locations.
    """
    if overall.config.get("qtrfile") is not None:
        overall_grid = overall.quadtree_grid
    else:
        overall_grid = overall.grid
    overall_crs = overall_grid.crs

    if detail.config.get("qtrfile") is not None:
        detail_crs = detail.quadtree_grid.crs
    else:
        detail_crs = detail.grid.crs

    transformer = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)

    # Detect whether the overall grid uses 0-360 longitudes.
    overall_degrees_west = False
    if overall_crs.is_geographic:
        try:
            grid_data = overall_grid.data
            if grid_data is None:
                overall_grid.read()
                grid_data = overall_grid.data
            x_max = grid_data.grid.face_coordinates[:, 0].max()
            if x_max > 180.0:
                overall_degrees_west = True
        except Exception as e:
            print(f"Could not determine if overall model is in degrees west: {e}")

    overall_names = overall.observation_points.list_names
    boundary_gdf = detail.water_level.gdf

    # Extract coordinates and transform them in a vectorised manner for efficiency.
    coords = np.array([geom.coords[0] for geom in boundary_gdf.geometry.values])
    xs, ys = coords[:, 0], coords[:, 1]

    # Transform all points at once
    xs_t, ys_t = transformer.transform(xs, ys)

    # Apply 0–360 correction vectorised
    if overall_degrees_west:
        xs_t = np.where(xs_t < 0, xs_t + 360.0, xs_t)

    # Create names vectorised
    names = [f"{detail.name}_{i+1:04d}" for i in range(len(xs))]

    # Build GeoDataFrame in one shot
    gdf = gpd.GeoDataFrame(
        {"name": names},
        geometry=gpd.points_from_xy(xs_t, ys_t),
        crs=overall_crs,
    )

    # Filter existing
    gdf = gdf[~gdf["name"].isin(overall_names)]

    if not gdf.empty:
        try:
            overall.observation_points.set(gdf, merge=True, skip_validation=True)
        except Exception as e:
            print(f"Error adding nesting points batch: {e}")