"""Nest 1 script for nesting XBeach within hydromt_sfincs SfincsModel.

Adds observation points to the overall SfincsModel at the flow boundary point
locations of the detail XBeach model.
"""

from typing import Any

import geopandas as gpd
from pyproj import Transformer


def nest1_xbeach_in_sfincs(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SfincsModel at XBeach flow boundary locations.

    Parameters
    ----------
    overall : hydromt_sfincs.SfincsModel
        The coarse SfincsModel that receives the new observation points.
    detail : XBeach
        The fine XBeach model whose flow boundary points are used.
    """
    if overall.config.get("qtrfile") is not None:
        overall_crs = overall.quadtree_grid.crs
    else:
        overall_crs = overall.grid.crs

    transformer = Transformer.from_crs(detail.crs, overall_crs, always_xy=True)

    xs = [p.geometry.x for p in detail.flow_boundary_point]
    ys = [p.geometry.y for p in detail.flow_boundary_point]
    names = [f"{detail.name}_{i+1:04d}" for i in range(len(xs))]

    xs, ys = transformer.transform(xs, ys)

    gdf = gpd.GeoDataFrame(
        {"name": names},
        geometry=gpd.points_from_xy(xs, ys),
        crs=overall_crs,
    )

    overall.observation_points.set(gdf, merge=True, skip_validation=True)