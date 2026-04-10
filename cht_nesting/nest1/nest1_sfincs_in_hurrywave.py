"""Nest 1 script for nesting hydromt_sfincs SfincsModel within hydromt_hurrywave HurrywaveModel.

Adds observation points to the overall HurrywaveModel at the SnapWave boundary
point locations of the detail SfincsModel.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincs_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurrywaveModel at SfincsModel SnapWave boundary locations.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel that receives the new observation points.
    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel whose SnapWave boundary points are used.
    """
    overall_crs = overall.crs

    if detail.config.get("qtrfile") is not None:
        detail_crs = detail.quadtree_grid.crs
    else:
        detail_crs = detail.grid.crs

    transformer = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)

    boundary_gdf = detail.snapwave_boundary_conditions.gdf

    for ind, point in boundary_gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x = point["geometry"].coords[0][0]
        y = point["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points.add_point(x, y, name)
