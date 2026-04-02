"""Nest 1 script for nesting cht_sfincs SFINCS within hydromt_hurrywave HurrywaveModel.

Adds observation points to the overall HurrywaveModel at the SnapWave boundary
point locations of the detail SFINCS model (cht_sfincs API).
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincscht_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurrywaveModel at SFINCS SnapWave boundary locations.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel that receives the new observation points.
    detail : cht_sfincs.SFINCS
        The fine SFINCS model whose SnapWave boundary points are used.
    """
    overall_crs = overall.crs
    detail_crs = detail.crs

    transformer = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)

    boundary_gdf = detail.snapwave.boundary_conditions.gdf

    for ind, point in boundary_gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x = point["geometry"].coords[0][0]
        y = point["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points.add_point(x, y, name)
