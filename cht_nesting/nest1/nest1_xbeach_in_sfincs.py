"""Nest 1 script for nesting XBeach within hydromt_sfincs SfincsModel.

Adds observation points to the overall SfincsModel at the flow boundary point
locations of the detail XBeach model.
"""

from typing import Any

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

    for ind, point in enumerate(detail.flow_boundary_point):
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points.add_point(x, y, name)
