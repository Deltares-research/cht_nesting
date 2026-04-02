"""Nest 1 script for nesting XBeach within cht_sfincs SFINCS.

Adds observation points to the overall SFINCS model (cht_sfincs API) at the
flow boundary point locations of the detail XBeach model.
"""

from typing import Any

from pyproj import Transformer


def nest1_xbeach_in_sfincscht(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SFINCS model at XBeach flow boundary locations.

    Parameters
    ----------
    overall : cht_sfincs.SFINCS
        The coarse SFINCS model that receives the new observation points.
    detail : XBeach
        The fine XBeach model whose flow boundary points are used.
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.flow_boundary_point):
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points.add_point(x, y, name)
