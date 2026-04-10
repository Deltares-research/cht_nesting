"""Nest 1 script for nesting BEWARE within hydromt_hurrywave HurrywaveModel.

Adds observation points to the overall HurrywaveModel at the wave boundary
point locations of the detail BEWARE model.
"""

from typing import Any

from pyproj import Transformer


def nest1_beware_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurrywaveModel at BEWARE wave boundary locations.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel that receives the new observation points.
    detail : BEWARE
        The fine BEWARE model whose wave boundary points are used.
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.wave_boundary_point):
        name = f"{detail.name}_{point.name}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points.add_point(x, y, name)
