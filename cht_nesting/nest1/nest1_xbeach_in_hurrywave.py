"""Nest 1 script for nesting XBeach within hydromt_hurrywave HurrywaveModel.

Adds observation points to the overall HurrywaveModel at the wave boundary
point locations of the detail XBeach model. Points are added to both the
spectral and regular observation point collections.
"""

from typing import Any

from pyproj import Transformer


def nest1_xbeach_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurrywaveModel at XBeach wave boundary locations.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse HurrywaveModel that receives the new observation points.
    detail : XBeach
        The fine XBeach model whose wave boundary points are used.
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, point in enumerate(detail.wave_boundary_point):
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x, y = transformer.transform(point.geometry.x, point.geometry.y)
        overall.observation_points_spectra.add_point(x, y, name)
        overall.observation_points.add_point(x, y, name)
