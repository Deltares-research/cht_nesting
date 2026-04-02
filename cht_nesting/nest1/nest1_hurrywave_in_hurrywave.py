"""Nest 1 script for nesting hydromt_hurrywave HurrywaveModel within HurrywaveModel.

Adds spectral observation points to the overall HurrywaveModel at the
boundary point locations of the detail HurrywaveModel.
"""

from typing import Any

from pyproj import Transformer


def nest1_hurrywave_in_hurrywave(overall: Any, detail: Any) -> None:
    """Add observation points to the overall HurrywaveModel at detail boundary locations.

    Parameters
    ----------
    overall : hydromt_hurrywave.HurrywaveModel
        The coarse model that receives the new observation points.
    detail : hydromt_hurrywave.HurrywaveModel
        The fine model whose boundary condition points are used.
    """
    transformer = Transformer.from_crs(detail.crs, overall.crs, always_xy=True)

    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points_spectra.add_point(x, y, name)
