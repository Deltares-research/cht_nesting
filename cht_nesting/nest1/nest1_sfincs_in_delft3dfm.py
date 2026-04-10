"""Nest 1 script for nesting hydromt_sfincs SfincsModel within Delft3D-FM.

Adds observation points to the overall Delft3D-FM model at the water level
boundary point locations of the detail SfincsModel.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincs_in_delft3dfm(overall: Any, detail: Any) -> None:
    """Add observation points to the overall Delft3D-FM model at SfincsModel boundary locations.

    Parameters
    ----------
    overall : Delft3DFM
        The coarse Delft3D-FM model that receives the new observation points.
    detail : hydromt_sfincs.SfincsModel
        The fine SfincsModel whose boundary points are used.
    """
    overall_crs = overall.crs

    if detail.config.get("qtrfile") is not None:
        detail_crs = detail.quadtree_grid.crs
    else:
        detail_crs = detail.grid.crs

    transformer = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)

    # Detect whether the overall grid uses 0-360 longitudes.
    overall_degrees_west = False
    if overall_crs.is_geographic:
        if overall.grid.data is None:
            overall.grid.read()
        x_max = overall.grid.data.grid.face_coordinates[:, 0].max()
        if x_max > 180.0:
            overall_degrees_west = True

    overall_names = overall.list_observation_names()
    boundary_gdf = detail.water_level.gdf

    for ind, row in boundary_gdf.iterrows():
        name = f"{detail.name}_{str(ind + 1).zfill(4)}"
        if name in overall_names:
            print(f"Observation point {name} already exists in the model")
            continue
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        if x < 0 and overall_degrees_west:
            x += 360.0
        overall.add_observation_point_gdf(x, y, name)
