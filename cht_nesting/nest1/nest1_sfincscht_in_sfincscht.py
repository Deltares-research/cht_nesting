"""Nest 1 script for nesting cht_sfincs SFINCS within cht_sfincs SFINCS.

Adds observation points to the overall SFINCS model at the water level boundary
point locations of the detail SFINCS model.  Both models use the cht_sfincs API.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincscht_in_sfincscht(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SFINCS model at detail boundary locations.

    Parameters
    ----------
    overall : cht_sfincs.SFINCS
        The coarse model that receives the new observation points.
    detail : cht_sfincs.SFINCS
        The fine model whose boundary points are used as observation locations.
    """
    overall_crs = overall.crs
    detail_crs = detail.crs

    transformer = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)

    # Detect whether the overall grid uses 0-360 longitudes.
    overall_degrees_west = False
    if overall_crs.is_geographic:
        try:
            grid_data = overall.grid.data
            if grid_data is None:
                overall.grid.read()
                grid_data = overall.grid.data
            x_max = grid_data.grid.face_coordinates[:, 0].max()
            if x_max > 180.0:
                overall_degrees_west = True
        except Exception as e:
            print(f"Could not determine if overall model is in degrees west: {e}")

    overall_names = overall.observation_points.list_names()
    boundary_gdf = detail.boundary_conditions.gdf

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
        overall.observation_points.add_point(x, y, name)
