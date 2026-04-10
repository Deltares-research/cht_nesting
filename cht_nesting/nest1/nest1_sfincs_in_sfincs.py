"""Nest 1 script for nesting hydromt_sfincs SfincsModel within SfincsModel.

Adds observation points to the overall SfincsModel at the water level boundary
point locations of the detail SfincsModel.
"""

from typing import Any

from pyproj import Transformer


def nest1_sfincs_in_sfincs(overall: Any, detail: Any) -> None:
    """Add observation points to the overall SfincsModel at detail boundary locations.

    Iterates over the water level boundary points of *detail*, transforms their
    coordinates to the CRS of *overall*, and adds them as observation points.
    Points that already exist in *overall* are silently skipped.

    Each point is named ``<detail.name>_<index>`` where *index* is the
    one-based, zero-padded GDF row index.

    Parameters
    ----------
    overall : hydromt_sfincs.SfincsModel
        The coarse model that receives the new observation points.
    detail : hydromt_sfincs.SfincsModel
        The fine model whose boundary points are used as observation locations.
    """
    if overall.config.get("qtrfile") is not None:
        overall_grid = overall.quadtree_grid
    else:
        overall_grid = overall.grid
    overall_crs = overall_grid.crs

    if detail.config.get("qtrfile") is not None:
        detail_crs = detail.quadtree_grid.crs
    else:
        detail_crs = detail.grid.crs

    transformer = Transformer.from_crs(detail_crs, overall_crs, always_xy=True)

    # Detect whether the overall grid uses 0-360 longitudes.
    overall_degrees_west = False
    if overall_crs.is_geographic:
        try:
            grid_data = overall_grid.data
            if grid_data is None:
                overall_grid.read()
                grid_data = overall_grid.data
            x_max = grid_data.grid.face_coordinates[:, 0].max()
            if x_max > 180.0:
                overall_degrees_west = True
        except Exception as e:
            print(f"Could not determine if overall model is in degrees west: {e}")

    overall_names = overall.observation_points.list_names
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
        overall.observation_points.add_point(x, y, name)
