# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting SFINCS in Delft3D-FM
"""

from pyproj import Transformer


def nest1_sfincs_in_delft3dfm(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    # Check if overall model is in geographic coordinates. If so, make sure that longitude is correct. Either 0-360 or -180-180.
    overall_degrees_west = False
    if overall.crs.is_geographic:
        # Get max longitude of the overall model
        if overall.grid.data is None:
            overall.grid.read()
        x_max = overall.grid.data.grid.face_coordinates[:, 0].max()
        if x_max > 180.0:
            overall_degrees_west = True

    # Get list of names of the observation points
    overall_names = overall.list_observation_names()
    # Loop over all points in the gdf of the boundary conditions
    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = detail.name + "_" + str(ind + 1).zfill(4)
        if name in overall_names:
            print(f"Observation point {name} already exists in the model")
            continue
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        if x < 0 and overall_degrees_west:
            x += 360
        overall.add_observation_point_gdf(x, y, name)
