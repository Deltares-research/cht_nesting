# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting SFINCS in HurryWave
"""

from pyproj import Transformer


def nest1_sfincs_in_hurrywave(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)

    for ind, point in detail.snapwave.boundary_conditions.gdf.iterrows():

        name = detail.name + "_" + str(ind + 1).zfill(4)
        x = point["geometry"].coords[0][0]
        y = point["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points_regular.add_point(x, y, name)    
