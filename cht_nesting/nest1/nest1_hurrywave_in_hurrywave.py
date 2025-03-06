# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting HurryWave in HurryWave
"""

from pyproj import Transformer


def nest1_hurrywave_in_hurrywave(overall, detail):

    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)

    for ind, row in detail.boundary_conditions.gdf.iterrows():
        name = detail.name + "_" + str(ind + 1).zfill(4)
        x = row["geometry"].coords[0][0]
        y = row["geometry"].coords[0][1]
        x, y = transformer.transform(x, y)
        overall.observation_points_sp2.add_point(x, y, name)    
