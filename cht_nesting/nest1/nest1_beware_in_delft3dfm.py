# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting BEWARE in Delft3D-FM
"""

from pyproj import Transformer


def nest1_beware_in_delft3dfm(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.flow_boundary_point):

        name = detail.name + "_" + point.name
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
        overall.add_observation_point_gdf(x, y, name)
