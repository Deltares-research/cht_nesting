# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting Delft3D-FM in Delft3D-FM
"""

from pyproj import Transformer


def nest1_delft3dfm_in_delft3dfm(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, bnd in enumerate(detail.boundary):        
        for ip, point in enumerate(bnd.point):
            x, y = transformer.transform(point.geometry.x,
                                         point.geometry.y)
            overall.add_observation_point_gdf(x, y, detail.name + "_" + point.name)
