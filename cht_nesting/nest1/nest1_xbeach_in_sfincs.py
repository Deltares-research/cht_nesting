# -*- coding: utf-8 -*-
"""
Nest 1 script for nesting XBeach in SFINCS
"""

from pyproj import Transformer


def nest1_xbeach_in_sfincs(overall, detail):
    
    transformer = Transformer.from_crs(detail.crs,
                                       overall.crs,
                                       always_xy=True)
    
    for ind, point in enumerate(detail.flow_boundary_point):

        name = detail.name + "_" + str(ind + 1).zfill(4)
        x, y = transformer.transform(point.geometry.x,
                                     point.geometry.y)
        overall.observation_points.add_point(x, y, name)
