"""
Handle this with more tried-and-true tom code
"""
import os
import xarray as xr
import six
from matplotlib import cm
import scipy.optimize as opt

import numpy as np
import matplotlib.pyplot as plt

from shapely import geometry
from stompy import utils
from stompy.spatial import wkb2shp, field
from stompy.grid import unstructured_grid, shadow_cdt

from stompy.spatial import linestring_utils
# from stompy.grid import front
from stompy.grid import paver
from stompy.grid import shadow_cdt, exact_delaunay

## 
bounds=wkb2shp.shp2geom('region-bounds-v00.shp')

xyz=np.loadtxt('apollo-xyz.txt') # run gen_scale to make this
apollo=field.PyApolloniusField(X=xyz[:,:2],F=xyz[:,2])

## 

six.moves.reload_module(exact_delaunay)
six.moves.reload_module(shadow_cdt)
#six.moves.reload_module(front)

##

# The actually grid generation step
g=unstructured_grid.UnstructuredGrid()

for geo in bounds['geom']:
    pnts=np.array(geo)

    A=pnts[:-1]
    B=pnts[1:]

    for a,b in zip(A,B):
        na=g.add_or_find_node(x=a)
        nb=g.add_or_find_node(x=b)
        try:
            j=g.add_edge(nodes=[na,nb])
        except g.GridException:
            pass

cycles=g.find_cycles(max_cycle_len=g.Nnodes())
polys=[geometry.Polygon( g.nodes['x'][cycle] )
       for cycle in cycles ]

from shapely import ops
full_poly=ops.cascaded_union( polys )


## 
# scale=field.ConstantField(100)
scale=apollo # out of order, from below

p=paver.Paving(geom=full_poly,density=scale)

p.verbose=1
p.pave_all()

output_dir='grid_v01'
os.path.exists(output_dir) or os.makedirs(output_dir)

p.write_suntans(output_dir)

##

g=unstructured_grid.SuntansGrid(output_dir)
g.write_ugrid(os.path.join(output_dir,'grid-v01.nc'))
