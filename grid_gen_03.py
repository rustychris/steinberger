"""
Use quads from Janet, and generate compatible triangles
Handle this with more tried-and-true tom code
"""
import os
import xarray as xr
import six
from matplotlib import cm
import scipy.optimize as opt

import numpy as np
import matplotlib.pyplot as plt

from stompy.model.delft import dfm_grid
from shapely import geometry
from stompy import utils
from stompy.spatial import wkb2shp, field
from stompy.grid import unstructured_grid, shadow_cdt
from stompy.grid import front
from stompy.plot import plot_wkb

from stompy.spatial import linestring_utils
from stompy.grid import paver
from stompy.grid import shadow_cdt, exact_delaunay

## 

six.moves.reload_module(exact_delaunay)
six.moves.reload_module(shadow_cdt)
six.moves.reload_module(front)

## 
g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/janet-out-edit04.nc')

g_janet.make_cells_from_edges(max_sides=4)
g_janet.edge_to_cells(recalc=True)

##

plt.figure(1).clf()
g_janet.plot_edges(color='k',lw=0.5)
# g_janet.plot_nodes(mask=[4299,4300],labeler='id')
g_janet.plot_nodes(clip=plt.axis(),labeler='id')

plt.axis('equal')

##

# HERE - some issue with doubling up edges at 4299/4300 ?
# yep, lots of edges duplicated.


# Try handing this to front
internal_pnt=[568572, 4153996]
node_string = g_janet.enclosing_nodestring(internal_pnt,max_nodes=-1)

has_cell_nbrs=[len(g_janet.node_to_cells(n))>0
               for n in node_string]

##

scale=field.PyApolloniusField(X=g_janet.nodes['x'],
                              F=50*np.ones(g_janet.Nnodes()))

## 

af=front.AdvancingTriangles(grid=g_janet,scale=scale)
af.cost_method='base' # trying that again..

af.add_curve(nodes=node_string)
af.initialize_boundaries()

##

af.loop(10)

## 

plt.cla()
af.grid.plot_edges(lw=0.5,
                   values=np.log10(1+np.arange(af.grid.Nedges())[::-1]),
                   cmap='jet')
plt.axis('equal')
