"""
Use quads from Janet, and generate compatible triangles
"""
from __future__ import print_function
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

# g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/janet-out-edit08.nc')
# g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/janet-out-edit10.nc')
# g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/janet-out-edit12.nc')
# g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/janet-out-stompy13.nc')
# g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/janet-out-edit16.nc')
g_janet=unstructured_grid.UnstructuredGrid.from_ugrid('janet/edit_26.nc')

g_janet.make_cells_from_edges(max_sides=4)
g_janet.edge_to_cells(recalc=True)


# Try handing this to front
internal_pnts=[
    [568572, 4153996], # largest area
    [568293, 4152059], # second largest
    [567487, 4151055], # southwestern chunk
    [572050, 4151616], # eastern shoal
    [568704, 4150650],
    [568808, 4150507],
    [568769, 4150312],
    [569065, 4150143],
    [571104, 4152122],
    [571454, 4152500],
    [571926, 4151130],
    [569398, 4151137],
    [572595, 4149360],
    [567825, 4154910],
]

# Point [568704, 4150722] returned no node string

# Had been 50, but maybe 65 will get us up to scale faster?
#scale=field.PyApolloniusField(X=g_janet.nodes['x'],
#                              F=65*np.ones(g_janet.Nnodes()))

af=front.AdvancingTriangles(grid=g_janet) # ,scale=scale)
af.cost_method='base' # trying that again..

# superceded by all_missing below


# furthermore, only the inward facing edges in the selected cycle should be flagged
# UNMESHED and subject to gridding.  Set all missing cells to UNDEFINED, then come
# back to set just the cycle we want.
all_missing=af.grid.edges['cells']<0
af.grid.edges['cells'][all_missing]=af.grid.UNDEFINED

if 1:
    # Select enclosing polygons for specific points:
    for internal_pnt in internal_pnts:
        node_string = g_janet.enclosing_nodestring(internal_pnt,max_nodes=-1)
        if node_string is None:
            print("Point %s returned no node string"%str(internal_pnt))
        af.add_curve(nodes=node_string)

        # Rather than call initialize boundaries, we want to just fix the existing
        # edges
        for a,b in utils.circular_pairs(node_string):
            he=af.grid.nodes_to_halfedge(a,b)
            af.grid.edges['cells'][he.j,he.orient]=af.grid.UNMESHED

        # Node rigid status:
        for n in node_string:
            if len(g_janet.node_to_cells(n))>0:
                af.grid.nodes['fixed'][n]=af.RIGID
else:
    # Set outside boundary to UNDEFINED, and all internal, unmeshed
    # areas to UNMESHED
    # Need to set the very outside to UNDEFINED, rather than UNMESHED.
    boundary=af.grid.boundary_cycle()
    for a,b in utils.circular_pairs(boundary):
        he=af.grid.nodes_to_halfedge(a,b)
        af.grid.edges['cells'][he.j,1-he.orient]=af.grid.UNDEFINED

    for n in node_string:
        if len(g_janet.node_to_cells(n))>0:
            af.grid.nodes['fixed'][n]=af.RIGID
        
    
# Edge rigid status
rigid_edges=g_janet.edges['cells'].max(axis=1)>=0
af.grid.edges['fixed'][rigid_edges]=af.RIGID

# Update scale to be edge lengths for existing boundary edges
scale_edges=(af.grid.edges['cells'].min(axis=1)==af.grid.UNMESHED)&(af.grid.edges['cells'].max(axis=1)>=0)

scale=field.PyApolloniusField(X=g_janet.edges_center()[scale_edges],
                              F=g_janet.edges_length()[scale_edges])
af.set_edge_scale(scale)

##


plt.figure(1).clf()
g_janet.plot_edges(color='k',lw=0.5)
P=np.array(internal_pnts)
plt.plot(P[:,0],P[:,1],'go')
plt.axis('equal')

##

# This often has scale issues at the end
# af.loop(1000)

## 

# try invoking the search tree approach
import logging
af.log.setLevel(logging.INFO)

af.current=af.root=front.DTChooseSite(af)

##

# can this be rewritten so that af.current holds enough state to always
# know the next thing to do?
def advance_greedy(af):
    """
    based on the state of af.current, take the next step.
    returns True indicating more work to do, or False to 
    indicate completion or stuck
    """
    # active_child refers to the last child of this node which was
    # attempted.  it is None/-1 if no children have been attempted.
    # it is len(children) if all children have been attempted and
    # failed (?)
    if af.current.active_child is None:
        af.current.active_child=-1
    af.current.active_child+=1

    if len(af.current.children)==0:
        return False
    
    if af.current.active_child>=len(af.current.children):
        print("Children at this node have been exhausted")
        af.current.revert_to_parent()
        return True
    else:
        # on success, af.current will have been updated
        success=af.current.try_child(af.current.active_child)
        print("Child success? ",success)
        return True

##

while 1:
    if not advance_greedy(af):
        break
##

# Something is wrong with the topology around [ 569242.9437, 4151191.7679]
# That is a point being fed to point_to_f, but it looks like it doesn't
# sit on a curve
# That's node 1129, which sure enough claims oring=0
# node_ring_f doesn't require that the n fall on the requested ring.
# could it require that the node fall on *a* ring?
# The request is for ring=7.
# this is during find_slide_conflicts(n=3978, delta_f=2285.397233532136)
# that delta_f appears to be a full run around the entire curve -- just shy of:
# 2285.502557288453
# Okay - one problem is that with such a large requested delta, *all* nodes
# will appear to be within the window, which means we can't even tell forward
# from backwards.
# So a first fix is to see why relax_slide_node() is choosing the long way
# around
# including why is find_slide_limits returning a range, [2285.502557288453, 2657.4545534723484]
# which doesn't include the node, f0=288.7246430585429?

## 

plt.cla()
af.grid.plot_edges(lw=0.25,
                   values=np.log10(1+np.arange(af.grid.Nedges())[::-1]),
                   cmap='jet')
plt.axis('equal')

##
af.grid.renumber()
af.grid.delete_node_field('vh')

af.grid.write_ugrid('janet/stompy27.nc',overwrite=True)


##

# stein_dfm.py is getting some nan cell centers.
# also, lots of messages that circumcenter is not inside cell.

g=unstructured_grid.UnstructuredGrid.from_ugrid('janet/stompy29.nc')

##
from shapely import geometry

# this is all finite, no nodes are nan.
cc=g.cells_center()

bad=[]
for c in range(g.Ncells()):
    poly=g.cell_polygon(c)
    pnt=geometry.Point(cc[c])
    if not poly.contains(pnt):
        bad.append(c)

##


plt.clf()
g.plot_edges(lw=0.2,color='k')
g.plot_cells(mask=bad)
plt.axis('equal')
