import xarray as xr
import six
from matplotlib import cm

import numpy as np
import matplotlib.pyplot as plt

from shapely import geometry
from stompy import utils
from stompy.spatial import wkb2shp, field
from stompy.grid import unstructured_grid, shadow_cdt

from stompy.spatial import linestring_utils
from stompy.grid import front
from stompy.grid import shadow_cdt, exact_delaunay

## 
bounds=wkb2shp.shp2geom('region-bounds-v00.shp')

## 

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

## 

# the base class for ShadowCDT provides more of the
# stuff we want.
cdt=shadow_cdt.ShadowCDT(g)

## 
plt.clf()
cdt.plot_edges(lw=0.5,color='k')

## 
import sys
# the new way - calculated voronoi points directly from the triangles
# in the delaunay triangulation, then match with edges with a hash
# on edge [a,b] node pairs
def subdivide():
    # what does cdt need to provide for this to work?
    vcenters = cdt.cells_center(refresh=True)

    n_edges = cdt.Nedges()
    to_subdivide = []

    min_edge_length=10.0

    for j_g in g.valid_edge_iter(): 
        a,b = g.edges['nodes'][j_g] 
        j=cdt.nodes_to_edge([a,b])
        cells= cdt.edges['cells'][j]

        assert cells.max()>=0

        for ci,c in enumerate(cells):
            if c<0:
                continue

            # Need the signed distance here:
            pntV = vcenters[c]
            pntA = cdt.nodes['x'][a]
            pntB = cdt.nodes['x'][b]
            AB=pntB-pntA
            AV=pntV-pntA
            # just switch sign on this
            left=utils.to_unit(np.array([-AB[1],AB[0]]))
            if ci==1:
                left*=-1
            line_clearance=np.dot(left,AV)
            v_radius = utils.mag(AV)

            if utils.mag(AB)<min_edge_length:
                continue
            # line_clearance=utils.point_line_distance(vcenters[c],
            #                                         cdt.nodes['x'][ [a,b] ] )
            # v_radius=utils.dist( vcenters[c], cdt.nodes['x'][a] )

            if (v_radius > 1.2*line_clearance) and (v_radius > min_edge_length):
                # second check - make sure that neither AC nor BC are also on the
                # boundary
                c_j=cdt.cell_to_edges(c)
                count=cdt.edges['constrained'][c_j].sum()

                if count == 1:
                    to_subdivide.append(j_g)
                    break
                elif count == 0:
                    print("While looking at edge %d=(%d,%d)"%(j_g,a,b))
                    raise Exception("We should have found at least 1 boundary edge")
                elif count == 3:
                    print("WARNING: Unexpected count of boundary edges in one element: ",count)

    for j in to_subdivide:
        sys.stdout.write(str(j)+",")
        sys.stdout.flush()
        g.split_edge(j)

    return len(to_subdivide)
## 
while 1:
    n_new = subdivide()
    print("Subdivide made %d new nodes"%n_new)
    if n_new == 0:
        break
## 


# Limit that to cdt cells for which the centroid is inside the domain
from shapely import ops, geometry

bounds=wkb2shp.shp2geom('region-bounds-v00.shp')

domain_poly=ops.cascaded_union(polys)


select_cells=[domain_poly.contains( geometry.Point(cxy) )
              for cxy in cdt.cells_centroid() ]
select_cells=np.array(select_cells)

## 

# come up with scales:
vc=cdt.cells_center(refresh=True)
diams=np.zeros(cdt.Ncells(),np.float64)
diams[:]=0

for c in cdt.valid_cell_iter():
    n=cdt.cells['nodes'][c,0]
    radius=utils.dist( vc[c], cdt.nodes['x'][n] )
    diams[c]=2*radius

plt.clf()
ccoll=cdt.plot_cells(values=diams,mask=select_cells)
ccoll.set_clim([1,500])

g.plot_edges(color='k')

## 
X=vc[select_cells]
# F=diams[select_cells]/2.0 # scale down to get 2-3 cells across
F=diams[select_cells] * 1.15 # aim for 1 cell across

apollo=field.PyApolloniusField(X=X,F=F)


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

# # 
cycles=g.find_cycles(max_cycle_len=g.Nnodes())

polys=[geometry.Polygon( g.nodes['x'][cycle] )
       for cycle in cycles ]

# # 

six.moves.reload_module(exact_delaunay)
six.moves.reload_module(shadow_cdt)
six.moves.reload_module(front)

# scale=field.ConstantField(100)
scale=apollo # out of order, from below
af=front.AdvancingTriangles(grid=g,scale=scale)
af.initialize_boundaries()

# May need to doctor up some things

# Nodes of degree>2 are fixed:
for n in g.valid_node_iter():
    if g.node_degree(n)>2:
        g.nodes['fixed'][n]=af.RIGID

# # 
g.nodes['oring']=0
Curve=front.Curve

node_rigid=af.grid.nodes['fixed']==af.RIGID

# All of these nodes go on Curves
for cycle in cycles:
    # Break cycles at fixed nodes:
    rigids=np.nonzero(node_rigid[cycle])[0]

    if len(rigids):
        # start on a rigid node
        cycle=np.roll(cycle,-rigids[0])
        rigids=np.nonzero(node_rigid[cycle])[0]

        for a,b in utils.circular_pairs(rigids):
            if a<b:
                seg_nodes=cycle[a:b+1] # include b
            else:
                seg_nodes=np.concatenate( ( cycle[a:],cycle[:b+1]) )

            # only non-end points of the curve
            seg_oring=af.grid.nodes['oring'][seg_nodes[1:-1]]
            if (len(seg_oring) < 1) or (seg_oring[0]==0):
                curve_idx=af.add_curve(nodes=seg_nodes,closed=False)
    else:
        print "Got a loop" # never been tested
        curve_idx=af.add_curve(nodes=cycle,closed=True)

# Before we can do any paving, have to mark edges as UNMESHED and UNDEFINED
g.edges['cells']=g.UNMESHED
# And then set any outward facing edges to UNDEFINED.
outside=g.boundary_cycle()

for a,b in utils.circular_pairs(outside):
    j=g.nodes_to_edge([a,b])
    side= int(g.edges['nodes'][j,0]==a)
    g.edges['cells'][j,side]=g.UNDEFINED


af.loop() 

## 

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

g.plot_edges(ax=ax,lw=0.5)

for curve in af.curves:
    curve.plot(ax=ax,alpha=0.5,zorder=-1,lw=2)
# ax.axis(zoom)

## 

grid_dir='grid_v00'
os.path.exists(grid_dir) or os.mkdir(grid_dir)

g_safe=unstructured_grid.UnstructuredGrid(grid=g)
g_safe.renumber()
g_safe.delete_node_field('vh')
g_safe.write_ugrid(os.path.join(grid_dir,'grid-v00.nc'),overwrite=True)

g_safe.write_edges_shp(os.path.join(grid_dir,'grid-v00.shp'))
