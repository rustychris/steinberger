import sys
import xarray as xr
import six
from matplotlib import cm

import numpy as np
import matplotlib.pyplot as plt

from shapely import ops, geometry

from stompy import utils
from stompy.spatial import wkb2shp, field
from stompy.grid import unstructured_grid, shadow_cdt

from stompy.spatial import linestring_utils
from stompy.grid import front
from stompy.grid import shadow_cdt, exact_delaunay

## 

bounds=wkb2shp.shp2geom('region-bounds-v00.shp')

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

while 1:
    n_new = subdivide()
    print("Subdivide made %d new nodes"%n_new)
    if n_new == 0:
        break
## 


# Limit that to cdt cells for which the centroid is inside the domain

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

X=vc[select_cells]
# F=diams[select_cells]/2.0 # scale down to get 2-3 cells across
F=diams[select_cells] * 1.15 # aim for 1 cell across

xyz=np.c_[X,F]
np.savetxt('apollo-xyz.txt',xyz)
