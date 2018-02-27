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

xyz=np.loadtxt('apollo-xyz.txt') # run gen_scale to make this
apollo=field.PyApolloniusField(X=xyz[:,:2],F=xyz[:,2])


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


af.loop(2695) 


## 
 
# grid_dir='grid_v00'
# os.path.exists(grid_dir) or os.mkdir(grid_dir)
# 
# g_safe=unstructured_grid.UnstructuredGrid(grid=g)
# g_safe.renumber()
# g_safe.delete_node_field('vh')
# g_safe.write_ugrid(os.path.join(grid_dir,'grid-v00.nc'),overwrite=True)
# 
# g_safe.write_edges_shp(os.path.join(grid_dir,'grid-v00.shp'))
# 
# ## 

## 

# with the check_edits() machinery, it's failing pretty fast.
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)
# zoom=(568019.3257816283, 568370.5310039715, 4152623.361854631, 4152885.0663912804)
af.grid.plot_edges(lw=0.5,color='k')

# af.choose_site().plot()
# af.grid.plot_nodes(clip=zoom,labeler=lambda i,r:str(i))
# af.grid.plot_edges(clip=zoom,labeler=lambda i,r:r['oring'])
ax.axis(zoom)

# up to 2694, looks good.
# 1147 is still just HINT
# resampled_success is True, though no changes.
# arguably 1147 should just be dropped.
# the Join is attempted and looks okay.
# seems like it doesn't have the same code as at home
# Following step is looking okay...
# wants to Resample.  That makes no changes, no improvement.
# possible that rolling that back is a problem?
# It's the join - it doesn't set oring for one of the
#  'extra' edges.

# when both are on a ring, it already requires that they
# be on the same ring, and that ring matches the edge's ring
# j_ac_ring

## 
self=af
cp=self.grid.checkpoint()

## 

site=self.choose_site()
resampled_success = self.resample_neighbors(site)
actions=site.actions()
metrics=[a.metric(site) for a in actions]
bests=np.argsort(metrics)
best=bests[1]
self.log.info("Chose strategy %s"%( actions[best] ) )
edits=actions[best].execute(site)

## 
opt_edits=self.optimize_edits(edits)
failures=self.check_edits(opt_edits)

## 
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)
# zoom=(570089.9221043529, 570258.4213675127, 4149369.4725253694, 4149495.031653724)
af.grid.plot_edges(lw=0.5,color='k')
ax.axis(zoom)

## 
from shapely import geometry
cc=af.grid.cells_center()

bad=[]
for ci in af.grid.valid_cell_iter():
    poly=af.grid.cell_polygon(ci)
    if not poly.contains( geometry.Point(cc[ci]) ):
        print "."
        bad.append(ci)

# 43 of em...


## 

# Debugging the cost function

cost_method='cc_py'
def cost_function(self,n):
    """
    Return a function which takes an x,y pair, and evaluates
    a geometric cost function for node n based on the shape and
    scale of triangle cells containing n
    """
    local_length = self.scale( self.grid.nodes['x'][n] )
    my_cells = self.grid.node_to_cells(n)

    if len(my_cells) == 0:
        return None

    cell_nodes = [self.grid.cell_to_nodes(c)
                  for c in my_cells ]

    # for the moment, can only deal with triangles
    cell_nodes=np.array(cell_nodes)

    # pack our neighbors from the cell list into an edge
    # list that respects the CCW condition that pnt must be on the
    # left of each segment
    for j in range(len(cell_nodes)):
        if cell_nodes[j,0] == n:
            cell_nodes[j,:2] = cell_nodes[j,1:]
        elif cell_nodes[j,1] == n:
            cell_nodes[j,1] = cell_nodes[j,0]
            cell_nodes[j,0] = cell_nodes[j,2] # otherwise, already set

    edges = cell_nodes[:,:2]
    edge_points = self.grid.nodes['x'][edges]

    def cost(x,edge_points=edge_points,local_length=local_length):
        return front.one_point_cost(x,edge_points,target_length=local_length)

    Alist=[ [ e[0],e[1] ]
            for e in edge_points[:,0,:] ]
    Blist=[ [ e[0],e[1] ]
            for e in edge_points[:,1,:] ]
    EPS=1e-5*local_length

    def cost_cc_and_scale_py(x0):
        C=list(x0)
        cc_cost=0
        scale_cost=0

        for A,B,cell_n in zip(Alist,Blist,cell_nodes):
            tri_cc=front.circumcenter_py(A,B,C)

            deltaAB=[ tri_cc[0] - A[0],
                      tri_cc[1] - A[1]]
            ABs=[B[0]-A[0],B[1]-A[1]]
            magABs=math.sqrt( ABs[0]*ABs[0] + ABs[1]*ABs[1])
            vecAB=[ABs[0]/magABs, ABs[1]/magABs]
            leftAB=vecAB[0]*deltaAB[1] - vecAB[1]*deltaAB[0] 

            deltaBC=[tri_cc[0] - B[0],
                     tri_cc[1] - B[1]]
            BCs=[C[0]-B[0], C[1]-B[1]]
            magBCs=math.sqrt( BCs[0]*BCs[0] + BCs[1]*BCs[1] )
            vecBC=[BCs[0]/magBCs, BCs[1]/magBCs]
            leftBC=vecBC[0]*deltaBC[1] - vecBC[1]*deltaBC[0]

            deltaCA=[tri_cc[0] - C[0],
                     tri_cc[1] - C[1]]
            CAs=[A[0]-C[0],A[1]-C[1]]
            magCAs=math.sqrt(CAs[0]*CAs[0] + CAs[1]*CAs[1])
            vecCA=[CAs[0]/magCAs, CAs[1]/magCAs]
            leftCA=vecCA[0]*deltaCA[1] - vecCA[1]*deltaCA[0]

            cc_fac=-4. # not bad
            # cc_fac=-2. # a little nicer shape
            # clip to 100, to avoid overflow in math.exp
            if 0:
                # this can favor isosceles too much
                this_cc_cost = ( math.exp(min(100,cc_fac*leftAB/local_length)) +
                                  math.exp(min(100,cc_fac*leftBC/local_length)) +
                                  math.exp(min(100,cc_fac*leftCA/local_length)) )
            else:
                # maybe?
                this_cc_cost = ( math.exp(min(100,cc_fac*leftAB/magABs)) +
                                 math.exp(min(100,cc_fac*leftBC/magBCs)) +
                                 math.exp(min(100,cc_fac*leftCA/magCAs)) )
                
            this_scale_cost=( (magABs-local_length)**2 
                              + (magBCs-local_length)**2 
                              + (magCAs-local_length)**2 )
            this_scale_cost/=local_length*local_length
            print "(%5d,%5d,%5d) => %8.4f cc  %8.4f scale"%(cell_n[0],cell_n[1],cell_n[2],
                                                            this_cc_cost,this_scale_cost)
            cc_cost+=this_cc_cost
            scale_cost+=this_scale_cost

        # With even weighting between these, some edges are pushed long rather than
        # having nice angles.
        # 3 is a shot in the dark.
        return 3*cc_cost+scale_cost

    if self.cost_method=='base':
        return cost
    elif self.cost_method=='cc_py':
        return cost_cc_and_scale_py
    else:
        assert False

print cost_function(af,1377)(af.grid.nodes['x'][1377])

## 

# on desktop:

# when trying to slide 1150, fails because it only 
# finds one neighbor.
# the edge 1150-1147 (j=6393) has 0 oring, even though 
# both nodes have oring=15.
# this at loop_count=2695

# at loop_count=1000:
# the nodes 1150---1147 are all ring=15, then 1146 ring=11.
# And all of those edges have oring=15.  
af.grid.plot_edges( clip=zoom, labeler=lambda i,r: r['oring'])

## 

# When does one 1147 edges lose its ring?
while 1:
    af.loop(1)
    print af.loop_count
##     
for j in af.grid.node_to_edges(1147):
    a,b=af.grid.edges['nodes'][j]
    print "   %5d -- %5d   oring=%3d"%( a,b, af.grid.edges['oring'][j] )

