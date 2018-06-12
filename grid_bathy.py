import os

from stompy.grid import unstructured_grid
from stompy.grid import depth_connectivity
from stompy.spatial import field
from stompy.model.delft import dfm_grid

## 

dem=field.GdalGrid('master_dem_v00.tif')

g=unstructured_grid.UnstructuredGrid.from_ugrid('janet/stompy29.nc')

## 

edge_depths=depth_connectivity.edge_connection_depth(g,dem,edge_mask=None,centers='lowest')

basic_node_depths=dem( g.nodes['x'] )

node_depths=depth_connectivity.greedy_edgemin_to_node(g,basic_node_depths,edge_depths)

## 


# def greedy_edgemean_to_node(g,orig_node_depth,edge_min_depth):
#     """
#     Try to adjust node depths such that the average of the nodes per edge
#     is at least as deep as edge_min_depth.  This errs on the side of
#     making edges deeper, useful for barely resolved conveyances.
#     """

orig_node_depth=g.nodes['depth'].copy()

target_edge_depth=edge_depths

conn_depth=orig_node_depth.copy()

node_mean=conn_depth[g.edges['nodes']].mean(axis=1)
errors=node_mean - target_edge_depth
errors[ np.isnan(errors) ] = 0.0


##
from scipy.optimize import fmin

potential=np.ones(g.Nedges())

## 

for loop in range(1000):
    verbose= (loop%100==0)
        
    # Find an offender
    j_bad=np.argmax(potential*errors)
    if potential[j_bad]==0:
        print "DONE"
        break
    
    potential[j_bad]=0 # only visit each edge once.

    
    # Get the neighborhood of nodes:
    # nodes=
    jj_nbrs=np.concatenate( [ g.node_to_edges(n)
                              for n in g.edges['nodes'][j_bad] ] )
    jj_nbrs=np.unique(jj_nbrs)
    jj_nbrs = jj_nbrs[ np.isfinite(target_edge_depth[jj_nbrs]) ]

    n_bad=g.edges['nodes'][j_bad]
    
    def cost(ds):
        # Cost function over the two depths of the ends of j_bad:
        conn_depth[n_bad]=ds
        new_errors=conn_depth[g.edges['nodes'][jj_nbrs]].mean(axis=1) - target_edge_depth[jj_nbrs]
        # weight high edges 10x more than low edges:
        cost=new_errors.clip(0,np.inf).sum() - 0.5 * new_errors.clip(-np.inf,0).sum()
        return cost
    ds0=conn_depth[n_bad]
    cost0=cost(ds0)
    
    ds=fmin(cost,ds0,disp=False)
    costn=cost(ds)
    conn_depth[n_bad]=ds

    if verbose:
        print("Loop %d: %d/%d edges  starting error: j=%d => %.4f"%(loop,potential.sum(),len(potential),
                                                                    j_bad,errors[j_bad]))
    
    node_mean=conn_depth[g.edges['nodes']].mean(axis=1)
    errors=node_mean - target_edge_depth
    errors[ np.isnan(errors) ] = 0.0

    if verbose:
        print("    ending error: j=%d => %.4f"%(j_bad,errors[j_bad]))

## 


node_mean=conn_depth[g.edges['nodes']].mean(axis=1)
errors=node_mean - target_edge_depth
errors[ np.isnan(errors) ] = 0.0

# plot
plt.figure(10).clf()
ncoll=g.plot_nodes(values=g.nodes['depth'])
ncoll.set_clim([-1,2])
ncoll.set_cmap('jet')

# ecoll=g.plot_edges(values=node_mean)
# ecoll.set_clim([-1,2])
# ecoll.set_cmap('jet')

ecoll=g.plot_edges(values=errors)
ecoll.set_clim([-1,1])
ecoll.set_cmap('seismic')
# g.plot_edges(lw=0.4,color='k')


# plt.axis( (568259., 568614., 4150269., 4150641.) )
plt.axis( (569564.3238793822, 569761.5851150813, 4151337.141236072, 4151543.8487844667) )


##

g.add_node_field('depth',conn_depth,on_exists='overwrite')

# g.add_edge_field('depth_mean', g.nodes['depth'][g.edges['nodes']].mean(axis=1),on_exists='overwrite' )

##

grid_out_fn='dfm/stein_04_net.nc'

if os.path.exists(grid_out_fn):
    os.unlink(grid_out_fn)
    
dfm_grid.write_dfm(g,grid_out_fn)

