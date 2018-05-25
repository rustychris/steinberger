import os

from stompy.grid import unstructured_grid
from stompy.grid import depth_connectivity
from stompy.spatial import field
from stompy.model.delft import dfm_grid

## 

dem=field.GdalGrid('master_dem_v00.tif')

g=dfm_grid.DFMGrid('janet/janet-out_net.nc')

## 

edge_depths=depth_connectivity.edge_connection_depth(g,dem,edge_mask=None,centers='lowest')
basic_node_depths=dem( g.nodes['x'] )

node_depths=depth_connectivity.greedy_edgemin_to_node(g,basic_node_depths,edge_depths)

g.add_node_field('depth',node_depths)


grid_out_fn='dfm/stein_02_net.nc'

if os.path.exists(grid_out_fn):
    os.unlink(grid_out_fn)
    
dfm_grid.write_dfm(g,grid_out_fn)
