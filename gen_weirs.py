# Borrowed from lsb_dfm

# v04: Add code to process water control structures.
# v03: use linear_features-v02.shp, to bring in some manual features,
#      mostly these are the features we want, although gate-like features
#      are also in here.
# v02: pull max elevation within 5m -- at some places the levee line does not
# quite follow the crest of the levee.  
# v01: separate out by 'Class'

import os
import matplotlib.pyplot as plt
from matplotlib import cm

import numpy as np

from shapely import geometry
from scipy import ndimage

from stompy.spatial import (field,wkb2shp,linestring_utils)
from stompy import utils
import stompy.plot.cmap as scm
from stompy.grid import unstructured_grid
import stompy.model.delft.io as dio

data_root='.'
out_dir='derived'

## 

os.path.exists(out_dir) or os.makedirs(out_dir)
## 


dem=field.GdalGrid(os.path.join(data_root,"master_dem_v00.tif"))

## 

inv_fn="region-bounds-v00.shp"
inv=wkb2shp.shp2geom(inv_fn)

## 

# For the moment, focus on just the parts of the inventory which overlap with the DEM.
bounds_xyxy=np.array(dem.extents)[ [0,2,1,3] ]
box=geometry.box(*bounds_xyxy)

# For steinberger, use all of the features in region-bounds.
sel=np.ones(len(inv),'b1')

## 

# New in v02:
# Apply a max filter across the DEM.

# for a 5m radius and 2m pixels, not much hope in really resolving
# a disc, but here goes
footprint=np.array( [[0,1,1,1,0],
                     [1,1,1,1,1],
                     [1,1,1,1,1],
                     [1,1,1,1,1],
                     [0,1,1,1,0]] )

## 

# The real deal - update dem in place (in RAM, not on disk)
dem.F=ndimage.maximum_filter(dem.F,footprint=footprint)

## 
res=field.ConstantField(5.0) # target output linear resolution

# count total features so that files can be concatenated without issue.
total_count=0

# no longer split by class, all filtering already complete in the input
# shapefile

g=unstructured_grid.UnstructuredGrid(extra_node_fields=[ ('elev_m','f8')],
                                     extra_edge_fields=[ ('mean_elev_m','f8')] )

for ix,sel_i in enumerate(np.nonzero(sel)[0]):
    geom=inv['geom'][sel_i]
    coords=np.array(geom)
    # not 100% sure about whether it's necessary to test for closed_ring
    new_coords=linestring_utils.upsample_linearring(coords,res,closed_ring=0)

    nodes=[g.add_or_find_node(x=xy,tolerance=0.0,
                              elev_m=np.nan)
           for xy in new_coords]

    for a,b in zip(nodes[:-1],nodes[1:]):
        j=g.nodes_to_edge([a,b])
        if j is None: # go easy in case there are overlaps.
            j=g.add_edge(nodes=[a,b])

# pull out point elevations at the nodes:
g.nodes['elev_m'] = dem( g.nodes['x'] )

# drastic, but go ahead and delete any nodes which failed to get an elevation

missing=np.isnan( g.nodes['elev_m'] )

for n in np.nonzero(missing)[0]:
    g.delete_node_cascade(n)

## 

# to replicate the 5 fields, where last two are just 10.0
# not entirely sure what these /should/ be, but this is what
# I've seen in previous input files.
g.add_node_field('sill_left',10*np.ones_like(g.nodes['elev_m']))
g.add_node_field('sill_right',10*np.ones_like(g.nodes['elev_m']))

## 

pli_data=dio.grid_to_pli_data(g,node_fields=['elev_m','sill_left','sill_right'],
                             labeler=lambda i: "L%04d"%(total_count+i))
dio.write_pli(os.path.join(out_dir,'fixed_weirs-v00.pli'),pli_data)

## 


def write_node_shp(self,shpname,extra_fields=[]):
    """ Write a shapefile with each node.  Fields will attempt to mirror
    self.nodes.dtype

    extra_fields: goal is similar to write_cells_shp and write_edges_shp, 
    but not yet supported.
    """
    assert len(extra_fields)==0 # not yet supported!
    
    # zero-based index of node (why does write_edge_shp create 1-based ids?)
    base_dtype = [('node_id',np.int32)]

    node_geoms=[geometry.Point( self.nodes['x'][i] )
                for i in self.valid_node_iter() ]

    node_data=self.nodes[~self.nodes['deleted']].copy()

    # don't need to write all of the original fields out:
    node_data=utils.recarray_del_fields(node_data,['x','deleted'])
    

    wkb2shp.wkb2shp(shpname,input_wkbs=node_geoms,fields=node_data,
                    overwrite=True)


nodes_fn=os.path.join(out_dir,inv_fn.replace(".shp","-nodes.shp"))
write_node_shp(g,nodes_fn)


