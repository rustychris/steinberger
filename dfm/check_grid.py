# Getting an error that the network is not orthogonal??
import matplotlib.pyplot as plt
import xarray as xr
from stompy.model.delft import dfm_grid
from stompy.grid import front

from stompy.grid import unstructured_grid

## 
g=dfm_grid.DFMGrid('runs/test00/stein_00_net.nc')

## 

plt.clf()
g.plot_edges(lw=0.5,color='k')
g.plot_cells(mask=bad)

# Sure enough 34 cells are bad.

## 

six.moves.reload_module(unstructured_grid)
six.moves.reload_module(dfm_grid)
six.moves.reload_module(front)

## 
af=front.AdvancingTriangles(grid=g)

af.grid.nodes['fixed']=af.FREE

e2c=af.grid.edge_to_cells()

boundary_edge=np.any(e2c<0,axis=1)
boundary_nodes=np.unique( af.grid.edges['nodes'][ boundary_edge ] )

af.grid.nodes['fixed'][ boundary_nodes ] = af.RIGID

## 

cc=af.grid.cells_center()

bad=[]

for ci,xy in enumerate(cc):
    checked=af.grid.select_cells_nearest(xy,inside=True)
    if checked!=ci:
        bad.append(ci)

## 

bad_i=bad[2]

nodes=af.grid.cell_to_nodes(bad_i)
cpoly=af.grid.cell_polygon(bad_i)
xyxy=cpoly.bounds
L=xyxy[2] - xyxy[0]
zoom=[xyxy[0]-1.5*L,xyxy[2]+1.5*L,
      xyxy[1]-1.5*L,xyxy[3]+1.5*L]


plt.clf()
fig,axs=plt.subplots(2,1,sharex=True,sharey=True,num=1)

af.grid.plot_edges(lw=0.5,color='k',ax=axs[0])
af.grid.plot_cells(mask=[bad_i],ax=axs[0])
axs[0].axis(zoom)

## 

from stompy.spatial import field
af.scale=field.ConstantField(0.4*L)

af.optimize_nodes(nodes,max_levels=4)

axs[1].collections=[]
af.grid.plot_edges(lw=0.5,color='k',ax=axs[1])
af.grid.plot_cells(mask=[bad_i],ax=axs[1])
axs[1].axis(zoom)

## 

from collections import defaultdict
from shapely import geometry


print check_changes(af,dict(cells=[bad[3]]))

    
