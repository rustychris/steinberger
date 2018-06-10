import matplotlib.pyplot as plt
import xarray as xr
from stompy.model.delft import dfm_grid
from stompy.grid import front

from stompy.grid import unstructured_grid
##

six.moves.reload_module(unstructured_grid)
six.moves.reload_module(dfm_grid)
six.moves.reload_module(front)

## 
g=dfm_grid.DFMGrid('runs/test02/stein_03_net.nc')

## 

plt.clf()
ax=plt.gca()

g.plot_edges(lw=0.5,color='k',ax=ax)

scat=ax.scatter(g.nodes['x'][:,0],
                g.nodes['x'][:,1],
                30,g.nodes['depth'])

ax.axis('equal')

## 

cc=g.cells_center()

bad=[]

for ci,xy in enumerate(cc):
    checked=g.select_cells_nearest(xy,inside=True)
    if checked!=ci:
        bad.append(ci)

## 

# Nothing obvious bad with the horizontal status.

