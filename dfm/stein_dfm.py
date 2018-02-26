import numpy as np
import os
import six

import stompy.model.delft.io as dio
from stompy.spatial import wkb2shp
from stompy.io.local import noaa_coops
from stompy import filters
from stompy.model.delft import dfm_grid


## 

mdu=dio.MDUFile('template.mdu')

mdu['geometry','NetFile']='stein_00_net.nc'

grid=dfm_grid.DFMGrid('stein_00_net.nc')

run_base_dir='runs/test00'
mdu.set_time_range(start=np.datetime64('2010-01-01'),stop =np.datetime64('2010-01-05'))

os.path.exists(run_base_dir) or os.makedirs(run_base_dir)
mdu.set_filename(os.path.join(run_base_dir,'flowfm.mdu'))


# Set BCs from shapefile
import bcs
six.moves.reload_module(bcs)
from bcs import *

def factory(feat):
    return eval(feat['src'])

bc_shp='forcing.shp'
bcs=wkb2shp.shp2geom(bc_shp)

for bc in bcs:
    data_src=factory(bc)
    data_src.write(mdu,bc)

mdu.write()

mdu['geometry','NetFile']='stein_00_net.nc'

dfm_grid.write_dfm(grid,mdu.filepath(['geometry','NetFile']),
                   overwrite=True)
