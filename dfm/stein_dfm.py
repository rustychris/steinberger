import numpy as np
import os
import shutil
import six

import stompy.model.delft.io as dio
from stompy.spatial import wkb2shp
from stompy.io.local import noaa_coops
from stompy import filters
from stompy.model.delft import dfm_grid

## 

mdu=dio.MDUFile('template.mdu')

mdu['geometry','NetFile']='stein_03_net.nc'

grid=dfm_grid.DFMGrid(mdu['geometry','NetFile'])

run_base_dir='runs/test01'
if os.path.exists(run_base_dir):
    shutil.rmtree(run_base_dir) # Safer - blow it away

mdu.set_time_range(start=np.datetime64('2010-01-01'),stop =np.datetime64('2010-01-05'))

os.path.exists(run_base_dir) or os.makedirs(run_base_dir)
mdu.set_filename(os.path.join(run_base_dir,'flowfm.mdu'))

## -

ext_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

# Clear any pre-existing BC file:
os.path.exists(ext_fn) and os.unlink(ext_fn)

# Set BCs from shapefile
import bcs
six.moves.reload_module(bcs)
from bcs import *

def factory(feat):
    return eval(feat['src'])

bc_shp='forcing.shp'
bc_shp_data=wkb2shp.shp2geom(bc_shp)

for bc in bc_shp_data:
    data_src=factory(bc)
    data_src.write(mdu=mdu,feature=bc,grid=grid)

###

fixed_weir_out="../derived"
if 1: # fixed weir file is just referenced as static input
    shutil.copyfile( os.path.join(fixed_weir_out,'fixed_weirs-v00.pli'),
                     os.path.join(run_base_dir,'fixed_weirs-v00.pli') )
    mdu['geometry','FixedWeirFile'] = 'fixed_weirs-v00.pli'


## 
    
mdu.write()

dfm_grid.write_dfm(grid,mdu.filepath(['geometry','NetFile']),
                   overwrite=True)

## 

