import numpy as np
import os
import shutil
import six
import subprocess

import logging
log=logging.getLogger('stein_dfm')

import stompy.model.delft.io as dio
from stompy.spatial import wkb2shp
from stompy.io.local import noaa_coops
from stompy import filters
from stompy.model.delft import dfm_grid

## 

# Set BCs from shapefile
import bcs
six.moves.reload_module(bcs)
from bcs import *

##

dfm_bin_dir="/opt/software/delft/dfm/r53925-dbg/bin"
dfm_output_count=0
def dflowfm(mdu_fn,args=['--autostartstop'],nprocs=1):
    global dfm_output_count
    
    cmd=[os.path.join(dfm_bin_dir,"dflowfm")] + args
    if mdu_fn is not None:
        cmd.append(os.path.basename(mdu_fn))

    run_base_dir=os.path.dirname(mdu_fn)
    
    if nprocs>1:
        cmd=["%s/mpiexec"%dfm_bin_dir,"-n","%d"%nprocs] + cmd

    # This is more backwards compatible than 
    # passing cwd to subprocess()
    pwd=os.getcwd()
    dfm_output_count+=1
    log_file=os.path.join(run_base_dir,'dflowfm-log-%d.log'%dfm_output_count)
    with open(log_file, 'wt') as fp:
        log.info("Command '%s' logging to %s"%(cmd,log_file))
        try:
            os.chdir(run_base_dir)
            res=subprocess.call(cmd,stderr=subprocess.STDOUT,stdout=fp)
        finally:
            os.chdir(pwd)
            log.info("Command '%s' completed"%cmd)
    return res

def factory(feat):
    code='[' + feat['src'] + ']'
    return eval(code)

def run_all(run_base_dir,storm_start_h,storm_duration_h,storm_flow,sources=None,force=False):    
    mdu=dio.MDUFile('template.mdu')

    mdu['geometry','NetFile']='stein_03_net.nc'

    grid=dfm_grid.DFMGrid(mdu['geometry','NetFile'])
    
    if os.path.exists(run_base_dir):
        if force:
            shutil.rmtree(run_base_dir) # Safer - blow it away
        else:
            log.warning("Will not run %s -- already exists"%run_base_dir)
            return False

    mdu.set_time_range(start=np.datetime64('2010-01-01'),stop =np.datetime64('2010-01-05'))

    os.path.exists(run_base_dir) or os.makedirs(run_base_dir)
    mdu.set_filename(os.path.join(run_base_dir,'flowfm.mdu'))

    ext_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

    # Clear any pre-existing BC file:
    os.path.exists(ext_fn) and os.unlink(ext_fn)

    # Monkey patch the parameters:
    Storm.storm_flow=storm_flow
    Storm.storm_duration_h=storm_duration_h
    Storm.storm_start_h=storm_start_h

    bc_shp='forcing_with_q.shp'
    bc_shp_data=wkb2shp.shp2geom(bc_shp)

    for bc in bc_shp_data:
        if sources is not None and bc['name'] not in sources:
            print("Skipping %s"%bc['name'])
            continue
        for data_src in factory(bc):
            data_src.write(mdu=mdu,feature=bc,grid=grid)

    fixed_weir_out="../derived"
    if 1: # fixed weir file is just referenced as static input
        shutil.copyfile( os.path.join(fixed_weir_out,'fixed_weirs-v00.pli'),
                         os.path.join(run_base_dir,'fixed_weirs-v00.pli') )
        mdu['geometry','FixedWeirFile'] = 'fixed_weirs-v00.pli'

    mdu.write()

    dfm_grid.write_dfm(grid,mdu.filepath(['geometry','NetFile']),
                       overwrite=True)

    dflowfm(mdu.filename,['-t','1','--autostartstop'])

