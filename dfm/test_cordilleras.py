"""
Run each source for multiple durations.
"""
import six
import stein_dfm_script
import os
import shutil

from stompy.spatial import wkb2shp

six.moves.reload_module(stein_dfm_script)

# about 1.5GB per run.

sources=wkb2shp.shp2geom('forcing_with_q.shp')

# Set BCs from shapefile
from bcs import *

start_h=2
duration_h=10

flow=10.0
run_dir="runs/cordilleras_test"
if os.path.exists(run_dir):
    shutil.rmtree(run_dir)
    
stein_dfm_script.run_all(run_dir,start_h,duration_h,flow,sources=['redwood','steinberger','cordilleras'])
