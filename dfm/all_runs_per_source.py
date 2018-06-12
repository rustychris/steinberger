"""
Run each source for multiple durations.
"""
import six
import stein_dfm_script
from stompy.spatial import wkb2shp

six.moves.reload_module(stein_dfm_script)

# about 1.5GB per run.

sources=wkb2shp.shp2geom('forcing_with_q.shp')

# Set BCs from shapefile
from bcs import *

for start_h in [42,45,48]:
    for duration_h in [3,6,7.3]:
        for source in sources:
            name=source['name']
            
            objs=stein_dfm_script.factory(source)
            flow=None
            for obj in objs:
                if isinstance(obj,Storm):
                    flow=obj.storm_flow
                    break
            else:
                # This source is not a stormwater source
                continue
            
            run_dir="runs/S%.1f_D%.1f_Q%.2f_%s"%(start_h,duration_h,flow,name)
            stein_dfm_script.run_all(run_dir,start_h,duration_h,flow,sources=['redwood','steinberger',name])
