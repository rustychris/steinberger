import six
import stein_dfm_script

six.moves.reload_module(stein_dfm_script)

# about 1.5GB per run.

for flow in [2,10,25]:
    for start_h in [42,45,48]:
        for duration_h in [3,6]:
            name="runs/S%.1f_D%.1f_Q%04.1f"%(start_h,duration_h,flow)
            
            stein_dfm_script.run_all(name,start_h,duration_h,flow)
