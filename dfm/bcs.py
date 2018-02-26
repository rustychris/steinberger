import os

import numpy as np
import stompy.model.delft.io as dio
from stompy.io.local import noaa_coops
from stompy import utils, filters

cache_dir='cache'
os.path.exists(cache_dir) or os.mkdir(cache_dir)

class BC(object):
    pass

class NoaaTides(BC):
    var_names=['ssh']
    def __init__(self,station,datum='NAVD88',z_offset=0.0):
        self.station=station
        self.datum=datum
        self.z_offset=z_offset
    def write(self,mdu,feat):
        print "Feature: %s"%(feat['name'])

        name=feat['name']
        old_bc_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

        for var_name in self.var_names:
            if feat['geom'].type=='LineString':
                pli_data=[ (name, np.array(feat['geom'].coords)) ]
                base_fn=os.path.join(mdu.base_path,"%s_%s"%(name,var_name))
                pli_fn=base_fn+'.pli'
                dio.write_pli(pli_fn,pli_data)

                if var_name=='ssh':
                    quant='waterlevelbnd'
                else:
                    assert False

                with open(old_bc_fn,'at') as fp:
                    lines=["QUANTITY=%s"%quant,
                           "FILENAME=%s%s.pli"%(name,var_name),
                           "FILETYPE=9",
                           "METHOD=3",
                           "OPERAND=O",
                           ""]
                    fp.write("\n".join(lines))

                self.write_data(mdu,feat,var_name,base_fn)
            else:
                assert False

    def write_data(self,mdu,feat,var_name,base_fn):
        tides=noaa_coops.coops_dataset_product(self.station,'water_level',
                                               mdu.time_range()[1],mdu.time_range()[2],
                                               days_per_request='M',cache_dir=cache_dir)
        tide=tides.isel(station=0)
        water_level=utils.fill_tidal_data(tide.water_level) + self.z_offset
        # IIR butterworth.  Nicer than FIR, with minor artifacts at ends
        # 3 hours, defaults to 4th order.
        water_level[:] = filters.lowpass(water_level[:].values,
                                         utils.to_dnum(water_level.time),
                                         cutoff=3./24)

        ref_date=mdu.time_range()[0]
        elapsed_minutes=(tide.time.values - ref_date)/np.timedelta64(60,'s')

        # just write a single node
        tim_fn=base_fn + "_0001.tim"
        data=np.c_[elapsed_minutes,water_level]
        np.savetxt(tim_fn,data)


class Storm(BC):
    var_names=['q']
    def __init__(self,name):
        self.name=name
    def write(self,mdu,feat):
        # obvious copy and paste from above.
        # not quite ready to abstract, though
        print "Feature: %s"%(feat['name'])

        name=feat['name']
        old_bc_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

        for var_name in self.var_names:
            if feat['geom'].type=='LineString':
                pli_data=[ (name, np.array(feat['geom'].coords)) ]
                base_fn=os.path.join(mdu.base_path,"%s_%s"%(name,var_name))
                pli_fn=base_fn+'.pli'
                dio.write_pli(pli_fn,pli_data)

                if var_name=='q':
                    quant='dischargebnd'
                else:
                    assert False

                with open(old_bc_fn,'at') as fp:
                    lines=["QUANTITY=%s"%quant,
                           "FILENAME=%s%s.pli"%(name,var_name),
                           "FILETYPE=9",
                           "METHOD=3",
                           "OPERAND=O",
                           ""]
                    fp.write("\n".join(lines))

                self.write_data(mdu,feat,var_name,base_fn)
            else:
                assert False

    def write_data(self,mdu,feat,var_name,base_fn):
        ref_date,run_start,run_stop=mdu.time_range()

        # trapezoid hydrograph
        times=np.array( [run_start,
                         run_start+np.timedelta64(47,'h'),
                         run_start+np.timedelta64(48,'h'),
                         run_start+np.timedelta64(48+3,'h'),
                         run_start+np.timedelta64(48+4,'h'),
                         run_stop+np.timedelta64(1,'D')] )
        flows=np.array( [0.0,0.0,
                         10.0,10.0,0.0,0.0] )
        elapsed_minutes=(times - ref_date)/np.timedelta64(60,'s')

        # just write a single node
        tim_fn=base_fn + "_0001.tim"
        data=np.c_[elapsed_minutes,flows]
        np.savetxt(tim_fn,data)

