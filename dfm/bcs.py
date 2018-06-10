from __future__ import print_function

import os

import numpy as np
import stompy.model.delft.io as dio
from stompy.io.local import noaa_coops
from stompy import utils, filters
from stompy.model.delft import dfm_grid

cache_dir='cache'
os.path.exists(cache_dir) or os.mkdir(cache_dir)

class BC(object):
    pass

class Constant(BC):
    def __init__(self,var_name,value):
        """
        quantity: 'salinity','temperature'
        value: floating point
        """
        self.var_name=var_name
        self.value=value
    def write(self,mdu,feature,grid):
        print("Feature: %s"%(feature['name']))

        name=feature['name']
        old_bc_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

        assert feature['geom'].type=='LineString'
        pli_data=[ (name, np.array(feature['geom'].coords)) ]
        base_fn=os.path.join(mdu.base_path,"%s_%s"%(name,self.var_name))
        pli_fn=base_fn+'.pli'
        dio.write_pli(pli_fn,pli_data)

        if self.var_name=='salinity':
            quant='salinitybnd'
        elif self.var_name=='temperature':
            quant='temperaturebnd'
        else:
            assert False

        with open(old_bc_fn,'at') as fp:
            lines=["QUANTITY=%s"%quant,
                   "FILENAME=%s_%s.pli"%(name,self.var_name),
                   "FILETYPE=9",
                   "METHOD=3",
                   "OPERAND=O",
                   ""]
            fp.write("\n".join(lines))

        self.write_data(mdu,feature,self.var_name,base_fn)

    def write_data(self,mdu,feature,var_name,base_fn):
        ref_date,start_date,end_date = mdu.time_range()
        period=np.array([start_date,end_date])
        elapsed_minutes=(period - ref_date)/np.timedelta64(60,'s')

        # just write a single node
        tim_fn=base_fn + "_0001.tim"
        with open(tim_fn,'wt') as fp:
            for t in elapsed_minutes:
                fp.write("%g %g\n"%(t,self.value))
        
class NoaaTides(BC):
    var_names=['ssh']
    def __init__(self,station,datum='NAVD88',z_offset=0.0):
        self.station=station
        self.datum=datum
        self.z_offset=z_offset
    def write(self,mdu,feature,grid):
        print("Feature: %s"%(feature['name']))

        name=feature['name']
        old_bc_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

        for var_name in self.var_names:
            if feature['geom'].type=='LineString':
                pli_data=[ (name, np.array(feature['geom'].coords)) ]
                base_fn=os.path.join(mdu.base_path,"%s_%s"%(name,var_name))
                pli_fn=base_fn+'.pli'
                dio.write_pli(pli_fn,pli_data)

                if var_name=='ssh':
                    quant='waterlevelbnd'
                else:
                    assert False

                with open(old_bc_fn,'at') as fp:
                    lines=["QUANTITY=%s"%quant,
                           "FILENAME=%s_%s.pli"%(name,var_name),
                           "FILETYPE=9",
                           "METHOD=3",
                           "OPERAND=O",
                           ""]
                    fp.write("\n".join(lines))

                self.write_data(mdu,feature,var_name,base_fn)
            else:
                assert False

    def write_data(self,mdu,feature,var_name,base_fn):
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
    dredge_depth=-1.0
    storm_flow=10.0
    storm_duration_h=3.0
    storm_start_h=48.0
    
    def __init__(self,name):
        self.name=name
    def write(self,mdu,feature,grid):
        # obvious copy and paste from above.
        # not quite ready to abstract, though
        print("Feature: %s"%(feature['name']))

        name=feature['name']
        old_bc_fn=mdu.filepath( ['external forcing','ExtForceFile'] )

        for var_name in self.var_names:
            if feature['geom'].type=='LineString':
                pli_data=[ (name, np.array(feature['geom'].coords)) ]
                base_fn=os.path.join(mdu.base_path,"%s_%s"%(name,var_name))
                pli_fn=base_fn+'.pli'
                dio.write_pli(pli_fn,pli_data)

                if var_name=='q':
                    quant='dischargebnd'
                else:
                    assert False

                with open(old_bc_fn,'at') as fp:
                    lines=["QUANTITY=%s"%quant,
                           "FILENAME=%s_%s.pli"%(name,var_name),
                           "FILETYPE=9",
                           "METHOD=3",
                           "OPERAND=O",
                           ""]
                    fp.write("\n".join(lines))

                self.write_data(mdu,feature,var_name,base_fn)

                dfm_grid.dredge_boundary(grid,pli_data[0][1],self.dredge_depth)
            else:
                assert False

    def write_data(self,mdu,feature,var_name,base_fn):
        ref_date,run_start,run_stop=mdu.time_range()

        # trapezoid hydrograph
        times=np.array( [run_start,
                         run_start+np.timedelta64(self.storm_start_h-1,'h'),
                         run_start+np.timedelta64(self.storm_start_h,'h'),
                         run_start+np.timedelta64(self.storm_start_h+self.storm_duration_h,'h'),
                         run_start+np.timedelta64(self.storm_start_h+self.storm_duration_h+1,'h'),
                         run_stop+np.timedelta64(1,'D')] )
        flows=np.array( [0.0,0.0,
                         self.storm_flow,self.storm_flow,0.0,0.0] )
        elapsed_minutes=(times - ref_date)/np.timedelta64(60,'s')

        # just write a single node
        tim_fn=base_fn + "_0001.tim"
        data=np.c_[elapsed_minutes,flows]
        np.savetxt(tim_fn,data)

