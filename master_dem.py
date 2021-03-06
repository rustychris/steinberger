"""
Importable field, for multiprocessing while producing tiles
"""
import os
import numpy as np

from stompy.spatial import field
from stompy.spatial import interp_coverage

def factory(attrs):
    geo_bounds=attrs['geom'].bounds

    if attrs['src_name']=='usgs_2m_topobathy_steinberger':
        fns=['steinberger-dem-crop.tif',
             ('/media/idrive/BASELAYERS/Elevation_DerivedProducts/'
              'LiDAR 2005-2012 entire Bay Area from AECOM/USGS_TopoBathy/'
              'San_Francisco_TopoBathy_Elevation_2m.tif') ]
        for fn in fns:
            if os.path.exists(fn):
                return field.GdalGrid(fn,geo_bounds=geo_bounds)
        else:
            raise Exception("Couldn't find the USGS 2m data")
    if attrs['src_name'].startswith('py:'):
        expr=attrs['src_name'][3:]
        # something like 'ConstantField(-1.0)'
        # a little sneaky... make it look like it's running
        # after a "from stompy.spatial.field import *"
        # and also it gets fields of the shapefile
        field_hash=dict(field.__dict__)
        # convert the attrs into a dict suitable for passing to eval
        attrs_dict={}
        for name in attrs.dtype.names:
            attrs_dict[name]=attrs[name]

        return eval(expr,field_hash,attrs_dict)
        
    assert False

src_shp='dem_sources-v00.shp'

mbf=field.CompositeField(shp_fn=src_shp,
                         factory=factory,
                         priority_field='priority',
                         data_mode='data_mode',
                         alpha_mode='alpha_mode')

dem=mbf.to_grid(dx=2,dy=2,
                bounds=[566300,573500,4149200,4155600])
fn='master_dem_v00.tif'
os.path.exists(fn) and os.unlink(fn)
dem.write_gdal(fn)
