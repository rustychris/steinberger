import xarray as xr
from shapely import geometry
from stompy import utils
from stompy.spatial import wkb2shp, field
from stompy.grid import unstructured_grid
from stompy.spatial import linestring_utils

## 
bounds=wkb2shp.shp2geom('region-bounds-v00.shp')

## 

g=unstructured_grid.UnstructuredGrid()

for geo in bounds['geom']:
    pnts=np.array(geo)

    A=pnts[:-1]
    B=pnts[1:]

    for a,b in zip(A,B):
        na=g.add_or_find_node(x=a)
        nb=g.add_or_find_node(x=b)
        try:
            j=g.add_edge(nodes=[na,nb])
        except g.GridException:
            pass

## 

cycles=g.find_cycles(max_cycle_len=g.Nnodes())

## 

polys=[geometry.Polygon( g.nodes['x'][cycle] )
       for cycle in cycles ]

wkb2shp.wkb2shp('regions_from_bounds.shp',
                polys,overwrite=True)

## 
from shapely import ops
one_poly=ops.cascaded_union(polys)
bounds_xyxy=one_poly.bounds

## 

dem=field.GdalGrid( ('/media/idrive/BASELAYERS/Elevation_DerivedProducts/'
                     'LiDAR 2005-2012 entire Bay Area from AECOM/USGS_TopoBathy/'
                     'San_Francisco_TopoBathy_Elevation_2m.tif'),
                    geo_bounds=[bounds_xyxy[0],bounds_xyxy[2],
                                bounds_xyxy[1],bounds_xyxy[3]])
## 


plt.figure(1).clf()
g.plot_edges()
dem.plot(vmin=-2,vmax=2.0,cmap='jet',interpolation='nearest')

## 

# 2. Loop over polygons, calculating wet area and volume for each tidal datum

datums=[('MHHW',2.155),
        ('MHW',1.962),
        ('MSL',0.996),
        ('MLW',0.020),
        ('MLLW',-0.346)]

# 4s
poly_masks=[dem.polygon_mask(poly)
            for poly in polys]

## 
import xarray as xr

dem_pix_A=dem.dx*dem.dy

ds=xr.Dataset()
ds['poly']=('poly',),polys
ds['datum']=('datum',),[d[0] for d in datums]
ds['z']=('datum',),[d[1] for d in datums]
ds['wet_area']=('poly','datum'),np.zeros( (len(ds.poly),len(ds.datum)) )
ds['volume']=('poly','datum'),ds.wet_area.values.copy()


for datum_i,(datum_name,datum_z) in enumerate(datums):
    for poly_i,poly_mask in enumerate(poly_masks):
        
        area=dem_pix_A * (dem.F[poly_mask]<=datum_z).sum()
        volume=dem_pix_A * (datum_z-dem.F[poly_mask]).clip(0,np.inf).sum()
        ds.wet_area.values[poly_i,datum_i]=area
        ds.volume.values[poly_i,datum_i]=volume

## 


one_poly_mask=dem.polygon_mask(one_poly)

dem_mask=dem.copy()
dem_mask.F=np.where(one_poly_mask,dem.F,np.nan)

## 

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)
fig.set_size_inches([11,8],forward=True)

ecoll=g.plot_edges(ax=ax,color='k')

img=dem.downsample(4).plot(vmin=-8,vmax=5.0,cmap='gray',
                           interpolation='nearest',ax=ax)

boundaries=np.sort(ds.z.values)
cset=dem_mask.contourf(boundaries,ax=ax,cmap='winter')
cbar=plt.colorbar(cset,spacing='proportional')
ax.xaxis.set_visible(0)
ax.yaxis.set_visible(0)
ax.axis('scaled')

labels=["% .2fm %s"%(datum_z,datum_name)
        for datum_name,datum_z in datums[::-1]]
cbar.set_ticklabels(labels)
fig.tight_layout()
fig.savefig('regions_and_datums.png',dpi=200)

## 
recs=[]
for poly_i in range(len(polys)):
    rec={}
    for datum in ds.datum:
        datum_name=datum.item()
        rec['A_%s'%datum_name]= float(ds.isel(poly=poly_i).sel(datum=datum).wet_area.values)
        rec['V_%s'%datum_name]= float(ds.isel(poly=poly_i).sel(datum=datum).volume.values)
    rec['id']=poly_i
    recs.append(rec)

wkb2shp.wkb2shp('results_v00.shp',
                polys,fields=recs,
                overwrite=True)


## 

# Profiles:
redwood=[ [570627,4153666],
          [570994,4153282] ]
steinberger=[ [568364,4155573],
              [568644,4155351] ]

def pnts_to_xA(pnts,datum='MSL'):
    samples=linestring_utils.upsample_linearring(pnts,5,closed_ring=False)
    z=dem(samples)
    dists=utils.dist_along(samples)
    z0=ds.z.sel(datum=datum).item()
    depths=(z0-z).clip(0,np.inf)

    xA=np.trapz(depths,dists)
    return xA

# hand-calc of prism for the combined Redwood polygons
redwood_prism=7.97e6 # m3
steinberger_prism=8.39e6 # m3
# m3 -- approx. for observed drainage divides
offset=3.9e6 
redwood_prism+=offset
steinberger_prism-=offset
redwood_xA=pnts_to_xA(redwood)
steinberger_xA=pnts_to_xA(steinberger)

redwood_Utidal=redwood_prism / redwood_xA / (6*3600.)
redwood_Utidal_max=np.pi/2 * redwood_Utidal

steinberger_Utidal=steinberger_prism / steinberger_xA / (6*3600.)
steinberger_Utidal_max=np.pi/2 * steinberger_Utidal

m3_to_ft3=35.314667
print "Redwood:     %.2f mean, %.2f max"%(redwood_Utidal,redwood_Utidal_max)
print "    mean tidal flux: %.1f cfs"%( redwood_prism / (6*3600.) * m3_to_ft3 )
print "Steinberger: %.2f mean, %.2f max"%(steinberger_Utidal,steinberger_Utidal_max)
print "    mean tidal flux: %.1f cfs"%( steinberger_prism / (6*3600.) * m3_to_ft3 )



## 

