import xarray as xr
import numpy as np

## 
his_path='/media/hpc/opt/data/delft/sfb_dfm_v2/runs/wy2013a/DFM_OUTPUT_wy2013a/wy2013a_0000_20120801_000000_his.nc'

his_nc=xr.open_dataset(his_path)

## 

redwood_idx=np.nonzero( (his_nc.station_name=='NOAA_Redwood').values )[0][0]
alameda_idx=np.nonzero( (his_nc.station_name=='ALA').values )[0][0]

## 

redwood_s1=his_nc['waterlevel'].isel(stations=redwood_idx)
alameda_s1=his_nc['waterlevel'].isel(stations=alameda_idx)

## 

fig,ax=plt.subplots(1,1,num=1)
ax.cla()

redwood_s1.plot(ax=ax)
alameda_s1.plot(ax=ax)

## 

msl_delta=redwood_s1.mean() - alameda_s1.mean()

print "MSL(redwood) - MSL(alameda) = %.3fm"%(msl_delta)

