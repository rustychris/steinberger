import xarray as xr
import numpy as np
from stompy import utils

## 
# looks like SFB1304 is here in the Slough -
adcp=xr.open_dataset('/home/rusty/data/noaa/ports/SFB1304-2013.nc')


## 

U=np.array( [adcp.u_davg,adcp.v_davg] ).T

## 

theta=utils.principal_theta(U)

print "Rotating to positive is %d compass"%((90-(theta * 180/np.pi)) % 360)

# theta gives the direction of flood-positive currents.  To get that
# component into the 'x' component, rotate by -theta
Urot=utils.rot(-theta,U)

## 

fig=plt.figure(2)
fig.clf()
fig.set_size_inches([14,8],forward=True)

fig,(ax_spring,ax_neap)=plt.subplots(2,1,sharey=True,num=2)

for ax in [ax_spring,ax_neap]:
    ax.plot( adcp.time,Urot[:,0],label='Depth avg vel.')
    ax.set_ylabel('m/s floodward')
    ax.grid(1)

# starts with ebb, ends with flood
spring=[735041.84472477518,
        735044.96242746431]
neap=[735033.47076959629,
      735036.62203205191]
ax_spring.axis(xmin=spring[0],xmax=spring[1])
ax_neap.axis( xmin=neap[0],xmax=neap[1])

fig.tight_layout()

adcp_dn=utils.to_dnum(adcp.time)

spring_sel=(adcp_dn>=spring[0])&(adcp_dn<spring[1])
neap_sel  =(adcp_dn>=neap[0])&(adcp_dn<neap[1])
flood_sel=Urot[:,0]>0
ebb_sel=Urot[:,0]<0

spring_flood=np.mean( Urot[spring_sel&flood_sel,0] )
spring_ebb=np.mean( Urot[spring_sel&ebb_sel,0] )
neap_flood=np.mean( Urot[neap_sel&flood_sel,0] )
neap_ebb=np.mean( Urot[neap_sel&ebb_sel,0] )

# Pick a greater spring:
big_spring=[735043.09911698522, 735043.6522579958]
big_spring_sel=(adcp_dn>=big_spring[0])&(adcp_dn<big_spring[1])

small_neap=[735035.31628801185, 735035.80342522147]
small_neap_sel=(adcp_dn>=small_neap[0])&(adcp_dn<small_neap[1])

big_spring_flood=np.mean( Urot[big_spring_sel&flood_sel,0] )
big_spring_ebb=np.mean( Urot[big_spring_sel&ebb_sel,0] )

small_neap_flood=np.mean( Urot[small_neap_sel&flood_sel,0] )
small_neap_ebb=np.mean( Urot[small_neap_sel&ebb_sel,0] )

print "Spring mean velocity: %.2f %.2f"%(spring_flood,spring_ebb)
print "Big spring mean velocity: %.2f %.2f"%(big_spring_flood,big_spring_ebb)
print "Neap mean velocity: %.2f %.2f"%(neap_flood,neap_ebb)
print "Small neap mean velocity: %.2f %.2f"%(small_neap_flood,small_neap_ebb)



## 

fig.savefig('sfb1304-flood_directed_current-spring_neap.png')
