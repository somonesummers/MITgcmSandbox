# %%
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import utm
import cmocean
import gsw

 # summary: Temperature and salinity profiles from 7 eXpendable Conductivity Temperature Depth (XCTD)-1 probes collected in Sermilik Fjord, East Greenland, on September 14 2011.
 #   title: XCTD profiles of temperature and salinity from Sermilik Fjord during September 2012
 #   dimensions(sizes): z(450), profile(7)
 #   dimensions(sizes): z(450), profile(7)
 #   variables(dimensions): 
 # float64 profile(profile), 
 # float64 time(profile), 
 # float64 lat(profile), 
 # float64 lon(profile), 
 # float64 depth(z), 
 # float64 pres(z), 
 # float64 sal(z, profile), 
 # float64 temp(z, profile)


file2read = 'Sermilik2012_XCTD.nc'
file_id = Dataset(file2read)
time = file_id.variables['time']
lat = np.asarray(file_id.variables['lat'])
lon = np.asarray(file_id.variables['lon'])
xy = utm.from_latlon(lat,lon)
utm_x = xy[0]
utm_y = xy[1]
depth = file_id.variables['depth']
temp = np.asarray(file_id.variables['temp'])
salt = np.asarray(file_id.variables['sal'])
pressure = file_id.variables['pres']

pressure = np.zeros_like(salt)
for i in range(7):
    pressure[:,i] = depth[:]  * 1020 * 9.81 /1e4
# print(pressure[:,1])
density = np.zeros_like(salt)
CT = gsw.CT_from_t(salt,temp,pressure)  
# print(CT)
density = gsw.rho(salt,CT,0) - 1000
densityLevels = np.linspace(22,28,31)


Trange = np.linspace(-2,4,31)
Srange = np.linspace(27,35,31)

# %%
#toss random very neg values

lat = lat[lat[:] > -90]
lon = lon[lon[:] > -90]
print(np.shape(lat))
print(np.shape(lon))

# %%
nMix = 25
mixingT = np.zeros([2,nMix])
mixingS = np.zeros([2,nMix])
mixingT[1,:] = np.linspace(-10,10,nMix)
mixingS[1,:] = np.ones([1,nMix])*50

nMelt = 25
meltT = np.ones([2,nMelt])*-90
meltS = np.zeros([2,nMelt])
meltT[1,:] = np.ones([1,nMelt])*10
meltS[1,:] = np.linspace(-10,10,nMelt)+32

#(s,t)
freezeS = [0,50]
freezeT = [0,-2.809]
freezeS100 = [0,50]
freezeT100 = [-.011,-2.919]


# for i in range(7):
#     plt.figure(i)
#     cp = plt.scatter(salt[:,i],temp[:,i],c=depth[:],linestyle='-',marker='o')
#     plt.colorbar(cp) 
#     plt.plot(mixingS,mixingT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
#     plt.plot(meltS,meltT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
#     plt.plot(freezeS,freezeT,linewidth=.5,color='red',alpha=.5,linestyle='--')
#     plt.plot(freezeS100,freezeT100,linewidth=.5,color='red',alpha=.5,linestyle='--')
#     ax = plt.gca()
#     ax.set_xlim([26,35.5])
#     ax.set_ylim([-2.5,3.75])      
#     plt.show()
#     plt.close()

# %%
#Sub 1, Along fjord
start = 0
stop = 6
sTime = datetime.datetime.fromtimestamp(int(time[start]))
eTime = datetime.datetime.fromtimestamp(int(time[stop]))


picks = [0,1,2,3,4,5,6]
distAlong = np.zeros(len(picks))
for i in range(len(picks)):
    if(i == 0):
        distAlong[i] = 0
    else:
        distTemp = np.sqrt((utm_x[picks[i]] - utm_x[picks[i-1]])**2 + (utm_y[picks[i]] - utm_y[picks[i-1]])**2)
        distAlong[i] = distAlong[i-1] + distTemp /1e3

plt.figure(1)
f=plt.scatter(utm_x,utm_y,c=picks)
plt.colorbar(f)
plt.title('Our Most Exclusive CTD casts\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))


fig= plt.figure(2,figsize=(12, 8))
ax1 = plt.subplot(211)
cf = ax1.contourf(
    distAlong,
    depth[:],
    temp[:,picks],
    Trange,
    cmap = 'cmo.thermal')
cc = plt.contour(
    distAlong,
    depth[:],
    density[:,picks],
    densityLevels,
    colors='black',
    linewidths=0.5,
    alpha=0.5
)
plt.clabel(cc, inline=3, fontsize=8)
cbar = plt.colorbar(cf)
cbar.set_label('Temp')
# ax1.set_title('Temp\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
ax1.set_xlabel('Along Profile [km]')
ax1.set_ylabel('Depth')
ax1.invert_yaxis()

ax2 = plt.subplot(212)
cf = plt.contourf(
    distAlong,
    depth[:],
    salt[:,picks],
    Srange,
    cmap='cmo.haline')
cc = plt.contour(
    distAlong,
    depth[:],
    density[:,picks],
    densityLevels,
    colors='black',
    linewidths=0.5,
    alpha=0.5
)
ax2.scatter(distAlong,800*np.ones_like(distAlong),marker='*',color='black')
ax2.clabel(cc, inline=3, fontsize=8)
cbar = plt.colorbar(cf)
cbar.set_label('Salt')
# ax2.set_title('Salt\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
ax2.set_xlabel('Along Profile [km]')
ax2.set_ylabel('Depth')
ax2.invert_yaxis()
plt.tight_layout()
plt.show()


# %%


# %%



