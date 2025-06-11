# %%
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import cmocean
import gsw
import utm

import geopandas as gpd
## 10.1029/2018GL077000 is paper 
## https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0171277 has data
#
# summary: 31 Temperature and salinity profiles from Sermilik Fjord in East Greenland.  
# Collected with a Seabird SBE25plus CTD (#0251108) during August 2015 from the R/V Adolf Jensen.  
# The profiles are 1-m bin averages.
#    acknowledgement: Funding from NSF OCE 1434041 and NSF OCE 1536856
#    references: One manuscript that uses and describes this data is: Beaird, N., F. Straneo, and W. Jenkins (2018): "Export of strongly diluted Greenland meltwater from a major glacial fjord."
#    keywords: Oceans > Ocean Temperature > Water Temperature, Oceans > Salinity/Density > Salinity
#    title: Temperature and salinity profiles from Sermilik Fjord, East Greenland, collected with a Seabird SBE25plus CTD (#0251108) during August 2015
#    dimensions(sizes): z(900), profile(31)
#    variables(dimensions): 
#       float64 profile(profile), 
#       float64 time(profile), 
#       float64 lat(profile), 
#       float64 lon(profile), 
#       float64 z(z), 
#       float64 sal(z, profile), 
#       float64 sal_qc(z, profile), 
#       float64 temp(z, profile), 
#       float64 temp_qc(z, profile)

# Open the coast shapefile
shapefile_path = "Greenland_coast/Greenland_coast.shp"

gdf = gpd.read_file(shapefile_path)
newProj = gdf.to_crs(epsg=32624) #UTM N 24


file2read = 'SF2015CTD.nc'
file_id = Dataset(file2read)
time = file_id.variables['time']
lat = np.asarray(file_id.variables['lat'])
lon = np.asarray(file_id.variables['lon'])
xy = utm.from_latlon(lat,lon)
utm_x = xy[0]
utm_y = xy[1]
depth = np.asarray(file_id.variables['z'])
temp = np.asarray(file_id.variables['temp'])
temp_qc = file_id.variables['temp']
salt = np.asarray(file_id.variables['sal'])
profile = np.asarray(file_id.variables['profile'])

pressure = np.zeros_like(salt)
for i in range(31):
    pressure[:,i] = depth[:]  * 1020 * 9.81 /1e4
# print(pressure[:,1])
density = np.zeros_like(salt)
CT = gsw.CT_from_t(salt,temp,pressure)  
# print(CT)
density = gsw.rho(salt,CT,0) - 1000
densityLevels = np.linspace(25,28,16)


# %%
# print(np.max(depth))
# print(np.shape(lat))
# print(np.shape(lon))
# print(np.shape(profile))
# print(lat,lon,profile)

Trange = np.linspace(-2,4,31)
Srange = np.linspace(27,35,31)

## show the entire domain
#Aug 2014 is 450 - 625 or so
newProj.plot()
ax = plt.gca()
f=plt.scatter(utm_x,utm_y,c=profile)
picks = [8,9,10]
plt.scatter(utm_x[picks],utm_y[picks],marker='*',color='black')
plt.colorbar(f)
plt.title('All Beautiful 2015 CTD casts')
plt.xlabel('Easting UTM [m]')
plt.ylabel('Northing UTM [m]')
ax.set_xlim([505000,585500])
ax.set_ylim([7.26e6,7.37e6])
xRange = plt.xlim()
yRange = plt.ylim()
strname = 'Map2015All.png'
plt.savefig(strname, format='png', dpi=400)
plt.show()
plt.close()
# %%
nMix = 25 #number mixing lines
mixingT = np.zeros([2,nMix])
mixingS = np.zeros([2,nMix])
mixingT[1,:] = np.linspace(-10,10,nMix)
mixingS[1,:] = np.ones([1,nMix])*50

nMelt = 25 #number melt lines
meltT = np.ones([2,nMelt])*-90
meltS = np.zeros([2,nMelt])
meltT[1,:] = np.ones([1,nMelt])*10
meltS[1,:] = np.linspace(-10,10,nMelt)+32

#(s,t) #freezing lines
freezeS = [0,50]
freezeT = [0,-2.809]
freezeS100 = [0,50]
freezeT100 = [-.011,-2.919]

## Pick a few fun TS and density profiles
for i in [29]:
    plt.figure(i)
    ax1 = plt.subplot(121)
    cp = ax1.scatter(salt[:,i],temp[:,i],c=depth[:],linestyle='-',marker='o')
    plt.colorbar(cp) 
    ax1.plot(mixingS,mixingT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
    ax1.plot(meltS,meltT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
    ax1.plot(freezeS,freezeT,linewidth=.5,color='red',alpha=.5,linestyle='--')
    ax1.plot(freezeS100,freezeT100,linewidth=.5,color='red',alpha=.5,linestyle='--')
    ax1.set_title(i)
    ax1.set_xlim([26,35.5])
    ax1.set_ylim([-2.5,4.75]) 
    ax1.set_ylabel('Temp')
    ax1.set_xlabel('Salinity')
    
    ax2 = plt.subplot(122)
    ax2.plot(density[:,i],depth)
    ax2.set_ylabel('Depth')
    ax2.set_xlabel('Density')
    ax2.set_ylim([-10,900])
    ax2.set_xlim([25,28])
    ax2.invert_yaxis()
    plt.tight_layout()
    plt.show()
    plt.close()
    np.savez('shelfProfile2015',T=temp[:,i],S=salt[:,i],z=depth[:],density=density[:,i])
# %%

#Sub 0, on shelf
# start = 28
# stop = 31
# order = [1,2,3]
# sTime = datetime.datetime.fromtimestamp(int(time[start]))
# eTime = datetime.datetime.fromtimestamp(int(time[stop-1]))


# plt.figure(1)
# f = plt.scatter(lon[start:stop],lat[start:stop],c=order)
# plt.plot(lon[start:stop],lat[start:stop])
# plt.colorbar(f)
# plt.xlim(xRange)
# plt.ylim(yRange)
# plt.title('Our Most Exclusive CTD casts\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))

# fig = plt.figure(2)
# cf = plt.contourf(
#     order,
#     depth[:],
#     temp[:,start:stop],
#     Trange,
#     cmap = 'cmo.thermal')
# plt.colorbar(cf)
# cc = plt.contour(
#     order,
#     depth[:],
#     density[:,start:stop],
#     densityLevels,
#     colors='black',
#     linewidths=0.5,
#     alpha=0.5
# )
# plt.clabel(cc, inline=3, fontsize=8)
# plt.title('Temp\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
# plt.xlabel('profile')
# plt.ylabel('Depth')
# ax = fig.gca()
# ax.invert_yaxis()

# fig = plt.figure(3)
# cf = plt.contourf(
#     order,
#     depth[:],
#     salt[:,start:stop],
#     Srange,
#     cmap='cmo.haline')
# cc = plt.contour(
#     order,
#     depth[:],
#     density[:,start:stop],
#     densityLevels,
#     colors='black',
#     linewidths=0.5,
#     alpha=0.5
# )
# plt.clabel(cc, inline=3, fontsize=8)
# plt.colorbar(cf)
# plt.title('Salt\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
# plt.xlabel('profile')
# plt.ylabel('Depth')
# ax = fig.gca()
# ax.invert_yaxis()
# plt.show()
# plt.close()
# %%



#Sub 1, Along fjord
picks = [10,9,8,12,14,15,16,19,21,22,29]
order = [1,2,3,4,5,6,7,8,9,10,12]
distAlong = np.zeros(len(picks))
for i in range(len(picks)):
    if(i == 0):
        distAlong[i] = 0
    else:
        distTemp = np.sqrt((utm_x[picks[i]] - utm_x[picks[i-1]])**2 + (utm_y[picks[i]] - utm_y[picks[i-1]])**2)
        distAlong[i] = distAlong[i-1] + distTemp /1e3
sTime = datetime.datetime.fromtimestamp(int(time[picks[0]]))
eTime = datetime.datetime.fromtimestamp(int(time[picks[-1]]))
print(distAlong)
newProj.plot()
f=plt.scatter(utm_x[picks],utm_y[picks],c=distAlong,cmap='jet')
plt.plot(utm_x[picks],utm_y[picks])
plt.xlabel('Easting UTM [m]')
plt.ylabel('Northing UTM [m]')
plt.colorbar(f)
plt.xlim(xRange)
plt.ylim(yRange)
plt.title('Our Most Exclusive xCTDs \n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
strname = 'Map2015AlongFjord.png'
plt.savefig(strname, format='png', dpi=400)

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

strname = '2015AlongFjordCustom.png'
plt.savefig(strname, format='png', dpi=400)
plt.show()
plt.close()
# %%
#Sub 1, midfjord cross section
# start = 18
# stop = 21
# order = [1,2,3]
# sTime = datetime.datetime.fromtimestamp(int(time[start]))
# eTime = datetime.datetime.fromtimestamp(int(time[stop-1]))

# plt.figure()
# f=plt.scatter(utm_x[start:stop],utm_y[start:stop],c=order)
# plt.plot(utm_x[start:stop],utm_y[start:stop])
# plt.colorbar(f)
# plt.xlim(xRange)
# plt.ylim(yRange)
# plt.title('Our Most Exclusive CTD casts\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))


# fig= plt.figure()
# cf = plt.contourf(
#     order,
#     depth[:],
#     temp[:,start:stop],
#     Trange,
#     cmap = 'cmo.thermal')
# cc = plt.contour(
#     order,
#     depth[:],
#     density[:,start:stop],
#     densityLevels,
#     colors='black',
#     linewidths=0.5,
#     alpha=0.5
# )
# plt.clabel(cc, inline=3, fontsize=8)
# plt.colorbar(cf)
# plt.title('Temp\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
# plt.xlabel('profile')
# plt.ylabel('Depth')
# ax = fig.gca()
# ax.invert_yaxis()

# fig= plt.figure()
# cf = plt.contourf(
#     order,
#     depth[:],
#     salt[:,start:stop],
#     Srange,
#     cmap='cmo.haline')
# cc = plt.contour(
#     order,
#     depth[:],
#     density[:,start:stop],
#     densityLevels,
#     colors='black',
#     linewidths=0.5,
#     alpha=0.5
# )
# plt.clabel(cc, inline=3, fontsize=8)
# plt.colorbar(cf)
# plt.title('Salt\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
# plt.xlabel('profile')
# plt.ylabel('Depth')
# ax = fig.gca()
# ax.invert_yaxis()
# plt.show()
# # %%



