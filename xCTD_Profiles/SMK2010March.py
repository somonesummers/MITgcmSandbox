# %%
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import utm
import cmocean
import gsw

import geopandas as gpd

 # Paul believes these are IN the melange
# Access DOI: DATASET | Published 2021 | doi:10.18739/A2M03XZ2K
#  summary: Temperature and salinity data from 5 XCTDs and 2 XBTs collected in Sermilik Fjord, March 15th-16th, 2010. 
#    title: Profiles of temperature and salinity from Sermilik Fjord during March 2010
#    dimensions(sizes): z(901), profile(7)
#    variables(dimensions): 
# float64 profile(profile),
# float64 time(profile),
# float64 lat(profile),
# float64 lon(profile),
# float64 probe_type(profile),
# float64 depth(z),
# float64 sal(z, profile),
# float64 temp(z, profile)

#Ref for mapping
file2read = 'SF2015CTD.nc'
zfile_id = Dataset(file2read)
zlat = np.asarray(zfile_id.variables['lat'])
zlon = np.asarray(zfile_id.variables['lon'])
zxy = utm.from_latlon(zlat,zlon)
zutm_x = zxy[0]
zutm_y = zxy[1]



# Open the coast shapefile
shapefile_path = "Greenland_coast/Greenland_coast.shp"

gdf = gpd.read_file(shapefile_path)
newProj = gdf.to_crs(epsg=32624) #UTM N 24




file2read = 'SF2010xctd.nc'
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

pressure = np.zeros_like(salt)
for i in range(7):
    pressure[:,i] = depth[:]  * 1020 * 9.81 /1e4
# print(pressure[:,1])
density = np.zeros_like(salt)
CT = gsw.CT_from_t(salt,temp,pressure)  
# print(CT)
density = gsw.rho(salt,CT,0) - 1000
densityLevels = np.linspace(22,28,31)


Trange = np.linspace(-2,5,31)
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

picks = [3,4,5,2,1,0,6]
distAlong = np.zeros(len(picks))
for i in range(len(picks)):
    if(i == 0):
        distAlong[picks[i]] = 0
    else:
        distTemp = np.sqrt((utm_x[picks[i]] - utm_x[picks[i-1]])**2 + (utm_y[picks[i]] - utm_y[picks[i-1]])**2)
        distAlong[picks[i]] = distAlong[picks[i-1]] + distTemp /1e3

# gate for all TS profiles
if(True):
    for i in picks:
        plt.figure(figsize=(12, 6))
        
        plt.subplot(131)
        cp = plt.scatter(salt[:,i],temp[:,i],c=depth[:],linestyle='-',marker='o')
        cbar = plt.colorbar(cp) 
        cbar.set_label('Depth [m]')
        plt.plot(mixingS,mixingT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
        plt.plot(meltS,meltT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
        plt.plot(freezeS,freezeT,linewidth=.5,color='red',alpha=.5,linestyle='--')
        plt.plot(freezeS100,freezeT100,linewidth=.5,color='red',alpha=.5,linestyle='--')
        plt.title(f"xCTD {i}, {distAlong[i]:.2f} km along fjord")
        ax = plt.gca()
        ax.set_xlabel('Salinity')
        ax.set_ylabel('Temp [C]')
        ax.set_xlim([26,35.5])
        ax.set_ylim([-2,5])      

        plt.subplot(132)
        cp = plt.scatter(salt[:,i],depth[:],linestyle='-',marker='o')
        ax = plt.gca()
        ax.set_xlabel('Salinity')
        ax.set_ylabel('Depth [m]')
        plt.grid(alpha = .5)
        ax.invert_yaxis()

        plt.subplot(133)
        cp = plt.scatter(temp[:,i],depth[:],linestyle='-',marker='o')
        ax = plt.gca()
        ax.set_xlabel('Temp [C]')
        ax.set_ylabel('Depth [m]')
        ax.invert_yaxis()
        plt.grid(alpha = .5)

        plt.tight_layout()
        plt.show()
        plt.close()
        if(i == 6):
            np.savez('shelfProfile2015',T=temp[:,i],S=salt[:,i],z=depth[:],density=density[:,i])

# %%
#Sub 1, Along fjord
start = 0
stop = 6
sTime = datetime.datetime.fromtimestamp(int(time[start]))
eTime = datetime.datetime.fromtimestamp(int(time[stop]))


newProj.plot()
ax = plt.gca()
# plt.scatter(zutm_x,zutm_y,alpha=.5,color='black')
plt.plot(utm_x[picks],utm_y[picks])
f = plt.scatter(utm_x,utm_y,c=picks)
ax.set_xlim([505000,585500])
ax.set_ylim([7.26e6,7.37e6])
strname = 'Map2010AlongFjord.png'
plt.savefig(strname, format='png', dpi=400)

plt.colorbar(f)
plt.title('Our Most Exclusive CTD casts\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
picks2 = [3,5,2,1,6] #only picks with salinity



fig= plt.figure(2,figsize=(12, 8))
ax1 = plt.subplot(211)
cf = ax1.contourf(
    distAlong[picks2],
    depth[:],
    temp[:,picks2],
    Trange,
    cmap = 'cmo.thermal')
cc = plt.contour(
    distAlong[picks2],
    depth[:],
    density[:,picks2],
    densityLevels,
    colors='black',
    linewidths=0.5,
    alpha=0.5
)
plt.clabel(cc, inline=3, fontsize=8)
ax1.scatter(distAlong[picks2],800*np.ones_like(distAlong[picks2]),marker='*',color='black')
cbar = plt.colorbar(cf)
cbar.set_label('Temp')
# ax1.set_title('Temp\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
ax1.set_xlabel('Along Profile [km]')
ax1.set_ylabel('Depth')
ax1.invert_yaxis()

ax2 = plt.subplot(212)
cf = plt.contourf(
    distAlong[picks2],
    depth[:],
    salt[:,picks2],
    Srange,
    cmap='cmo.haline')
cc = plt.contour(
    distAlong[picks2],
    depth[:],
    density[:,picks2],
    densityLevels,
    colors='black',
    linewidths=0.5,
    alpha=0.5
)
ax2.scatter(distAlong[picks2],800*np.ones_like(distAlong[picks2]),marker='*',color='black')
ax2.clabel(cc, inline=3, fontsize=8)
cbar = plt.colorbar(cf)
cbar.set_label('Salt')
# ax2.set_title('Salt\n' + sTime.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + eTime.strftime('%Y-%m-%d %H:%M:%S'))
ax2.set_xlabel('Along Profile [km]')
ax2.set_ylabel('Depth')
ax2.invert_yaxis()
plt.tight_layout()

strname = '2010AlongFjord.png'
plt.savefig(strname, format='png', dpi=400)
plt.show()
plt.close()


# %%


# %%



