# %%
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import datetime
import utm
import cmocean
import gsw
import argparse
import geopandas as gpd

parser = argparse.ArgumentParser(description='Plot xCTDs from Sermilik')
parser.add_argument('-y','--year', nargs='?', type=str ,default = None,
                    help='Specify year to plot [default = 2010]')
parser.add_argument('-p','--profile', action='count', default=0,
                    help='option to show all profiles')
args = parser.parse_args()

 
 #DATASET | All found on this site | doi:10.18739/A2348GJ3X
 # contributor_name: Fiamma Straneo, 
 #    summary: Temperature and salinity profiles from  eXpendable Conductivity Temperature Depth (XCTD)-1 probes collected in Sermilik Fjord, East Greenland,
 #    title: XCTD profiles of temperatu
 #    dimensions(sizes): z(#), profile(#)
 #    variables(dimensions): 
 # float64 profile(profile),
 # float64 time(profile),
 # float64 lat(profile),
 # float64 lon(profile),
 # float64 depth(z),
 # float64 pres(z),
 # float64 sal(z, profile),
 # float64 temp(z, profile)

# Open the coast shapefile
shapefile_path = "Greenland_coast/Greenland_coast.shp"

gdf = gpd.read_file(shapefile_path)
newProj = gdf.to_crs(epsg=32624) #UTM N 24
if(args.year == None):
    year = '2012'
else:
    year = args.year
file2read = f'Sermilik{year}_XCTD.nc'
file_id = Dataset(file2read)
time = file_id.variables['time']
lat = np.asarray(file_id.variables['lat'])
lon = np.asarray(file_id.variables['lon'])
xy = utm.from_latlon(lat,lon)
utm_x = xy[0]
utm_y = xy[1]
if(year=='2015SF'):
    depth = file_id.variables['z']
else:
    depth = file_id.variables['depth']
temp = np.asarray(file_id.variables['temp'])
salt = np.asarray(file_id.variables['sal'])
# pressure = file_id.variables['pres']

pressure = np.zeros_like(salt)
for i in range(len(lat)):
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

if(args.profile > 0):
    for i in range(len(lat)):
        plt.figure()
        cp = plt.scatter(salt[:,i],temp[:,i],c=depth[:],linestyle='-',marker='o')
        plt.colorbar(cp) 
        plt.plot(mixingS,mixingT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
        plt.plot(meltS,meltT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
        plt.plot(freezeS,freezeT,linewidth=.5,color='red',alpha=.5,linestyle='--')
        plt.plot(freezeS100,freezeT100,linewidth=.5,color='red',alpha=.5,linestyle='--')
        plt.title(i)
        ax = plt.gca()
        ax.set_xlim([26,35.5])
        ax.set_ylim([-2.5,5])
        plt.show()
        plt.close()
if(args.profile > 0):
    for i in range(len(lat)):
        plt.figure()
        cp = plt.plot(temp[:,i],-depth[:],linestyle='-',marker='o')
        plt.title(i)
        plt.show()
        plt.close()

# %%
#Sub 1, Along fjord
start = 0
stop = len(lat)-1
sTime = datetime.datetime.fromtimestamp(int(time[start]))
eTime = datetime.datetime.fromtimestamp(int(time[stop]))


#Picks grab sub set, throwing out bad data, and also sets order such that they go down fjord as we go
if(year == '2010marSF'):
    picks = [3,5,2,1,6]
elif(year == '2010m'):
    picks = [6,1,0,2,5,4,3]
    picks = picks[::-1]
elif(year == '2015SF'):
    picks = [10,9,8,12,14,15,16,19,21,22,29]
elif(year == '2016'):
    # picks = [3,1,4]
    picks = [3,1,4,5,6,7,8,9,10,11]
elif(year == '2019'):
    picks = [1,2,3,4,0] #6,5 is the plume I think, so thats cool
elif(year == '2021'):
    picks = [2,3,4,5]
elif(year == '2023'):
    picks = [0,1,2,3,4,5,9,10,11]
else:
    picks = np.arange(0,len(lat))



distAlong = np.zeros(len(picks))
for i in range(len(picks)):
    if(i == 0):
        distAlong[i] = 0
    else:
        distTemp = np.sqrt((utm_x[picks[i]] - utm_x[picks[i-1]])**2 + (utm_y[picks[i]] - utm_y[picks[i-1]])**2)
        distAlong[i] = distAlong[i-1] + distTemp /1e3

newProj.plot()
# plt.scatter(zutm_x,zutm_y,alpha=.5,color='black')
f=plt.scatter(utm_x[picks],utm_y[picks],c=distAlong)
plt.plot(utm_x[picks],utm_y[picks])
plt.colorbar(f)
plt.title(f'Our Most Exclusive CTD casts {year}')
ax = plt.gca()
ax.set_xlim([505000,585500])
ax.set_ylim([7.26e6,7.37e6])
strname = f'{year}Map.png'
plt.savefig(strname, format='png', dpi=400)

# distAlong = [1,2,3,4,5,6,7]
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
plt.suptitle(f'Along Fjord Profile from {year}')
plt.clabel(cc, inline=3, fontsize=8)
ax1.scatter(distAlong,800*np.ones_like(distAlong),marker='*',color='black')
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

strname = f'{year}AlongFjord.png'
plt.savefig(strname, format='png', dpi=400)

# fig= plt.figure(3,figsize=(12, 8))
# f=plt.scatter(lat[picks],lon[picks],c=picks)
# plt.plot(lat[picks],lon[picks])
# plt.colorbar(f)

plt.show()
plt.close()


# %%


# %%



