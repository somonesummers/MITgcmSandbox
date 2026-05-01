from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import xarray as xr
import argparse
import gsw

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for fresh water flux over time')
parser.add_argument('-s','--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-n','--numFrames', nargs=1, default=[250],type=int,
                    help='Max number of samples from time series [defaut = 250]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()


# Pick cross section to view from file or default
plotDPI = 400
manualMax = None

folders = ['.']
# folders = ['hotel_fringe_sgd1300','foxtrot_s_sgd1300']
# labels = ['fringe','none']
# folders = ['hotel_Umin03_sgd1300','foxtrot_s_sgd1300','hotel_Umin05_sgd1300']
# labels = ['0.03','0.04','0.05']
# folders = ['hotel_zhaoMelt_sgd1300','hotel_4x_sgd1300','foxtrot_s_sgd1300']
# labels = ['Zhao','4x','1x']
colors = ['xkcd:blue','xkcd:gray','xkcd:red']
resultFolder = '/results'
fileEnding = ""

if(os.path.isfile('input/plotHelperLocal.py')):
    sys.path.append('input')
    from plotHelperLocal import *
    print('Found experiment plotting settings')
elif(os.path.isfile('../plotHelper.py')):
    sys.path.append('../')
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')

plotDPI = 400
print('Plot DPI:',plotDPI)

figLabels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
fig = plt.figure(figsize=(15, 6),layout="tight")
# plt.figure(figsize=(16, 6))
ax1 = plt.subplot(2,5,(1,4))
ax2 = plt.subplot(2,5,(5,10))
ax3 = plt.subplot(2,5,(6,9))

axes = [ax1,ax2,ax3]
#Open figure before looping
## Plotting water production as function of depth
for i in range(len(folders)):
    folder = folders[i]


    y = mds.rdmds("results/YC")/1e3
    x = mds.rdmds("results/XC")/1e3
    z = np.squeeze(mds.rdmds("%s%s/RC" %(folder,resultFolder))) #midpoints of cells
    dz = np.load("%s/input/dz.npy" %(folder)) # we now save this upon grid generation, much better

    if(os.path.isfile('input/bathymetry.bin')):
        topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
        topo = topo.reshape(np.shape(x))
    else:
        topo = np.zeros(np.shape(x))

    ic_t = np.fromfile('input/T.init', dtype='>f8')
    ic_t = ic_t.reshape([np.shape(z)[0],np.shape(x)[0],np.shape(x)[1]])
    ic_s = np.fromfile('input/S.init', dtype='>f8')
    ic_s = ic_s.reshape([np.shape(z)[0],np.shape(x)[0],np.shape(x)[1]])

    ic_v = np.fromfile('input/NsBCv.bin', dtype='>f8')
    ic_v = ic_v.reshape([480,32,255]) #length of time vector, depth, length of x
    ic_v[ic_v ==0] = np.nan

    initDir = '/Users/psummers8/Documents/MITgcm/MITgcm/analysis/uploadFiles/Oscar_init'
    # initDir = 'input'
    openFrac = np.fromfile(f'{initDir}/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape([np.shape(z)[0],np.shape(x)[0],np.shape(x)[1]])
    openFrac[:,topo==0] = np.nan
    openFrac = 1- openFrac
    openFrac[openFrac == 0] = np.nan
cp = ax1.pcolormesh(x[0,:],y[:,0],topo,cmap='cmo.deep')
cv = ax1.contourf(x[0,:],y[:,0],np.transpose(ic_v[0,0,:, None] + 0*y[:,0]),np.linspace(-.03,.00,31),cmap='Blues_r',extend='neither',alpha=.5)
cc = ax1.pcolormesh(x[0,:],y[:,0],openFrac[0,:,:],vmin=.2,vmax=.8,cmap='gray_r',alpha=0.5)
ax1.set_aspect('equal')
ax1.set_xlabel('Along Fjord [km]')
ax1.set_ylabel('Across Fjord [km]')

data = mds.rdmds("results/dynDiag", 5760)
data[:,data[2,:,:,:] == 0] = np.nan

# cv = ax3.contourf(x[0,:],z,np.transpose(ic_v[0,0,:, None] + 0*z),np.linspace(-.03,.00,31),cmap='Blues_r',extend='neither',alpha=.5)
ct = ax3.contourf(x[0,:],z,np.nanmean(data[0,:,:,:],axis=1),tempRange,cmap=tempCmap)#,extend='both')
cc = ax3.pcolormesh(x[0,:],z,np.nanmean(openFrac[:,12:13,:],axis=1),vmin=.2,vmax=.8,cmap='gray_r',alpha=0.5)

ax3.set_xlabel('Along Fjord [km]')
ax3.set_ylabel('Depth [m]')
ax3.set_aspect(1/50)
cbaxes0 = fig.add_axes([0.05, 0.5, 0.2, 0.03])
cbaxes1 = fig.add_axes([0.3 , 0.5, 0.2, 0.03]) #left, bottom,width height
cbaxes2 = fig.add_axes([0.55, 0.5, 0.2, 0.03])

cbar = plt.colorbar(cc,cax=cbaxes0,orientation='horizontal',fraction=0.06,format='%.1f')
cbar.set_label('Ice Fraction [ ]')
cbar = plt.colorbar(ct,cax=cbaxes1,orientation='horizontal',fraction=0.06,format='%.1f')
cbar.set_label('Temperature [C]')
cbar = plt.colorbar(cv,cax=cbaxes2,orientation='horizontal',fraction=0.06,format='%.2f')
cbar.set_label('Across Fjord Speed [m/s]')

pressure = -1 * np.ones(np.shape(ic_s[:,0,0])) * 1020 * 9.81 * z /(1e4)
CT = gsw.CT_from_t(ic_s[:,0,0], ic_t[:,0,0], 0)
density = gsw.rho(ic_s[:,0,0], CT, 0) - 1025 #in-stu density less 1000

ax2.plot(ic_t[:,0,0],z,label='Temperature [C]',color='xkcd:red')
ax2.plot(ic_s[:,0,0] - 30,z,label='Salinty - 30 [PSU]',color='xkcd:green')
ax2.plot(density,z,label='Density - 1025 [kg/m^3]',color='xkcd:squash')
ax2.set_ylabel('Depth [m]')
ax2.grid(alpha = .25)
ax2.legend()

for label,ax in zip(figLabels,axes):
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')

plt.savefig('figs/initial_conditions%s.png' %fileEnding, format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()
   