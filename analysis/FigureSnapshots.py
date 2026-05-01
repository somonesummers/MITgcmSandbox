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
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = [645, 860],
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-nb','--nobergs', action='count', default=0,
                    help='optional to turn off bergs]')
args = parser.parse_args()

plt.rcParams.update({'font.size': 18})
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
fig = plt.figure(figsize=(18, 6),layout="tight")
# plt.figure(figsize=(16, 6))
ax1 = plt.subplot(2,1,1)
ax3 = plt.subplot(2,1,2)

axes = [ax1,ax3]
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

data = mds.rdmds("results/dynDiag", args.timeRange[0]*86400/15)
data[:,data[2,:,:,:] == 0] = np.nan
if(args.nobergs == 0):
    dataBRG = mds.rdmds("results/BRGFlx", args.timeRange[0]*86400/15)
    openFrac = dataBRG[5,:,:,:]
    openFrac[data[2,:,:,:] == 0] = np.nan
else:
    openFrac = np.ones_like(data[2,:,:,:])
openFrac[openFrac == 1] = np.nan

# cv = ax3.contourf(x[0,:],z,np.transpose(ic_v[0,0,:, None] + 0*z),np.linspace(-.03,.00,31),cmap='Blues_r',extend='neither',alpha=.5)
ct = ax1.contourf(x[0,:],z,np.nanmean(data[0,:,:,:],axis=1),tempRange,cmap=tempCmap)#,extend='both')
cc = ax1.pcolormesh(x[0,:],z,np.nanmean(1-openFrac[:,12:13,:],axis=1),vmin=.2,vmax=.8,cmap='gray_r',alpha=0.5)
zSkip = 2
ySkip = 10
u = np.squeeze(np.nanmean(data[2, ::zSkip, :, ::ySkip],axis=1))
w = np.squeeze(np.nanmean(data[3, ::zSkip, :, ::ySkip],axis=1))
qv = ax1.quiver(
    x[0,::ySkip],
    np.squeeze(z[::zSkip]),
    u,#/np.sqrt(u**2 + w**2 + 1e-12),
    w,#/np.sqrt(u**2 + w**2 + 1e-12),
    alpha=.4,
    angles='xy',
    pivot='tail',
    width = .002, #width of line
    scale = 3.0
    )
ax1.set_title('Winter')
# plt.quiverkey(qv,.75,1.30,.01,'1 cm/s')

cbar = plt.colorbar(ct,orientation='vertical',fraction=0.05,format='%.1f')
cbar.set_label('Temperature [C]')
ax1.set_xlabel('Along Fjord [km]')
ax1.set_ylabel('Depth [m]')
ax1.set_aspect(1/50)


data = mds.rdmds("results/dynDiag", args.timeRange[1]*86400/15)
data[:,data[2,:,:,:] == 0] = np.nan
if(args.nobergs == 0):
    dataBRG = mds.rdmds("results/BRGFlx", args.timeRange[1]*86400/15)
    openFrac = dataBRG[5,:,:,:]
    openFrac[data[2,:,:,:] == 0] = np.nan
else:
    openFrac = np.ones_like(data[2,:,:,:])
openFrac[openFrac == 1] = np.nan

# cv = ax3.contourf(x[0,:],z,np.transpose(ic_v[0,0,:, None] + 0*z),np.linspace(-.03,.00,31),cmap='Blues_r',extend='neither',alpha=.5)
ct = ax3.contourf(x[0,:],z,np.nanmean(data[0,:,:,:],axis=1),tempRange,cmap=tempCmap)#,extend='both')
cc = ax3.pcolormesh(x[0,:],z,np.nanmean(1-openFrac[:,12:13,:],axis=1),vmin=.05,vmax=.8,cmap='gray_r',alpha=0.5)
zSkip = 2
ySkip = 10
u = np.squeeze(np.nanmean(data[2, ::zSkip, :, ::ySkip],axis=1))
w = np.squeeze(np.nanmean(data[3, ::zSkip, :, ::ySkip],axis=1))
qv = ax3.quiver(
    x[0,::ySkip],
    np.squeeze(z[::zSkip]),
    u,#/np.sqrt(u**2 + w**2 + 1e-12),
    w,#/np.sqrt(u**2 + w**2 + 1e-12),
    alpha=.4,
    angles='xy',
    pivot='tail',
    width = .002, #width of line
    scale = 3.0
    )
plt.quiverkey(qv,1.05,1.25,.10,'10 cm/s')
ax3.set_title('Summer')
ax3.set_xlabel('Along Fjord [km]')
ax3.set_ylabel('Depth[m]')
ax3.set_aspect(1/50)

cbar = plt.colorbar(cc,orientation='vertical',fraction=0.05,format='%.1f')
cbar.set_label('Ice Fraction [ ]')
# cbar = plt.colorbar(cv,cax=cbaxes2,orientation='horizontal',fraction=0.06,format='%.2f')
# cbar.set_label('Across Fjord Speed [m/s]')


# for label,ax in zip(figLabels,axes):
#     ax.text(
#         ax.get_xlim()[0], ax.get_ylim()[1], label,
#         fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')

timeString = ''
if(args.timeRange != None):
    timeString = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'
plt.savefig(f'figs/snapshot{fileEnding}{timeString}.png', format='png',dpi=plotDPI)
if(args.silent == 0):
    plt.show()
plt.close()
   