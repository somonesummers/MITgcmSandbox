from MITgcmutils import mds
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import sys
import cmocean
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -10
plotDPI = 300
cleanPNGs = True

def red_yellow_white_cyan_blue():
    N = 256
    cols = ['#000055', '#0000f5', '#008cff',
            # '#7affff', '#ffffff', '#ffff83',
            '#9ee6e6', '#ececec', '#eaea99',
            '#ff9d0c', '#ff0500', '#5f0000']
    return LinearSegmentedColormap.from_list('custom', cols, N)


folders = ['wiskey_b1','wiskey_b0']
labels = ['Iceberg Drag', 'No Drag']
# colors = ['xkcd:blue','xkcd:green']
resultFolder = '/results_0025'
fileEnding = ''
maxStep = 40320 #7 days @ dt 15
# maxStep = 24192 #7 days @ dt 25
# maxStep = 5976 #7 days @ dt 100
# maxStep = 43128 #50 days


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
print('Plot DPI:',plotDPI,'; clean PNGs?',cleanPNGs)


zSlice = 2
x_zSlice = 3

ySection = -2100
# ySlice = 2
# x_ySlice = 400
# Hughes plotting

dataset = xr.open_dataset('~/Documents/MITgcm/MITgcm/Sandbox_Berg/10366091/fields.deltaS_3.gamma_multiple_4.T0_2.nc')
# downscaled_dataset = dataset.coarsen(X=40, Y=40,boundary='pad').mean()
# dataset = downscaled_dataset.compute()
x_t = np.array(dataset['theta'])
x_u = np.array(dataset['U'])
# x_v = dataset['V']
x_x = dataset['X']+8000
x_y = dataset['Y']
x_z = dataset['Z']
## U vel  now ================================================

# fig = plt.figure(figsize=(12, 6))
# ax1 = plt.subplot(2,3,1)
# ax2 = plt.subplot(2,3,2)
# ax3 = plt.subplot(2,3,3)
# ax4 = plt.subplot(2,3,4)
# ax5 = plt.subplot(2,3,5)
# ax6 = plt.subplot(2,3,6)
# axes = [ax1, ax2, ax3, ax4, ax5, ax6]
labels = ["A","B","C","D","E","F"]
fig, axes = plt.subplots(3, 2, figsize=(14, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[2,0]
ax4 = axes[0,1]
ax5 = axes[1,1]
ax6 = axes[2,1]

lvl = np.linspace(-.07, .07, 127)
# plt.subplot(131)

x_u[x_u == 0] = np.nan
cmapTemp = red_yellow_white_cyan_blue()
cmapTemp.set_bad(color='gray')
cp = ax1.pcolormesh(
            x_x,
            x_y,
            x_u[4,x_zSlice,:,:],
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            #extend="both",
            cmap=cmapTemp,
            )
# cbar = plt.colorbar(cp)
# cbar.set_label('[C]')
ax1.set_title('Hughes 2024')
ax1.set_aspect('equal')
ax1.set_xlim([0, 25000])
# ax1.set_xlabel('Along Fjord [m]')
ax1.set_xticks([])
ax1.set_ylabel('Across Fjord [m]')

# plt.show()


#Blocking Plotting
# plt.subplot(132)

y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[0],resultFolder))[:,0]-2550-150)
x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[0],resultFolder))[0,:])
z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[0],resultFolder)))

print(" Hughes Z: %i Summers Z: %s" %(x_z[x_zSlice],z[zSlice]))
data = mds.rdmds("%s%s/dynDiag"%(folders[0],resultFolder),maxStep)
cp = ax2.pcolormesh(
            x[:],
            y[1:-1],
            data[2,zSlice,1:-1,:],
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            #extend="both",
            cmap=red_yellow_white_cyan_blue(),
            )
# cbar = plt.colorbar(cp)
# cbar.set_label('[C]')
ax2.set_title('Blocking + Drag')
ax2.set_aspect('equal')
ax2.set_xlim([0, 25000])
# ax2.set_xlabel('Along Fjord [m]')
ax2.set_xticks([])
ax2.set_ylabel('Across Fjord [m]')
# plt.ylabel('Across Fjord [m]')
# plt.show()

#Non Blocking Plotting
# plt.subplot(133)
# maxStep = 5976
y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[1],resultFolder))[:,0]-2550-150)
x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[1],resultFolder))[0,:])
z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[1],resultFolder)))

data = mds.rdmds("%s%s/dynDiag"%(folders[1],resultFolder),maxStep)
cp = ax3.pcolormesh(
            x[:],
            y[1:-1],
            data[2,zSlice,1:-1,:],
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            #extend="both",
            cmap=red_yellow_white_cyan_blue(),
            )
cbar = plt.colorbar(cp,ticks=[-.06,-.04,-.02,0,.02,.04,.06],orientation="horizontal")
cbar.set_label('[m/s]')
ax3.set_title('No Blocking + No Drag')
ax3.set_aspect('equal')
ax3.set_xlim([0, 25000])
ax3.set_xlabel('Along Fjord [m]')
ax3.set_ylabel('Across Fjord [m]')
# plt.show()

# plt.savefig('figs/mapU%s.png' % fileEnding, format='png',dpi=plotDPI)
# plt.close()

## Temp  now ================================================

# fig, axs = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")

lvl = np.linspace(0, 2, 128)
# plt.subplot(131)
x_t[x_t == 0] = np.nan

cmapTemp = mpl.colormaps.get_cmap('cmo.thermal')
cmapTemp.set_bad(color='gray')
cp = ax4.pcolormesh(
            x_x,
            x_y,
            x_t[4,x_zSlice,:,:],
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            # extend="both",
            cmap=cmapTemp,
            )
# cbar = plt.colorbar(cp)
# cbar.set_label('[C]')
ax4.set_title('Hughes 2024')
ax4.set_aspect('equal')
ax4.set_xlim([0, 25000])

# ax4.set_xlabel('Along Fjord [m]')
# ax4.set_ylabel('Across Fjord [m]')
ax4.set_xticks([])
ax4.set_yticks([])
# plt.show()


#Blocking Plotting
# plt.subplot(132)
# maxStep = 5976 #7 days
# maxStep = 43128 #50 days
y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[0],resultFolder))[:,0]-2500)
x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[0],resultFolder))[0,:])
z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[0],resultFolder)))

data = mds.rdmds("%s%s/dynDiag"%(folders[0],resultFolder),maxStep)
cp = ax5.pcolormesh(
            x[:],
            y[1:-1],
            data[0,zSlice,1:-1,:],
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            # extend="both",
            cmap='cmo.thermal',
            )
# cbar = plt.colorbar(cp)
# cbar.set_label('[C]')
ax5.set_title('Blocking + Drag')
ax5.set_aspect('equal')
ax5.set_xlim([0, 25000])
# ax5.set_xlabel('Along Fjord [m]')
ax5.set_yticks([])
ax5.set_xticks([])
# plt.ylabel('Across Fjord [m]')
# plt.show()

#Non Blocking Plotting
# plt.subplot(133)
# maxStep = 5976
y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[1],resultFolder))[:,0]-2500)
x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[1],resultFolder))[0,:])
z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[1],resultFolder)))

data = mds.rdmds("%s%s/dynDiag"%(folders[1],resultFolder),maxStep)
cp = ax6.pcolormesh(
            x[:],
            y[1:-1],
            data[0,zSlice,1:-1,:],
            vmin=np.min(lvl),
            vmax=np.max(lvl),
            # extend="both",
            cmap='cmo.thermal',
            )
cbar = plt.colorbar(cp,ticks=[0,.5,1,1.5,2],orientation="horizontal")
cbar.set_label('[C]')
ax6.set_title('No Blocking + No Drag')
ax6.set_aspect('equal')
ax6.set_xlim([0, 25000])
ax6.set_xlabel('Along Fjord [m]')
ax6.set_yticks([])
# plt.ylabel('Across Fjord [m]')
# plt.show()

# plt.tight_layout()

for label,ax in zip(labels,axes.reshape(6,1)):
    ax=ax[0]
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')
plt.savefig('figs/mapHughes%s.png' %fileEnding, format='png',dpi=plotDPI)
plt.show()

plt.close()

if(False):
    ## Side view speed  now ================================================
    # fig = plt.figure(figsize=(12, 6))
    # ax1 = plt.subplot(2,3,1)
    # ax2 = plt.subplot(2,3,2)
    # ax3 = plt.subplot(2,3,3)
    # ax4 = plt.subplot(2,3,4)
    # ax5 = plt.subplot(2,3,5)
    # ax6 = plt.subplot(2,3,6)
    # axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    labels = ["A","B","C","D","E","F"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), layout="constrained")
    # print(axes)
    ax1 = axes[0,0]
    ax2 = axes[0,1]
    ax3 = axes[0,2]
    ax4 = axes[1,0]
    ax5 = axes[1,1]
    ax6 = axes[1,2]

    adjustedDim = np.array(x_y[:])
    x_ySlice = np.argmin(np.abs(adjustedDim - ySection))

    # fig, axs = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")

    lvl = np.linspace(-.07, .07, 127)
    # plt.subplot(131)
    cp = ax1.contourf(
                x_x,
                x_z,
                x_u[4,:,x_ySlice,:],
                lvl,
                extend="both",
                cmap=red_yellow_white_cyan_blue(),
                )
    # cc = plt.contour(
    #     x_x,
    #     x_z,
    #     x_u[4,:,x_ySlice,:],
    #     [0],
    #     colors='gray',
    #     linewidths=0.5
    #     )
    # plt.clabel(cc, inline=3, fontsize=8)
    # cbar = plt.colorbar(cp)
    # cbar.set_label('[C]')
    ax1.set_title('Hughes 2024')
    ax1.set_xlim([-7750, 25000])
    ax1.set_ylim([-300, 0])
    ax1.set_xlabel('Along Fjord [m]')
    ax1.set_ylabel('Depth [m]')
    # plt.show()


    ySlice = np.argmin(np.abs(y - 3000 - ySection))
    print(" Hughes Y: %i Summers Y: %s" %(x_y[x_ySlice],y[ySlice]-3000))
    #Blocking Plotting
    # plt.subplot(132)
    # maxStep = 5976 #7 days
    # maxStep = 43128 #50 days
    y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[0],resultFolder))[:,0])
    x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[0],resultFolder))[0,:])
    z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[0],resultFolder)))

    data = mds.rdmds("%s%s/dynDiag"%(folders[0],resultFolder),maxStep)
    cp = ax2.contourf(
                x,
                z,
                # np.nanmean(data[2,:,1:-1,:],axis=1),
                data[2,:,ySlice,:],
                lvl,
                extend="both",
                cmap=red_yellow_white_cyan_blue(),
                )
    # cc = plt.contour(
    #     x[:],
    #     z,
    #     # np.mean(data[0,:,1:-1,:],axis=1),
    #     data[2,:,ySlice,:],
    #     [0],
    #     colors='gray',
    #     linewidths=0.5
    #     )
    # plt.clabel(cc, inline=3, fontsize=8)
    # cbar = plt.colorbar(cp)
    # cbar.set_label('[C]')
    ax2.set_title('Iceberg Drag')
    ax2.set_xlim([-7750, 25000])
    ax2.set_ylim([-300, 0])
    ax2.set_xlabel('Along Fjord [m]')
    ax2.set_yticks([])
    # plt.ylabel('Across Fjord [m]')
    # plt.show()

    #Non Blocking Plotting
    # plt.subplot(133)
    # maxStep = 5976
    y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[1],resultFolder))[:,0])
    x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[1],resultFolder))[0,:])
    z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[1],resultFolder)))

    data = mds.rdmds("%s%s/dynDiag"%(folders[1],resultFolder),maxStep)
    cp = ax3.contourf(
                x[:],
                z,
                # np.mean(data[2,:,1:-1,:],axis=1),
                data[2,:,ySlice,:],
                lvl,
                extend="both",
                cmap=red_yellow_white_cyan_blue(),
                )
    # cc = plt.contour(
    #     x[:],
    #     z,
    #     # np.mean(data[0,:,1:-1,:],axis=1),
    #     data[2,:,ySlice,:],
    #     [0],
    #     colors='gray',
    #     linewidths=0.5
    #     )
    # plt.clabel(cc, inline=3, fontsize=8)
    cbar = plt.colorbar(cp,ticks=[-.07,-.06,-.04,-.02,0,.02,.04,.06,.07])
    cbar.set_label('[m/s]')
    ax3.set_title('No Drag')
    ax3.set_xlim([-7750, 25000])
    ax3.set_ylim([-300, 0])
    ax3.set_xlabel('Along Fjord [m]')
    ax3.set_yticks([])
    # plt.ylabel('Across Fjord [m]')
    # plt.show()
    # plt.savefig('figs/mapUSide%s.png' %fileEnding, format='png',dpi=plotDPI)
    # plt.close()

    ## Side view Temp  now ================================================

    # fig, axs = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")

    lvl = np.linspace(0, 2, 128)
    # plt.subplot(131)
    cp = ax4.contourf(
                x_x,
                x_z,
                # np.nanmean(x_t[4,:,:,:],axis=1),
                x_t[4,:,x_ySlice,:],
                lvl,
                extend="both",
                cmap='cmo.thermal',
                )
    # cbar = plt.colorbar(cp)
    # cbar.set_label('[C]')
    ax4.set_title('Hughes')
    ax4.set_xlim([-7750, 25000])
    ax4.set_ylim([-300, 0])
    ax4.set_xlabel('Along Fjord [m]')
    ax4.set_ylabel('Depth [m]')
    # plt.show()


    #Blocking Plotting
    # plt.subplot(132)
    # maxStep = 5976 #7 days
    # maxStep = 43128 #50 days
    y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[0],resultFolder))[:,0])
    x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[0],resultFolder))[0,:])
    z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[0],resultFolder)))

    data = mds.rdmds("%s%s/dynDiag"%(folders[0],resultFolder),maxStep)
    cp = ax5.contourf(
                x,
                z,
                # np.nanmean(data[0,:,1:-1,:],axis=1),
                data[0,:,ySlice,:],
                lvl,
                extend="both",
                cmap='cmo.thermal',
                )
    # cbar = plt.colorbar(cp)
    # cbar.set_label('[C]')
    ax5.set_title('Iceberg Drag')
    ax5.set_xlim([-7750, 25000])
    ax5.set_ylim([-300, 0])
    ax5.set_xlabel('Along Fjord [m]')
    ax5.set_yticks([])
    # plt.ylabel('Across Fjord [m]')
    # plt.show()

    #Non Blocking Plotting
    # plt.subplot(133)
    # maxStep = 5976
    y = np.squeeze(mds.rdmds("%s%s/YC"%(folders[1],resultFolder))[:,0])
    x = np.squeeze(mds.rdmds("%s%s/XC"%(folders[1],resultFolder))[0,:])
    z = np.squeeze(mds.rdmds("%s%s/RC"%(folders[1],resultFolder)))



    data = mds.rdmds("%s%s/dynDiag"%(folders[1],resultFolder),maxStep)
    cp = ax6.contourf(
                x[:],
                z,
                # np.mean(data[0,:,1:-1,:],axis=1),
                data[0,:,ySlice,:],
                lvl,
                extend="both",
                cmap='cmo.thermal',
                )
    cbar = plt.colorbar(cp,ticks=[0,.5,1,1.5,2])
    cbar.set_label('[C]')
    ax6.set_title('No Drag')
    ax6.set_xlim([-7750, 25000])
    ax6.set_ylim([-300, 0])
    ax6.set_xlabel('Along Fjord [m]')
    ax6.set_yticks([])
    # plt.ylabel('Across Fjord [m]')
    # plt.show()

    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # cbar = plt.colorbar(cp,ticks=[0,.5,1,1.5,2],cax=cbar_ax)
    # cbar.set_label('[C]')

    # plt.tight_layout()
    for label,ax in zip(labels,axes.reshape(6,1)):
        ax=ax[0]
        ax.text(
            ax.get_xlim()[0], ax.get_ylim()[1], label,
            fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')

    plt.savefig('figs/mapSide%s.png' %fileEnding, format='png',dpi=plotDPI)
    plt.show()
    plt.close()

