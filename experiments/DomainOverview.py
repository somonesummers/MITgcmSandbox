from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
# import argparse

# parser = argparse.ArgumentParser(description='Plot dynamics at ySlice')
# parser.add_argument('yCrossSection', nargs='?', const=0.0, type=float,
#                     help='optional slice location [m]')
# args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 400
cleanPNGs = True
showQuiver = False
showZeros = True
showDensity = True
makeMovie = False   

# if(args.yCrossSection != None):
#     print('** Manual ySlice detected **')
#     yCrossSection = args.yCrossSection
# #Overwrite local settings here if desired
# # usePcolor = True

print('Plot DPI:',plotDPI)

labels = ["A","B","C","D","E","F"]
fig, axes = plt.subplots(3, 2, figsize=(14, 6), layout="constrained",width_ratios=[1, 1])
ax1 = axes[0,0]
ax2 = axes[1,0]
ax3 = axes[2,0]
ax4 = axes[0,1]
ax5 = axes[1,1]
ax6 = axes[2,1]

# Forced Flow
folder = 'rbcs_12_20/'
x = mds.rdmds("%sresults/XC" %folder)[0,:]
y = mds.rdmds("%sresults/YC" %folder)[:,0] - 1400
z = np.squeeze(mds.rdmds("%sresults/RC" %folder))
hFacC = mds.rdmds("%sresults/hFacC" %folder)

topo = np.fromfile('rbcs_12_20/input/bathymetry.bin', dtype='>f8')
topo = topo.reshape([len(y),len(x)])

T = np.fromfile('rbcs_12_20/input/EBCt.bin', dtype='>f8')
T = T.reshape([len(z),len(y)])
S = np.fromfile('rbcs_12_20/input/EBCs.bin', dtype='>f8')
S = S.reshape([len(z),len(y)])
U = np.fromfile('rbcs_12_20/input/EBCu.bin', dtype='>f8')
U = U.reshape([len(z),len(y)])

ax1.set_title('Forced Flow Domain')
cp = ax1.pcolormesh(
            x,
            y,
            topo,
            # extend="both",
            cmap='cmo.deep',
            )
ax1.set_aspect('equal')
# cbar = plt.colorbar(cp)
# cbar.set_label('Bathymetry [m]')
ax1.set_ylabel('Across Fjord [m]')
ax1.set_xlabel('Along Fjord [m]')

ax4.plot(T[:,0],z,label='Temp [˙C]')
ax4.plot(U[:,0]/.12,z,label='Along Fjord Velocity [u/U]')
ax4.plot(S[:,0]-34,z,label='Salt Anomonly [PSU - 34]')
ax4.legend()
ax4.set_ylabel('Depth [m]')
ax4.set_xlabel('')
ax4.grid(alpha=.5)

# Hughes 2024

folder = 'wiskey_b1/'
x = mds.rdmds("%sresults_0025/XC" %folder)[0,:]
y = mds.rdmds("%sresults_0025/YC" %folder)[:,0] - 2700
z = np.squeeze(mds.rdmds("%sresults_0025/RC" %folder))
hFacC = mds.rdmds("%sresults_0025/hFacC" %folder)

if(os.path.isfile('%sinput/bathymetry.bin'%folder)):
    topo = np.fromfile('%sinput/bathymetry.bin'%folder, dtype='>f8')
    topo = topo.reshape([len(y),len(x)])
else:
    topo = np.zeros(np.shape(x))

T = np.fromfile('%sinput/EBCt.bin' %folder, dtype='>f8')
T = T.reshape([len(z),len(y)])
S = np.fromfile('%sinput/EBCs.bin' %folder, dtype='>f8')
S = S.reshape([len(z),len(y)])
V = np.fromfile('%sinput/EBCv.bin' %folder, dtype='>f8')
V = V.reshape([len(z),len(y)])

ax2.set_title('Mélange Melt Domain')
cp = ax2.pcolormesh(
            x,
            y,
            topo,
            # extend="both",
            cmap='cmo.deep',
            )
ax2.set_aspect('equal')
# cbar = plt.colorbar(cp)
# cbar.set_label('Bathymetry [m]')
ax2.set_ylabel('Across Fjord [m]')
ax2.set_xlabel('Along Fjord [m]')

ax5.plot(T[:,0],z,label='Temp [˙C]')
ax5.plot(V[:,0],z,label='Coastal Current V [m/s]')
ax5.plot(S[:,0]-34,z,label='Salt Anomonly [PSU - 34]')
ax5.legend()
ax5.set_ylabel('Depth [m]')
ax5.set_xlabel('')
ax5.grid(alpha=.5)

# Fjord

ax3.set_title('Fjord Scale Domain')
folder = '../analysis/fjord_b1/'
x = mds.rdmds("%sresults/XC" %folder)[0,:]
y = mds.rdmds("%sresults/YC" %folder)[:,0] - 3400
z = np.squeeze(mds.rdmds("%sresults/RC" %folder))
hFacC = mds.rdmds("%sresults/hFacC" %folder)

if(os.path.isfile('%sinput/bathymetry.bin'%folder)):
    topo = np.fromfile('%sinput/bathymetry.bin'%folder, dtype='>f8')
    topo = topo.reshape([len(y),len(x)])
else:
    topo = np.zeros(np.shape(x))

nt = 25 #time forcing option
T = np.fromfile('%sinput/EBCt.bin' %folder, dtype='>f8')
T = T.reshape([nt,len(z),len(y)])
S = np.fromfile('%sinput/EBCs.bin' %folder, dtype='>f8')
S = S.reshape([nt,len(z),len(y)])
V = np.fromfile('%sinput/EBCv.bin' %folder, dtype='>f8')
V = V.reshape([nt,len(z),len(y)])

cp = ax3.pcolormesh(
            x,
            y,
            topo,
            # extend="both",
            cmap='cmo.deep',
            )
ax3.set_aspect('equal')
cbar = plt.colorbar(cp,orientation='horizontal')
cbar.set_label('Bathymetry [m]')
ax3.set_ylabel('Across Fjord [m]')
ax3.set_xlabel('Along Fjord [m]')

ax6.plot(T[0,:,0],z,label='Temp [˙C]')
ax6.plot(V[0,:,0],z,label='Coastal Current V [m/s]')
ax6.plot(S[0,:,0]-34,z,label='Salt Anomonly [PSU - 34]')
ax6.legend()
ax6.set_ylabel('Depth [m]')
ax6.set_xlabel('')
ax6.grid(alpha=.5)

str = "figs/modelDomainS.png"
for label,ax in zip(labels,axes.reshape(6,1)):
    ax=ax[0]
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')
plt.savefig(str, format='png', dpi=plotDPI)
plt.show()
plt.close()



        