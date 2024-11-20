from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput

# Plot settings from local helper
zDepth = -50
plotDPI = 150
cleanPNGs = True

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

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

hfactorName = ['hFacC', 'hFacS', 'hFacW']

zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

for k in range(len(hfactorName)):
    data = mds.rdmds("results/%s" % hfactorName[k])
    lvl = [0,.01, 0.04, 0.05, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 1.00]
    cm = "cmo.gray"
    plt.figure()
    cp = plt.pcolor(
        np.squeeze(x),
        np.squeeze(y),
        np.squeeze(data[zSlice, :, :]),
        cmap=cm,
    )
    cbar = plt.colorbar(cp)
    # cbar.set_label(units[k])
    plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[zSlice, :, :]),np.nanmax(data[zSlice, :, :]),np.max(np.isnan(data[zSlice, :, :]))))
    plt.ylabel('Depth [m]')
    plt.title("%s depth %f" % (hfactorName[k], z[zSlice,0,0]))
    
    str = "figs/%s.png" % (hfactorName[k])
    
    plt.savefig(str, format='png',dpi=plotDPI)
    plt.close()
