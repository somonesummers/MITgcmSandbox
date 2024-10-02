from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os


maxStep = 21600
sizeStep = 2160
startStep = 2160

name = ["Temp", "Sal", "U", "W"]
for k in [0,1,2,3]:
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/dynDiag", i)
        if k == 0:
            lvl = np.linspace(-3, 1, 32)
            cm = "hot_r"
        elif k == 1:
            lvl = np.linspace(33, 35, 128)
            cm = "viridis"
        elif k == 2:
            lvl = np.linspace(-0.01, 0.01, 32)
            cm = "bwr"
        elif k == 3:
            lvl = np.linspace(-0.001, 0.001, 32)
            cm = "bwr"
        x = mds.rdmds("results/XC")
        z = mds.rdmds("results/RC")

        cp = plt.contourf(
            np.squeeze(x[0,:]),
            np.squeeze(z),
            np.squeeze(data[k, :, 5, :]),
            lvl,
            extend="both",
            cmap=cm,
        )
        plt.colorbar(cp)
        plt.title("%s at %i" % (name[k], i))
        j = i/startStep
        
        str = "figs/%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay 5 figs/%s*.png -colors 256 -depth 256 figs/auto_%s.gif' %(name[k], name[k]))

# BCT = np.fromfile("T.bound", dtype=">f8")
# plt.plot(BCT)
# plt.show()
#
# BCS = np.fromfile("S.bound", dtype=">f8")
# plt.plot(BCS)
# plt.show()
#
# BCU = np.fromfile("U.bound", dtype=">f8")
# plt.plot(BCU)
# plt.show()
#
#
