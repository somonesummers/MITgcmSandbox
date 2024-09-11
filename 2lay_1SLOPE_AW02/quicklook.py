from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np

snapshotTimeStep = 43200
modelTimeStep = 20
maxStep = 864000
startStep = snapshotTimeStep/modelTimeStep
for i in np.arange(startStep, maxStep + 1, startStep*10):
    data = mds.rdmds("dynDiag", i)

    k = 2

    if k == 0:
        lvl = 16
        cm = "viridis"
        x = mds.rdmds("XC")
        z = mds.rdmds("RC")
    elif k == 1:
        lvl = np.linspace(34, 35, 20)
        cm = "viridis"
        x = mds.rdmds("XC")
        z = mds.rdmds("RC")
    elif k == 2:
        lvl = np.linspace(-0.1, 0.1, 32)
        cm = "RdBu"
        x = mds.rdmds("XG")
        z = mds.rdmds("RC")
    elif k == 3:
        lvl = np.linspace(-0.01, 0.01, 32)
        cm = "RdBu"
        x = mds.rdmds("XC")
        z = mds.rdmds("RC")

    cp = plt.contourf(
        np.squeeze(x),
        np.squeeze(z),
        np.squeeze(data[k, :, :, :]),
        lvl,
        extend="both",
        cmap=cm,
    )
    plt.colorbar(cp)
    if k == 0:
        plt.title("Temp at {}".format(i))
    elif k == 1:
        plt.title("Sal at {}".format(i))
    elif k == 2:
        plt.title("hor Vel at {}".format(i))
    elif k == 3:
        plt.title("Ver Vel at {}".format(i))
    j = i/startStep
    str = "quick%05i.png" % j
    # plt.savefig(str, format='png')
    # plt.close()
    plt.show()
    
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
