from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean

os.system('rm -f figs/map*.png')
os.system('rm -f figs/autoMap*.gif')

depth = -100 #default value
try:
    with open('input/plotPoint.txt', 'r') as file:
        lines = file.readlines()
        depth = float(lines[1]) #reads the 1st line in the doc
        print('Depth read from file', depth)
except FileNotFoundError:
    print('plot point file does not exist, using default')
#Pick Depth you want to see


maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('results'):
    # print(file)
    if "dynDiag" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            sizeStep = int(words[1])
            startStep = int(words[1])

if(maxStep/sizeStep > 50):   #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print(startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile('input/data.iceberg')):
    isBerg = True
else:
    isBerg = False

y = mds.rdmds("results/YC")
x = mds.rdmds("results/XC")
z = mds.rdmds("results/RC")



if(os.path.isfile('input/topog.slope')):
    topo = np.fromfile('input/topog.slope', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

if(isBerg):
    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1


zSlice = np.argmin(np.abs(z[:,0,0]- depth))
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

# iceEdge = np.interp(z[zSlice,0,0],ice[0,:],x[0,:])

name = ["Temp", "Sal", "U", "W", "V"]
cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

if(maxStep/sizeStep > 50):   #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

if(isBerg):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx']
    name = ["Temp", "Sal", "U", "W", "V","BRGmltRt"]
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]"]
else:
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
    name = ["Temp", "Sal", "U", "W", "V"]
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

for k in range(len(name)):
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = np.linspace(-1.5, 1.5, 128)
            cm = "cmo.thermal"
        elif k == 1:
            lvl = np.linspace(32, 34, 128)
            cm = "cmo.haline"
        elif k == 2 or k == 4:
            lvl = np.linspace(-.5, .5, 128)
            cm = "cmo.balance"
        elif k == 3:
            lvl = np.linspace(-0.08, 0.08, 128)
            cm = "cmo.curl"
        elif k == 5:
            lvl = np.linspace(0, .5, 128)
            cm = "cmo.rain"
        if(k == 5):
            kk = 2
        else:
            kk = k
        cp = plt.contourf(
            np.squeeze(x),
            np.squeeze(y),
            np.squeeze(data[kk, zSlice, :, :]),
            lvl,
            extend="both",
            cmap=cm,
        )
        # plt.plot(iceEdge*np.ones(np.shape(x)),np.squeeze(y[:,0]),color='gray')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[kk, zSlice, :, :]),np.nanmax(data[kk, zSlice, :, :]),np.max(np.isnan(data[kk, zSlice, :, :]))))
        plt.ylabel('Across Fjord [m]')
        plt.title("%s depth %f at %i" % (name[k], z[zSlice,0,0] ,i))
        if(isBerg):
            cp2 = plt.contourf(np.squeeze(x),
                np.squeeze(y),
                np.squeeze(openFrac[zSlice, :, :]),
                [.1,.5,.9],
                extend="min",
                alpha=.2,
                cmap='cmo.gray')
            cbar2 = plt.colorbar(cp2)
            cbar2.set_label('Ocean Fraction')
        j = i/startStep
        str = "figs/map%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay 5 figs/map%s*.png -colors 256 -depth 256 figs/autoMap%s.gif' %(name[k], name[k]))

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
