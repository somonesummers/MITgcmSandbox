from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean
import fileinput


depth = -100 #default value
try:
    with open('input/plotPoint.txt', 'r') as file:
        lines = file.readlines()
        depth = float(lines[1]) #reads the 2nd line in the doc
        print('Depth read from file', depth)
except FileNotFoundError:
    print('plot point file does not exist, using default')

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('results'):
    # print(file)
    if "dynDiag.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            startStep = int(words[1])
        if abs(int(words[1]) - startStep) < sizeStep and abs(int(words[1]) - startStep) > 0:
            sizeStep = abs(int(words[1]) - startStep)


if((maxStep-startStep)/sizeStep > 60):   #if more than 50 frames, downscale to be less than 50
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/60)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

#Clean up old gifs and pngs
os.system('rm -f figs/map*.png')
os.system('rm -f figs/autoMap*.gif')

y = mds.rdmds("results/YC")
x = mds.rdmds("results/XC")
z = mds.rdmds("results/RC")



if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
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
        if(isBerg and os.path.isfile('results/BRGFlx.%010i.001.001.data' % i)):
            localBergs = True
        else:
            localBergs = False
        if((not localBergs) and k==5):
            #fill melt image with 0s if bergs in run but not frame
            data = np.zeros(np.shape(mds.rdmds("results/%s"%(dynName[k-1]), i)))
        else:
            data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = np.linspace(-0.5, 3, 128)
            cm = "cmo.thermal"
        elif k == 1:
            lvl = np.linspace(32, 35, 128)
            cm = "cmo.haline"
        elif k == 2 or k == 4:
            lvl = np.linspace(-.5, .5, 127)
            cm = "cmo.balance"
        elif k == 3:
            lvl = np.linspace(-0.005, 0.005, 127)
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
        plt.title("%s depth %f at %.02f days" % (name[k], z[zSlice,0,0] ,i/86400.0*dt))
        if(localBergs):
            cp2 = plt.contourf(np.squeeze(x),
                np.squeeze(y),
                np.squeeze(openFrac[zSlice, :, :]),
                [.1,.5,.9],
                extend="min",
                alpha=.1,
                cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')
        j = i/sizeStep + startStep
        str = "figs/map%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay %f figs/map%s*.png -colors 256 -depth 256 figs/autoMap%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
os.system('rm -f figs/map*.png')
