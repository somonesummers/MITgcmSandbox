from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean

#Pick cross section
crossSection = 10000

maxStep = 0
sizeStep = 1e10
startStep = 1e10

#Decide if iceBerg data files exist
if(os.path.isfile('input/data.iceberg')):
    isBerg = True
else:
    isBerg = False

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

if(maxStep/sizeStep > 50):  #if more than 50 frames, downscale to be less than 50
    dwnScale = round((maxStep/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print(startStep,sizeStep,maxStep)

os.system('rm -f figs/sideX*.png')
os.system('rm -f figs/autosideX*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

topo = np.fromfile('input/topog.slope', dtype='>f8')
ice = np.fromfile('input/icetopo.exp1', dtype='>f8')
topo = topo.reshape(np.shape(x))
ice = ice.reshape(np.shape(x))

if(isBerg):
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1


xSlice = np.argmin(np.abs(x[0,:] - crossSection))
print('cross section is x =', x[0,xSlice])


name = ["Temp", "Sal", "U", "W", "V"]
cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

for k in [0,1,2,3,4]:
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/dynDiag", i)
        if k == 0:
            lvl = np.linspace(-1.5, 1.5, 128)
            cm = "cmo.thermal"
        elif k == 1:
            lvl = np.linspace(32, 34, 128)
            cm = "cmo.haline"
        elif k == 2:
            lvl = np.linspace(-0.5, 0.5, 128)
            cm = "cmo.balance"
        elif k == 3:
            lvl = np.linspace(-0.02, 0.02, 128)
            cm = "cmo.curl"
        elif k == 4:
            lvl = np.linspace(-0.5, 0.5, 128)
            cm = "cmo.balance"

        cp = plt.contourf(
            np.squeeze(y[:,xSlice]),
            np.squeeze(z),
            np.squeeze(data[k, :, :, xSlice]),
            lvl,
            extend="both",
            cmap=cm,
        )

        plt.plot(y[:,xSlice],topo[:,xSlice],color='black')
        plt.plot(y[:,xSlice],ice[:,xSlice],color='gray')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        if(isBerg):
            cp2 = plt.contourf(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(openFrac[:, :, xSlice]),
                [.1,.5,.9],
                extend="min",
                alpha=.2,
                cmap='cmo.gray')
            cbar2 = plt.colorbar(cp2)
            cbar2.set_label('Ocean Fraction')
        plt.xlabel('Across Fjord [m]')
        plt.ylabel('Depth [m]')
        plt.title("%s y = %i at %i" % (name[k], x[0,xSlice], i))
        j = i/startStep
        
        str = "figs/sideX%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png')
        plt.close()
        #plt.show()

    os.system('magick -delay 5 figs/sideX%s*.png -colors 256 -depth 256 figs/autosideX_%s.gif' %(name[k], name[k]))

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
