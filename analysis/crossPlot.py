#!/usr/bin/env python
# coding: utf-8
from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import cmocean

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = 50
try:
    with open('input/plotPoint.txt', 'r') as file:
        lines = file.readlines()
        zDepth = float(lines[1]) #reads the 2nd line in the doc
        yCrossSection = float(lines[2]) #reads the 3rd line in the doc
        xCrossSection = float(lines[3]) #reads the 4th line in the doc
        print('cross sections read from file\nx:', xCrossSection,'\ny:',yCrossSection,'\nz:',zDepth)
except FileNotFoundError:
    print('plot point file does not exist, using default')

#Find Diagnostic file, iterate through them
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
if((maxStep-startStep)/sizeStep > 50):   #if more than 50 frames, downscale to be less than 50
    dwnScale = round(((maxStep-startStep)/sizeStep)/50)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale
print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile('input/data.iceberg')):
    isBerg = True
else:
    isBerg = False

#Clean up old gifs and pngs
os.system('rm -f figs/cross_*.png')
os.system('rm -f figs/autoCross_*.gif')

#Import grid
x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

#Import bathymetry
if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))
                
# Print actual cross section values
xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
ySlice = np.argmin(np.abs(y[:,0] - yCrossSection))
zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('cross section is x =', x[0,xSlice], 'index', xSlice)
print('cross section is y =', y[ySlice,0], 'index', ySlice)
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

# if(isBerg):
#     dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx']
#     name = ["Temp", "Sal", "U", "W", "V","BRGmltRt"]
#     cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]"]
# else:
dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
name = ["Temp", "Sal", "U", "W", "V"]
cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

#In KM for x,y for better axis labeling
x = x/1000
y = y/1000

#NOTE matplotlib x and y and MITgcm x,y are FLIPPED below. Be careful.
for k in range(len(name)):
    #print('k,',k)
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d',computed_zorder=False)
        # ax = fig.add_subplot(111)
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = np.linspace(-0.5, 3.0, 128)
            cm = "cmo.thermal"
        elif k == 1:
            lvl = np.linspace(32, 34, 128)
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

        #Long profile
        XX,ZZ = np.meshgrid(np.squeeze(x[ySlice,:]),np.squeeze(z))
        cp = ax.contourf(
            np.squeeze(data[kk, :, ySlice, :]),
            XX,
            ZZ,
            levels=lvl,
            extend="both",
            cmap=cm,
            zdir='x',offset=y[ySlice,0],zorder=2
        )

        #viewers top half profile
        YY,ZZ = np.meshgrid(np.squeeze(y[:,xSlice]),np.squeeze(z))    
        cp = ax.contourf(
            YY[:,0:ySlice+1],
            np.squeeze(data[kk, :, 0:ySlice+1, xSlice]),
            ZZ[:,0:ySlice+1],
            levels=lvl,
            extend="both",
            cmap=cm,
            zdir='y',offset=x[0,xSlice],zorder=3
        )
       
        #viewer's bottom width half profile
        YY,ZZ = np.meshgrid(np.squeeze(y[:,xSlice]),np.squeeze(z))    
        cp = ax.contourf(
            YY[:,ySlice:],
            np.squeeze(data[kk, :, ySlice:, xSlice]),
            ZZ[:,ySlice:],
            levels=lvl,
            extend="both",
            cmap=cm,
            zdir='y',offset=x[0,xSlice],zorder=0
        )
        
        #map view ghost
        #Front
        ghost = ax.contourf(
            np.squeeze(y[0:ySlice+1,0:xSlice+1]),
            np.squeeze(x[0:ySlice+1,0:xSlice+1]),
            np.ones(np.shape(x[0:ySlice+1,0:xSlice+1])),
            levels=[0, 1, 2],
            alpha=.1,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=4
        )
        #Back
        ghost = ax.contourf(
            np.squeeze(y[ySlice:,xSlice:]),
            np.squeeze(x[ySlice:,xSlice:]),
            np.ones(np.shape(x[ySlice:,xSlice:])),
            levels=[0, 1, 2],
            alpha=.1,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=-1
        )
        #Left
        ghost = ax.contourf(
            np.squeeze(y[ySlice:,0:xSlice]),
            np.squeeze(x[ySlice:,0:xSlice]),
            np.ones(np.shape(x[ySlice:,0:xSlice])),
            levels=[0, 1, 2],
            alpha=.1,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=1
        )
        #Right
        ghost = ax.contourf(
            np.squeeze(y[0:ySlice+1,xSlice:]),
            np.squeeze(x[0:ySlice+1,xSlice:]),
            np.ones(np.shape(x[0:ySlice+1:,xSlice:])),
            levels=[0, 1, 2],
            alpha=.1,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=2
        )
        #map view, projected onto bottom of frame
        cp = ax.contourf(
            np.squeeze(y),
            np.squeeze(x),
            np.squeeze(data[kk, zSlice, :, :]),
            levels=lvl,
            extend="both",
            cmap=cm,
            zdir='z',offset=z.min(),zorder=0
        )

        cbar = fig.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        
        #Topography along fjord
        ax.plot(x[ySlice,:xSlice],topo[ySlice,:xSlice],color='black',zdir='x',zs=y[ySlice,0],zorder=2)        
        ax.plot(x[ySlice,xSlice:],topo[ySlice,xSlice:],color='black',zdir='x',zs=y[ySlice,0],zorder=2)
        
        #Topography across fjord
        ax.plot(y[ySlice:,xSlice],topo[ySlice:,xSlice],color='black',zdir='y',zs=x[0,xSlice],zorder=0)
        ax.plot(y[:ySlice+1,xSlice],topo[:ySlice+1,xSlice],color='black',zdir='y',zs=x[0,xSlice],zorder=3)
        ax.invert_xaxis()
        plt.title("%s at x,y,z (%i,%i,%i) at %i" % (name[k], x[0,xSlice]*1000,y[ySlice,0]*1000,z[zSlice,0,0], i))
        ax.set_xlabel('Width [km]')
        ax.set_ylabel('Along [km]')
        ax.set_zlabel('Depth [m]')
        ax.axes.set_xlim3d(left=y.max(), right=y.min())
        ax.axes.set_ylim3d(bottom=x.min(), top=x.max()) 
        ax.axes.set_zlim3d(bottom=z.min(), top=z.max()) 

        ax.set_box_aspect([1,2,1])
        j = i/sizeStep + startStep

        str = "figs/cross_%s%05i.png" % (name[k],j)
        plt.savefig(str, format='png')
        plt.close()
        # plt.show()

    os.system('magick -delay %f figs/cross_%s*.png -colors 256 -depth 256 figs/autoCross_%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))




