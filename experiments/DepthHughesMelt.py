from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
cleanPNGs = True
hFacWeighted = True

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

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

#Find time steps to take
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

os.system('rm -f figs/depthPlotX*.png')
# os.system('rm -f figs/depthPlotX*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
print('cross section is x =', x[0,xSlice],'index', xSlice)

dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag','dynDiag']
name = ["Temp", "Sal", "U", "W", "V", "N2"]
units = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[rad^2 s^-2 ]"]

iceEnd = 16
print('Mélange ends at index',iceEnd)

if(hFacWeighted):
    hFacC = mds.rdmds("results/hFacC") #must weight by ocean fraction
    hFacC[hFacC==0] = 1 #this makes the averaging happy, and we no 'dry' cells from bergs alone
for k in range(len(name)):
    print('\t' + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = [np.min(tempRange),np.max(tempRange)]
            cm = "xkcd:raspberry"
            plotData = np.squeeze(np.average(data[k,:,:,1:iceEnd],weights=hFacC[:,:,1:iceEnd],axis=2)) 
        elif k == 1:
            lvl = [np.min(saltRange),np.max(saltRange)]
            cm = "xkcd:green"
            plotData = np.squeeze(np.average(data[k,:,:,1:iceEnd],weights=hFacC[:,:,1:iceEnd],axis=2)) 
        elif k == 2 or k == 4:
            lvl = [np.min(uRange),np.max(uRange)]  
            cm = "xkcd:rose"
            plotData = np.squeeze(np.average(data[k,:,:,1:iceEnd],weights=hFacC[:,:,1:iceEnd],axis=2)) 
        elif k == 3:
            lvl = [np.min(wRange),np.max(wRange)]
            cm = "xkcd:violet"
            plotData = np.squeeze(np.average(data[k,:,:,1:iceEnd],weights=hFacC[:,:,1:iceEnd],axis=2))  
        elif k == 5:
            lvl = [0,1e-3]
            cm = "xkcd:olive green"
            salt = np.squeeze(data[1,:,:,xSlice])
            if(i == startStep): #only calc pressure once
                pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
            CT = gsw.CT_from_t(salt, data[0,:,:,xSlice], pressure)
            tmp = gsw.Nsquared(salt, CT, pressure) 
            plotData = tmp[0]
            if(i == startStep):
                z = -1 * tmp[1] * 10e3 / (1020 * 9.81) # defined at midpoints, so new z here on first i loop
        plt.figure()
        for j in range(np.shape(y[1:-1,:])[0]):
            plt.plot(plotData[:,j+1],z,linewidth=.5,alpha=.5,color=cm)
        plt.plot(np.mean(plotData[:,1:-1],1),z,linewidth=2,color=cm)
        ax = plt.gca()
        ax.set_xlim(lvl)

        plt.xlabel(name[k] + " " + units[k] + ' %.3e %.3e nan: %i' %(np.nanmin(plotData[:,1:-1]),np.nanmax(plotData[:,1:-1]),np.max(np.isnan(plotData[:,1:-1]))))
        plt.ylabel('Depth [m]')
        plt.title("%s x = %i at %.02f days" % (name[k], x[0,xSlice], i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/depthAvgX%s%05i.png" % (name[k],j)
        plt.grid(alpha=.5)
        plt.savefig(str, format='png',dpi=plotDPI)
        plt.close()
        # plt.show()

    os.system('magick -delay %f figs/depthAvgX%s*.png -colors 256 -depth 256 figs/depthAvgX%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/depthAvgX*.png')
