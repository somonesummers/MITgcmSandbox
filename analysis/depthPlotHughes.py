from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot Properties as function of depth at one xCrossSection')
parser.add_argument('xCrossSection', nargs='?', const=0.0, type=float,
                    help='optional cross section location [m]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
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

if(args.xCrossSection != None):
    print('** Manual xCrossSection detected **')
    xCrossSection = args.xCrossSection

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

os.system('rm -f figs/depth*.png')
# os.system('rm -f figs/depthPlotX*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

# xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
# print('cross section is x =', x[0,xSlice],'index', xSlice)

dynName = ['dynDiag', 'dynDiag', 'dynDiag']
# dynName = ['dynDiag', 'dynDiag', 'dynDiag']
name = ["U", "W", "V"]
units = ["[m/s]", "[m/s]", "[m/s]"]

iceStart = np.argmin(np.abs(x[0,:] - 13100)) 
iceEnd = np.argmin(np.abs(x[0,:] - 15400))

print('Average is is x =', x[0,iceStart],',',x[0,iceEnd],'index', iceStart,',',iceEnd)

hFacC = mds.rdmds("results/hFacC")

for k in range(len(name)):
    print('\t' + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0 or k == 2:
            lvl = [np.min(uRange),np.max(uRange)]
            cm = "xkcd:rose"
            plotData = np.squeeze(data[k+2,:,:,:]) * 1 
        elif k == 1:
            lvl = [np.min(wRange),np.max(wRange)]
            cm = "xkcd:violet"
            plotData = np.squeeze(data[k+2,:,:,:]) * 1
        plt.figure()
        for ii in range(iceEnd - iceStart):
            for j in range(np.shape(y[1:-1,:])[0]):
                plt.plot(plotData[:,j+1,iceStart + ii],z,linewidth=.5,alpha=.1,color=cm)

        plt.plot(np.mean(plotData[:,1:-1,iceStart:iceEnd],axis=(1,2)),z,linewidth=2,color=cm)
        ax = plt.gca()
        ax.set_xlim(lvl)
        plt.plot(np.cos(z * np.pi /600)* 0.12, z,linewidth=1,color='gray',linestyle='--')
        plt.xlabel(name[k] + " " + units[k] + ' %.3f %.3f nan: %i' %(np.nanmin(plotData[:,1:-1]),np.nanmax(plotData[:,1:-1]),np.max(np.isnan(plotData[:,1:-1]))))
        plt.ylabel('Depth [m]')
        plt.title("%s over melange at %.02f days" % (name[k], i/86400.0*dt))
        if (k == 0 ):
            plt.xlim([0, .25])
        plt.ylim([-300, 0])
        j = i/sizeStep + startStep
        
        str = "figs/depth%s%05i.png" % (name[k],j)
        plt.savefig(str, format='png',dpi=plotDPI)
        plt.close()
        plt.show()
    if(args.xCrossSection != None):
        os.system('magick -delay %f figs/depth%s*.png -colors 256 -depth 256 figs/depth%i%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], args.xCrossSection, name[k]))
    else:
        os.system('magick -delay %f figs/depth%s*.png -colors 256 -depth 256 figs/depth%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/depth*.png')
