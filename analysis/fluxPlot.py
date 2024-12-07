from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Plot fluxes as function of depth at CrossSections')
parser.add_argument('xCrossSection', nargs='?', const=0.0, type=float,
                    help='optional x cross section location [m]')

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

os.system('rm -f figs/fluxPlot*.png')
os.system('rm -f figs/fluxPlot*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

# Print actual cross section values
xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
ySlice = np.argmin(np.abs(y[:,0] - yCrossSection))
zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('cross section is x =', x[0,xSlice], 'index', xSlice)
print('cross section is y =', y[ySlice,0], 'index', ySlice)
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

dynName = ['fluxMassDiag', 'fluxMassDiag', 'fluxMassDiag', 'fluxMassDiag','fluxMassDiag','fluxMassDiag']
name = ['UTHMASS','USLTMASS','VTHMASS','VSLTMASS','WTHMASS','WSLTMASS']
units = ["[˙C m/s]", "[PSU m/s]", "[˙C m/s]", "[PSU m/s]", "[˙C m/s]", "[PSU m/s]"]

for k in range(6):
    print('\t' + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if k == 0:
            lvl = [-.05,.05]
            cm = "xkcd:raspberry"
        elif k == 1:
            lvl = [-2,2]
            cm = "xkcd:green"
        elif k == 2:
            lvl = [-.05,.05]
            cm = "xkcd:raspberry"
        elif k == 3:
            lvl = [-2,2]
            cm = "xkcd:green"
        elif k == 4:
            lvl = [-.05,.05]
            cm = "xkcd:raspberry"
        elif k == 5:
            lvl = [-2,2]
            cm = "xkcd:green"
        plt.figure()
        if(k == 0 or k ==1): # U values
            dataPlot = data[k,:,1:-1,xSlice]
            yPlot = np.squeeze(z)
            yStr = 'Depth [m]'
            avgDirection = 1
            dirStr = 'x='
            sliceVal = x[0,xSlice]
        elif(k==2 or k==3): # V fluxes
            dataPlot = data[k,:,ySlice,1:]
            yPlot = np.squeeze(z)
            yStr = 'Depth [m]'
            avgDirection = 1
            dirStr = 'y='
            sliceVal = y[ySlice,0]
        elif(k==4 or k==5): # W fluxes
            dataPlot = data[k,zSlice,1:-1,:]
            yPlot = np.squeeze(x[0,:])
            yStr = 'Along Fjord [m]'
            avgDirection = 0
            dirStr = 'z='
            sliceVal = z[zSlice,0,0]
        for j in range(np.shape(dataPlot)[avgDirection]):
            if(avgDirection == 1):
                plt.plot(dataPlot[:,j],yPlot,linewidth=.5,alpha=.5,color=cm)
            else:
                plt.plot(dataPlot[j,:],yPlot,linewidth=.5,alpha=.5,color=cm)
        plt.plot(np.mean(dataPlot,avgDirection),yPlot,linewidth=1,color=cm)
        fullMean = np.nanmean(dataPlot)
        plt.plot([fullMean,fullMean],[yPlot[0],yPlot[-1]],linewidth=1,color='black',linestyle='--')
        ax = plt.gca()
        ax.set_xlim(lvl)

        plt.xlabel(name[k] + " " + units[k] + ' %.3f %.3f nan: %i' %(np.nanmin(dataPlot),np.nanmax(dataPlot),np.max(np.isnan(dataPlot))))
        plt.ylabel(yStr)
        plt.title("%s %s %i at %.02f days" % (name[k], dirStr, sliceVal, i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/fluxPlot%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png',dpi=plotDPI)
        plt.close()
        plt.show()
    if(args.xCrossSection != None):
        os.system('magick -delay %f figs/fluxPlot%s*.png -colors 256 -depth 256 figs/fluxPlot%i%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], args.xCrossSection, name[k]))
    else:
        os.system('magick -delay %f figs/fluxPlot%s*.png -colors 256 -depth 256 figs/fluxPlot%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/fluxPlot*.png')
