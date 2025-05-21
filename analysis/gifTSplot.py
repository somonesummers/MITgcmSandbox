from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Plot TS at one xCrossSection')
parser.add_argument('-x','--xCrossSection', nargs=1, type=float,default = None,
                    help='optional cross section location [m]')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last frame only, double to show plot(s)')
parser.add_argument('-n','--numFrames', nargs=1, type=int, default = [60],
                    help='optional specification of numFrames [default = 60]')
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

if((maxStep-startStep)/sizeStep > args.numFrames):   #if more than numFrames, downscale to be less than numFrames
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/args.numFrames)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

# os.system('rm -f figs/TSPlot*.png')
# os.system('rm -f figs/TSPlot*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")
 
xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
print('cross section is x =', x[0,xSlice],'index', xSlice)

# dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
# name = ["Temp", "Sal", "U", "W", "V"]
# units = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

nMix = 25
mixingT = np.zeros([2,nMix])
mixingS = np.zeros([2,nMix])
mixingT[1,:] = np.linspace(-10,10,nMix)
mixingS[1,:] = np.ones([1,nMix])*50

nMelt = 25
meltT = np.ones([2,nMelt])*-90
meltS = np.zeros([2,nMelt])
meltT[1,:] = np.ones([1,nMelt])*10
meltS[1,:] = np.linspace(-10,10,nMelt)+32

freezeS = [0,50]
freezeT = [0,-2.809]
freezeS100 = [0,50]
freezeT100 = [-.011,-2.919]

if(args.quick > 0):
    startStep = maxStep
    cleanPNGs = False

for i in np.arange(startStep, maxStep + 1, sizeStep):
    data = mds.rdmds("results/dynDiag", i)
    plt.figure()
    # data_old = mds.rdmds("results/dynDiag", 20*86400/dt)  #breaks if not a 20 day run, fix later
    # sc=plt.scatter(np.mean(data_old[1,:,1:-1,xSlice],1),np.mean(data_old[0,:,1:-1,xSlice],1),
    #                alpha=.5,s=25,color='black',edgecolor='none')
    for j in range(np.shape(y[1:-1,:])[0]):
        plt.scatter(data[1,:,j+1,xSlice],data[0,:,j+1,xSlice],c=np.squeeze(z),
                   alpha=.25,s=10,cmap='viridis')
    sc=plt.scatter(np.mean(data[1,:,1:-1,xSlice],1),np.mean(data[0,:,1:-1,xSlice],1),c=np.squeeze(z),
                   alpha=1.,s=25,cmap='viridis')
    plt.plot(mixingS,mixingT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
    plt.plot(meltS,meltT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
    plt.plot(freezeS,freezeT,linewidth=.5,color='red',alpha=.5,linestyle='--')
    plt.plot(freezeS100,freezeT100,linewidth=.5,color='red',alpha=.5,linestyle='--')
    cbar = plt.colorbar(sc)
    cbar.set_label('Depth [m]')
    ax = plt.gca()
    ax.set_xlim([np.min(saltRange), np.max(saltRange)])
    ax.set_ylim([np.min(tempRange), np.max(tempRange)])
    
    plt.xlabel('Salt [PSU]')
    plt.ylabel('Temperature [C]')
    plt.title("TS x = %i at %.02f days" % (x[0,xSlice], i/86400.0*dt))
    j = i/sizeStep + startStep
    

    str = "figs/tmpTSPlot%05i.png" % (j)
    
    plt.savefig(str, format='png', dpi=plotDPI)
    if(args.quick > 1):
        plt.show()
    plt.close()
if(args.quick == 0):
    if(args.xCrossSection != None):
        os.system('magick -delay %f figs/tmpTSPlot*.png -colors 256 -depth 256 figs/TSPlot%i.gif' %(500/((maxStep-startStep)/sizeStep),args.xCrossSection[0]))
    else:
        os.system('magick -delay %f figs/tmpTSPlot*.png -colors 256 -depth 256 figs/TSPlot.gif' %(500/((maxStep-startStep)/sizeStep)))

# #Clean up intermediate pngs
if(cleanPNGs):
    os.system('rm -f figs/tmpTSPlot*.png')
