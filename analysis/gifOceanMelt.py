from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot mélange melt rates from ocean only')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last frame only, double to show plot(s)')
parser.add_argument('-xl','--xlimit', nargs='?', type=int, default = 40000,
                    help='optional max x [default = entire frame]')
parser.add_argument('-dpi','--dpi', nargs='?', default=175,type=int,
                    help='dpi to print [default = 175]')
args = parser.parse_args()


plotDPI = args.dpi
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

#Overwrite local settings here if desired
# usePcolor = True

print('Plot DPI:',plotDPI)

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
    if "BRGFlx.0" in file:
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

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))
dx = 2*x[0,0]
dy = 2*y[0,0]


if(args.quick > 0):
    startStep = maxStep
    cleanPNGs = False

shiftBack = 2 #we don't plot the glacier front, GLACIOME smooths this out, very high
for i in np.arange(startStep, maxStep + 1, sizeStep):
    data = mds.rdmds("results/BRGFlx", i)
    lambdaHelper = 1-data[5,0,:,shiftBack:] #only surface filled-ness
    lambdaHelper[lambdaHelper == 0] = 0.01 #avoid divide by 0
    b_mitgcm = 1*np.nanmean(np.nansum(data[0,:,:,shiftBack:],axis=0)/lambdaHelper,axis=0)/(dx*dy)*(24*3600)
    
    b_mitgcmFit = b_mitgcm[b_mitgcm != 0]
    xFit = x[0,shiftBack:]
    g = np.ones([len(b_mitgcmFit),2])
    g[:,1] = xFit[b_mitgcm != 0]
    m = np.linalg.pinv(g) @ b_mitgcmFit

    plt.plot(x[0,shiftBack:],b_mitgcm,color='red',marker='o')
    plt.plot(x[0,shiftBack:],m[0] + m[1]*x[0,shiftBack:],color='red',linestyle='--',alpha=.5)
    plt.title("Mélange melt rate at %.02f days" % (i/86400.0*dt))
    j = i/sizeStep
    plt.xlabel('Distance along fjord [m]')
    plt.ylabel('Avg Mélange Melt Rate [m/day]')
    str = "figs/melangeMelt%05i.png" % (j)
    if(args.xlimit != None):
        plt.xlim([0,args.xlimit])        
    plt.ylim([0,2])
    plt.savefig(str, format='png', dpi=plotDPI)
    if(args.quick > 1):
        plt.show()
    plt.close()

if(args.quick == 0):
    os.system('magick -delay %f figs/melangeMelt*.png -colors 256 -depth 256 figs/melangeMelt.gif' %(500/((maxStep-startStep)/sizeStep)))
    if(cleanPNGs):
        os.system('rm -f figs/melangeMelt*.png')