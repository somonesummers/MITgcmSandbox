from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot Plume variables')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last frame only')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = None,
                    help='option specification of views to plot [default = all]')
args = parser.parse_args()

# Plot settings from local helper
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
    if "plumeDiag.0" in file:
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

os.system('rm -f figs/plumePlot*.png')
# os.system('rm -f figs/plumePlot*.gif')

x = mds.rdmds("results/XC")
dx = x[0,0]*2
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")



dynName = 'plumeDiag'
name = ["W", "Temp", "Sal", "CellMeltRate", "RadiusThickness","Flux","DensityABS"]
units = ["[m/s]", "[C]", "[PSU]", "[m/day]", "[m]","[m^3/s]","[kg/m^3]"]

 #amb is T,S,U,W,V

if(os.path.isfile('input/plumeMask.bin')):
    plumeMask = np.fromfile('input/plumeMask.bin', dtype='>f8')
    plumeMask = plumeMask.reshape(np.shape(x))
    plumeLocations = np.where(plumeMask > 1)

else:
    print('No Plumes Found, exiting PlumePlot')
    exit()

plumeLoc = [plumeLocations[0][0],plumeLocations[1][0]]
plumeType = np.max(plumeMask)
print('Plume Locaion is grid', plumeLoc )

peakDownFjord = 2
farFjord = -20 #try to get out of boundary layer

if(args.quick > 0):
    startStep = maxStep
    cleanPNGs = False

if(args.kValues == None):
    kList = range(len(name))
else:
    kList = args.kValues
for k in kList:
    print("\t %s" %name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        j = plumeLoc[0]
        data = mds.rdmds("results/%s"%(dynName), i)
        if k == 0:
            lvl = [-.5,3] 
            cm = "xkcd:raspberry"
            plotData = data[k,:,j,plumeLoc[1]]
            ambData = mds.rdmds("results/dynDiag", i)
            plotData2 = np.squeeze(ambData[3,:,j,plumeLoc[1]])
            plotData22 = np.squeeze(ambData[3,:,j,plumeLoc[1]+peakDownFjord])
            plotData3 = np.squeeze((ambData[2,:,j,plumeLoc[1]]**2 + ambData[3,:,j,plumeLoc[1]]**2 + ambData[4,:,j,plumeLoc[1]]**2)**(0.5))
            plotData4 = np.squeeze(ambData[3,:,j,farFjord])
        elif k == 1:
            lvl = tempRange
            cm = "xkcd:rose"
            plotData = data[k,:,j,plumeLoc[1]]
            ambData = mds.rdmds("results/dynDiag", i)
            plotData2 = np.squeeze(ambData[0,:,j,plumeLoc[1]])
            plotData22 = np.squeeze(ambData[0,:,j,plumeLoc[1]+peakDownFjord])
            plotData3 = np.squeeze((ambData[2,:,j,plumeLoc[1]]**2 + ambData[3,:,j,plumeLoc[1]]**2 + ambData[4,:,j,plumeLoc[1]]**2)**(0.5))
            plotData4 = np.squeeze(ambData[0,:,j,farFjord])
        elif k == 2:
            lvl = saltRange
            cm = "xkcd:green"
            plotData = data[k,:,j,plumeLoc[1]]
            ambData = mds.rdmds("results/dynDiag", i)
            plotData2 = np.squeeze(ambData[1,:,j,plumeLoc[1]])
            plotData22 = np.squeeze(ambData[1,:,j,plumeLoc[1]+peakDownFjord])
            plotData3 = np.squeeze((ambData[2,:,j,plumeLoc[1]]**2 + ambData[3,:,j,plumeLoc[1]]**2 + ambData[4,:,j,plumeLoc[1]]**2)**(0.5))
            plotData4 = np.squeeze(ambData[1,:,j,farFjord])
        elif k == 3:
            lvl = [-.1, 9]
            cm = "xkcd:violet"
            plotData = data[k,:,j,plumeLoc[1]]
        elif k == 4:
            lvl = [-1, 45]
            cm = "xkcd:lavender"
            plotData = data[k,:,j,plumeLoc[1]]
        elif k == 5:
            lvl = [-10,30000]
            cm = "xkcd:twilight"
            if(plumeType == 2):
                plotData = data[0,:,j,plumeLoc[1]]*dx*data[4,:,j,plumeLoc[1]]
        elif k == 6:
            rho0 = 1020
            lvl = [0,10]
            ambData = mds.rdmds("results/dynDiag", i)
            cm = "xkcd:pumpkin"
            salt = np.squeeze(data[2,:,j,plumeLoc[1]])
            if(i == startStep): #only calc pressure once
                pressure = -1 * np.ones(salt.shape) * 1028 * 9.81 /10e3
            CT = gsw.CT_from_t(salt, data[1,:,j,plumeLoc[1]], pressure)
            plotData = gsw.rho(salt, CT, pressure) - rho0 #in-stu density less rho0
            ambS = np.squeeze(ambData[1,:,j,plumeLoc[1]])
            ambCT = gsw.CT_from_t(salt, ambData[0,:,j,plumeLoc[1]], pressure)
            plotData2 = gsw.rho(ambS, ambCT, pressure) - rho0 #in-stu density less rho0

            ambS = np.squeeze(ambData[1,:,j,plumeLoc[1]+peakDownFjord])
            ambCT = gsw.CT_from_t(salt, ambData[0,:,j,plumeLoc[1]+peakDownFjord], pressure)
            plotData22 = gsw.rho(ambS, ambCT, pressure) - rho0 #in-stu density less rho0
            
            plotData3 = np.squeeze((ambData[2,:,j,plumeLoc[1]]**2 + ambData[3,:,j,plumeLoc[1]]**2 + ambData[4,:,j,plumeLoc[1]]**2)**(0.5))

            ambS = np.squeeze(ambData[1,:,j,farFjord])
            ambCT = gsw.CT_from_t(salt, ambData[0,:,j,farFjord], pressure)
            plotData4 = gsw.rho(ambS, ambCT, pressure) - rho0 #in-stu density less rho0

        plt.figure()
        plt.plot(plotData,np.squeeze(z),linewidth=1,color=cm)
        ax = plt.gca()
        ax.grid(alpha=.5)
        ax.set_xlim([np.min(lvl),np.max(lvl)])
        ax.set_xlabel(name[k] + " " + units[k] + ' %.3f %.3f nan: %i' %(np.nanmin(plotData),np.nanmax(plotData),np.max(np.isnan(plotData))))
        ax.set_ylabel('Depth [m]')
        ax.set_title("%s x = %i at %.02f days" % (name[k], x[0,plumeLoc[1]], i/86400.0*dt))
        if(k == 6 or k < 3):
            ax.plot(plotData2,np.squeeze(z),linewidth=1,color=cm,linestyle='--',label="Ambient")
            ax.plot(plotData22,np.squeeze(z),linewidth=1,color=cm,linestyle=':',label="Ambient %i cells Down Fjord" %peakDownFjord)
            ax.plot([],[],linewidth=1,color='xkcd:periwinkle',linestyle='--',label='Amb Speed [m/s]') # add to legend on ax
            ax.plot(plotData4,np.squeeze(z),linewidth=1,color='xkcd:gray',linestyle='--',label='fjord mouth',alpha=.5)
            ax.legend(loc='lower left')
            ax_2 = ax.twiny()
            ax_2.plot(plotData3,np.squeeze(z),linewidth=1,color='xkcd:periwinkle',linestyle='--',label='Amb Speed [m/s]')
            ax_2.set_xlabel('Amb Speed [m/s]',color='xkcd:periwinkle')
            ax_2.tick_params(axis='x',labelcolor='xkcd:periwinkle')
            ax_2.set_xlim([-0.02,0.42])

        j = i/sizeStep + startStep
        
        str = "figs/plumePlot%s%05i.png" % (name[k],j)
        plt.tight_layout()
        plt.savefig(str, format='png',dpi=plotDPI)
        if(args.quick > 1):
            plt.show()
        plt.close()
    if(args.quick == 0):
        os.system('magick -delay %f figs/plumePlot%s*.png -colors 256 -depth 256 figs/plumePlot%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/plumePlot*.png')
