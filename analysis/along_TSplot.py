from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse
from matplotlib import colormaps

distColorMap = 'Spectral_r'

def getLengthColor(value, cbar):
    vMin = 0
    vMax = cbar.ax.get_ylim()[1]
    # print(f'asking for {value}')
    colorMapForLine = colormaps[distColorMap]
    n_colors = 128
    widthColors = colorMapForLine(np.linspace(0,1,n_colors+1))
    return widthColors[int((value-vMin)/(vMax-vMin)*n_colors)]

parser = argparse.ArgumentParser(description='Plot TS at many xCrossSections downfjord')

parser.add_argument('-t','--timeRange', nargs=1, type=int, default = [685],
                    help='optional specification of time to plot in DAYS [default = 100]')
args = parser.parse_args()

# Pick cross section to view from file or default
xCrossSection = 15000
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

# if(args.xCrossSection != None):
#     print('** Manual xCrossSection detected **')
#     xCrossSection = args.xCrossSection

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

if(args.timeRange != None):
    startStep = args.timeRange[0] * 86400 / dt
    # maxStep = args.timeRange[1] * 86400 / dt

# os.system('rm -f figs/TSPlot*.png')
# os.system('rm -f figs/TSPlot*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")
 
if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.ones(np.shape(x))

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
# freezeS100 = [0,50]
# freezeT100 = [-.011,-2.919]

i = startStep
data = mds.rdmds("results/dynDiag", i)
plt.figure(figsize=(12, 6))
ax3 = plt.subplot(133)
ax2 = plt.subplot(132)
ax1 = plt.subplot(131)

locations_to_plot = [6e3,12e3,25e3,50e3,70e3,90e3]
for xCrossSection in locations_to_plot:
    xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
    print('cross section is x =', x[0,xSlice],'index', xSlice)

    for j in range(np.shape(y)[0]):
        if(topo[j,xSlice] == 0): #this is a wall
            data[:,:,j,xSlice] = np.nan

    sc_color = ax3.scatter([],[],[],[],vmin=0,vmax=100,cmap=distColorMap)
    if(xCrossSection == locations_to_plot[0]):
        cbar_color = plt.colorbar(sc_color)
    for j in range(np.shape(y)[0]):
            ax3.plot(data[0,:,j,xSlice],np.squeeze(z),alpha = .1, color=getLengthColor(xCrossSection/1e3,cbar_color))
    cp = ax3.plot(np.nanmean(data[0,:,:,xSlice],1),np.squeeze(z),c=getLengthColor(xCrossSection/1e3,cbar_color),linestyle='-',marker='.')
    if(xCrossSection == locations_to_plot[0]):
        cbar_color.set_label('Distance Along Fjord [m]')
        ax3.set_title('Temperature')
        ax3.set_xlabel('Temperature [C]')
        ax3.set_ylabel('Depth [m]')
        ax3.grid(alpha = .5)

    for j in range(np.shape(y)[0]):
            ax2.plot(data[1,:,j,xSlice],np.squeeze(z),alpha = .1, color=getLengthColor(xCrossSection/1e3,cbar_color))
    cp = ax2.plot(np.nanmean(data[1,:,:,xSlice],1),np.squeeze(z),c=getLengthColor(xCrossSection/1e3,cbar_color),linestyle='-',marker='.')
    if(xCrossSection == locations_to_plot[0]):
        ax2.set_title('Salt')
        ax2.set_xlabel('Salinity')
        ax2.set_ylabel('Depth [m]')
        ax2.grid(alpha = .5)

    for j in range(np.shape(y)[0]):
        ax1.plot(data[1,:,j,xSlice],data[0,:,j,xSlice],alpha = .1, color=getLengthColor(xCrossSection/1e3,cbar_color))
    sc=ax1.plot(np.nanmean(data[1,:,:,xSlice],1),np.nanmean(data[0,:,:,xSlice],1),c=getLengthColor(xCrossSection/1e3,cbar_color),
                   alpha=1.,linestyle='-',marker='.')
    if(xCrossSection == locations_to_plot[0]):
        ax1.plot(mixingS,mixingT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
        ax1.plot(meltS,meltT,linewidth=.5,color='gray',alpha=.5,linestyle='--')
        ax1.plot(freezeS,freezeT,linewidth=.5,color='red',alpha=.5,linestyle='--')
        # ax1.plot(freezeS100,freezeT100,linewidth=.5,color='red',alpha=.5,linestyle='--')
        # cbar = plt.colorbar(sc)
        # cbar.set_label('Depth [m]')
        ax1.set_xlim([np.min(saltRange), np.max(saltRange)])
        ax1.set_ylim([np.min(tempRange), np.max(tempRange)])
        ax1.set_xlabel('Salt [PSU]')
        ax1.set_ylabel('Temperature [C]')
        ax1.set_title('Temp/Salt')
        plt.suptitle("TS along fjord at %.02f days" % (i/86400.0*dt))


str = f"figs/along_TSPlot_{args.timeRange}.png"
plt.tight_layout()
plt.savefig(str, format='png', dpi=plotDPI)

plt.show()
plt.close()

