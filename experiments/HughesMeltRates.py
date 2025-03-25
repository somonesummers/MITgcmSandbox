from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import xarray as xr

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 150
cleanPNGs = True

folders = ['wiskey_b1','wiskey_b0','wiskey_b1_ra','wiskey_b1_ra2']
labels = ['Iceberg Drag', 'No Drag','ra','ra2']
colors = ['xkcd:blue','xkcd:green','xkcd:red','xkcd:teal']
resultFolder = '/results_0025'
fileName = "_ra"

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

#Open figure before looping

dataset = xr.open_dataset('~/Documents/MITgcm/MITgcm/Sandbox_Berg/10366091/fields.deltaS_3.gamma_multiple_4.T0_2.nc')
x_fw = dataset['MELfwFlx']* -1/998  #kg/s -> m^3/s
x_z = dataset['Z']
x_t =  np.array(dataset['t'])
x_dz =  np.array(dataset['dz'])


#Total melt over time
plt.figure()
plt.plot(x_t/86400,np.nansum(x_fw[:,:,:,:],axis=(1,2,3)),color='xkcd:gray',label='Hughes 2024')
for j in range(len(folders)):
    folder = folders[j]

    dt = 0.0   
    for line in fileinput.input('%s/input/data' %folder):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    print('dt is loaded as', dt)

    #Find time steps to take
    maxStep = 0
    sizeStep = 1e10
    startStep = 1e10

    for file in os.listdir('%s%s' %(folder,resultFolder)):
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

    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

    dynName = ['BRGFlx']
    name = ['BRGfwFlx']
    units = ["[m^3/s]"]

    timeSteps = np.arange(startStep, maxStep + 1, sizeStep)
    fwOverTime = np.zeros(np.shape(timeSteps))
    for i in range(len(timeSteps)):
        data = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), timeSteps[i])
        fwOverTime[i] = np.nansum(data[0,:,:,:])
           
   

    plt.plot(timeSteps*dt/86400,fwOverTime,label=labels[j],color=colors[j])
# Wrap up after plotting everything

str = "figs/fwProduction%s.png" %fileName
plt.grid(alpha=.5)
plt.legend()
plt.ylabel('Total water flux [m^3/s]')
plt.xlabel('Time [days]')
plt.savefig(str, format='png',dpi=plotDPI)
plt.show()
plt.close()

#Open figure before looping
# Meltrate as a function of depth at start and end
plt.figure()

# hMax = np.nanmax(np.nansum(x_fw[0,:,:,:],axis=(1,2))/x_dz) # for plotting relative rates
hMax = 1
plt.plot(np.nansum(x_fw[0,:,:,:],axis=(1,2))/x_dz/hMax,x_z,\
    label='Hughes 12 hr',color='xkcd:gray',linestyle='--',alpha=.5)
plt.plot(np.nansum(x_fw[-1,:,:,:],axis=(1,2))/x_dz/hMax,x_z,\
    label='Hughes 7 day',color='xkcd:gray')
for i in range(len(folders)):
    folder = folders[i]

    dt = 0.0   
    for line in fileinput.input('%s/input/data' %folder):
            if "deltaT=" in line:
                dt = float(line[8:-2])
    print('dt is loaded as', dt)

    #Find time steps to take
    maxStep = 0
    sizeStep = 1e10
    startStep = 1e10

    for file in os.listdir('%s%s' %(folder,resultFolder)):
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

    print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

    z = np.squeeze(mds.rdmds("%s%s/RC" %(folder,resultFolder))) #midpoints of cells
    dz = np.append(-z[0]*2,z[:-1]-z[1:]) #z is negative values, but we want pos dz

    dynName = ['BRGFlx']
    name = ['BRGfwFlx']
    units = ["[m^3/s]"]

    dataStart = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), 12*3600/dt)
    dataEnd = mds.rdmds("%s%s/%s"%(folder,resultFolder, dynName[0]), maxStep-sizeStep)
    
    # sMax = np.nanmax(np.nansum(dataStart[0,:,:,:],axis=(1,2))/dz) #for plotting relative rates
    sMax = 1
    if(i == 0):
        plt.plot(np.nansum(dataStart[0,:,:,:],axis=(1,2))/dz/sMax,z,linestyle='--',color=colors[0],\
            label='12 hour',alpha=.5)
    plt.plot(np.nansum(dataEnd[0,:,:,:],axis=(1,2))/dz/sMax,z,label=labels[i],color=colors[i])
# Wrap up after plotting everything

str = "figs/fwDepth%s.png" %fileName
plt.grid(alpha=.5)
plt.ylim([-300, 0])
plt.ylabel('Depth [m]')
plt.xlabel('Fresh water flux per unit depth [m^2/s]')
plt.legend()
plt.savefig(str, format='png',dpi=plotDPI)
plt.show()
plt.close()


# plt.figure()
# x_u = dataset['U']
# x_v = dataset['V']
# x_x = np.array(dataset['X'])
# xSlice = np.argmin(np.abs(x_x))
# x_spd = (x_u[-1,:,:,:xSlice]**2 + x_v[-1,:,:,:xSlice]**2)**(0.5)
# x_spdBoolean = x_spd > 0.04
# plt.plot(np.mean(x_spdBoolean,axis=(1,2)),x_z,label='Hughes',color='xkcd:gray')

hFacC = mds.rdmds("%s%s/hFacC" %(folders[0],resultFolder))
for j in range(len(folders)):
    folder = folders[j]

    x = mds.rdmds("%s%s/XC" %(folder,resultFolder))-8000
    xSlice = np.argmin(np.abs(x[0,:]))
    
    data = mds.rdmds("%s%s/%s"%(folder,resultFolder, 'dynDiag'), maxStep)
    spd = (data[2,:,:,:xSlice]**2 + data[4,:,:,:xSlice]**2)**(0.5)
    spdBoolean = (spd > 0.04) * (hFacC[:,:,:xSlice] < .95)
    plt.plot(np.mean(spdBoolean,axis=(1,2)),z,label=labels[j],color=colors[j])

str = "figs/fwMeltSpd%s.png" %fileName
plt.grid(alpha=.5)
plt.ylim([-300, 0])
plt.ylabel('Depth [m]')
plt.xlabel('Fraction ocean cells above min melt speed [ ]')
plt.legend()
plt.savefig(str, format='png',dpi=plotDPI)
plt.show()
plt.close()


   