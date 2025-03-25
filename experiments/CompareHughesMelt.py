from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
import sys
import cmocean
import fileinput
import gsw
import netCDF4

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 300
cleanPNGs = True

print('Plot DPI:',plotDPI,'; clean PNGs?',cleanPNGs)

resultFolder = '/results'
hFacWeighted = True

# dirNames = ['rbcs_12_20','rbcs_12_20_b1_d0','rbcs_12_20_b0_d0']
dirNames = ['rbcs_12_20','rbcs_12_20_b1_d0','rbcs_12_20_b0_d0','rbcs_melt_b1_d1','rbcs_melt_b1_d0','rbcs_melt_b0_d0']
legendNames = ['Block + Drag','Block + No Drag','No Block + No Drag']
fileName = 'meltRates'
titleText = '$\\lambda = .2$ U = 0.12 m/s'
mainColor = 'lavender'
spds = [.12]*len(dirNames)


dynName = ['dynDiag', 'dynDiag', 'BRGFlx']

name = ["U","W","Freshwater Flux"]
units = ["[u(z)/U]","[m/s]","[m^2/s]"]


figLabels = ["A","B"]
fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
ax1 = axes[0]
ax2 = axes[1]

#Calculate freezing depth
lam1 = -5.75e-2
lam2 = -7.61e-2
lam3 = 9.01e-4
s_b = 34
t_b = -8.3
z_freeze = -(lam1 * 34 + lam2) / (t_b/600 - lam3)


ax_i = 0
for k in [0,2]: #controling which values we plot in each frame
    print('\t' + name[k])
    # plt.figure()
    for l in range(len(dirNames)):
        #Find time steps to take
        maxStep = 0
        sizeStep = 1e10
        startStep = 1e10

        dt = 0.0   
        for line in fileinput.input(dirNames[l] + '/input/data'):
                if "deltaT=" in line:
                    dt = float(line[8:-2])
        # print('dt is loaded as', dt)

        for file in os.listdir(dirNames[l] + resultFolder):
            # print(file)
            if "dynMassDiag.0" in file:
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

        # print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)


        dirName = dirNames[l]
        #Comparing Diff Runs
        x = mds.rdmds(dirName + resultFolder +"/XC")
        y = mds.rdmds(dirName + resultFolder + "/YC")
        z = np.squeeze(mds.rdmds(dirName + resultFolder+ "/RC"))

        dz = np.append(-z[0]*2,z[:-1]-z[1:]) #this is strictly wrong, but very close. 

        for i in range(len(z)):
            if i == 0:
                dz[i] = -2*z[i]
            else:
                dz[i] = dz[i]*2 - dz[i-1]
        # print(z)
        # print(dz)
        iceStart = np.argmin(np.abs(x[0,:] - 13000)) 
        iceEnd = np.argmin(np.abs(x[0,:] - 15500))
        # print('Average is is x =', x[0,iceStart],',',x[0,iceEnd],'index', iceStart,',',iceEnd)
        
        # i = maxStep * 4/9
        i = maxStep

        data = mds.rdmds(dirName + resultFolder + "/%s"%(dynName[k]), i)
        if k == 0:
            lvl = [-0,2]
            seedColor = mainColor
            plotData = np.squeeze(data[2,:,:,:]) * 1 / spds[l]
        elif k == 1:
            lvl = [-0.001, 0.001]
            seedColor = 'red'
            plotData = np.squeeze(data[3,:,:,:]) * 1
        elif k == 2:
            lvl = [-0.01,0.6]
            seedColor = 'orange'
            plotData = np.squeeze(data[0,:,:,:]) * 1

        n = 3
        color1 = 'xkcd:light ' + seedColor
        color2 = 'xkcd:dark ' + seedColor
        rlin = np.linspace(mcolors.to_rgb(color1)[0],mcolors.to_rgb(color2)[0],n)
        glin = np.linspace(mcolors.to_rgb(color1)[1],mcolors.to_rgb(color2)[1],n)
        blin = np.linspace(mcolors.to_rgb(color1)[2],mcolors.to_rgb(color2)[2],n)
        colors = np.stack((rlin,glin,blin), axis=0)
        # print(colors)
        # for ii in range(iceEnd - iceStart):
        #     for j in range(np.shape(y[1:-1,:])[0]):
        #         plt.plot(plotData[:,j+1,iceStart + ii],z,linewidth=.5,alpha=.1,color=cm)
        
        hFacC = mds.rdmds(dirName + resultFolder + "/hFacC") #must weight by ocean fraction
        if(l < 3):
            if(k == 0 or k == 1):
                axes[ax_i].plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=legendNames[l]+"+ No Melt",linestyle='--',color=colors[:,l])
            else:
                pass
                # axes[ax_i].plot(np.nansum(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2))/dz,z,linewidth=2,label=legendNames[l]+"+ No Melt",linestyle='--',color=colors[:,l])
        else:
            if(k == 0 or k == 1):
                axes[ax_i].plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=legendNames[l-3]+"+ Melt",color=colors[:,l-3])
            else:
                axes[ax_i].plot(np.nansum(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2))/dz,z,linewidth=2,label=legendNames[l-3]+"+ Melt",color=colors[:,l-3])
        
        if(l == 0):
            bergDepth = -1*np.fromfile(dirName + '/input/icebergs_depths.bin', dtype='>f8')
            bergDepth[bergDepth == 0] = np.nan #dont avg ones that don't exist
            axes[ax_i].plot([-1,3],[np.nanpercentile(bergDepth,10),np.nanpercentile(bergDepth,10)],color='red',linestyle='--',alpha=.5)
            axes[ax_i].plot([-1,3],[np.nanmedian(bergDepth),np.nanmedian(bergDepth)],color='red',linestyle='--',alpha=.5)

            if(k == 2):
                axes[ax_i].plot([-1,3],[z_freeze,z_freeze],color='xkcd:sky blue',linestyle='--',alpha=.5,label='Liquidus Line')
                axes[ax_i].plot([-1,3],[np.nanmin(bergDepth),np.nanmin(bergDepth)],color='black',linestyle='--',alpha=.5,label="Deepest Draft")

    axes[ax_i].set_xlim(lvl)
    if k == 0:
        axes[ax_i].plot(np.cos(z * np.pi /600)* 1, z,linewidth=1,color='gray',linestyle=':')
    axes[ax_i].set_xlabel(name[k] + " " + units[k])
    axes[ax_i].set_ylabel('Depth [m]')
    # plt.title("%s over melange at %.02f days" % (name[k], i/86400.0*dt))
    axes[ax_i].set_title(titleText)
    axes[ax_i].set_ylim([-300, 0])
    axes[ax_i].grid(alpha=.5)
    j = i/sizeStep + startStep
    # plt.show()
    axes[ax_i].legend()
    ax_i += 1
for label,ax in zip(figLabels,axes.flatten()):
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')
str = "figs/windMelt%s.png" % (fileName)
plt.savefig(str, format='png',dpi=plotDPI)
plt.show() 
plt.close()


