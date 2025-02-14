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

resultFolder = '/results'
hFacWeighted = True
compareHughes = True

plt.figure(figsize=(12, 8))
ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,3)
ax3 = plt.subplot(2,2,2)
ax4 = plt.subplot(2,2,4)
axes = [ax1, ax2, ax3, ax4]
labels = ["A","C","B","D"]

for m in [0,1]:
    if(m == 0):
        dirNames = ['rbcs_12_02','rbcs_12_05','rbcs_12_10','rbcs_12_20','rbcs_12_40']
        legendNames = ['$\\lambda = 0.02$','$\\lambda = 0.05$','$\\lambda = 0.10$','$\\lambda = 0.20$','$\\lambda = 0.40$']
        fileName = 'Lambda_0025_ra'
        titleText = 'U = 0.12 m/s'
        mainColor = 'green'
        spds = [.12]*len(dirNames)
    elif(m == 1):
        dirNames = ['rbcs_02_20','rbcs_05_20','rbcs_12_20','rbcs_25_20','rbcs_40_20']
        legendNames = ['$u = 0.02$ m/s','$u = 0.05$ m/s','$u = 0.12$ m/s','$u = 0.24$ m/s','$u = 0.40$ m/s',]
        fileName = 'Speed_0025'
        titleText = '$\\lambda = 0.20$'
        mainColor = 'blue'
        spds = [.02,.05,.12,.25,.40]

    #Normalized Errors fields
    u1 = np.zeros((3,len(dirNames),150))
    z1 = np.ones((len(dirNames),150))*-600
    u2 = np.zeros((3,len(dirNames),46))
    z2 = np.zeros((len(dirNames),46))

    if(hFacWeighted):
        dynName = ['dynDiag', 'dynDiag', 'dynDiag']
    else:
        dynName = ['dynMassDiag', 'dynMassDiag', 'dynMassDiag']

    name = ["U"]
    units = ["[u(z)/U]"]

    for k in range(len(name)):
        print('\t' + name[k])
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
            iceStart = np.argmin(np.abs(x[0,:] - 13000)) 
            iceEnd = np.argmin(np.abs(x[0,:] - 15500))
            # print('Average is is x =', x[0,iceStart],',',x[0,iceEnd],'index', iceStart,',',iceEnd)
            
            i = maxStep

            data = mds.rdmds(dirName + resultFolder + "/%s"%(dynName[k]), i)
            if k == 0:
                lvl = [-0,1.2]
                seedColor = mainColor
                plotData = np.squeeze(data[k+2,:,:,:]) * 1 / spds[l]
            elif k == 1:
                lvl = [-0.001,0.001]
                seedColor = 'red'
                plotData = np.squeeze(data[k+2,:,:,:]) * 1
            elif k == 2:
                lvl = [-.01,0.01]
                seedColor = 'peach'
                plotData = np.squeeze(data[k+2,:,:,:]) * 1 
            
            n = len(dirNames) - 1 

            color1 = 'xkcd:light ' + seedColor
            color2 = 'xkcd:dark ' + seedColor
            rlin = np.linspace(mcolors.to_rgb(color1)[0],mcolors.to_rgb(color2)[0],n)
            glin = np.linspace(mcolors.to_rgb(color1)[1],mcolors.to_rgb(color2)[1],n)
            blin = np.linspace(mcolors.to_rgb(color1)[2],mcolors.to_rgb(color2)[2],n)
            colors = np.stack((rlin,glin,blin), axis=0)

            
            if(hFacWeighted):
                hFacC = mds.rdmds(dirName + resultFolder + "/hFacC") #must weight by ocean fraction

            u1[k,l,:len(z)] = np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2))
            z1[l,:len(z)] = z
            if(l != 0 ):
                if(hFacWeighted):
                    axes[m*2].plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=legendNames[l],color=colors[:,l-1])
                else:
                    axes[m*2].plot(np.mean(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color=colors[:,l-1])
            else:
                if(hFacWeighted):
                    axes[m*2].plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=legendNames[l],color='xkcd:gray')
                else:
                    axes[m*2].plot(np.mean(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color='xkcd:gray')
            if(l == 3):
                bergDepth = -1*np.fromfile(dirName + '/input/icebergs_depths.bin', dtype='>f8')
                bergDepth[bergDepth == 0] = np.nan #dont avg ones that don't exist
                axes[m*2].plot([-1,3],[np.nanpercentile(bergDepth,10),np.nanpercentile(bergDepth,10)],color='red',linestyle='--',alpha=.5)
                axes[m*2].plot([-1,3],[np.nanmedian(bergDepth),np.nanmedian(bergDepth)],color='red',linestyle='--',alpha=.5)
            # Compare against Hughes Directly
        if(compareHughes):
            if (m == 0):
                pack = ['020','050','100','200','400']
            else:
                pack = ['020','050','120','250','400']
            for l in range(len(pack)):
                if(m == 0):
                    nc = netCDF4.Dataset('/Users/psummers8/Documents/MITgcm/MITgcm/Sandbox_Berg/6555200/run.lambda_%s.U_120.nc' % pack[l])
                else:
                    nc = netCDF4.Dataset    ('/Users/psummers8/Documents/MITgcm/MITgcm/Sandbox_Berg/6555200/run.lambda_200.U_%s.nc' % pack[l])
                Z = nc['Z'][:]
                if(k == 0):
                    plotData = nc['U'][:,1:-1,100:350]/spds[l]
                elif(k == 1):
                    plotData = nc['W'][:,1:-1,100:350]
                elif(k == 2):
                    plotData = nc['V'][:,1:-1,100:350]
                plotData[plotData == 0] = np.nan
                u2[k,l,:] = np.nanmean(plotData,axis=(1,2))
                z2[l,:] = Z
                if(l == 0):
                    axes[m*2].plot(np.nanmean(plotData,axis=(1,2)),Z,linewidth=2,color='xkcd:gray',linestyle='--',alpha=.5)
                else:
                    axes[m*2].plot(np.nanmean(plotData,axis=(1,2)),Z,linewidth=2,color=colors[:,l-1],linestyle='--',alpha=.5)

        axes[m*2].set_xlim(lvl)
        axes[m*2].plot(np.cos(z * np.pi /600)* 1, z,linewidth=1,color='gray',linestyle=':')
        axes[m*2].set_xlabel(name[k] + " " + units[k])
        axes[m*2].set_ylabel('Depth [m]')
        # plt.title("%s over melange at %.02f days" % (name[k], i/86400.0*dt))
        axes[m*2].set_title(titleText)
        axes[m*2].set_ylim([-300, 0])
        axes[m*2].grid(alpha=.5)
        j = i/sizeStep + startStep
        # plt.show()
        axes[m*2].legend()

        totalError = np.zeros(np.shape(z1)[1])
        for l in range(len(dirNames)):
            zSlice = np.argmin(np.abs(z1[l,:] + 100))
            zSlice_Deep = np.argmin(np.abs(z1[l,:] + 275))
            u2_interp = np.interp(-z1[l,:],-z2[l,:],u2[k,l,:])
            forcingSpeed = np.cos(z1[l,:] * np.pi /600)
            error = (u1[k,l,:] - u2_interp)/forcingSpeed
            localError = error**2
            totalError = totalError+ localError/len(dirNames)
            # if l !=0:
            #     axes[m*2+1].plot(error,z1[l,:],label=legendNames[l],color=colors[:,l-1])
            # else:
            #     axes[m*2+1].plot(error,z1[l,:],label=legendNames[l],color='xkcd:light gray')
            top100error = np.sqrt(np.nanmean(localError[:zSlice]))*100
            top250error = np.sqrt(np.nanmean(localError[:zSlice_Deep]))*100
            if l !=0:
                axes[m*2+1].plot(error,z1[l,:],label=("%.1f%%, %.1f%%") 
                    %(top100error,top250error),color=colors[:,l-1])
            else:
                axes[m*2+1].plot(error,z1[l,:],label=("%.1f%%, %.1f%%") 
                    %(top100error,top250error),color='xkcd:gray')
        axes[m*2+1].plot([-1,3],[np.nanpercentile(bergDepth,10),np.nanpercentile(bergDepth,10)],color='red',linestyle='--',alpha=.5)
        axes[m*2+1].plot([-1,3],[np.nanmedian(bergDepth),np.nanmedian(bergDepth)],color='red',linestyle='--',alpha=.5)
        axes[m*2+1].set_title("Residual for "+ titleText)
        if k == 0:
            axes[m*2+1].set_xlim([-.5, 0.5])
        axes[m*2+1].set_ylim([-300, 0])
        # axes[m*2+1].set_xlabel("∆"+name[k] + " " + units[k] +"\n $\\sum\\epsilon^2$: %.2e" %(np.nanmean(totalError[:zSlice])))
        # axes[m*2+1].set_xlabel("∆"+name[k] + " " + units[k])
        axes[m*2+1].set_xlabel("Relative Residual [ ]")
        axes[m*2+1].set_ylabel('Depth [m]')
        axes[m*2+1].grid(alpha=.5)
        axes[m*2+1].legend(title="RMSE: 100 m, 275 m")    
        # str = "figs/depthError%s%s.png" % (name[k],fileName)

plt.tight_layout()

for label,ax in zip(labels,axes):
    ax.text(
        ax.get_xlim()[0], ax.get_ylim()[1], label,
        fontsize='x-large', va='bottom',ha='center', fontfamily='sans serif')

str = "figs/HughesLamdbdaSpd_relU.png"
plt.savefig(str, format='png', dpi=plotDPI)
plt.show()
plt.close()