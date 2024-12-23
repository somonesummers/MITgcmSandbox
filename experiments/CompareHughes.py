from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
import sys
import cmocean
import fileinput
import gsw

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
elif(os.path.isfile('plotHelper.py')):
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')
print('Plot DPI:',plotDPI,'; clean PNGs?',cleanPNGs)

resultFolder = '/results'
hFacWeighted = True
# dirNames = ['Hughes_smag_LR_30','Hughes_smag_30','Hughes_smag_LR_50','Hughes_smag_50','Hughes_smag_100']
# dirNames = ['Hotel_40','Hotel_20','Hotel_10','Hotel_05','Hotel_02']
# dirNames = ['india_40','india_25','india_16','india_05','india_02']
# dirNames = ['Hotel_80','Hotel_60','Hotel_40','Hotel_20','Hotel_10','Hotel_05','Hotel_02']
# dirNames = ['juliet_001','juliet_002','juliet_005','juliet_01','juliet_02','juliet_05','juliet_10']
dirNames = ['kilo_15','kilo_30','kilo_60','kilo_90','kilo_120']

# spds = [.4,.25,.16,.05,.02]
spds = [.12]*len(dirNames)

if(hFacWeighted):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag']
else:
    dynName = ['dynMassDiag', 'dynMassDiag', 'dynMassDiag']

name = ["U", "W", "V"]
units = ["[u(z)/U]", "[m/s]", "[m/s]"]


for k in range(len(name)):
    print('\t' + name[k])
    plt.figure()
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

        #Comparing Diff resolutions
        x = mds.rdmds(dirName + resultFolder +"/XC")
        y = mds.rdmds(dirName + resultFolder + "/YC")
        z = np.squeeze(mds.rdmds(dirName + resultFolder+ "/RC"))
        iceStart = np.argmin(np.abs(x[0,:] - 13100)) 
        iceEnd = np.argmin(np.abs(x[0,:] - 15400))
        # print('Average is is x =', x[0,iceStart],',',x[0,iceEnd],'index', iceStart,',',iceEnd)
        
        i = maxStep
        data = mds.rdmds(dirName + resultFolder + "/%s"%(dynName[k]), i)
        if k == 0:
            lvl = [-0,2.5]
            seedColor = 'orange'
            plotData = np.squeeze(data[k+2,:,:,:]) * 1 / spds[l]
        elif k == 1:
            lvl = [-0.01,0.01]
            seedColor = 'red'
            plotData = np.squeeze(data[k+2,:,:,:]) * 1
        elif k == 2:
            lvl = [-.01,0.01]
            seedColor = 'green'
            plotData = np.squeeze(data[k+2,:,:,:]) * 1 
        
        n = len(dirNames) - 1 
        
        color1 = 'xkcd:dark ' + seedColor
        color2 = 'xkcd:light ' + seedColor
        rlin = np.linspace(mcolors.to_rgb(color1)[0],mcolors.to_rgb(color2)[0],n)
        glin = np.linspace(mcolors.to_rgb(color1)[1],mcolors.to_rgb(color2)[1],n)
        blin = np.linspace(mcolors.to_rgb(color1)[2],mcolors.to_rgb(color2)[2],n)
        colors = np.stack((rlin,glin,blin), axis=0)
        # print(colors)
        # for ii in range(iceEnd - iceStart):
        #     for j in range(np.shape(y[1:-1,:])[0]):
        #         plt.plot(plotData[:,j+1,iceStart + ii],z,linewidth=.5,alpha=.1,color=cm)
        if(hFacWeighted):
            hFacC = mds.rdmds(dirName + resultFolder + "/hFacC") #must weight by ocean fraction
        if(l != n ):
            if(hFacWeighted):
                plt.plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color=colors[:,l])
            else:
                plt.plot(np.mean(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color=colors[:,l])
        else:
            if(hFacWeighted):
                plt.plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color='xkcd:light gray')
            else:
                plt.plot(np.mean(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color='xkcd:light gray')
        ax = plt.gca()
        if(l == 3):
            bergDepth = -1*np.fromfile(dirName + '/input/icebergs_depths.bin', dtype='>f8')
            bergDepth[bergDepth == 0] = np.nan #dont avg ones that don't exist
            plt.plot([-1,3],[np.nanpercentile(bergDepth,10),np.nanpercentile(bergDepth,10)],color='red',linestyle='--',alpha=.5)
            plt.plot([-1,3],[np.nanmedian(bergDepth),np.nanmedian(bergDepth)],color='red',linestyle='--',alpha=.5)
        ax.set_xlim(lvl)
        plt.plot(np.cos(z * np.pi /600)* 1, z,linewidth=1,color='gray',linestyle='--')
        plt.xlabel(name[k] + " " + units[k] + ' %.3f %.3f nan: %i' %(np.nanmin(plotData[:,1:-1]),np.nanmax(plotData[:,1:-1]),np.max(np.isnan(plotData[:,1:-1]))))
        plt.ylabel('Depth [m]')
        plt.title("%s over melange at %.02f days" % (name[k], i/86400.0*dt))
        plt.ylim([-300, 0])
        plt.grid(alpha=.5)
        j = i/sizeStep + startStep
        # plt.show()
    plt.legend()
    # plt.show()    
    str = "figs/depth%s%05i.png" % (name[k],j)
    plt.savefig(str, format='png',dpi=plotDPI)
    plt.close()
