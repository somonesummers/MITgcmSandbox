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
compareHughes = True

# dirNames = ['Example_L20_U12','Example_L20_U12_NLFS']
# legendNames = ['Linear FS','Non-Linear FS']
# fileName = 'NLFS'
# titleText = 'Non-linear Free Surface'
# mainColor = 'red'
# spds = [.12]*len(dirNames)

# dirNames = ['rbcs_12_20','rbcs_12_20_b1_d0','rbcs_12_20_b0_d0']
# legendNames = ['block + drag','block + no drag','no block + no drag']
# fileName = 'controlCases'
# titleText = '$\\lambda = .2$ U = 0.12 m/s'
# mainColor = 'lavender'
# spds = [.12]*len(dirNames)

dirNames = ['rbcs_12_02','rbcs_12_05','rbcs_12_10','rbcs_12_20','rbcs_12_40']
legendNames = ['$\\lambda = 0.02$','$\\lambda = 0.05$','$\\lambda = 0.10$','$\\lambda = 0.20$','$\\lambda = 0.40$']
fileName = 'Lambda_0025_temp'
titleText = 'U = 0.12 m/s'
mainColor = 'green'
spds = [.12]*len(dirNames)

# dirNames = ['rbcs_02_20','rbcs_05_20','rbcs_12_20','rbcs_25_20','rbcs_40_20']
# legendNames = ['$u = 0.02$ m/s','$u = 0.05$ m/s','$u = 0.12$ m/s','$u = 0.24$ m/s','$u = 0.40$ m/s',]
# fileName = 'Speed_0025'
# titleText = '$\\lambda = 0.20$'
# mainColor = 'blue'
# spds = [.02,.05,.12,.25,.40]

# dirNames = ['rbcs_dx_100','rbcs_dx_200','rbcs_dx_400','rbcs_dx_800','rbcs_dx_1200','rbcs_dx_2400',]
# legendNames = ['dx = 100 m','dx = 200 m','dx = 400 m','dx = 800 m','dx = 1200 m','dx = 2400 m']
# fileName = 'Dx_0025'
# titleText = 'U = 0.12 m/s $\\lambda = 0.20$'
# mainColor = 'teal'
# spds = [.12]*len(dirNames)

# dirNames = ['rbcs_dz_150','rbcs_dz_100','rbcs_dz_75','rbcs_dz_50','rbcs_dz_25','rbcs_dz_12','rbcs_dz_10',]
# legendNames = ['Nr = 150','Nr = 100','Nr = 75','Nr = 50','Nr = 25','Nr = 12','Nr = 10']
# # dirNames = ['rbcs_dz_12','rbcs_dz_10',]
# # legendNames = ['Nr = 12','Nr = 10']
# fileName = 'Dz_0025'
# titleText = 'U = 0.12 m/s $\\lambda = 0.20$'
# mainColor = 'pink'
# spds = [.12]*len(dirNames)

#Normalized Errors fields
u1 = np.zeros((3,len(dirNames),150))
z1 = np.ones((len(dirNames),150))*-600
u2 = np.zeros((3,len(dirNames),46))
z2 = np.zeros((len(dirNames),46))




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
        #Comparing Diff Runs
        x = mds.rdmds(dirName + resultFolder +"/XC")
        y = mds.rdmds(dirName + resultFolder + "/YC")
        z = np.squeeze(mds.rdmds(dirName + resultFolder+ "/RC"))
        iceStart = np.argmin(np.abs(x[0,:] - 13000)) 
        iceEnd = np.argmin(np.abs(x[0,:] - 15500))
        # print('Average is is x =', x[0,iceStart],',',x[0,iceEnd],'index', iceStart,',',iceEnd)
        
        # i = maxStep * 4/9
        i = maxStep

        data = mds.rdmds(dirName + resultFolder + "/%s"%(dynName[k]), i)
        if k == 0:
            lvl = [-0,1.6]
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
        # print(colors)
        # for ii in range(iceEnd - iceStart):
        #     for j in range(np.shape(y[1:-1,:])[0]):
        #         plt.plot(plotData[:,j+1,iceStart + ii],z,linewidth=.5,alpha=.1,color=cm)
        
        if(hFacWeighted):
            hFacC = mds.rdmds(dirName + resultFolder + "/hFacC") #must weight by ocean fraction

        u1[k,l,:len(z)] = np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2))
        z1[l,:len(z)] = z
        if(l != 0 ):
            if(hFacWeighted):
                plt.plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=legendNames[l],color=colors[:,l-1])
            else:
                plt.plot(np.mean(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color=colors[:,l-1])
        else:
            if(hFacWeighted):
                plt.plot(np.average(plotData[:,1:-1,iceStart:iceEnd], weights=hFacC[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=legendNames[l],color='xkcd:light gray')
            else:
                plt.plot(np.mean(plotData[:,1:-1,iceStart:iceEnd], axis=(1,2)),z,linewidth=2,label=dirName,color='xkcd:light gray')
        ax = plt.gca()
        if(l == 0):
            bergDepth = -1*np.fromfile(dirName + '/input/icebergs_depths.bin', dtype='>f8')
            bergDepth[bergDepth == 0] = np.nan #dont avg ones that don't exist
            plt.plot([-1,3],[np.nanpercentile(bergDepth,10),np.nanpercentile(bergDepth,10)],color='red',linestyle='--',alpha=.5)
            plt.plot([-1,3],[np.nanmedian(bergDepth),np.nanmedian(bergDepth)],color='red',linestyle='--',alpha=.5)
        # Compare against Hughes Directly
    if(compareHughes):
        # pack = ['020','050','100','200','400']
        pack = ['200']
        # spds=[.12]
        # pack = ['020','050','120','250','400']
        # spds = [.02,.05,.12,.25,.40]
        for l in range(len(pack)):
            nc = netCDF4.Dataset('/Users/psummers8/Documents/MITgcm/MITgcm/Sandbox_Berg/6555200/run.lambda_%s.U_120.nc' % pack[l])
            # nc = netCDF4.Dataset    ('/Users/psummers8/Documents/MITgcm/MITgcm/Sandbox_Berg/6555200/run.lambda_200.U_%s.nc' % pack[l])
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
                plt.plot(np.nanmean(plotData,axis=(1,2)),Z,linewidth=2,color='xkcd:light gray',linestyle='--',alpha=.5)
            else:
                plt.plot(np.nanmean(plotData,axis=(1,2)),Z,linewidth=2,color=colors[:,l-1],linestyle='--',alpha=.5)

    ax.set_xlim(lvl)
    plt.plot(np.cos(z * np.pi /600)* 1, z,linewidth=1,color='gray',linestyle='--')
    plt.xlabel(name[k] + " " + units[k])
    plt.ylabel('Depth [m]')
    # plt.title("%s over melange at %.02f days" % (name[k], i/86400.0*dt))
    plt.title(titleText)
    plt.ylim([-300, 0])
    plt.grid(alpha=.5)
    j = i/sizeStep + startStep
    # plt.show()
    plt.legend()
      
    str = "figs/depth%s%s.png" % (name[k],fileName)
    plt.savefig(str, format='png',dpi=plotDPI)
    plt.show() 
    plt.close()

    plt.figure()
    totalError = np.zeros(np.shape(z1)[1])
    for l in range(len(dirNames)):
        # u2_interp = np.interp(-z1[l,:],-z2[l,:],u2[k,l,:])
        u2_interp = np.interp(-z1[l,:],-z2[0,:],u2[k,0,:])
        error = u1[k,l,:] - u2_interp
        totalError = totalError + error**2
        # plt.plot(u1[k,l,:],z1[l,:],label='Summers')
        # plt.plot(u2[k,l,:],z2[l,:],label='Hughes')
        if l !=0:
            plt.plot(error,z1[l,:],label=legendNames[l],color=colors[:,l-1])
        else:
            plt.plot(error,z1[l,:],label=legendNames[l],color='xkcd:light gray')
    plt.plot([-1,3],[np.nanpercentile(bergDepth,10),np.nanpercentile(bergDepth,10)],color='red',linestyle='--',alpha=.5)
    plt.plot([-1,3],[np.nanmedian(bergDepth),np.nanmedian(bergDepth)],color='red',linestyle='--',alpha=.5)
    plt.title("Residual for "+ titleText)
    if k == 0:
        plt.xlim([-.5, .5])
    elif k == 1:
        plt.xlim([-.0005, .0005])
    elif k == 2:
        plt.xlim([-.005, .005])
    plt.ylim([-300, 0])
    zSlice = np.argmin(np.abs(z1[l,:] + 100))
    plt.xlabel("∆"+name[k] + " " + units[k] +" MSE: %.2f,%.2f" %(np.nansum(totalError[:zSlice]),np.nansum(totalError)))
    plt.ylabel('Depth [m]')
    plt.grid(alpha=.5)
    plt.legend()    
    str = "figs/depthError%s%s.png" % (name[k],fileName)
    plt.savefig(str, format='png',dpi=plotDPI)
    plt.show()
    plt.close()

