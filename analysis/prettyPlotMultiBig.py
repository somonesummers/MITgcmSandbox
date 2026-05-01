#!/usr/bin/env python3

# Paul Summers November 2025
# Made for making figures for melange collapse paper


import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from matplotlib import colormaps
import sys
from scipy.integrate import simpson
import argparse
import os
import fileinput
# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome
from MITgcmutils import mds
import cmocean
import glob
import pickle
import gsw
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
from matplotlib import patches

saveProfiles = False

constant = constants()

plt.rcParams.update({'font.size': 12})

# Take input options. Some defaults are set here, so be aware
parser = argparse.ArgumentParser(description='Plot options for melange changed over time')
parser.add_argument('-f','--files', nargs='*', default=['couplingResults/MITgcmRun_'],
                    help='file string used to name pickle files FILES00000.pickle [default = couplingResults/MITgcmRun_]')
parser.add_argument('-l','--labels', nargs='*', default=None,
                    help='file string label plots [default = None]')
parser.add_argument('-sgd','--sgd', nargs='*', type=float, default = None,
                    help='SGD, any value loads SGD from input, needs results RY')
parser.add_argument('--silent', action='count', default=0,
                    help='Option to silence showing of plots')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = [80, 80+365], #94, 459 
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-w','--windowMean', nargs='?', type=int, default=10,
                    help='window width for time averaging melt rates [default = 10]')
parser.add_argument('-ff','--fixForce', action='count', default=0,
                    help='Option to fix forcing left scale to +/- 5%%')
parser.add_argument('-ll','--legLoc', nargs='?', default=0,type=int,
                    help='legend location [default = 0]')
parser.add_argument('-dpi','--dpi', nargs='?', default=200,type=int,
                    help='dpi to print [default = 200]')
parser.add_argument('-y','--yCrossSection', nargs=1, type=float,default = 5000,
                    help='Y slice location [default 5000 m]')
parser.add_argument('-s','--shadow', action='count', default=1,
                    help='option of shadow for mélange [default (off)]')
parser.add_argument('-xl','--xlimit', nargs='?', type=int, default = 30000,
                    help='optional max x [default = 30 km]')
parser.add_argument('-Y','--Years', nargs='*', type=int, default = [1,2],
                    help='years to plot, default = [1,2]')
args = parser.parse_args()

showDensity = False
showZero = False
showQuiver = True

# print(args)
figLabels = ["(a)","(b)","(c)","(d)","(e)","(f)"] 
# These are out of order as I shuffled them after making the script. Sorry it is terrible
plt.figure(figsize=(12, 10),layout="tight")
# plt.figure(figsize=(16, 6))
ax5 = plt.subplot(3,2,(1,2))
ax6 = None #plt.subplot(2,4,(3,4))
ax1 = plt.subplot(3,2,3)
ax2 = plt.subplot(3,2,4) # I am so sorry. 
ax3 = plt.subplot(3,2,5) # It is really horrible. 
ax4 = plt.subplot(3,2,6)
if(ax6 != None):
    axes = [ax5, ax6, ax1, ax2, ax3, ax4]
else:
    axes = [ax5, ax1, ax2, ax3, ax4]
lStyle = ['-',':','--','-.',(0, (3, 2, 1, 2, 1, 2)),(0, (2, 3, 1, 2, 1, 2, 1, 2,))]
meltColors = ['xkcd:tomato red','xkcd:rose','xkcd:wine','xkcd:grape purple','xkcd:grape','xkcd:violet']
lWidth = 2
DOY_to_plot = 155
#colors matching phase
colorMapForLine = colormaps['cmo.phase']
n_colors = 5
colorTmp = colorMapForLine(np.linspace(0,1,n_colors+1))
colorList = colorTmp[:-1]
# Paul's Custom colors
# colorList = ['xkcd:light blue',
#                 'xkcd:bright orange','xkcd:tomato',
#                 'xkcd:faded pink','xkcd:dark blue']

if(args.silent > 0): #print status for script
    print(f"rampPlot: {args}")

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

x = mds.rdmds("results/XC") 
x = x - x[0,1]
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

years = args.Years
for i in range(len(years)):
    print(f'Year {years[i]}, index {i}')
    fileStr = 'couplingResults/MITgcmRun_'
    files = sorted(glob.glob('%s*.pickle'%fileStr))
    jShift = 0
    if(args.timeRange[0] == 0):
        toIterate = np.arange(len(files))
    else:
        if(args.timeRange[1] == 0):
            toIterate = np.arange(len(files[args.timeRange[0]:]))
        else:
            toIterate = np.arange(len(files[args.timeRange[0]:args.timeRange[1]]))
        jShift = args.timeRange[0]

    H0Time = np.zeros(np.shape(toIterate))
    VTime = np.zeros(np.shape(toIterate))
    UcTime = np.zeros(np.shape(toIterate))
    UtTime = np.zeros(np.shape(toIterate))
    lengthTime = np.zeros(np.shape(toIterate))
    iterationNumber = np.zeros(np.shape(toIterate))
    pressTime = np.zeros(np.shape(toIterate))
    BTime = np.zeros(np.shape(toIterate))
    BShapeTime = np.zeros([np.shape(toIterate)[0],20])
    H_fullTime =  np.zeros([np.shape(toIterate)[0],22])
    X_fullTime =  np.zeros([np.shape(toIterate)[0],22])
    gTime = np.zeros(np.shape(toIterate))
    spdTime = np.zeros(np.shape(toIterate))
    timeTime = np.zeros(np.shape(toIterate))
    muSTime = np.zeros(np.shape(toIterate))
    xTermTime = np.zeros(np.shape(toIterate))
    bFluxTime = np.zeros(np.shape(toIterate))
    endFluxTime = np.zeros(np.shape(toIterate))
    inFluxTime = np.zeros(np.shape(toIterate))
    for j in toIterate:
        if(j + jShift + years[i]*365 >= len(files)):
            print('Truncating final year, seem to have ended early')
            break #stop building data if end of file
        file = files[j + jShift + years[i]*365]
        with open(files[j + jShift + years[i]*365], 'rb') as fileName:
            data = pickle.load(fileName)
            fileName.close()
        it = file.replace(fileStr,'').replace('./', '').replace('.pickle', '')
        H0Time[j]=data.H0
        xTermTime[j] = data.X[0]
        X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
        H = np.concatenate(([data.H0], data.H, [data.HL]))
        W = np.concatenate(([data.W0], data.W, [data.WL]))
        H_fullTime[j,:] = H
        X_fullTime[j,:] = X_
        VTime[j] = simpson(H*W, x=X_)*1e-9
        lengthTime[j]=data.X[-1] - data.X[0]  #difference between these
        UcTime[j] = data.Uc
        UtTime[j] = data.Ut
        BTime[j] = -1 * np.mean(data.B) # we flip this for plotting purposes
        BShapeTime[j,:] = -1 * data.B /365.25
        timeTime[j] = data.t*365 
        pressTime[j] = data.force()#data.H0*data.pressure(data.H0) #should I use force?
        spdTime[j]=np.mean(data.U)
        gTime[j]=np.mean(data.gg)
        muSTime[j] = data.param.muS
        iterationNumber[j]=int(it)
        endFluxTime[j] = data.HL*data.U[-1] #m^2/y
        if(j != toIterate[0]): #include mass loss to length loss
            endFluxTime[j] += (data.H[-2]) * (lengthTime[j-1] - lengthTime[j]) / (data.dt) #m^2/y
        inFluxTime[j] = data.Ht*data.Uc #m^2/y
        bFluxTime[j] = -1*simpson(data.B, x=data.X_)  #m^3/s
    timeTime[timeTime == 0] = np.nan
    timeTime = timeTime - timeTime[0] #times of 0 should stay 0
    # if(not args.timeRange == None):
    #     timeTime = timeTime + args.timeRange[0]
    timeLabel = 'Day of year'
    
    from scipy import interpolate
    from MITgcmutils import mds
    remoteResult = args.files[0].rsplit('/',1)[0]
    # print(f'remoteResult is ({remoteResult})')
    dy = 2*y[0,0]
    for line in fileinput.input(f'{remoteResult}/../input/data'):
        if "ExternForcingPeriod=" in line:
            period = float(line[21:-2])
        elif "ExternForcingCycle=" in line:
            cycle = float(line[20:-2])
    nt = int(cycle/period)
    plumeMask = np.fromfile(f'{remoteResult}/../input/plumeMask.bin', dtype='>f8') #(2 is line, 3 is semi-cone)
    rad = np.fromfile(f'{remoteResult}/../input/runoffRad.bin', dtype='>f8')
    rad = rad.reshape((nt,len(plumeMask)))
    vel = np.fromfile(f'{remoteResult}/../input/runoffVel.bin', dtype='>f8')
    vel = vel.reshape((nt,len(plumeMask)))
    T_bc = np.fromfile(f'{remoteResult}/../input/EBCt.bin', dtype='>f8')
    extraDim = int(np.size(T_bc) / 32 / 24) #480 for sierra, 120 for Quebec
    T_bc = T_bc.reshape((extraDim, 32, 24))  

    plumeType = np.max(plumeMask)
    cleanRad = rad[:,plumeMask == plumeType]
    cleanVel = vel[:,plumeMask == plumeType]
    seasonTime = np.linspace(0,cycle/86400,nt)
    if(plumeType == 2):
        runoff = dy*cleanRad[:,0]*cleanVel[:,0]
    else: #semi-rad
        runoff = np.pi*cleanRad**2*cleanVel
    # print(f'seasonTime {seasonTime}, nt {nt}')
    # print(f'runoff {runoff}')
    # print(f'plume type {plumeType}')
    # print(np.shape(seasonTime))
    # print(np.shape(runoff))
    f = interpolate.interp1d(seasonTime, runoff,fill_value='0')
    f_temp = interpolate.interp1d(seasonTime, np.mean(T_bc,axis=(1,2)),fill_value='0')
    if(args.sgd != None):
        TempTime = f_temp(timeTime % (cycle/86400))
        timeTime = f(timeTime % (cycle/86400))
        timeLabel = 'SGD [$\\mathrm{m^3/s}$]'
    N = args.windowMean #window over to smooth in time points (days)
    #Pad edges of BTime to have smooth edges and work for short time series
    padL = int(N/2)
    padR = int((N-1)/2)  #This ensures odds work OK
    BTimeSmooth = np.convolve(np.concatenate((BTime[0]*np.ones(padL),BTime,BTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    # timeTime = (timeTime - timeTime[0])
    forceTimeRough = BTime
    forceTime = BTimeSmooth
    bFluxSmooth = np.convolve(np.concatenate((bFluxTime[0]*np.ones(padL),bFluxTime,bFluxTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    endFluxSmooth = np.convolve(np.concatenate((endFluxTime[0]*np.ones(padL),endFluxTime,endFluxTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')

    forceLabel = 'Melt Rate [m/d]'
    forceString = 'melt'
    # plt.suptitle('Mélange over Time')
    tempColor = 'xkcd:magenta'
    
    # ax1.plot(timeTime + years[i]*365,f((timeTime  + args.timeRange[0] + years[i]*365) % (cycle/86400)),color='xkcd:azure',linestyle=lStyle[i],alpha=.5)
    ax1.plot(timeTime,f((timeTime  + args.timeRange[0] + years[i]*365) % (cycle/86400)),
                color='xkcd:azure',
                linestyle=lStyle[i],
                alpha=.5,
                label=f'Year {years[i]}',
                linewidth=lWidth)
    ax1.set_ylabel('SGD [$\\mathrm{m^3/s}$]',color='xkcd:azure')
    ax1.tick_params(axis='y',labelcolor='xkcd:azure')
    ax1.set_xlabel(timeLabel)
    ax1.grid(alpha=.5)
    ax1.set_xticks([0,60,121,182,244,305,365])
    ax1.set_xticklabels(['Jan','Mar','May','Jul','Sep','Nov','Jan'])
    if(i == 1):
        prev_ylim = ax1.get_ylim()
        ax1.plot([DOY_to_plot,DOY_to_plot],prev_ylim,'--',color='black',alpha=.2,label = f'Day {DOY_to_plot}')
        ax1.set_ylim(prev_ylim)

    if(False): #plot temp or not
        if(i == 0):
            ax1_2 = ax1.twinx()
        ax1_2.plot(timeTime,f_temp(((timeTime  + args.timeRange[0] + years[i]*365) % (cycle/86400))),color=tempColor,linestyle=lStyle[i])
        # ax1.scatter(timeTime[0],forceTime[0],s=50,marker='*',color='black')
        
        ax1_2.set_ylabel('Depth avg shelf temp [C]',color=tempColor)
        ax1_2.tick_params(axis='y',labelcolor=tempColor)

    if(args.labels != None):
        ax1.legend(loc=args.legLoc)
    xString =''
    ax1.legend()

    drivingVariable = timeTime
    drivingLabel = timeLabel
    xString = 'Time'
    if(args.sgd != None):
        xString = 'SGD'
    roughLT = lengthTime
    roughHT = H0Time
    lengthTime = np.convolve(np.concatenate((lengthTime[0]*np.ones(padL),lengthTime,lengthTime[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    H0Time = np.convolve(np.concatenate((H0Time[0]*np.ones(padL),H0Time,H0Time[-1]*np.ones(padR))), np.ones(N)/N, mode='valid')
    ax2.plot(timeTime,forceTime/365.25,color=meltColors[0],linestyle=lStyle[i],linewidth=lWidth)
    # ax2.plot(timeTime,forceTimeRough/365.25,color=meltColors[0],linestyle='-',alpha=.25)
    ax2.set_ylabel('Melt rate, $\\bar b$ [m/d]',color=meltColors[0])
    ax2.tick_params(axis='y',labelcolor=meltColors[0])
    print(f'Year {years[i]} avg melt rate is {np.mean(forceTimeRough)/365.25} [m/d]')

    ax2.grid(alpha=.5)
    # ax2.set_xticks([0,73,73*2,73*3,73*4,73*5])
    ax2.set_xticks([0,60,121,182,244,305,365])
    ax2.set_xticklabels(['Jan','Mar','May','Jul','Sep','Nov','Jan'])
    # sca=ax2.scatter(BTime,lengthTime,s=None,c=timeTime,cmap='cividis')
    # cbar=plt.colorbar(sca)
    # cbar.set_label('Iteration')
    if(i == (len(years)-1)): #gate patch addition
        init_xlim = ax2.get_xlim()
        init_ylim = ax2.get_ylim()
        p1 = patches.Rectangle([0   ,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[0],alpha=.2,edgecolor=None,zorder=-1)
        p2 = patches.Rectangle([73  ,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[1],alpha=.2,edgecolor=None,zorder=-1)
        p3 = patches.Rectangle([73*2,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[2],alpha=.2,edgecolor=None,zorder=-1)
        p4 = patches.Rectangle([73*3,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[3],alpha=.2,edgecolor=None,zorder=-1)
        p5 = patches.Rectangle([73*4,init_ylim[0]],73,init_ylim[1]-init_ylim[0],facecolor=colorList[4],alpha=.2,edgecolor=None,zorder=-1)
        ax2.add_patch(p1)
        ax2.add_patch(p2)
        ax2.add_patch(p3)
        ax2.add_patch(p4)
        ax2.add_patch(p5)

    if(True): #gate for 2nd plot
        if(i == 0):
            ax2_2 = ax2.twinx()
        # ax2_2.plot(drivingVariable,roughLT,color='xkcd:blueberry',alpha=0.25,linewidth=lWidth)
        ax2_2.plot(drivingVariable,lengthTime,color='xkcd:blueberry',linestyle=lStyle[i],linewidth=lWidth)
        ax2_2.set_ylabel('Mélange length, $L$    [m]',color='xkcd:blueberry')
        ax2_2.tick_params(axis='y',labelcolor='xkcd:blueberry')
    ax2.set_xlabel(drivingLabel)
    if(i == 0):
        ax2.plot([DOY_to_plot,DOY_to_plot],[0,1.6],'--',color='black',alpha=.2)
    ax2.set_ylim([0,1.6])
    # ax2.set_title('Size')

    if(i == 0): #Plot blanks for legend 
        ax3.plot([],[],color='xkcd:sand',linestyle='-',label='$Q_{in}$')   
        ax3.plot([],[],color='xkcd:apple',linestyle='-',label='$Q_{out}$')  
        ax3.plot([],[],color='xkcd:red',linestyle='-',label='$Q_{melt}$')  
        ax3.plot([],[],color='xkcd:sky blue',linestyle='-',label='$Q_{net}$',alpha=.5)
        # ax3.plot([],[],color='xkcd:indigo',linestyle='-',label='Volume')
    ax3.plot(timeTime,inFluxTime,color='xkcd:sand',linestyle=lStyle[i],linewidth=lWidth)   
    # ax3.plot(timeTime,endFluxTime,color='xkcd:apple',linestyle=lStyle[i],alpha = .25)  
    # ax3.plot(timeTime,bFluxTime,color='xkcd:red',linestyle=lStyle[i],alpha = .25)   
    ax3.plot(timeTime,endFluxSmooth,color='xkcd:apple',linestyle=lStyle[i],linewidth=lWidth)  
    ax3.plot(timeTime,bFluxSmooth,color='xkcd:red',linestyle=lStyle[i],linewidth=lWidth)
    ax3.plot(timeTime,inFluxTime-endFluxSmooth-bFluxSmooth,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.75,linewidth=lWidth)
    # ax3.plot(timeTime,inFluxTime-endFluxTime-bFluxTime,color='xkcd:sky blue',linestyle=lStyle[i],alpha=.25,linewidth=lWidth)
    # ax3.plot([timeTime[0],timeTime[-1]],[0,0],color = 'gray',alpha=.5,linestyle='--')    
    ax3.set_ylabel('Volume flux [$\\mathrm{m^2/y}$]')
    ax3.set_xlabel(timeLabel)
    if(i == 0):
        ax3.plot([DOY_to_plot,DOY_to_plot],[-6.2e6, 7.9e6],'--',color='black',alpha=.2)
    ax3.set_ylim([-6.2e6, 7.9e6])
    ax3.grid(alpha=.5)
    # ax3.set_title('')
    # ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=4)
    ax3.legend(fontsize='medium',ncols = 4,framealpha=.2,loc=9)
    ax3.set_xticks([0,60,121,182,244,305,365])
    ax3.set_xticklabels(['Jan','Mar','May','Jul','Sep','Nov','Jan'])


    ## Melt shapes
    if(i == 0): #Plot blanks for legend 
        ax4.plot([],[],color=colorList[0],linestyle=lStyle[0],label='Late Winter')   
        ax4.plot([],[],color=colorList[1],linestyle=lStyle[0],label='Spring')  
        ax4.plot([],[],color=colorList[2],linestyle=lStyle[0],label='Summer')  
        ax4.plot([],[],color=colorList[3],linestyle=lStyle[0],label='Autumn')
        ax4.plot([],[],color=colorList[4],linestyle=lStyle[0],label='Early Winter')
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[  0: 73,:],axis=0),alpha = 1,linestyle=lStyle[i],color=colorList[0],linewidth=lWidth)
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[ 73:164,:],axis=0),alpha = 1,linestyle=lStyle[i],color=colorList[1],linewidth=lWidth)
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[164:219,:],axis=0),alpha = 1,linestyle=lStyle[i],color=colorList[2],linewidth=lWidth)
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[219:292,:],axis=0),alpha = 1,linestyle=lStyle[i],color=colorList[3],linewidth=lWidth)
    ax4.plot(np.linspace(0,1,20),np.mean(BShapeTime[292:365,:],axis=0),alpha = 1,linestyle=lStyle[i],color=colorList[4],linewidth=lWidth)
    ax4.grid(alpha=.5)
    ax4.legend(loc=9,fontsize='small')
    # if(i == 0):
    #     ax4_2=ax4.twinx()
    # ax4_2.plot(drivingVariable,roughG,color='xkcd:lilac',alpha=0.25)
    # ax4_2.plot(drivingVariable,gTime,color='xkcd:lilac',linestyle=lStyle[i])
    # ax4_2.scatter(drivingVariable[0],gTime[0],s=50,marker='*',color='black')
    # ax4_2.set_ylabel('fluidity $g^{\\prime}$ [$yr^{-1}$]   ',color='xkcd:lilac')
    # ax4_2.tick_params(axis='y',labelcolor='xkcd:lilac')
    # ax4.set_title('Avg Speed and Fluidity $g^{\\prime}$')
    ax4.set_xlabel('Fraction along mélange $\\xi$ [ ]')
    ax4.set_ylabel('Melt rate, $b$ [m/d]')
    ax4.set_ylim([0,1.6])
    if(i == 0):    
        ## Do temp, U plotting first, only for first year
        ySlice = np.argmin(np.abs(y[:,0] - args.yCrossSection))
        print('cross section is y =', y[ySlice,0], 'index', ySlice)

        shift = DOY_to_plot # add this many to day 0 of viewed year
        data = mds.rdmds("results/dynDiag", int((np.max([args.timeRange[0],1]) + shift + 365*years[i])*86400.0/dt))
        data[:,data[1,:,:,:] == 0] = np.nan #if salt = 0 its a wall, nan it out for averging later

        if(args.shadow > 0): #enable berg shadows here
            dataBergs = mds.rdmds("results/BRGFlx",int((np.max([args.timeRange[0],1]) + shift + 365*years[i])*86400.0/dt))
            if(dataBergs.shape[0] < 6):
                dataBergPlot = dataBergs[0,:,:,:].copy()
                dataBergPlot[dataBergPlot != 0] = .2
                dataBergPlot[data[0,:,:,:] < -1] = .2 #attemp to recapture freezing limit ice
                dataBergPlot[:,np.nansum(dataBergs[0,:,:,:],axis=0) == 0] = 0 #attemp to recapture freezing limit ice
            else:
                dataBergPlot = dataBergs[5,:,:,:]
            dataBergPlot[dataBergPlot == 0] = 1 #np.nan
        ## Temp
        lvl = np.linspace(-2, 4.0, 31)
        cm = 'RdYlBu_r'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            plotData = np.squeeze(np.nanmean(data[0, :, :, :],axis=1))
        cp = ax5.contourf(
                        np.squeeze(x[ySlice,:]),
                        np.squeeze(z),
                        plotData,
                        lvl,
                        extend="both",
                        cmap=cm,
                        )
        if(True and showQuiver):
            zSkip = 1
            ySkip = 5
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                u = np.squeeze(np.nanmean(data[2, ::zSkip, :, ::ySkip],axis=1))
                w = np.squeeze(np.nanmean(data[3, ::zSkip, :, ::ySkip],axis=1))
            qv = ax5.quiver(
                x[ySlice,::ySkip],
                np.squeeze(z[::zSkip]),
                u,#/np.sqrt(u**2 + w**2 + 1e-12),
                w,#/np.sqrt(u**2 + w**2 + 1e-12),
                alpha=.4,
                angles='xy',
                pivot='tail',
                #width = .004, #width of line
                #scale = 40
                )
            plt.quiverkey(qv,1.08,.95,.1,'10 cm/s')
        if(args.shadow > 0):
                    # plt.plot(x[ySlice,:],-np.max(maxDepth,axis=0),color='gray',linestyle='dotted')
                    cp2 = ax5.contourf(
                        x[ySlice,:],
                        np.squeeze(z),
                        np.squeeze(1-dataBergPlot[:, ySlice, :])*100,
                        [5,20,40,60,80],
                        extend="max",
                        alpha=.35,
                        cmap='cmo.gray_r')
        if(showDensity):
            salt = np.squeeze(data[1,:,ySlice,:])
            pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
            CT = gsw.CT_from_t(salt, data[0,:,ySlice,:], pressure)
            density = gsw.rho(salt, CT, 0) - 1000 #in-stu density less 1000
            densityLevels = np.linspace(22,28,31)
            cc = ax5.contour(
                np.squeeze(x[ySlice,:]),
                np.squeeze(z),
                np.squeeze(density),
                densityLevels,
                colors='black',
                linewidths=0.5,
                alpha=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        if(args.xlimit != None):
            ax5.set_xlim([0, args.xlimit])
        if(i == 0):
            divider = make_axes_locatable(ax5)
            ax5_cb = divider.append_axes("right", size="2%", pad="15%")
            ax5_cb2 = divider.append_axes("right", size="2%", pad="10%")
            fig = ax5.get_figure()
            fig.add_axes(ax5_cb)
            fig.add_axes(ax5_cb2)
        ax5.set_xticks([0,5e3,10e3,15e3,20e3,25e3,30e3])
        ax5.set_xticklabels([0,5,10,15,20,25,30])
        ax5.yaxis.tick_right()
        ax5.yaxis.set_label_position("right")
        ax5.set_xlabel('Along fjord [km]')
        ax5.set_ylabel('Depth [m]')
        ax5.set_title(f"Fjord on day {shift:.0f} of year {years[0]}")
        cbar = plt.colorbar(cp,cax=ax5_cb,ticks=[-1, 0, 1,2,3,4], orientation="vertical",fraction=0.06,format='%.0f')
        cbar.set_label('Temperature [C]')
        cbar2 = plt.colorbar(cp2,cax=ax5_cb2, orientation="vertical",fraction=0.06,format='%.0f')
        cbar2.set_label('Ice fraction $\\varphi$ [%]')
        init_xlim = ax5.get_xlim()
        init_ylim = ax5.get_ylim()
        pat = patches.Rectangle([-3e3,init_ylim[0]],3e3,init_ylim[1]-init_ylim[0],facecolor='xkcd:ice blue',alpha=.9,clip_on=False)
        ax5.add_patch(pat)
        ax5.text(-1.5e3,(init_ylim[1]+init_ylim[0])/2,'Glacier',clip_on = False,rotation='vertical',va='center',ha='center',size='large')
    #melange model cross section
        mapName = 'cmo.phase'
        colorMapForMelange= colormaps[mapName]
        n_colors = len(toIterate)
        melange_Colors = colorMapForMelange(np.linspace(0,1,n_colors+1))
        
        for j in [shift]:
            p = ax5.plot(np.append(X_fullTime[j,:],X_fullTime[j,::-1]),np.append(-917/1020*H_fullTime[j,:],(1-917/1020)*H_fullTime[j,::-1]),
                color='xkcd:black', # colorList[int(j/73)] or melange_Colors[j]
                linewidth=2,
                alpha=.6,
                linestyle = '--',
                label=f'Continuum Mélange Profile')
        ax5.legend(loc=4)

    if(ax6 != None):
        if(i == 0):
            divider = make_axes_locatable(ax6)
            ax6_cb = divider.append_axes("right", size="2%", pad="2%")
            fig = ax6.get_figure()
            fig.add_axes(ax6_cb)
        ax6.plot([],[],color='gray',linestyle=lStyle[i],label=f'Year {years[i]}')  
        #Thick lines
        for j in toIterate[::73]:
            p = ax6.plot(np.append(X_fullTime[j,:],X_fullTime[j,::-1]),np.append(-917/1020*H_fullTime[j,:],(1-917/1020)*H_fullTime[j,::-1]),
                color=melange_Colors[j], # colorList[int(j/73)] or melange_Colors[j]
                linewidth=2, linestyle = lStyle[i],
                )
        ## thin lines
        # for j in toIterate[::28]:
        #     p = ax6.plot(np.append(X_fullTime[j,:],X_fullTime[j,::-1]),np.append(-917/1020*H_fullTime[j,:],(1-917/1020)*H_fullTime[j,::-1]),
        #         color=melange_Colors[j], # colorList[int(j/73)] or melange_Colors[j]
        #         linewidth=1, linestyle = lStyle[i],alpha=.3
        #         )
        if(i == len(years)-1):
            sc = ax6.scatter([],[],c=[],vmin=0,vmax=365,cmap=mapName,alpha=1)
            # cbar = plt.colorbar(sc,cax=ax6_cb,ticks=[0,31,59,90,120,151,181,212,243,273,304,334,365])
            # cbar.ax.set_yticklabels(['','Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
            cbar = plt.colorbar(sc,cax=ax6_cb,ticks=[0,73,73*2,73*3,73*4,73*5])
            # cbar.ax.set_yticklabels(['0\n(Jan 1)','73\n(Mar 14)', '146\n(May 26)', '219\n(Aug 7)', '292\n(Oct 19)', '365\n(Dec 31)'])
            cbar.set_label('Day of Year')
            ax6.legend(loc=4)
        ax6.set_ylabel('Depth [m]')
        ax6.set_xlabel('Length, $L$ [m]')

# plt.tight_layout()

## label and save
for label,ax in zip(figLabels,axes):
    if(ax != None):
        ax.text(
            ax.get_xlim()[0], ax.get_ylim()[1], label,
            fontsize='x-large', va='bottom',ha='right', fontfamily='sans serif')

dirStr = ''
if(os.path.isdir('figs')):
    dirStr = 'figs/'
multiStr = ''
timeString = ''
if(i > 0):
    multiStr = f'Multi{i+1}'
if(args.timeRange != None):
    multiStr = f'{multiStr}_T_{args.timeRange[0]}_{args.timeRange[1]}'

plt.savefig('%sprettyPlotBig%s%s%s.png' %(dirStr,forceString,xString,multiStr),format='png',dpi=args.dpi)
if(args.silent == 0):
    plt.show()
plt.close()
