import os
import sys
import numpy as np
import fileinput
from MITgcmutils import mds
import matplotlib.pyplot as plt
import cmocean
import glob,pickle
import warnings
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
from scipy.interpolate import interp1d

def write_bin(fname, data):
    # print(fname + " " + str(np.shape(data)))
    data.astype(">f8").tofile('input/'+fname)

#Get most recent melt file
maxStep = 0
for file in os.listdir('results'):
    # print(file)
    if "BRGFlx.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])

dataMITgcm = mds.rdmds("results/BRGFlx", maxStep)

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = -1*np.squeeze(mds.rdmds("results/RC"))

ny = np.squeeze(np.shape(x)[0]) #x,y are loaded as [ny,nx] arrays from MITgcm
nx = np.squeeze(np.shape(x)[1])
nz = np.squeeze(np.shape(z)[0])
dz = np.load('input/dz.npy')
deltaX = x[0,0]*2
deltaY = y[0,0]*2
maxBergs = 500

#icebergCoverLamba is enforced to remain constant left of here
icebergCoverLambda = .80
icebergRefreshingGate = 3

#icebergs are destroyed when they move past here
icebergRightHandGate = int(nx*.8)
couplingTimeStep = 1*24*3600 #[s]

files = sorted(glob.glob('couplingResults/MITgcmRun_[0-9][0-9][0-9][0-9][0-9].pickle'))
# print(files)
with open(files[-1], 'rb') as file:
    data = pickle.load(file)
    file.close()
    index = int(str(file).split(".")[1].split("_")[1])

U_glaciome1d = (data.U+data.Ut-data.Uc)/(24*3600*365) #raw is m/yr, convert to m/s
H = np.concatenate(([data.H0], data.H, [data.HL]))
X_glaciome1d = data.X + deltaX #meters, shift to align with front of MITgcm melange
X_ = np.concatenate(([data.X[0]],data.X_,[data.X[-1]])) + deltaX
glaciome_L = data.L + deltaX
speedLookup = interp1d(X_glaciome1d,U_glaciome1d,bounds_error=False,fill_value=[0])

# plt.figure
# plt.plot(X_glaciome1d,U_glaciome1d)
# plt.show()

#2D Fields
plumeMask  = np.fromfile('input/plumeMask.bin', dtype='>f8')
plumeMask  = plumeMask.reshape((ny,nx))
bergWidths = np.fromfile('input/icebergs_widths.bin', dtype='>f8')
bergWidths = bergWidths.reshape((maxBergs,ny,nx))
bergLength = np.fromfile('input/icebergs_length.bin', dtype='>f8')
bergLength = bergLength.reshape((maxBergs,ny,nx))
bergDepths = np.fromfile('input/icebergs_depths.bin', dtype='>f8')
bergDepths = bergDepths.reshape((maxBergs,ny,nx))
bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
bergMask = bergMask.reshape((ny,nx))
bergMaskNums = np.fromfile('input/bergMaskNums.bin', dtype='>f8')
bergMaskNums = bergMaskNums.reshape((ny,nx))
bergsPerCell = np.fromfile('input/numBergsPerCell.bin', dtype='>f8')
bergsPerCell = bergsPerCell.reshape((ny,nx)).astype(np.int64) 
driftMask = np.fromfile('input/driftMask.bin', dtype='>f8')
driftMask = driftMask.reshape((ny,nx))
barrierMask = np.fromfile('input/barrierMask.bin', dtype='>f8')
barrierMask = barrierMask.reshape((ny,nx))
meltMask = np.fromfile('input/meltMask.bin', dtype='>f8')
meltMask = meltMask.reshape((ny,nx))
#3D Fields
totalBergArea = np.fromfile('input/totalBergArea.bin', dtype='>f8')
totalBergArea = totalBergArea.reshape((nz,ny,nx))
openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
openFrac = openFrac.reshape((nz,ny,nx))

if(np.min([np.min(bergDepths),np.min(bergLength),np.min(bergWidths)]) != 0):
    print('negative iceberg dimentions values found')

for i in range(nx):
    for j in range(ny):
        if len(bergDepths[bergDepths[:,j,i] > 0,j,i]) != bergsPerCell[j,i]:
            print('depth error',j,i)
        if len(bergLength[bergLength[:,j,i] > 0,j,i]) != bergsPerCell[j,i]:
            print('length error',j,i)
        if len(bergWidths[bergWidths[:,j,i] > 0,j,i]) != bergsPerCell[j,i]:
            print('width error',j,i)
        if bergMask[j,i] and (barrierMask[j,i] + meltMask[j,i]) != 2:
            print('mask error',j,i)
        if np.nanmax(openFrac[:,j,i]) > 1 or np.nanmin(openFrac[:,j,i]) < 0:
            print('open frac error',j,i)



varphi = 1-openFrac
varphi[varphi == 1] = np.nan
depthHelper = bergDepths.copy()
depthHelper[depthHelper == 0] = np.nan
effectiveDepth = np.nanmean(np.nansum((1-openFrac) * dz[:,None,None],axis=0),axis=0)
# effectiveDepth[1] = np.nanmean(np.nansum(1-openFrac[:,,]))
print('Max lambda is:',np.nanmax(varphi[0,:,:]))

print(openFrac[:,6,1])
print(effectiveDepth[1:2])

# print(bergDepths[:,6,1],bergWidths[:,6,1],bergLength[:,6,1])

plt.figure(3)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    pc = plt.pcolormesh(x[0,:],-z,np.nanmean(varphi[:,1:-1,:],axis=1),cmap='gray_r')
cbar = plt.colorbar(pc)
depthHelper = bergDepths.copy()
depthHelper[depthHelper == 0] = np.nan
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    plt.plot(x[0,:],-effectiveDepth,color='xkcd:red',linestyle='--')
    plt.plot(x[0,:],-np.nanmean(depthHelper[:,1:-1,:],axis=(0,1)),color='xkcd:orange',linestyle='--')
plt.plot(X_,-H*917/1030,color='xkcd:green',linestyle='--')
plt.plot([glaciome_L,glaciome_L],[20,-600],linestyle='--',color='cyan')
plt.suptitle('$\\varphi$')
plt.ylabel('Depth [m]')
plt.xlabel('Along fjord [m]')
plt.ylim([-600, 0])
cbar.set_label('ice fraction')


plt.figure(1)
plt.subplot(221)
plt.plot(x[0,:],np.sum(bergsPerCell[1:-1,:],axis=0))
plt.title('berg count')
plt.ylabel('Count [ ]')
plt.xlabel('Along fjord [m]')

plt.subplot(222)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    pc = plt.pcolormesh(x[0,:],y[:,0],varphi[0,:,:],cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('$\\lambda$')
plt.ylabel('Across fjord [m]')
plt.xlabel('Along fjord [m]')
cbar.set_label('ice fraction')
plt.clim([0,1])

plt.subplot(223)
# pc = plt.pcolormesh(x[0,:],y[:,0],meltMask,alpha=.2,cmap='bone_r')
# pc = plt.pcolormesh(x[0,:],y[:,0],bergMask,alpha=.2,cmap='bone_r')
pc = plt.pcolormesh(x[0,:],y[:,0],barrierMask,alpha=.2,cmap='bone_r')
# pc = plt.pcolormesh(x[0,:],y[:,0],driftMask,alpha=.2,cmap='bone_r')
# plt.colorbar(pc)
plt.suptitle('Masks')
plt.ylabel('Depth [m]')
plt.xlabel('Along fjord [m]')
# cbar.set_label('ice fraction')

plt.subplot(224)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    pc = plt.pcolormesh(x[0,:],y[:,0],varphi[0,:,:],cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('$\\lambda$')
plt.ylabel('Across fjord [m]')
plt.xlabel('Along fjord [m]')
cbar.set_label('ice fraction')
plt.clim([0,1])
plt.tight_layout() 

bergDepths[bergDepths == 0] = np.nan
bergWidths[bergWidths == 0] = np.nan
bergLength[bergLength == 0] = np.nan

plt.figure(2)
plt.subplot(311)
plt.hist(bergWidths.ravel(),bins=100)
plt.xlabel('Widths [m]')
plt.yscale('log')

plt.subplot(312)
plt.hist(bergLength.ravel(),bins=100)
plt.xlabel('Lengths [m]')
plt.yscale('log')

plt.subplot(313)
plt.hist(bergDepths.ravel(),bins=100)
plt.xlabel('Depths [m]')

plt.tight_layout() 
plt.show()
# plt.savefig('figs/advectBergs.png',format='png',dpi=150)
plt.close()

