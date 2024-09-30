#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import matplotlib.pylab as plt
import pickle
import importlib
import shutil

from scipy.interpolate import make_interp_spline

import sys
sys.path.append('/Users/psummers8/Documents/MITgcm/MITgcm/elizaScripts/main_scripts')

import build_domain_funcs as build_domain 
import run_config_funcs as rcf # import helpter functions


# ## Main run configuration

# In[2]:


# set high level run configurations

run_config = {}
run_config['ncpus_xy'] = [1, 1] # cpu distribution in the x and y directions
run_config['run_name'] = 'BergSolid500'
run_config['ndays'] = 20 # simulaton time (days)
run_config['test'] = True # if True, run_config['nyrs'] will be shortened to a few time steps

run_config['horiz_res_m'] = 500 # horizontal grid spacing (km)
run_config['Lx_m'] = 25000 # domain size in x (m)
run_config['Ly_m'] = 6000 # domain size in y (m)
# NOTE: the number of grid points in x and y should be multiples of the number of cpus.

run_config['evolve_salt'] = False
run_config['use_GMRedi'] = False # should be set to false for eddy permitting resolutions
run_config['periodic_forcing'] = False # note: the code is not yet set up to handle time-dependent forcing

MITgcm_release = 'MITgcm-checkpoint69a' #Sept 2024 release
#MITgcm_code_dir = os.path.join(group_home_dir, 'shared/mitgcm_releases', MITgcm_release)

# you probably don't need to touch this
run_config['use_MPI'] = True # for multi-processing
run_config['lf'] = '\r\n' # linebreak characters 
run_config['exps_dir'] = os.path.join('/Users/psummers8/Documents/MITgcm/MITgcm/experiments') 
run_config['run_dir'] = os.path.join(run_config['exps_dir'], run_config['run_name'])


# In[3]:


run_config


# ## Generate new experiment directory and copy over defaults

# In[4]:


run_config


# In[5]:


# create experimentary directory on SCRATCH and copy over default configuration
# NOTE: this step does not overwrite existing directories. 
run_subdir_list = ['build', 'code', 'input', 'results']
for subdir in run_subdir_list:
    run_config['%s_dir'% subdir] = os.path.join(run_config['run_dir'], subdir)
    os.makedirs(run_config['%s_dir'% subdir], exist_ok=True)
     
# copy over defaults
default_dirs = os.listdir('/Users/psummers8/Documents/MITgcm/MITgcm/DEFAULT/')
for dir00 in default_dirs:
    print(dir00)
    if dir00.startswith('.'):
        continue
        
    default_dir = '/Users/psummers8/Documents/MITgcm/MITgcm/DEFAULT/%s/'%dir00
    default_files = os.listdir(default_dir)
    dst_dir = os.path.join(run_config['run_dir'], dir00)
    
    for file in default_files:

        if file.startswith('.'):
            continue
        else:
            src_fpath = os.path.join(default_dir, file)
            shutil.copy2(src_fpath, dst_dir)
            #print(src_fpath, '>', dst_dir)


# In[6]:


# just to see what the default files
os.listdir('/Users/psummers8/Documents/MITgcm/MITgcm/DEFAULT/code')


# In[7]:


print(run_config['run_dir'])
print(os.listdir(run_config['run_dir']))


# In[8]:


# create new analysis sub-dir in your home directory
analysis_dir = '../analysis/%s'%run_config['run_name']
os.makedirs(analysis_dir, exist_ok=True)
print(analysis_dir)
os.getcwd()


# ## Domain  and grid parameters

# In[9]:


secsInDay = 24*60*60
secsInYear = 365*secsInDay


# set domain size
domain_params = {}
domain_params['Lx'] = run_config['Lx_m'] # domain size in x (m)
domain_params['Ly'] = run_config['Ly_m'] # domain size in y (m)
domain_params['L_sponge'] = 5000 # width of eastern sponge layer (m)
domain_params['H'] = 200 # max domain depth (m)


# In[10]:


#---grid parameters ---# 

# NOTE: the only thing you may need to change here is the number of z-grid points)
grid_params = {}

domain_params['Lx']/(run_config['horiz_res_m'])
domain_params['Ly']/(run_config['horiz_res_m'])

grid_params['nSx'] = 1 # num of tiles per processor in x-direction
grid_params['nSy'] = 1 # num of tiles per processor in y-direction
grid_params['nTx'] = 1 # num of threads per processor in x-direction
grid_params['nTy'] = 1 # num of threads per processor in y-direction
grid_params['OLx'] = 3 # num of overlapping x-gridpoints per tile
grid_params['OLy'] = 3 # num of overlapping y-gridpoints per tile
grid_params['Nr'] = 20 # num of z-grid points

grid_params['nPx'] = run_config['ncpus_xy'][0] #num of processors in x-direction
grid_params['nPy'] = run_config['ncpus_xy'][1] #num of processors in x-direction

# grid_params['nSx'] = domain_params['Lx']/(run_config['horiz_res_m']) # num of x points in sub grid
# grid_params['nSy'] = domain_params['Ly']/(run_config['horiz_res_m']) # num of y points in sub grid

# grid_params['Nx'] = grid_params['sNx'] * grid_params['nSx'] * grid_params['nPx']
# grid_params['Ny'] = grid_params['sNy'] * grid_params['nSy'] * grid_params['nPy']

grid_params['Nx'] = domain_params['Lx']/(run_config['horiz_res_m']) # num of x points
grid_params['Ny'] = domain_params['Ly']/(run_config['horiz_res_m']) # num of y points

print("Nx: %s" %grid_params['Nx'])
print("Ny: %s" %grid_params['Ny'])

grid_params['sNx'] = grid_params['Nx']/grid_params['nPx']#num of x-gridpoints per tile
grid_params['sNy'] = grid_params['Ny']/grid_params['nPy'] #num of y-gridpoints per tile

print("sNx: %s" %grid_params['sNx'])
print("sNy: %s" %grid_params['sNy'])

# NOTE: sNx and sNy should be whole numbers/integers. As long we keep the horizontal resolution,
# domain dimesions, and number of cpus to be multiples of five, we should be ok. 

for key, param  in grid_params.items():
    assert param%1==0, "grid parameter needs to be an integer"
    grid_params[key] = int(param)
    
print(grid_params)
#run_config['grid_params'] = grid_params


# In[11]:


# grid_params cont'd
grid_params['usingCartesianGrid'] = True
grid_params['usingSphericalPolarGrid'] = False 

# horizontal grid spacing
grid_params['delX'] = (domain_params['Lx']/grid_params['Nx'])*np.ones(grid_params['Nx'])
grid_params['delY'] = (domain_params['Ly']/grid_params['Ny'])*np.ones(grid_params['Ny'])


# vertical grid spacing 
# spacing increases with depth---can be modified
zidx = np.arange(1, grid_params['Nr']+1)
aa = 10
bb = 10
dz1 = 2*domain_params['H']/grid_params['Nr']/(aa+1)
dz2 = aa*dz1
dz = dz1 + ((dz2-dz1)/2)*(1+np.tanh((zidx-((grid_params['Nr']+1)/2))/aa))
zz1 = np.append([0], np.cumsum(dz))
zz = -(zz1[:-1] + np.diff(zz1)/2) # layer midpoints

grid_params['delZ'] = dz
grid_params['hFacMinDr'] = dz.min()


# In[ ]:





# 
# ## Physical parameters

# In[25]:


#---physical params---#

params01 = {} 

# physical constants
g = 9.81 # acc. due to gravity (m/s**2)
Omega = 2*np.pi*366/365/86400 # planetary rotation rate 
Rp = 6400*1000 # planetary radius (m)
lat_min = -70 # latitude at southern boundary (degrees)
#f0 = 2*Omega*np.sin(np.deg2rad(lat_min)) # coriolis param (1/s)
#beta = (2*Omega*np.cos(np.deg2rad(lat_min))/Rp) # beta param


# momentum scheme
params01['vectorInvariantMomentum'] = True

#Note: here and elsewhere, we need to be explicit about floats vs ints. E.g., use 12.0 to represent float and
# 12 for int

 #Tref = 40*0.,
 #Sref = 40*34.,
 #viscAz=0.00001,
 #viscC2smag=2.2,
 #diffKzT=0.00001,
 #diffKhT=30,
 #diffKzS=0.00001,
 #diffKhS=30,
 #no_slip_sides=.FALSE.,
 #no_slip_bottom=.TRUE.,
 #rigidLid=.FALSE.,
 #implicitFreeSurface=.TRUE.,
 #staggerTimeStep=.TRUE.,
readBinaryPrec=64
 #eosType='JMD95Z', 
 #tempAdvScheme=33,
 #saltAdvScheme=33,
f0=0.
beta=0.E-11
 #selectAddFluid=1,


# viscosity parameters
params01['viscAz'] = 0.00001 # Vertical viscosity
params01['viscC2smag'] = 2.2 # ??? viscosity

# advection and time stepping
params01['tempAdvScheme'] = 33 # needs to be int
params01['saltAdvScheme'] = 33 # needs to be int
#params01['tempStepping'] = True
#params01['saltStepping'] = run_config['evolve_salt']
params01['staggerTimeStep'] = True

# diffusivity
params01['diffKzT'] = 0.00001 # Vertical temp diffusion 
params01['diffKhS'] = 30 # ???diffusion
params01['diffKhT'] = 30 # Horizontal temp diffusion

# equation of state
params01['eosType'] = 'JMD95Z'
params01['Tref'] = np.ones(grid_params['Nr'])*0. #ref temp
params01['Sref'] = np.ones(grid_params['Nr'])*34. #ref salt

# boundary conditions
params01['no_slip_sides'] = False
params01['no_slip_bottom'] = True
params01['rigidLid'] = False
params01['implicitFreeSurface'] = True

# physical parameters
params01['f0'] = f0
params01['beta'] = beta
params01['gravity'] = g

# misc
params01['selectAddFluid'] = 1


# AS: Don't allow partial cell height to fall below min grid spacing
#params01['hFacMinDr'] = np.min(dz)


# ## Check for numericl stability?

# ## Numeric solvers and I/O controls

# In[27]:


# numeric solver parameters 

params02 = {}
params02['cg2dMaxIters'] = 300
params02['cg2dTargetResidual'] = 1e-13
params02['cg3dMaxIters'] = 20
params02['cg3dTargetResidual'] = 1e-8

# time stepping parameters 
params03 = {}
params03['nIter0'] = 0
#params03['endTime'] = 864000.0
deltaT = 100.0
params03['abEps'] = 0.1

#if run_config['testing']:
    
params03['chkptFreq'] = 0.0
params03['pChkptFreq'] = 0.0
params03['taveFreq'] = 0.0
params03['dumpFreq'] = 86400.0
params03['taveFreq'] = 0.0
params03['monitorFreq'] = 86400.0
params03['monitorSelect'] = 2


params03['periodicExternalForcing'] = False
params03['ExternForcingPeriod'] = 100.0
params03['ExternForcingCycle'] = 1000.0 


# In[28]:


if run_config['test']:
    nTimeSteps = 10
else:
    nTimeSteps = np.ceil(run_config['ndays']*secsInDay/detlaT)

simTimeAct = nTimeSteps*deltaT

params03['endTime'] = int(params03['nIter0']*deltaT+simTimeAct)
params03['deltaT'] = np.round(deltaT)
grid_params['Nt'] = nTimeSteps


# ## Create 'data' files

# In[29]:


# NOTE: These steps generate the data text files located int the input directory


# In[30]:


# gather params for data file 
params04 = {} #<-- using params04 to be consistent with ordering in Andrew's code
params04['usingCartesianGrid'] = grid_params['usingCartesianGrid']
params04['usingSphericalPolarGrid'] = grid_params['usingSphericalPolarGrid']
params04['delX']  = grid_params['delX']
params04['delY'] = grid_params['delY']
params04['delZ'] = dz


# In[31]:


# get data fnames param
params05 = {}
params05['bathyFile'] ='bathymetry.bin'
params05['hydrogThetaFile'] = 'tempini.bin'
params05['hydrogSaltFile'] = 'saltini.bin'


# In[32]:


data_params = [params01, params02, params03, params04, params05]
rcf.write_data(run_config, data_params, group_name='data', lf=run_config['lf'])


# In[ ]:





# ## Make SIZE.h File

# In[19]:


# generate size.h file
rcf.createSIZEh(run_config, grid_params)


# ## Specifiy relaxation of temperature and salinity 

# In[20]:


OBCS = {}


# ## Set boundary conditions

# In[ ]:





# ## Specify Diagnostics

# In[21]:


# Here we specify variables that should saved (i.e., written to disk) at various time intervals

# adjust output frequency
if run_config['test']:
    run_config['inst_freq'] = 1 # multiples of timestep
    run_config['tavg_freq'] = 5 # multiples of timestep
    
else:
    run_config['inst_freq'] = 1 # multiples of year
    run_config['tavg_freq'] = 5 # multiples of year


#---------specify time averaged fields------#
# NOTE: many more options available see mitgcm docs
diag_fields_avg = ['UVEL', 'VVEL', 'WVEL','UVELTH', 'VVELTH','WVELTH','THETA']

# diag_fields_avg = ['UVEL', 'VVEL', 'WVEL', 'UVELSQ', 'VVELSQ', 'WVELSQ',
#                   'UVELTH', 'VVELTH', 'WVELTH', 'THETA', 'THETASQ',
#                   'PHIHYD', 'LaUH1TH', 'LaVH1TH', 'LaHw1TH','LaHs1TH']

numdiags_avg = len(diag_fields_avg)
diag_phase_avg = 0.0

if run_config['test'] == True:
    diag_freq_inst = -run_config['inst_freq']*deltaT # negative values indicate snapshots at given interval
    diag_freq_avg = run_config['tavg_freq']*deltaT # positive values indicate time average over specified interval
else:
    diag_freq_inst = -run_config['inst_freq']*secsInYear 
    diag_freq_avg = run_config['tavg_freq']*secsInYear
    
    
diag_params01 = {}
diag_params01['diag_mnc'] = False #<---you would need to modify this if you want netcdf output

for ii in range(numdiags_avg):  
    n = ii+1
    diag_params01['fields(1,%s)'%n] = diag_fields_avg[ii]
    diag_params01['fileName(%s)'%n] = diag_fields_avg[ii]
    diag_params01['frequency(%s)'%n] = diag_freq_avg
    diag_params01['timePhase(%s)'%n] = diag_phase_avg

    
#--------specify instanteous fields (i.e. snapshots)--------#
diag_fields_inst = ['UVEL', 'VVEL', 'THETA']
numdiags_inst = len(diag_fields_inst)
diag_phase_inst = 0.0

for ii in range(numdiags_inst):
    n = numdiags_avg+ii+1
    diag_params01['fields(1,%s)'%n] = diag_fields_inst[ii]
    diag_params01['fileName(%s)'%n] = diag_fields_inst[ii] + '_inst'
    diag_params01['frequency(%s)'%n] = diag_freq_inst
    diag_params01['timePhase(%s)'%n] = diag_phase_inst
    
print(diag_params01)
Ndiags = n


# In[22]:


diag_params02={}
diag_params = [diag_params01, diag_params02]
rcf.write_data(run_config, diag_params, group_name='diagnostics')


# In[23]:


## create DIAGNOSTICS_SIZE.h
Nlevels = np.max([grid_params['Nr'], Nlayers])
rcf.createDIAGSIZEh(run_config, Ndiags, Nlevels)


# ## Generate eedata 

# In[33]:


# create eedata
rcf.create_eedata(run_config, grid_params['nTx'], grid_params['nTy'])


# ## Estimate wall clock time

# ## Next steps
# 
# Once you've successfully set up the model experiment, you will need to do the following:
# 
# + Navigate to the build directory and compile the code (i.e., `sh build.sh`). This may take a couple minutes and will generates lots of new files in the build directory. The key file is the `mitgcmuv` executable.
# 
# 
# 
# + Next, go to the results directory and submit the job to Sherlock queue (`sh run.sh`). It is recommended that you do this from a login node. This should generate a unique output_xxx file, which logs error messages for each run. You will get emails when the job is released from the queue and when it is completed. Test runs normally take a few minutes. A successful run will generate dozens of output files, including data files like `THETA_inst_000000.data`. 
# 
# 
# To check if the job is running:   `squeue | grep ejdawson`
# 
