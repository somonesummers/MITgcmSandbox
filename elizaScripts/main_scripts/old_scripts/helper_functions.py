import numpy as np
import os
import imp
import matplotlib.pylab as plt
import MITgcmutils
from IPython.core.debugger import set_trace
import pickle
import run_utils as ru

imp.reload(ru)

debug_here = set_trace

#home_dir = '/home/earlew/research/scripts/MITgcm_py'

#home_dir = '/home/users/earlew/research/MITgcm_python/idealized_channel' #<--would be nice to not have this hard coded
#home_dir = os.environ['HOME']

#output_dir = os.path.join(home_dir, 'output/')

# TODO: import these directories
# analysis_dir = '/scratch/users/earlew/research/modeling/MITgcm_PG/analysis/'
# output_dir = os.path.join(analysis_dir, 'output/')
# plot_dir = os.path.join(analysis_dir, 'plots/')

# def set_exp_dir(exp_dir):
    
#     with open(os.path.join(home_dir, 'exp_dir_loc.txt'), 'w') as f:
#         f.write(exp_dir)



def get_exp_params_py(exp_path, iter0=0, add_hFac=True, add_forcing=True, 
                      vname='THETA', periodic_forcing_tlen=365):
    
    """
    Function to load run configurations and parameters. This script has the same
    functionality as get_exp_params(), which imports run parameters from a file generated by
    Andrew Stewart's matlab script.
    
    
    NOTE: this probably should be in run_utils
    
    
    """
    
    #get exp_dir
#     with open(os.path.join(home_dir, 'exp_dir_loc.txt'), 'r') as f:
#         exp_home_dir = f.read()

    #exp_path = os.path.join(exp_root_dir, exp_name)
    exp_name = exp_path.split('/')[-1]
    results_path = os.path.join(exp_path, 'results')
    input_path = os.path.join(exp_path, 'input')
    
    # open run config file
    with open(os.path.join(input_path, 'params.p'), 'rb') as f:
        run_config = pickle.load(f)
       
    # import parameters (Note: using params structure from previous code for convenience)
    params = {}
    params['Nx'] = run_config['grid_params']['Nx']
    params['Ny'] = run_config['grid_params']['Ny']
    params['Nr'] = run_config['grid_params']['Nr']

    params['delX'] = run_config['grid_params']['delX']
    params['delY'] = run_config['grid_params']['delY']
    params['delR'] = run_config['grid_params']['delR']


    params['Lx'] = run_config['domain_params']['Lx']
    params['Ly'] = run_config['domain_params']['Ly']
    params['H'] = run_config['domain_params']['H']

    params['xx'] = np.cumsum(np.append([0], params['delX'][:-1])) + params['delX']/2 - params['Lx']/2
    params['yy'] = np.cumsum(np.append([0], params['delY'][:-1])) + params['delY']/2 
    params['zz'] = -(np.cumsum(np.append([0], params['delR'][:-1])) + params['delR']/2) 

    params_IO = run_config['mitgcm_params']['data'][2]
    params['deltaT'] = params_IO['deltaT']

    # TODO: use available data in results directory
#     params['startTime'] = params_IO['nIter0']*params_IO['deltaT']
#     params['nTimeSteps'] = params_IO['nIter0'] + np.ceil((params_IO['endTime'] - 
#                                                           params['startTime'])/params['deltaT'])
    
    params['dumpIters'] = ru.getSavedIters(exp_name, vname, results_path)
    params['dumpIters'] = params['dumpIters'][params['dumpIters'] > iter0]
    #print(iters)
    #iter0 = iters[0]
    params['startTime'] = params['dumpIters']*params_IO['deltaT']
    params['nDumps'] = len(params['dumpIters'])
    params['nTimeSteps'] = np.ceil((params_IO['endTime'] - params['startTime'])/params['deltaT'])
    

    params['rho0'] = 1000 # kg/m3
    params['Cp'] = 4e3  # J/kg/K

    
    # extract dump freq (varies based on variable)
    diags = run_config['mitgcm_params']['data_diagnostics'][0]
    
    #print(diags)

    # map field name to dump freq
    vble_dump_freq = {}
    for key in diags:

        if key.startswith('fileName'):
            fname = diags[key]

            if fname not in vble_dump_freq:
                fnum = int(key[len('fileName')+1: -1])
                #print(fnum)
                vble_dump_freq[fname] = abs(diags['frequency(%i)' %fnum])
    
    params['dumpFreq'] = vble_dump_freq[vname]
#     params['dumpFreq'] = np.abs(params_IO['dumpFreq'])

#     params['nDumps'] = np.floor(params['nTimeSteps']*params['deltaT']/params['dumpFreq']) 
#     params['dumpIters'] = np.round(np.arange(1, params['nDumps']+1)*params['dumpFreq']/params['deltaT'])
#     params['dumpIters'] = params['dumpIters'][params['dumpIters'] > iter0]
    
    # pass mitgcm params and grid_params
    params['mitgcm_params'] = run_config['mitgcm_params']
    params['grid_params'] = run_config['grid_params']
    
    if add_hFac:
        params['DRF'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'DRF'))
        params['hFacS'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacS'))
        params['hFacW'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacW'))
        params['hFacC'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacC'))

    # tack on exp_name and directory paths
    params['exp_name'] = exp_name
    params['exp_path'] = exp_path
    params['input_path'] = input_path
    params['results_path'] = results_path
    
    
    # get forcing and bathy files names
    input_fnames = run_config['mitgcm_params']['data'][4]
    params['bathyFile'] = input_fnames['bathyFile']
    params['zonalWindFile'] = input_fnames['zonalWindFile']
    params['surfQfile'] = input_fnames['surfQfile']
   
    # read in forcing data
    forcing = {}
    
    bathy_path = os.path.join(input_path, params['bathyFile'])
    forcing['bathy'] = np.fromfile(bathy_path, dtype='>f8').reshape(params['Ny'], params['Nx'])
    
    params03 = params['mitgcm_params']['data'][2]
    
    if 'periodicExternalForcing' in params03 and params03['periodicExternalForcing']:
        periodic_forcing = True
        forcing_shape = (periodic_forcing_tlen, params['Ny'], params['Nx'])
    else:
        forcing_shape = (params['Ny'], params['Nx'])

    if add_forcing:
        wind_path = os.path.join(input_path, params['zonalWindFile'])
        forcing['zonalWind'] = np.fromfile(wind_path, dtype='>f8').reshape(forcing_shape)

        wind_path = os.path.join(input_path, params['surfQfile'])
        forcing['surfQ'] = np.fromfile(wind_path, dtype='>f8').reshape(forcing_shape)
        
        if 'EmPmRFile' in input_fnames:
            params['EmPmRFile'] = input_fnames['EmPmRFile']
            emp_path = os.path.join(input_path, params['EmPmRFile'])
            forcing['EmPmR'] = np.fromfile(emp_path, dtype='>f8').reshape(forcing_shape)

        
    else:
        forcing['zonalWind'] = np.zeros(forcing_shape)
        forcing['surfQ'] = np.zeros(forcing_shape)
        forcing['EmPmR'] = np.zeros(forcing_shape)
                                        
        

    
    return params, forcing
    
    


def get_exp_params(exp_name, vname='THETA', add_hFac=True, add_forcing=True):
    
    """
    
    Note: vname='THETA' is only used to get the dump frequency.
    """
    
    # set experiment directory
    #exp_name = 'gyre_ctrl_run_128n_50yr_SS_WAP_DP_NWR'

    # set plot directory
    #plot_dir = os.path.join('/home/earlew/research/plots/MITgcm/', exp_name)
    #os.makedirs(plot_dir, exist_ok=True)
    
    # get exp_dir
    with open(os.path.join(home_dir, 'exp_dir_loc.txt'), 'r') as f:
        exp_dir = f.read()
    
    exp_path = os.path.join(exp_dir, exp_name)
    input_path = os.path.join(exp_path, 'input')
    results_path = os.path.join(exp_path, 'results')
    
    if vname.lower()=='ssh':
        vname='PHIHYD'
        
    elif vname in ['BT_psi', 'SPD', 'SPD_inst']:
        vname='UVEL'
    
    # get model params
    params = scrape_params_m(input_path)
    
    # derive and set additional model parameters
    params['Nx'] = len(params['delX'])
    params['Ny'] = len(params['delY'])
    params['Nr'] = len(params['delR'])

    params['Lx'] = np.sum(params['delX'])
    params['Ly'] = np.sum(params['delY'])
    params['H'] = np.sum(params['delR'])

    params['xx'] = np.cumsum(np.append([0], params['delX'][:-1])) + params['delX']/2 - params['Lx']/2
    params['yy'] = np.cumsum(np.append([0], params['delY'][:-1])) + params['delY']/2 
    params['zz'] = -(np.cumsum(np.append([0], params['delR'][:-1])) + params['delR']/2) 

    params['startTime'] = params['nIter0']*params['deltaT']
    params['nTimeSteps'] = params['nIter0'] + np.ceil((params['endTime'] - params['startTime'])/params['deltaT'])

    params['rho0'] = 1000 # kg/m3
    params['Cp'] = 4e3  # J/kg/K

    params['dumpFreq'] = np.abs(params['diag_frequency'][vname])
    params['nDumps'] = np.floor(params['nTimeSteps']*params['deltaT']/params['dumpFreq']) 
    params['dumpIters'] = np.round(np.arange(1, params['nDumps']+1)*params['dumpFreq']/params['deltaT'])
    params['dumpIters'] = params['dumpIters'][params['dumpIters'] > params['nIter0']]
    
    #debug_here()

    if add_hFac:
        params['DRF'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'DRF'))
        params['hFacS'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacS'))
        params['hFacW'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacW'))
        params['hFacC'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacC'))

    # tack on exp_name and directory paths
    params['exp_name'] = exp_name
    params['exp_path'] = exp_path
    params['input_path'] = input_path
    params['results_path'] = results_path
    
    # read in forcing data
    forcing = {}
    
    bathy_path = os.path.join(input_path, params['bathyFile'])
    forcing['bathy'] = np.fromfile(bathy_path, dtype='>f8').reshape(params['Ny'], params['Nx'])

    if add_forcing:
        wind_path = os.path.join(input_path, params['zonalWindFile'])
        forcing['zonalWind'] = np.fromfile(wind_path, dtype='>f8').reshape(params['Ny'], params['Nx'])

        wind_path = os.path.join(input_path, params['surfQfile'])
        forcing['surfQ'] = np.fromfile(wind_path, dtype='>f8').reshape(params['Ny'], params['Nx'])
        
    else:
        forcing['zonalWind'] = np.zeros((params['Ny'], params['Nx']))
        forcing['surfQ'] = np.zeros((params['Ny'], params['Nx']))
                                        
        

    
    return params, forcing

def rdmdsWrapper(fpath,dumpIter):
    
    """
    Convenience wrapper for rdmds.m that tries to read the specified file
    'fname' and iteration 'dumpIter'. If that name/iteration can't be
    found, we try dumpIter+1 and dumpIter-1 in case of rounding errors
    between Matlab and MITgcm.
    """
    
    A = MITgcmutils.mds.rdmds(fpath, dumpIter)        
    if len(A)==0:
        A = MITgcmutils.mds.rdmds(fpath, dumpIter-1)
    if len(A)==0:
        A = MITgcmutils.mds.rdmds(fpath, dumpIter+1)
        
    return A

def readIters(vname, tr, params, mask_zeros=True, returnTimeAvg=True, printStatus=True):
    
    if printStatus:
        print("Loading %s for years %s-%s..." %(vname, tr[0], tr[1]))
    
    dumpIters = params['dumpIters']
    if vname.startswith("La"):
        zlen = len(params['layers_bounds'])-1
    else:
        zlen = params['Nr']
    
    #print(dumpIters)
    #vble_t = np.zeros((params['Nx'], params['Ny'], zlen))
    #avg = np.zeros((params['Nx'], params['Ny'], zlen))
    #navg = 0
    secsInYear = 86400*365
    dt_sec = params['deltaT']
    ii = 0
    tvec = []
    tyears = params['dumpIters']*dt_sec/secsInYear
    # loop through output iterations
    for nn in range(len(dumpIters)):
        
        dump_time = params['dumpIters'][nn]*dt_sec
        dump_time_in_range = dump_time>=(tr[0]*secsInYear-2*dt_sec) and dump_time<=(tr[-1]*secsInYear+2*dt_sec)
        
        if dump_time_in_range:
            # Read the next iteration and check the data were found  
            try:
                vble_tmp = rdmdsWrapper(os.path.join(params['results_path'], vname), dumpIters[nn])
            except OSError:
                #break # continue?
                continue
                
            if len(vble_tmp)==0:
                print('Ran out of data at ITER= %s/%s t= %s days.' %(dumpIters[nn], params['nDumps'], tyears))
                break
            else:
                vble_tmp = vble_tmp.T
                if ii==0:
                    vble_t = vble_tmp[..., np.newaxis]
                else:
                    vble_t = np.concatenate([vble_t, vble_tmp[..., np.newaxis]], axis=3)
                    
                #avg += avg_tmp.T
                ii += 1
                tvec.append(dump_time)
                
    if ii==0:
        debug_here()
        raise ValueError("No output files found.")

    vble_t = np.ma.masked_invalid(vble_t)
    if mask_zeros:
        vble_t = np.ma.masked_where(vble_t==0, vble_t)

    
    # calculate average
    if returnTimeAvg:
        return vble_t.mean(axis=3)
    else:
        return vble_t, np.array(tvec)
    

def check_temp_trend(exp_name, vname='THETA_inst', zr=[], tstep=1, load_previous=True, tw=30, dec_trend_thresh=0.005):
    
    from scipy import stats
    
    pdata = computeDomainInt(exp_name, vname=vname, zr=zr, tstep=tstep, load_previous=load_previous)
    
    tvec_yrs = np.array(pdata['tt'])
    ti = tvec_yrs>=tvec_yrs.max()-tw
    temp_avg_tt = np.array(pdata['vble_avg'])[ti]
    years_tt = tvec_yrs[ti]
    
    slope, intercept, r, p, se = stats.linregress(years_tt, temp_avg_tt)
    decadal_trend = slope*10
    print("Experiment: %s " %exp_name)
    print("%s trend for last %s years: %.4f C/dec" %(vname, tw, decadal_trend))
    print("Equilibrium threshold met? %s" %(decadal_trend<dec_trend_thresh))
    print("--------------------------------")
    
        

def plot_layer_avg_vals(vble, vname, zr, tr, params, clvls=[], cmap=plt.cm.rainbow, fz=14):
    
    if vname.lower() in ['ssh']:
        vble_zavg = vble[...,0]
    else:
        zz = np.abs(params['zz'])
        zi = np.logical_and(zz>=zr[0], zz<=zr[-1])
        vble_zavg = vble[:, :, zi].mean(axis=2)
    
    fig = plt.figure(figsize=(10, 8))
    YY,XX = np.meshgrid(params['yy']/1000, params['xx']/1000)
    if len(clvls)==0:
        im = plt.contourf(XX, YY, vble_zavg, 30, cmap=cmap)
        cb = plt.colorbar()
    else:
        im = plt.contourf(XX, YY, vble_zavg, clvls, cmap=cmap, extend='both')
        cb = plt.colorbar(im, extend='both')
    
    # add labels
    plt.xlabel("X (km)", fontsize=fz)
    plt.xlabel("Y (km)", fontsize=fz)
    plt.title("Mean %s-%sm %s for years %s-%s" %(zr[0], zr[1], vname, tr[0], tr[-1]), fontsize=fz)
    cb.set_label('%s (%s)' %(vname, get_units(vname.upper())), fontsize=fz)
    
    # set label size
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)
    cb.ax.tick_params(labelsize=fz) 
    
    plt.show()
    
def get_units(vname):

    if 'THETA' in vname:
        units = '$^{\circ}$C'
    elif vname in ['SALT']:
        units = 'PSU'
    elif vname in ['UVEL', 'VVEL', 'WVEL'] or 'VEL' in vname or vname.startswith('SPD'):
        units = 'm/s'
    elif vname in ['SSH']:
        units = 'm'
    elif vname in ['BT_psi']:
        units = 'Sv'
    
    elif vname in ['SALT', 'S', 'SALT_inst']:
        units = 'PSU'
    
    else:
        units = ''
    
    #debug_here()
    return units
    
def get_layer_mean(vble, zr, vname, params):
    
    if vname.lower() in ['ssh']:
        vble_zavg = vble[...,0]
    else:
        zz = np.abs(params['zz'])
        zi = np.logical_and(zz>=zr[0], zz<=zr[-1])
        vble_zavg = vble[:, :, zi].mean(axis=2)

    return vble_zavg

def get_zonal_mean(vble, xr, vname, params):
    
    if len(xr)==0:
        vble_xavg = vble.mean(axis=0) 
    else:
        xx = params['xx']
        xi = np.logical_and(xx>=xr[0], xx<=xr[-1])
        vble_xavg = vble[xi, :, :].mean(axis=0)
    
    return vble_xavg

def get_slice_mean(vble, rr, vname, params, coord='xx'):
    
    """
    Function that computes the mean of a slice in the x, y, or z direction.
    This function supercedes get_layer_mean() and get_zonal_mean().
    """
    
    assert coord in ['xx', 'yy', 'zz'], "invalid coordinates"

    ii = np.logical_and(params[coord]>=rr[0], params[coord]<=rr[-1])
    
    # using if statements here is a bit clumsy
    # having vble be an xarray with labelled dimensions/coordinates would make this a lot cleaner
    if coord=='xx':
        vble_avg = vble[ii, :, :].mean(axis=0)
    elif coord=='yy':
        vble_avg = vble[:, ii, :].mean(axis=1)
    elif coord=='zz':
        vble_avg = vble[:, :, ii].mean(axis=2)
        
    return vble_avg

def compute_MLD_fast(temp, params, dT=0.2, xr=[], get_xavg=False):
    
    z = np.array(params['zz'])
    temp = np.ma.masked_invalid(temp)
    temp_mld = temp[:, :, 0] - dT
    idx = (np.ma.abs(temp - temp_mld[:, :, np.newaxis])).argmin(axis=2)
    mld = z[idx]
    
    if get_xavg:
        assert len(xr)==2
        xx = params['xx']
        xi = np.logical_and(xx>=xr[0], xx<=xr[-1])
        return mld[xi,...].mean(axis=0) 
    else:
        return mld

def calc_BT_streamfunction(exp_name, tr):   
    """
    Function to compute barotropic streamfunction. Translated from Andrew Stewart 
    """
    
    # load control params
    params, _ = get_exp_params_py(exp_name)
    
    # load uvel
    uvel = readIters(vname='UVEL', tr=tr, params=params)
    #theta = readIters(vname='THETA', tr=tr, params=params)
    
    #uvel = np.ma.masked_where(theta.mask, uvel)
    
    # Calculate depth-integrated zonal velocity
    hFacW = params['hFacW'].T
    DZ_3d = params['delR'][np.newaxis, np.newaxis, :]
    uvel_zint = np.sum(uvel*DZ_3d*hFacW, axis=2)
    
    # Calculate barotropic streamfunction
    Nx = params['Nx']
    Ny = params['Ny']
    DY_2d = params['delY'][np.newaxis, :]
    
    Psi = np.zeros((Nx+1, Ny+1))
    Psi[:Nx, :Ny] = -np.cumsum(uvel_zint*DY_2d, 1)
    Psi[:, Ny] = 0
    Psi[Nx, :] = Psi[0, :]
    
    return Psi

def calc_streamfunction(exp_name, tr, zr):   
    """
    Function to compute barotropic streamfunction. Translated from Andrew Stewart 
    """
    
    # load control params
    params, _ = get_exp_params_py(exp_name)
    
    # load uvel
    uvel = readIters(vname='UVEL', tr=tr, params=params)
    #theta = readIters(vname='THETA', tr=tr, params=params)
    
    #uvel = np.ma.masked_where(theta.mask, uvel)
    zz = np.abs(params['zz'])
    zi = np.logical_and(zz>=zr[0], zz<=zr[-1])
    
    # Calculate depth-integrated zonal velocity
    hFacW = params['hFacW'].T
    hFacW_z = hFacW[:, :, zi]
    DZ_3d = params['delR'][np.newaxis, np.newaxis, :]
    DZ_3d_z = DZ_3d[..., zi]
    uvel_zint = np.sum(uvel[..., zi]*DZ_3d_z*hFacW_z, axis=2)
    
    #debug_here()
    
    # Calculate barotropic streamfunction
    Nx = params['Nx']
    Ny = params['Ny']
    DY_2d = params['delY'][np.newaxis, :]
    
    Psi = np.zeros((Nx+1, Ny+1))
    Psi[:Nx, :Ny] = -np.cumsum(uvel_zint*DY_2d, 1)
    Psi[:, Ny] = 0
    Psi[Nx, :] = Psi[0, :]
    
    return Psi

 
def computeDomainInt(exp_name, vname, zr=[], xr=[], yr=[], tstep=1, getKE=False, load_previous=True, 
                     suppress_err_msg=True):
    
    #home_dir = '/home/earlew/research/scripts/MITgcm_py/'
    #output_dir = os.path.join(home_dir, 'output/')
    
    if getKE:
        vname = 'UVEL'
        
    #load control params
    params, forcing = get_exp_params_py(exp_name, vname=vname)
    nDumps = int(params['nDumps'])-1
    dumpIters = params['dumpIters']
    
#     print(nDumps)
#     print(dumpIters)
    
    zz = params['zz']
    xx = params['xx']
    yy = params['yy']
    
    
    def get_domain_subset(sr, ss, dim):
        
        if len(sr)==2:
            assert sr[-1]>=sr[0]
            si = np.logical_and(ss>=sr[0], ss<sr[-1])
            if dim=='z':
                s_str = "_%s=%s-%sm" %(dim, np.abs(sr).min(), np.abs(sr).max())
            else:
                s_str = "_%s=%s-%skm" %(dim, np.abs(sr).min()/1e3, np.abs(sr).max()/1e3)
        else:
            si = np.logical_and(ss>=ss.min(), ss<=ss.max())
            s_str = ''
            
        return si, s_str
    
#     if len(zr)==2:
#         assert zr[-1]>=zr[0]
#         zi = np.logical_and(zz>=zr[0], zz<zr[-1])
#         z_str = "_%s-%sm" %(np.abs(zr).min(), np.abs(zr).max())
#     else:
#         zi = np.logical_and(zz>=zz.min(), zz<=zz.max())
#         z_str = ''

    zi, z_str = get_domain_subset(zr, zz, 'z')
    yi, y_str = get_domain_subset(yr, yy, 'y')
    xi, x_str = get_domain_subset(xr, xx, 'x')
        
    s_str = z_str+y_str+x_str
    # Get domain volume
    DX = params['delX'][xi, np.newaxis, np.newaxis]
    DY = params['delY'][np.newaxis, yi, np.newaxis]
    DZ = params['delR'][np.newaxis, np.newaxis, zi]
    DV = DX*DZ*DY*params['hFacC'][zi, :, :][:, yi, :][:, :, xi].T
    vol = np.sum(DV)
    
    #print(vol)
    
    rho0 = 1025 # kg/m3


            
    # pre-allocate containers and counters
    tt = []
    iters_list = [] # logs dumped files actually found
    ptlen = 0
    n0 = 0
    
    if getKE:
        TKE = []
        TEKE = []
        pdata = {'tt': [], 'TKE': [], 'TEKE':[], 'dumpIters': []}
        fname = '%s_KE_tseries_tstep%s.p'%(exp_name, tstep)
    else:
        vble_avg = []
        pdata = {'tt': [], 'vble_avg':[], 'dumpIters': []}
        fname = '%s_%s_avg%s_tseries_tstep%s.p'%(exp_name, vname, s_str, tstep)

    # load previous output
    output_fpath = os.path.join(output_dir, fname)
    
    if load_previous:
        try:
            pdata = pickle.load(open(output_fpath, 'rb'))
            n0 = np.argmin(np.abs(pdata['dumpIters'][-1]-dumpIters))+1 # set n0 so loop starts from the last saved dumpIter
            assert n0<=nDumps
        except FileNotFoundError:
            pass

    if suppress_err_msg:
        print_err = False
    else:
        print_err = True

    #debug_here()
    for n in range(n0, nDumps):          
        if n==0 or (n%tstep)==0 or n==(nDumps-1):

            try:
                if getKE:
                    uvel = np.ma.masked_invalid(rdmdsWrapper(os.path.join(params['exp_path'], 'results/UVEL'), dumpIters[n]))
                    vvel = np.ma.masked_invalid(rdmdsWrapper(os.path.join(params['exp_path'], 'results/VVEL'), dumpIters[n]))
                    wvel = np.ma.masked_invalid(rdmdsWrapper(os.path.join(params['exp_path'], 'results/WVEL'), dumpIters[n]))
                    uvelsq = np.ma.masked_invalid(rdmdsWrapper(os.path.join(params['exp_path'], 'results/UVELSQ'),
                                                                  dumpIters[n]))
                    vvelsq = np.ma.masked_invalid(rdmdsWrapper(os.path.join(params['exp_path'], 'results/VVELSQ'),
                                                                  dumpIters[n]))
                    wvelsq = np.ma.masked_invalid(rdmdsWrapper(os.path.join(params['exp_path'], 'results/WVELSQ'),
                                                                  dumpIters[n]))

                    TKE_n = 0.5*(uvelsq + vvelsq + wvelsq)*DV.T*rho0
                    TEKE_n = TKE_n - 0.5*(uvel**2 + vvelsq**2 + wvelsq**2)*DV.T*rho0

                    # Calculate domain-mean potential temperature
                    TKE.append(np.sum(TKE_n))
                    TEKE.append(np.sum(TEKE_n))
                    
                else:
                    vble = rdmdsWrapper(os.path.join(params['exp_path'], 'results/%s' %vname), dumpIters[n])
                    # Calculate domain-mean potential temperature
                    vble_avg.append(np.sum(vble[zi, :, :][:, yi, :][:, :, xi].T*DV)/vol)
                    
                    
                tt.append(dumpIters[n]*params['deltaT']/86400/365)
                iters_list.append(dumpIters[n])
                    
                # Increment counter
                ptlen += 1

            except OSError as e:
                if print_err:
                    print(str(e))
                    print("skipping over missing files...")
                    print_err = False
                continue

    if ptlen==0:
        print("Warning: no new data loaded")
    else:
        print("%s updated to year %s" %(exp_name, np.round(tt[-1])))

    # save output
    for ii in range(len(tt)):
        pdata['tt'].append(tt[ii])
        pdata['dumpIters'].append(iters_list[ii])
        
        if getKE:
            pdata['TKE'].append(TKE[ii])
            pdata['TEKE'].append(TEKE[ii])
        else:
            pdata['vble_avg'].append(vble_avg[ii])

    pickle.dump(pdata, open(output_fpath, 'wb')) 
    
    return pdata

def find_tracer_depth(vname, tr, exp_name, exp_name_alias, clvls, testing=False, plot_test=False):
    
    from scipy import interpolate
    import xarray as xr
    import time
    
    fz = 14
    
    # load coordinates
    params, forcing = get_exp_params_py(exp_name)
    xx = params['xx']/1e3
    yy = params['yy']/1e3
    zz = forcing['bathy']
    zz[-2:, :] = np.nan

    # load tracer field
    vble = readIters(vname=vname, tr=tr, params=params)
    
    xlen = len(params['xx'])
    ylen = len(params['yy'])
    zlen = len(params['zz'])
    
    surf_z = np.ma.masked_all((xlen, ylen, len(clvls)))
    
    print("Mapping tracer values to depth. This may take a few seconds...")
    
    exit_loop = False # for testing
    testing_str = ''
    t0 = time.time()
    for ii in range(xlen):
        if exit_loop:
            break
        for jj in range(ylen):
            
            vble_ij = np.ma.masked_invalid(vble[ii, jj, :])
            
            if len(vble_ij.compressed()) == 0:
                continue
            else:
                intp = interpolate.interp1d(vble_ij[~vble_ij.mask], params['zz'][~vble_ij.mask], bounds_error=False)
                surf_z[ii, jj, :] = intp(clvls)
 
        if ii==0:
            t1 = time.time()
            dt = t1-t0
            est_comp_time = dt*xlen/60
            print("Estimated completion time: %.1f minutes" %est_comp_time)
        
        if testing:
            print("testing complete.")
            testing_str = '_test'
            exit_loop = True
            break
     
    # make a test plot
    if plot_test:
        xx = [2, 4]
        fig, axes = plt.subplots(2, 1, figsize=(15, 10))
        for ii,idx in enumerate(xx):
            plt.sca(axes[ii])
            XX, YY = np.meshgrid(params['xx']/1e3, params['yy']/1e3)
            im = plt.contourf(XX, YY, surf_z[:, :, idx].T, 25, cmap=plt.cm.rainbow_r)
            cb = plt.colorbar(im)
            cb.set_label('Depth (m)', fontsize=fz) 
            plt.xlabel("X (km)", fontsize=fz)
            plt.ylabel("Y (km)", fontsize=fz)
            plt.title("Depth of %s=%s surface" %(vname, clvls[idx]), fontsize=fz)
            #axes[ii].set_facecolor('k')
            
        plt.subplots_adjust(hspace=0.4)
    # save data
    surf_z_ds = xr.Dataset({'xx': (['x'], params['xx']), 'yy': (['y'], params['yy']), 'clvls': (['p'], clvls), 
                          'surf_z': (['x', 'y', 'p'], surf_z)})
    
    if tr[0]==tr[1]:
        t_str = '_yr%i'%tr[0]
    else:
        t_str = '_yr%i-%i'%(tr[0], tr[1])
        
    fname = '%s_surf_z_%s%s%s.nc' %(vname, exp_name_alias, t_str, testing_str)
    fpath = os.path.join(output_dir, fname)
    surf_z_ds.to_netcdf(fpath)     
    
    return surf_z_ds


def getFinalYear(exp):
    params,_ = get_exp_params_py(exp)
    lastIter = ru.getLastIter(exp)
    pickupSec = int(lastIter)*params['deltaT']
    pickupYear = pickupSec/(3600*24*365)
    
    return np.round(pickupYear)


def scrape_params_m(input_path):
    
    params_dict = {}
    params_dict['diag_name'] = {}
    params_dict['diag_freq_tmp'] = {}
    params_dict['diag_tp_tmp'] = {}
    with open(os.path.join(input_path, 'params.m'), 'r') as f: 
        for line in f:
            if line.startswith("%") or len(line.strip())==0 or '=' not in line:
                continue
            else:
                line_strp = line.strip() # remove whitespace padding
                vname = line_strp.split('=')[0]
                val = line_strp.split('=')[1][:-1] 
                if val.startswith('[') and val.endswith(']'):
                    # parse list object
                    var_list = list(val[1:-1].strip().split(' '))
                    var_list = [v for v in var_list if len(v.strip())>0]
                    params_dict[vname] = np.array(var_list).astype(float)
                else:  
                    try:
                        # try to convert numbers to float
                        if vname.startswith('diag_frequency'):
                            # these variables are stored as 'diag_frequency(1)', 'diag_frequency(2)', etc. 
                            # this is an attempt to combine them. We use an unordered dict first then sort them later.
                            params_dict['diag_freq_tmp'][int(vname[15:-1])] = float(val)
                            
                        elif vname.startswith('diag_fileNames'):
                            # these variables are stored as 'diag_fileNames{1}', 'diag_fileNames{2}', etc. 
                            # this is an attempt to combine them. We use an unordered dict first then sort them later.
                            params_dict['diag_name'][int(vname[15:-1])] = str(val.strip("'")) # strip double quotes
                            
                        elif vname.startswith('diag_timePhase'):
                            # these variables are stored as 'diag_timePhase(1)', 'diag_timePhase(1)', etc. 
                            # this is an attempt to combine them. We use an unordered dict first then sort them later.
                            params_dict['diag_tp_tmp'][int(vname[15:-1])] = float(val) 
                            # side note: the fact that these variable names have the same length is an interesting coincidence
                            
                        else:  
                            params_dict[vname] = float(val)
                            
                    except ValueError:
                        # if val is not a number store as string
                        params_dict[vname] = str(val.strip("'")) # strip double quotes
            
    # condense diagnostic stats into two dicts
    params_dict['diag_frequency'] = {}
    params_dict['diag_timePhase'] = {}
    for n in range(len(params_dict['diag_name'])):
        params_dict['diag_frequency'][params_dict['diag_name'][n+1]] = params_dict['diag_freq_tmp'][n+1]
        params_dict['diag_timePhase'][params_dict['diag_name'][n+1]] = params_dict['diag_tp_tmp'][n+1]
    
    # clean up
    params_dict.pop('diag_name')
    params_dict.pop('diag_freq_tmp')
    params_dict.pop('diag_tp_tmp')
    
    #print(params_dict['diag_frequency'])
    
    return params_dict
# if __name__ == "__main__":
    
#     # execute only if run as a script
#     get_exp_params_py(*sys.argv[1:])
    
    
