#!/usr/bin/python
#
# Lacerda@Saco - 3/Sep/2014
import numpy as np
import matplotlib as mpl
from pycasso import fitsQ3DataCube
from matplotlib import pyplot as plt
from CALIFAUtils import paths 
from CALIFAUtils.scripts import sort_gals
from CALIFAUtils.scripts import loop_cubes

# All this part before the function definitions are just to format the plot and to
# configure which version of CALIFA is used to do the plot (basically to retrieve
# the filename).
mpl.rcParams['font.size']       = 14
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

outputImgSuffix = 'png'
debug = False
#debug = True
v_run = -1

CALIFAWorkDir = paths.califa_work_dir
#galaxiesListFile = CALIFAWorkDir + 'list15.txt'
galaxiesListFile = CALIFAWorkDir + 'listv20_q050.d15a.txt'
#galaxiesListFile = CALIFAWorkDir + 'listDR2.txt'
#DRVersion = 'q046'
DRVersion = 'DR3'
versionSuffix = 'v20_q053.d22a512.mE'

paths.set_v_run(v_run)

imgDir = CALIFAWorkDir + 'images/'

gals, _ = sort_gals(galaxiesListFile)
N_gals = len(gals)
maxGals = None
if debug:
    maxGals = 10
    if N_gals > maxGals:
        N_gals = maxGals

def CALIFAResidualSpectraStack(x_ini, x_fin, dx, x_label, y_ini, y_fin, dy, y_label, z, z_label, fileName):
    y, x = np.mgrid[slice(y_ini, y_fin + dy, dy),
                    slice(x_ini, x_fin + dx, dx)]
    f = plt.figure()
    f.set_size_inches(10, 4)
    plt.pcolormesh(x, y, z, cmap = mpl.cm.Blues_r, vmax = 10, vmin = -10)
    cb = plt.colorbar()
    plt.axis([x_ini, x_fin, y_ini, y_fin])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    cb.set_label(z_label)
    plt.tight_layout()
    f.savefig(fileName)
    
if __name__ == '__main__':
    print 'Plotting with %d galaxies' % N_gals
    
    redshift = np.ndarray((N_gals), dtype = ([('gal' , '|S10'), ('z', np.float64)]))
    f_res_perc = np.ndarray((N_gals, 1601))
    f_res_norm_perc = np.ndarray((N_gals, 1601))
    N_zones = 0 
    N_zonesOk = 0 
    N_zone1pix = 0  
    zonesOk = np.ones((N_gals), dtype = np.bool)

    for iGal, K in loop_cubes(gals.tolist(), imax = maxGals, v_run = v_run):
        if K is None or gals[iGal] == 'K0847':
            print '<<< miss file:', gals[iGal]
            zonesOk[iGal] = False
            continue
        zonesOk = K.filterResidual(w2=4600)
        N_zonesOk += np.asarray(zonesOk, dtype = int).sum()
        N_zones += K.N_zone
        N_zone1pix += np.where(K.zoneArea_pix == 1, 1, 0).sum()
        print K.califaID, K.N_zone, zonesOk.sum(), N_zones, np.where(K.zoneArea_pix == 1, 1, 0).sum(), N_zone1pix
        
        redshift['gal'][iGal] = K.califaID
        redshift['z'][iGal] = K.redshift
        
        # Building an array with all the nuclear spectra from DR2 galaxies. 
        f_res_perc[iGal, :] = (1 - K.f_syn[:, 0] / K.f_obs[:, 0]) * 100
        f_res_norm_perc[iGal, :] = ((K.f_obs[:, 0] - K.f_syn[:, 0]) / K.fobs_norm[0]) * 100
        
    # Align the spectra by redshift
    S = np.argsort(redshift['z'][zonesOk])
    f_res_perc_sorted = f_res_perc[zonesOk][S, :]
    f_res_norm_perc_sorted = f_res_norm_perc[zonesOk][S, :]
    
    print " TOTAL: %d %d %d" % (N_zones, N_zone1pix, N_zonesOk)
    
    #Stack spectra
    CALIFAResidualSpectraStack(x_ini = K.l_obs[0], x_fin = K.l_obs[-1], dx = 2., x_label = r'rest-frame wavelength $[\AA]$',  
                               y_ini = 0, y_fin = iGal, dy = 1., y_label = r'%s galaxy spectra' % DRVersion, 
                               z = f_res_perc_sorted, z_label = r'deviation from stellar model [%]',
                               fileName = '%sStackNuclearResidual_dev_perc.%s' % (DRVersion, outputImgSuffix))  

    CALIFAResidualSpectraStack(x_ini = K.l_obs[0], x_fin = K.l_obs[-1], dx = 2., x_label = r'rest-frame wavelength $[\AA]$',  
                               y_ini = 0, y_fin = iGal, dy = 1., y_label = r'%s galaxy spectra' % DRVersion, 
                               z = f_res_norm_perc_sorted, z_label = r'deviation from stellar model [%]',
                               fileName = '%sStackNuclearResidual_norm_perc.%s' % (DRVersion, outputImgSuffix))  
