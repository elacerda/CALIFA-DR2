#!/usr/bin/python
#
# Lacerda@Saco - 3/Sep/2014
import numpy as np
import matplotlib as mpl
from pycasso import fitsQ3DataCube
from matplotlib import pyplot as plt

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

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
galaxiesListFile = CALIFAWorkDir + 'list15.txt'
DRVersion = 'DR2'

versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'v20_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = CALIFAWorkDir + 'gal_fits/v20_q043.d14a/')
#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'v20_q046.d15a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = CALIFAWorkDir + 'gal_fits/v20_q046.d15a/')
#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'px1_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = CALIFAWorkDir + 'gal_fits/px1_q043.d14a/')
#versRun = dict(baseCode = 'Bgsd61', versionSuffix = 'v20_q036.d13c', othSuffix = '512.ps03.k2.mC.CCM.', SuperFitsDir = CALIFAWorkDir + 'gal_fits/v20_q036.d13c/')
#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'v20_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = '/Volumes/backupzeira/CALIFA/q043/v20/Bgsd6e/')
#versRun = dict(baseCode = 'Bgsd61', versionSuffix = 'v20_q036.d13c', othSuffix = '512.ps03.k2.mC.CCM.', SuperFitsDir = '/Volumes/backupzeira/CALIFA/q036/v20/Bgsd61/')
#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'px1_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = '/Volumes/backupzeira/CALIFA/q043/d14a/px1/')

imgDir = CALIFAWorkDir + 'images/'

f = open(galaxiesListFile, 'r')
listOfPrefixes = f.readlines()
f.close()

if debug:
    listOfPrefixes = listOfPrefixes[0:20]        # Q&D tests ...
    #listOfPrefixes = ['K0026\n']
    
N_gals = len(listOfPrefixes)

def CALIFAResidualSpectraStack(x_ini, x_fin, dx, x_label, y_ini, y_fin, dy, y_label, z, z_label, fileName):
    y, x = np.mgrid[slice(y_ini, y_fin + dy, dy),
                    slice(x_ini, x_fin + dx, dx)]
    f = plt.figure()
    f.set_size_inches(10, 4)
    plt.pcolormesh(x, y, z, cmap = plt.get_cmap('gray'), vmax = 10, vmin = -10)
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
    N_zone1pix = 0  

    for iGal in np.arange(N_gals):
        galName = listOfPrefixes[iGal][:-1]
        
        CALIFASuffix = '_synthesis_eBR_' + versRun['versionSuffix'] + versRun['othSuffix'] + versRun['baseCode'] + '.fits'
        CALIFAFitsFile = versRun['SuperFitsDir'] + galName + CALIFASuffix
        
        K = fitsQ3DataCube(CALIFAFitsFile)
        
        N_zones += K.N_zone
        N_zone1pix += np.where(K.zoneArea_pix == 1, 1, 0).sum()
        print K.califaID, K.N_zone, N_zones, np.where(K.zoneArea_pix == 1, 1, 0).sum(), N_zone1pix
        
        redshift['gal'][iGal] = K.califaID
        redshift['z'][iGal] = K.redshift
        
        # Building an array with all the nuclear spectra from DR2 galaxies. 
        f_res_perc[iGal, :] = (1 - K.f_syn[:, 0] / K.f_obs[:, 0]) * 100
        f_res_norm_perc[iGal, :] = ((K.f_obs[:, 0] - K.f_syn[:, 0]) / K.fobs_norm[0]) * 100
        
    # Align the spectra by redshift
    S = np.argsort(redshift['z'])
    f_res_perc_sorted = f_res_perc[S, :]
    f_res_norm_perc_sorted = f_res_norm_perc[S, :]
    
    print " TOTAL: %d %d" % (N_zones, N_zone1pix)
    
    #Stack spectra
    CALIFAResidualSpectraStack(x_ini = K.l_obs[0], x_fin = K.l_obs[-1], dx = 2., x_label = r'rest-frame wavelength $[\AA]$',  
                               y_ini = 0, y_fin = iGal, dy = 1., y_label = r'%s galaxy spectra' % DRVersion, 
                               z = f_res_perc_sorted, z_label = r'deviation from stellar model [%]',
                               fileName = '%sStackNuclearResidual_dev_perc.%s' % (DRVersion, outputImgSuffix))  

    CALIFAResidualSpectraStack(x_ini = K.l_obs[0], x_fin = K.l_obs[-1], dx = 2., x_label = r'rest-frame wavelength $[\AA]$',  
                               y_ini = 0, y_fin = iGal, dy = 1., y_label = r'%s galaxy spectra' % DRVersion, 
                               z = f_res_norm_perc_sorted, z_label = r'deviation from stellar model [%]',
                               fileName = '%sStackNuclearResidual_norm_perc.%s' % (DRVersion, outputImgSuffix))  
