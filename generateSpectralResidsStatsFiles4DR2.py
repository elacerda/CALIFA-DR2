#!/usr/bin/python
import numpy as np
from pycasso import fitsQ3DataCube
import os

debug = False
#debug = True
#excludeWei0 = True
excludeWei0 = False

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
    
#galaxiesListFile    = CALIFAWorkDir + 'listOf300GalPrefixes.txt'
galaxiesListFile = CALIFAWorkDir + 'listDR2.txt'

#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'v20_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = CALIFAWorkDir + 'gal_fits/v20_q043.d14a/')
versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'px1_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = CALIFAWorkDir + 'gal_fits/px1_q043.d14a/')
#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'v20_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = '/Volumes/backupzeira/CALIFA/q043/v20/Bgsd6e/')
#versRun = dict(baseCode = 'Bgsd61', versionSuffix = 'v20_q036.d13c', othSuffix = '512.ps03.k2.mC.CCM.', SuperFitsDir = '/Volumes/backupzeira/CALIFA/q036/v20/Bgsd61/')
#versRun = dict(baseCode = 'Bgsd6e', versionSuffix = 'px1_q043.d14a', othSuffix = '512.ps03.k1.mE.CCM.', SuperFitsDir = '/Volumes/backupzeira/CALIFA/q043/d14a/px1/')

imgDir = CALIFAWorkDir + 'images/'

f = open(galaxiesListFile, 'r')
listOfPrefixes = f.readlines()
f.close()

if debug:
    listOfPrefixes = listOfPrefixes[0:10]        # Q&D tests ...
    #listOfPrefixes = ['K0026\n']
    
N_gals = len(listOfPrefixes)
print 'Initial list contains %d galaxies' % N_gals

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    l_int = np.arange(3650, 7200, 2)
    Nl_int = len(l_int)
    Nl_obs = 1601
    N_spectra = 0
    N_pixels = 0
    
    excludeWei0Code = 'w1'
    
    if excludeWei0:
        excludeWei0Code = 'w0'

    rad_low = np.array([   0.0, 0.0, 0.0, 1.0 ])
    rad_upp = np.array([ 999.0, 0.01, 1.0, 999.0 ])
    radCode = ['Galaxy', 'Nucleus', 'Bulge', 'Disc']
            
    for iRad in range(len(radCode)):     
        # Now reset the cumulative resid-stat arrays, both in the rest frame (RF_*) and in the observed frame (OF_*)
        RF_SumR__l = np.ma.zeros((Nl_obs))
        RF_SumR2__l = np.ma.zeros((Nl_obs))
        RF_SumS__l = np.ma.zeros((Nl_obs))
        RF_SumS2__l = np.ma.zeros((Nl_obs))
        RF_SumU__l = np.ma.zeros((Nl_obs))
        RF_SumU2__l = np.ma.zeros((Nl_obs))
        RF_NOk__l = np.ma.zeros((Nl_obs))
        OF_SumR__o = np.ma.zeros((Nl_int))
        OF_SumR2__o = np.ma.zeros((Nl_int))
        OF_SumS__o = np.ma.zeros((Nl_int))
        OF_SumS2__o = np.ma.zeros((Nl_int))
        OF_SumU__o = np.ma.zeros((Nl_int))
        OF_SumU2__o = np.ma.zeros((Nl_int))
        OF_NOk__o = np.ma.zeros((Nl_int)) 
    
        # Reset total fluxes, in OtotR, OtotS & Otot0 incarnations, and the Mtot* versions
        RF_SumOtotR__l = np.ma.zeros((Nl_obs))
        RF_SumOtotS__l = np.ma.zeros((Nl_obs))
        RF_SumOtot0__l = np.ma.zeros((Nl_obs))
        RF_SumMtotR__l = np.ma.zeros((Nl_obs))
        RF_SumMtotS__l = np.ma.zeros((Nl_obs))
        RF_SumMtot0__l = np.ma.zeros((Nl_obs))
        OF_SumOtotR__o = np.ma.zeros((Nl_int))
        OF_SumOtotS__o = np.ma.zeros((Nl_int))
        OF_SumOtot0__o = np.ma.zeros((Nl_int))
        OF_SumMtotR__o = np.ma.zeros((Nl_int))
        OF_SumMtotS__o = np.ma.zeros((Nl_int))
        OF_SumMtot0__o = np.ma.zeros((Nl_int))
    
        # Setup stuff for histograms of R/S/U residuals 
        res_bin_low = -5.0
        res_bin_upp = 5.0
        res_bin_stp = 0.05
        res_bin_cen = np.arange(res_bin_low, res_bin_upp + res_bin_stp, res_bin_stp)
        res_bin_edges = np.append((res_bin_cen - res_bin_stp / 2.0), (res_bin_cen[-1] + res_bin_stp / 2.0))
        Nres_bins = len(res_bin_cen)
    
        histR = np.zeros((Nres_bins))     
        histS = np.zeros((Nres_bins))     
        histU = np.zeros((Nres_bins))
        
        Ng = N_gals

        for iGal in np.arange(N_gals):
            galName = listOfPrefixes[iGal][:-1]
            
            CALIFASuffix = '_synthesis_eBR_' + versRun['versionSuffix'] + versRun['othSuffix'] + versRun['baseCode'] + '.fits'
            CALIFAFitsFile = versRun['SuperFitsDir'] + galName + CALIFASuffix
            
            if not os.path.exists(CALIFAFitsFile):
                print '%s: file not found' % CALIFAFitsFile
                Ng = Ng - 1
                continue
            
            K = fitsQ3DataCube(CALIFAFitsFile)
    
            # Setup is* 1/0 flags to select (zone,lambda) pixels which are Ok in different senses; isOk__lz wraps them all together
            isFlagOk__lz = np.where(K.f_flag == 0, 1, 0)
            isErrOk__lz = np.where(K.f_err > 0, 1, 0)
            
            isWeiOk__lz = 0 * isFlagOk__lz + 1     
            
            # if input-given excludeWei0 is True than f_wei = 0 points are excluded form the stats
            # otherwise, em-lines and other masked lambdas (including the blue egde) DO enter the stats! 
            if excludeWei0:                         
                isWeiOk__lz = np.where(K.f_wei > 0, 1, 0)
            
            isOk__lz = isFlagOk__lz * isWeiOk__lz * isErrOk__lz
            
            # Define radius to zone center, in units of HLR, and a mask to pick the desired rad_low <= rad__z < rad_upp region
            rad__z = np.sqrt((K.zonePos['x'] - K.x0) ** 2. + (K.zonePos['y'] - K.y0) ** 2.) / K.HLR_pix
            isInRad__z = np.where(rad__z >= rad_low[iRad], 1, 0) * np.where(rad__z < rad_upp[iRad], 1, 0)
            
            print ' ---->', K.califaID, ' : ', radCode[iRad], ' : including ', isInRad__z.sum(), ' out of ', K.N_zone, ' zones'
            
            # Setup flag which filters zone in the irad radial bin and is Ok in lambda-space
            isOkAndInRad__lz = isOk__lz * isInRad__z
        
            # Define mean of Ok & em-line-free (continuum) fluxes for each zone. 
            # This gives an alternative normalization to fobs_norm, used in the definition of S__l
            avef_obs__z = (K.f_obs * isOk__lz).sum(axis = 0) / (isOk__lz.sum(axis = 0) + 1.e-100)
        
            # Define Otot total fluxes: OtotR = in units of fobs_norm; OtotS = in units of mean flux; Otot0 = in units of erg/s/cm2/Angs
            # Similarly, Mtot* are total Model fluxes
            RF_OtotR__lz = K.f_obs / K.fobs_norm
            RF_OtotS__lz = K.f_obs / (avef_obs__z + 1.e-100)
            RF_Otot0__lz = K.f_obs
    
            RF_MtotR__lz = K.f_syn / K.fobs_norm
            RF_MtotS__lz = K.f_syn / (avef_obs__z + 1.e-100)
            RF_Mtot0__lz = K.f_syn
    
            # Define spec resids: R = in units of fobs_norm; S = in units of mean flux; U = in units of error    
            R__lz = (K.f_obs - K.f_syn) / K.fobs_norm
            S__lz = (K.f_obs - K.f_syn) / (avef_obs__z + 1.e-100)
            U__lz = (K.f_obs - K.f_syn) / (K.f_err + 1.e-100)
    
            # Update histograms of R, S & U residuals
            histR += np.histogram(R__lz[isOkAndInRad__lz > 0], res_bin_edges)[0]
            histS += np.histogram(S__lz[isOkAndInRad__lz > 0], res_bin_edges)[0]
            histU += np.histogram(U__lz[isOkAndInRad__lz > 0], res_bin_edges)[0]
    
            # Update counters
            N_spectra += isInRad__z.sum()
            N_pixels += isOkAndInRad__lz.sum()
    
            # rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
            #         REST-FRAME stats: Update summations of R, S, U & their squares, & counter
            # rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
            RF_SumR__l += (R__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumR2__l += (R__lz ** 2 * isOkAndInRad__lz).sum(axis = 1)
            RF_SumS__l += (S__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumS2__l += (S__lz ** 2 * isOkAndInRad__lz).sum(axis = 1)
            RF_SumU__l += (U__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumU2__l += (U__lz ** 2 * isOkAndInRad__lz).sum(axis = 1)
            RF_NOk__l += isOkAndInRad__lz.sum(axis = 1)
    
            # Update summations of total fluxes (O & M)
            RF_SumOtotR__l += (RF_OtotR__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumOtotS__l += (RF_OtotS__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumOtot0__l += (RF_Otot0__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumMtotR__l += (RF_MtotR__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumMtotS__l += (RF_MtotS__lz * isOkAndInRad__lz).sum(axis = 1)
            RF_SumMtot0__l += (RF_Mtot0__lz * isOkAndInRad__lz).sum(axis = 1)
            # rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
    
            # 
            OF_l_obs = K.l_obs * (1. + K.header['REDSHIFT'])
        
            # ATT: The ugly explicit loop in z is because that's how I know how to use np.interp!
            for iz in np.arange(K.N_zone):
                OF_R__o = np.interp(l_int, OF_l_obs, R__lz[:, iz])
                OF_S__o = np.interp(l_int, OF_l_obs, S__lz[:, iz])
                OF_U__o = np.interp(l_int, OF_l_obs, U__lz[:, iz])
                OF_isOk__o = np.interp(l_int, OF_l_obs, isOk__lz[:, iz], left = 0, right = 0)
                OF_isOkAndInRad__o = np.where(OF_isOk__o >= 0.99, 1, 0) * isInRad__z[iz]
    
                OF_SumR__o += (OF_R__o * OF_isOkAndInRad__o)
                OF_SumR2__o += (OF_R__o ** 2 * OF_isOkAndInRad__o)
                OF_SumS__o += (OF_S__o * OF_isOkAndInRad__o)
                OF_SumS2__o += (OF_S__o ** 2 * OF_isOkAndInRad__o)
                OF_SumU__o += (OF_U__o * OF_isOkAndInRad__o)
                OF_SumU2__o += (OF_U__o ** 2 * OF_isOkAndInRad__o)
                OF_NOk__o += OF_isOkAndInRad__o
    
                # Update summations of total fluxes in Obs Frame (O & M)
                OF_OtotR__o = np.interp(l_int, OF_l_obs, RF_OtotR__lz[:, iz])
                OF_OtotS__o = np.interp(l_int, OF_l_obs, RF_OtotS__lz[:, iz])
                OF_Otot0__o = np.interp(l_int, OF_l_obs, RF_Otot0__lz[:, iz])
                OF_SumOtotR__o += (OF_OtotR__o * OF_isOkAndInRad__o)
                OF_SumOtotS__o += (OF_OtotS__o * OF_isOkAndInRad__o)
                OF_SumOtot0__o += (OF_Otot0__o * OF_isOkAndInRad__o)
    
                OF_MtotR__o = np.interp(l_int, OF_l_obs, RF_MtotR__lz[:, iz])
                OF_MtotS__o = np.interp(l_int, OF_l_obs, RF_MtotS__lz[:, iz])
                OF_Mtot0__o = np.interp(l_int, OF_l_obs, RF_Mtot0__lz[:, iz])
                OF_SumMtotR__o += (OF_MtotR__o * OF_isOkAndInRad__o)
                OF_SumMtotS__o += (OF_MtotS__o * OF_isOkAndInRad__o)
                OF_SumMtot0__o += (OF_Mtot0__o * OF_isOkAndInRad__o)
            #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            
            K.close()

        #===============================================================================================
        # REST-FRAME Final stats: ave & sig's of R, S & U
        RF_aveR__l = np.where(RF_NOk__l >= 1, RF_SumR__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveS__l = np.where(RF_NOk__l >= 1, RF_SumS__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveU__l = np.where(RF_NOk__l >= 1, RF_SumU__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveR2__l = np.where(RF_NOk__l >= 1, RF_SumR2__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveS2__l = np.where(RF_NOk__l >= 1, RF_SumS2__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveU2__l = np.where(RF_NOk__l >= 1, RF_SumU2__l / (RF_NOk__l + 1.e-100), 0)
        RF_sigR__l = np.sqrt(RF_aveR2__l - RF_aveR__l ** 2)
        RF_sigS__l = np.sqrt(RF_aveS2__l - RF_aveS__l ** 2)
        RF_sigU__l = np.sqrt(RF_aveU2__l - RF_aveU__l ** 2)
        
        # OBSERVED-FRAME Final stats: ave & sig's of R, S & U
        OF_aveR__o = np.where(OF_NOk__o >= 1, OF_SumR__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveS__o = np.where(OF_NOk__o >= 1, OF_SumS__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveU__o = np.where(OF_NOk__o >= 1, OF_SumU__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveR2__o = np.where(OF_NOk__o >= 1, OF_SumR2__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveS2__o = np.where(OF_NOk__o >= 1, OF_SumS2__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveU2__o = np.where(OF_NOk__o >= 1, OF_SumU2__o / (OF_NOk__o + 1.e-100), 0)
        OF_sigR__o = np.sqrt(OF_aveR2__o - OF_aveR__o ** 2)
        OF_sigS__o = np.sqrt(OF_aveS2__o - OF_aveS__o ** 2)
        OF_sigU__o = np.sqrt(OF_aveU2__o - OF_aveU__o ** 2)
        
        # Laundry
        del RF_SumR__l, RF_SumR2__l, RF_aveR2__l, RF_SumS__l, RF_SumS2__l, RF_aveS2__l, RF_SumU__l, RF_SumU2__l, RF_aveU2__l
        del OF_SumR__o, OF_SumR2__o, OF_aveR2__o, OF_SumS__o, OF_SumS2__o, OF_aveS2__o, OF_SumU__o, OF_SumU2__o, OF_aveU2__o
        
        # REST & OBS FRAME MEAN TOTAL FLUXES (O & M)
        RF_aveOtotR__l = np.where(RF_NOk__l >= 1, RF_SumOtotR__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveOtotS__l = np.where(RF_NOk__l >= 1, RF_SumOtotS__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveOtot0__l = np.where(RF_NOk__l >= 1, RF_SumOtot0__l / (RF_NOk__l + 1.e-100), 0)
        OF_aveOtotR__o = np.where(OF_NOk__o >= 1, OF_SumOtotR__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveOtotS__o = np.where(OF_NOk__o >= 1, OF_SumOtotS__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveOtot0__o = np.where(OF_NOk__o >= 1, OF_SumOtot0__o / (OF_NOk__o + 1.e-100), 0)
        
        RF_aveMtotR__l = np.where(RF_NOk__l >= 1, RF_SumMtotR__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveMtotS__l = np.where(RF_NOk__l >= 1, RF_SumMtotS__l / (RF_NOk__l + 1.e-100), 0)
        RF_aveMtot0__l = np.where(RF_NOk__l >= 1, RF_SumMtot0__l / (RF_NOk__l + 1.e-100), 0)
        OF_aveMtotR__o = np.where(OF_NOk__o >= 1, OF_SumMtotR__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveMtotS__o = np.where(OF_NOk__o >= 1, OF_SumMtotS__o / (OF_NOk__o + 1.e-100), 0)
        OF_aveMtot0__o = np.where(OF_NOk__o >= 1, OF_SumMtot0__o / (OF_NOk__o + 1.e-100), 0)
        
        del RF_SumOtotR__l, RF_SumOtotS__l, RF_SumOtot0__l
        del OF_SumOtotR__o, OF_SumOtotS__o, OF_SumOtot0__o
        #===============================================================================================
        
        print 'Found: %d galaxies' % Ng
        
        dicResidStats = {
            'RF_lambda' : K.l_obs,
            'RF_NOk__l' : RF_NOk__l,
            'RF_aveR__l' : RF_aveR__l,
            'RF_sigR__l' : RF_sigR__l,
            'RF_aveS__l' : RF_aveS__l,
            'RF_sigS__l' : RF_sigS__l,
            'RF_aveU__l' : RF_aveU__l,
            'RF_sigU__l' : RF_sigU__l,
            'OF_lambda' : l_int,
            'OF_NOk__o' : OF_NOk__o,
            'OF_aveR__o' : OF_aveR__o,
            'OF_sigR__o' : OF_sigR__o,
            'OF_aveS__o' : OF_aveS__o,
            'OF_sigS__o' : OF_sigS__o,
            'OF_aveU__o' : OF_aveU__o,
            'OF_sigU__o' : OF_sigU__o,
            'N_gals': Ng,
            'N_spectra' : N_spectra,
            'N_pixels' : N_pixels,
            'res_bin_cen' : res_bin_cen,
            'histR' : histR,
            'histS' : histS,
            'histU' : histU,
            'RF_aveOtotR__l' : RF_aveOtotR__l,
            'RF_aveOtotS__l' : RF_aveOtotS__l,
            'RF_aveOtot0__l' : RF_aveOtot0__l,
            'OF_aveOtotR__o' : OF_aveOtotR__o,
            'OF_aveOtotS__o' : OF_aveOtotS__o,
            'OF_aveOtot0__o' : OF_aveOtot0__o,
            'RF_aveMtotR__l' : RF_aveMtotR__l,
            'RF_aveMtotS__l':RF_aveMtotS__l,
            'RF_aveMtot0__l' : RF_aveMtot0__l,
            'OF_aveMtotR__o' : OF_aveMtotR__o,
            'OF_aveMtotS__o' : OF_aveMtotS__o,
            'OF_aveMtot0__o' : OF_aveMtot0__o,
        }
                
        suf = excludeWei0Code + '.' + radCode[iRad]
        D = dicResidStats

        fname1 = 'SpecResidStats4DR2_RestFrame.' + suf
        fname2 = 'SpecResidStats4DR2_ObsFrame.' + suf
        f1 = open(fname1, 'w')
        f2 = open(fname2, 'w')
        tabHeader = '# lambda   N   aveR         sigR         aveS         sigS         aveU         sigU         aveOtotR     aveOtotS     aveOtot0     aveMtotR     aveMtotS     aveMtot0\n'
        f1.write(tabHeader)
        f2.write(tabHeader)
        fmt = '%4i    %4i   ' + 12 * '%.5e  ' + ' \n'
        for i in range(K.Nl_obs):
            f1.write(fmt % (D['RF_lambda'][i], D['RF_NOk__l'][i], \
                             D['RF_aveR__l'][i], D['RF_sigR__l'][i], \
                             D['RF_aveS__l'][i], D['RF_sigS__l'][i], \
                             D['RF_aveU__l'][i], D['RF_sigU__l'][i], \
                             D['RF_aveOtotR__l'][i], D['RF_aveOtotS__l'][i], D['RF_aveOtot0__l'][i], \
                             D['RF_aveMtotR__l'][i], D['RF_aveMtotS__l'][i], D['RF_aveMtot0__l'][i]))            
        for i in range(Nl_int):
            f2.write(fmt % (D['OF_lambda'][i], D['OF_NOk__o'][i], \
                             D['OF_aveR__o'][i], D['OF_sigR__o'][i], \
                             D['OF_aveS__o'][i], D['OF_sigS__o'][i], \
                             D['OF_aveU__o'][i], D['OF_sigU__o'][i], \
                             D['OF_aveOtotR__o'][i], D['OF_aveOtotS__o'][i], D['OF_aveOtot0__o'][i], \
                             D['OF_aveMtotR__o'][i], D['OF_aveMtotS__o'][i], D['OF_aveMtot0__o'][i]))
        f1.close()
        f2.close()
    
        # R, S & U histograms
        fname3 = 'SpecResidStats4DR2_HistsRSU.' + suf
        f3 = open(fname3, 'w')
        tabHeader = '# res_bin_cen histR        histS        histU \n'
        f3.write(tabHeader)
        fmt = '%.5e  ' * 4 + ' \n'
        
        for i in range(Nres_bins):
            f3.write(fmt % (res_bin_cen[i], histR[i], histS[i], histU[i]))
        f3.close()
    
        # npz file with the full Dictionary
        fname4 = 'SpecResidStats4DR2_AllDict.' + suf
        np.savez_compressed(fname4, D)
