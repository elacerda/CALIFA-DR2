#!/usr/bin/python
#
# Lacerda@Saco - 20/Aug/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, MaxNLocator
from CALIFAUtils import paths
from CALIFAUtils.plots import add_subplot_axes

mpl.rcParams['font.size']       = 14
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

outputImgSuffix = 'pdf'

paths.set_v_run(0)
#versionSuffix = paths.get_config()['versionSuffix']
versionSuffix = 'v20_q053.d22a512.mE' 

radCode = ['Galaxy', 'Nucleus', 'Bulge', 'Disc']
#radCode = ['Galaxy']
radColor = dict(Galaxy = 'k', Nucleus = 'r', Bulge = 'y', Disc = 'b')
llow = 3800 
lup = 7200

###################################################################
###################################################################
###################################################################
def fixLocatorsAx(ax, prune):
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MaxNLocator(7, prune = prune))
    ax.xaxis.grid(which = 'minor')
    ax.yaxis.grid()

def plotObsNormAx(ax, data):
    for rad in radCode:
        color = radColor[rad]
        x = getAttribFromData(data, rad, 'l')
        m = (x >= llow) & (x <= lup)
        aveOtotR = getAttribFromData(data, rad, 'aveOtotR')
        N = getAttribFromData(data, rad, 'N')
        y = np.ma.masked_array(aveOtotR, mask = (N == 0))
        ax.set_xlim(3750, 6900) 
        ax.plot(x[m], y[m], color = color, label = rad, lw = 1.5)
        
    ax.legend(loc='upper center', bbox_to_anchor=(0.93, 0.74), prop={'size' : 10})
    ax.set_ylabel(r'$\langle O_\lambda / O_{\lambda5635} \rangle$')
    fixLocatorsAx(ax, 'lower')
    plt.setp(ax.get_xticklabels(), visible = False)        
        
def plotRAx(ax, data):
    for rad in radCode:
        color = radColor[rad]
        x = getAttribFromData(data, rad, 'l')
        m = (x >= llow) & (x <= lup)
        aveR = getAttribFromData(data, rad, 'aveR')
        N = getAttribFromData(data, rad, 'N')
        y = np.ma.masked_array(aveR, mask = (N == 0))
        ax.set_xlim(3750, 6900) 
        ax.plot(x[m], y[m], color = color, label = rad, lw = 1.5)
    
    ax.set_ylabel(r'$\langle R_\lambda \rangle$')
    fixLocatorsAx(ax, 'lower')
    plt.setp(ax.get_xticklabels(), visible = False)        
    
def plotUAx(ax, data):
    for rad in radCode:
        color = radColor[rad]
        x = getAttribFromData(data, rad, 'l')
        m = (x >= llow) & (x <= lup)
        aveU = getAttribFromData(data, rad, 'aveU')
        N = getAttribFromData(data, rad, 'N')
        y = np.ma.masked_array(aveU, mask = (N == 0))
        ax.set_xlim(3750, 6900) 
        ax.plot(x[m], y[m], color = color, label = rad, lw = 1.5)
    
    ax.set_ylabel(r'$\langle U_\lambda \rangle$')
    fixLocatorsAx(ax, 'lower')
    plt.setp(ax.get_xticklabels(), visible = False)        

def plotSigmaRAx(ax, data):
    #subax = add_subplot_axes(ax, [0.65,  0.25, 0.30, 0.65])
    for rad in radCode:
        color = radColor[rad]
        x = getAttribFromData(data, rad, 'l')
        m = (x >= llow) & (x <= lup)
        sigR = getAttribFromData(data, rad, 'sigR')
        N = getAttribFromData(data, rad, 'N')
        y = np.ma.masked_array(sigR, mask = (N == 0))
        #ax.set_ylim(0, 0.14)
        ax.set_xlim(3750, 6900)
        ax.plot(x[m], y[m], color = color, label = rad, lw = 1.5)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         subax.plot(x[m], y[m], color = color, label = rad, lw = 1.5)
#     subax.set_ylim(0, 0.14)
#     subax.set_xlim(3750, 4600)
#     subax.xaxis.set_major_locator(MultipleLocator(500))
#     subax.xaxis.set_minor_locator(MultipleLocator(100))
#     subax.yaxis.set_major_locator(MaxNLocator(7, prune = 'lower'))
#     subax.xaxis.grid(which = 'minor')
#     subax.yaxis.grid()
#     plt.setp(subax.get_xticklabels(), visible = True)
#     subax.tick_params(axis='both', which='major', labelsize=10)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ax.set_ylabel(r'$\sigma(R_\lambda)$')
    fixLocatorsAx(ax, 'lower')
    plt.setp(ax.get_xticklabels(), visible = False)        

def plotSigmaUAx(ax, data):
    for rad in radCode:
        color = radColor[rad]
        x = getAttribFromData(data, rad, 'l')
        m = (x >= llow) & (x <= lup)
        sigU = getAttribFromData(data, rad, 'sigU')
        N = getAttribFromData(data, rad, 'N')
        y = np.ma.masked_array(sigU, mask = (N == 0))
        ax.set_xlim(3750, 6900)
        ax.plot(x[m], y[m], color = color, label = rad, lw = 1.5) 
    
    ax.set_ylabel(r'$\sigma(U_\lambda)$')
    fixLocatorsAx(ax, 'lower')
    plt.setp(ax.get_xticklabels(), visible = False)        

def plotNRatioAx(ax, data):
    for rad in radCode:
        color = radColor[rad]
        x = getAttribFromData(data, rad, 'l')
        m = (x >= llow) & (x <= lup)
        y = getAttribFromData(data, rad, 'N_Ntot')
        ax.set_xlim(3750, 6900)
        #ax.plot(x[m], y[m] / (y[m]).sum(), color = color, label = rad, lw = 1.5)
        ax.plot(x[m], y[m], color = color, label = rad, lw = 1.5)
    
    ax.set_ylabel(r'$N^{Ok}_\lambda / N_{tot}$')
    ax.set_xlabel(r'wavelength')
    fixLocatorsAx(ax, None)

def getAttribFromData(data, rad, attrib):
    radDict = dict(Galaxy = 0, Nucleus = 1, Bulge = 2, Disc = 3)
    attribDict = dict(l = 0, N = 1, N_Ntot = 2, 
                      aveR = 3, sigR = 4,
                      aveS = 5, sigS = 6, 
                      aveU = 7, sigU = 8,
                      aveOtotR = 9, aveOtotS = 10, aveOtot0 = 11,
                      aveMtotR = 12, aveMtotS = 13, aveMtot0 = 14)
    iRad = radDict[rad]
    iAttrib = attribDict[attrib]

    if type(data) is np.ndarray:
        return data[iAttrib]
    else:
        return data[iRad][iAttrib] 

###################################################################
###################################################################
###################################################################
    
if __name__ == '__main__':
    for frame in [ 'ObsFrame', 'RestFrame' ]:
    #for frame in [ 'RestFrame' ]:
        #for wei in ['w0', 'w1']:
        for wei in ['w0']:
            data = []
            
            for rad in radCode:
                fileName = 'SpecResidStats4DR3_%s.%s.%s' % (frame, wei, rad)
                aux = np.loadtxt(fileName, unpack = True)
                data.append(aux)
                
            f, axArr = plt.subplots(6, 1)
            f.suptitle(r'Spec Resid Stats %s %s %s' % (frame, versionSuffix, wei))
            f.set_size_inches(10, 10)
            
            plotObsNormAx(axArr[0], data)
            plotRAx(axArr[1], data)
            plotSigmaRAx(axArr[2], data)
            plotUAx(axArr[3], data)
            plotSigmaUAx(axArr[4], data)
            plotNRatioAx(axArr[5], data)
            
            plt.setp(axArr[5].get_xticklabels(), visible = True)
            plt.setp([a.get_yticklabels() for a in f.axes], visible = True)
            f.savefig('SpecResidStats.%s.%s.%s.%s' % (versionSuffix, frame, wei, outputImgSuffix))
