#!/usr/bin/python
#
# Lacerda@Saco - 21/Aug/2014
#
# This script creates the SpecResidStats.{version}.RestFrame.{output} image.
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, MaxNLocator

mpl.rcParams['font.size']       = 14
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

outputImgSuffix = 'pdf'
#versionSuffix = 'v20_q036.d13c'
versionSuffix = 'v20_q043.d14a'
#versionSuffix = 'v20_q046.d15a'
#versionSuffix = 'px1_q043.d14a'

rad = 'Galaxy'
llow = 3700 
lup = 6850

def fixLocatorsAx(ax):
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(100))

def plotObsNormAx(ax, data):
    x = getAttribFromData(data, 'w1', 'l')
    m = (x >= llow) & (x <= lup)
    aveOtotR = getAttribFromData(data, 'w1', 'aveOtotR')
    aveMtotR = getAttribFromData(data, 'w1', 'aveMtotR')
    N = getAttribFromData(data, 'w1', 'N')
    y = np.ma.masked_array(aveOtotR, mask = (N == 0))
    y2 = np.ma.masked_array(aveMtotR, mask = (N == 0))
    R = y - y2
    ax.plot(x[m], y2[m] - 0.1, 'y-', lw = 1.5)
    ax.plot(x[m], y[m], color = 'k', lw = 1.5)
    ax.plot(x[m], R[m], color = 'r', lw = 1.0)
    ax.set_xlim([llow, lup])
    ax.set_ylim([-0.05, 1.5])
    ax.yaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
    ax.set_ylabel(r'$\langle O_\lambda / O_{\lambda5635} \rangle$')
    plt.setp(ax.get_xticklabels(), visible = False)        
    fixLocatorsAx(ax)

def plotRAx(ax, data):
    x = getAttribFromData(data, 'w0', 'l')
    m = (x >= llow) & (x <= lup)
    aveR = getAttribFromData(data, 'w0', 'aveR')
    N = getAttribFromData(data, 'w0', 'N')
    y = np.ma.masked_array(aveR, mask = (N == 0))
    ax.plot(x[m], y[m] * 100., lw = 1.5)
    ax.set_xlim([llow, lup])
    ax.set_ylim([-7, 7])
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.axhline(-3, color='k', ls = '--', lw = 0.1)
    ax.axhline(0, color='k', ls = '--', lw = 0.1)
    ax.axhline(3, color='k', ls = '--', lw = 0.1)
    ax.set_ylabel(r'$\langle R_\lambda  / O_{\lambda5635} \rangle\ [\%]$')
    plt.setp(ax.get_xticklabels(), visible = False)
    fixLocatorsAx(ax)        

def plotNRatioAx(ax, data):
    x = getAttribFromData(data, 'w0', 'l')
    m = (x >= llow) & (x <= lup)
    y = getAttribFromData(data, 'w0', 'N_Ntot')
    ax.plot(x[m], y[m], lw = 1, c = 'k')
    #ax.plot(x[m], y[m] / 169577., lw = 1, c = 'k')
    ax.set_ylabel(r'$N^{Ok}_\lambda / N_{tot}$')
    ax.set_xlabel(r'rest-frame wavelength $[\AA]$')
    ax.set_xlim([llow, lup])
    ax.set_ylim([-0.05, 1.15])
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    #fixLocatorsAx(ax)

def getAttribFromData(data, wei = None, attrib = None):
    attribDict = dict(l = 0, N = 1, N_Ntot = 2, 
                      aveR = 3, sigR = 4,
                      aveS = 5, sigS = 6, 
                      aveU = 7, sigU = 8,
                      aveOtotR = 9, aveOtotS = 10, aveOtot0 = 11,
                      aveMtotR = 12, aveMtotS = 13, aveMtot0 = 14)
    if wei:
        weiDict = dict(w0 = 0, w1 = 1)
        iWei = weiDict[wei]
        
    iAttrib = attribDict[attrib]

    if type(data) is np.ndarray:
        return data[iAttrib]
    elif type(data) is not np.ndarray and wei:
        return data[iWei][iAttrib] 
    else: 
        return None
    
if __name__ == '__main__':
    data = []
    for wei in ['w0', 'w1']:
        fileName = 'SpecResidStats4DR2_RestFrame.%s.Galaxy' % wei
        aux = np.loadtxt(fileName, unpack = True)
        data.append(aux) 
    
    f = plt.figure()
    gs = gridspec.GridSpec(4, 1)
    f.set_size_inches(10, 10)
    
    ax = plt.subplot(gs[0:2, 0])
    plotObsNormAx(ax, data)
    ax = plt.subplot(gs[2, 0])
    plotRAx(ax, data)
    ax = plt.subplot(gs[3, 0])
    plotNRatioAx(ax, data)
    plt.subplots_adjust(hspace=0.)
    plt.setp(ax.get_xticklabels(), visible = True)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = True)
    f.savefig('SpecResidStats.%s.RestFrame.%s' % (versionSuffix, outputImgSuffix))