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
from CALIFAUtils import paths

mpl.rcParams['font.size']       = 14
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

outputImgSuffix = 'pdf'

paths.set_v_run(-1)
#versionSuffix = paths.get_config()['versionSuffix']
versionSuffix = 'v20_q053.d22a512.me' 

rad = 'Galaxy'
llow = 3700 
lup = 6850

def fixLocatorsAx(ax):
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(100))

def plotObsNormAx(ax, data, c = ['y-', 'k', 'r'], lw = [1.5, 1.5, 1.], ver = 'q050', plus = 0):
    x = getAttribFromData(data, 'w1', 'l')
    m = (x >= llow) & (x <= lup)
    aveOtotR = getAttribFromData(data, 'w1', 'aveOtotR')
    aveMtotR = getAttribFromData(data, 'w1', 'aveMtotR')
    N = getAttribFromData(data, 'w1', 'N')
    y = np.ma.masked_array(aveOtotR, mask = (N == 0))
    y2 = np.ma.masked_array(aveMtotR, mask = (N == 0))
    R = y - y2
    #ax.plot(x[m], y2[m] - 0.1, c[0], lw = lw[0], label = r'$M_\lambda$ %s' % ver)
    ax.plot(x[m], y2[m] + plus, c[0], lw = lw[0], label = r'$M_\lambda$ %s' % ver)
    ax.plot(x[m], y[m] + plus, c[1], lw = lw[1], label = r'$O_\lambda$ %s' % ver)
    ax.plot(x[m], R[m], c[2], lw = lw[2], label = r'$R_\lambda$ %s' % ver)
    ax.set_xlim([llow, lup])
    ax.set_ylim([-0.05, 1.5])
    #ax.legend(loc = 'upper left', fontsize = 10, frameon = False)
    ax.yaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
    ax.set_ylabel(r'$\langle O_\lambda / O_{\lambda5635} \rangle$')
    plt.setp(ax.get_xticklabels(), visible = False)        
    fixLocatorsAx(ax)
    return m, x, y, y2, R

def plotRAx(ax, data, c = 'k', ls = '-', lw = 1.5, plot_hlines = True):
    x = getAttribFromData(data, 'w0', 'l')
    m = (x >= llow) & (x <= lup)
    aveR = getAttribFromData(data, 'w0', 'aveR')
    N = getAttribFromData(data, 'w0', 'N')
    if plot_hlines is True:
        ax.fill_between(x[m], -3, 3, color = '#7B8DBF')
        #ax.axhline(-3, c = 'w', ls = '--', lw = 1)
        ax.axhline(0, c = 'w', ls = '--', lw = 0.5)
        #ax.axhline(3, c ='w', ls = '--', lw = 1)
    y = np.ma.masked_array(aveR, mask = (N == 0))
    ax.plot(x[m], y[m] * 100., c, lw = lw, ls = ls)
    m2 = np.bitwise_and(m, np.greater_equal(y, 0.03))
    if (c == 'w'):
        c2 = '#7B8DBF'
    else:
        c2 = c
    ax.plot(np.ma.masked_array(x, mask = ~m2), np.ma.masked_array(y, mask = ~m2) * 100., c = c2, ls = ls, lw = lw)
    ax.set_xlim([llow, lup])
    ax.set_ylim([-7, 7])
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.set_ylabel(r'$\langle R_\lambda  / O_{\lambda5635} \rangle\ [\%]$')
    plt.setp(ax.get_xticklabels(), visible = False)
    fixLocatorsAx(ax)        

def plotNRatioAx(ax, data, c = 'k'):
    x = getAttribFromData(data, 'w0', 'l')
    m = (x >= llow) & (x <= lup)
    y = getAttribFromData(data, 'w0', 'N_Ntot')
    #ax.plot(x[m], y[m], lw = 1, c = c)
    ax.fill_between(x[m], 0, y[m], color = c)
    #ax.plot(x[m], y[m] / 169577., lw = 1, c = 'k')
    ax.set_ylabel(r'$N^{Ok}_\lambda / N_{tot}$')
    ax.set_xlabel(r'rest-frame wavelength $[\AA]$')
    ax.set_xlim([llow, lup])
    ax.set_ylim([-0.05, 1.15])
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    #fixLocatorsAx(ax)

def plotRR(ax, data, data_old):
    x = getAttribFromData(data, 'w1', 'l')
    aveOtotR = getAttribFromData(data, 'w1', 'aveOtotR')
    aveMtotR = getAttribFromData(data, 'w1', 'aveMtotR')
    N = getAttribFromData(data, 'w1', 'N')
    y = np.ma.masked_array(aveOtotR, mask = (N == 0))
    y2 = np.ma.masked_array(aveMtotR, mask = (N == 0))
    R = y - y2
    x = getAttribFromData(data_old, 'w1', 'l')
    m = (x >= llow) & (x <= lup)
    N = getAttribFromData(data_old, 'w1', 'N')
    ax.plot(x[m], (R - (np.ma.masked_array(getAttribFromData(data_old, 'w1', 'aveOtotR'), mask = (N == 0)) - np.ma.masked_array(getAttribFromData(data_old, 'w1', 'aveMtotR'), mask = (N == 0))))[m], 'b', lw = 1., label = r'$R_R$')
    
    
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
    #data_old = []
    #for wei in ['w0', 'w1']:
    for wei in ['w0']:
        fileName = 'SpecResidStats4DR3_RestFrame.%s.Galaxy' % wei
        aux = np.loadtxt(fileName, unpack = True)
        data.append(aux) 
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # fileName = '../v20_q046.d15a/SpecResidStats4DR2_RestFrame.%s.Galaxy' % wei
        # aux = np.loadtxt(fileName, unpack = True)
        # data_old.append(aux) 
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    f = plt.figure()
    gs = gridspec.GridSpec(4, 1)
    f.set_size_inches(10, 10)
    
    ax = plt.subplot(gs[0:2, 0])
    m, wl, O, M, R = plotObsNormAx(ax, data, ['g', '#7B8DBF', '#7B8DBF'], [1., 1., 0.3], 'q050')
    #m_old, wl_old, O_old, M_old, R_old = plotObsNormAx(ax, data_old, ['y', '#9FBAE2', '#9FBAE2'], [1., 1., 0.3], 'q046', 0.5)
    #ax.fill_between(wl[m], R[m], R_old[m], color = '#9FBAE2', edgecolor = 'none')
    #ax.fill_between(wl[m], R[m], R_old[m], color = 'k', edgecolor = 'none')
#    plotRR(ax, data, data_old)
    
    ax = plt.subplot(gs[2, 0])
    plotRAx(ax, data, 'w', '-', 0.8)
    #plotRAx(ax, data_old, '#9FBAE2', '-', 0.8, False)

    ax = plt.subplot(gs[3, 0])
    plotNRatioAx(ax, data, '#7B8DBF')
    #plotNRatioAx(ax, data_old, '#9FBAE2')
    plt.subplots_adjust(hspace=0.)
    plt.setp(ax.get_xticklabels(), visible = True)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = True)
    #f.suptitle(r'Spec Resid Stats Restframe %s' % versionSuffix)
    f.savefig('SpecResidStats.%s.RestFrame.%s' % (versionSuffix, outputImgSuffix))  