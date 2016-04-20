#!/usr/bin/python
#
# Lacerda@Saco - 19/Aug/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
#from astropy.modeling import models, fitting
from CALIFAUtils import paths

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

paths.set_v_run(-1)
#versionSuffix = paths.get_config()['versionSuffix'] 
versionSuffix = 'v20_q053.d22a512.me' 

outputImgSuffix = 'pdf'

#fit = False
fit = True

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

def plotFitHisto(binCenter, histo, title, fileName, fit = False):
    f = plt.figure()
    f.set_dpi(92)
    f.set_size_inches(2 * 3.45, 6)
    
    ax = plt.gca()
    #ax.bar(binCenter, histo, width = 0.05, edgecolor = 'none', color = mpl.cm.cool(0.7))
    #ax.bar(binCenter, histo, width = 0.05, linewidth = .5, edgecolor = mpl.cm.cool(0.4), facecolor = mpl.cm.cool(0.4))
    c = '#9FBAE2'
    #c = '#7B8DBF'
    ax.bar(binCenter, histo, width = 0.05, linewidth = .5, edgecolor = c, facecolor =c)
    #ax.set_title(title)
    
    print fileName, histo.sum()
    
    if fit:
        coeff, _ = curve_fit(gauss, binCenter, histo, p0 = [1000., 0., 1.])
        y = gauss(binCenter, *coeff)
        ax.plot(binCenter, y, color = '#BA3E04', lw = 1.5)
        txt = '$A=%.2f$\n$x_0=%.2f$\n$\sigma=%.2f$' % (coeff[0], coeff[1], np.abs(coeff[2]))
        print txt
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # ax.text(0.96, 0.96, txt,
        #         fontsize = 12, transform = ax.transAxes, va = 'top', ha = 'right',
        #         bbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.))
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        
    f.suptitle(r'A:%.2f  $x_0$:%.2f  $\sigma$:%.2f  $N_p$:%d' % (coeff[0], coeff[1], np.abs(coeff[2]), histo.sum()))

    ax.set_xlabel(r'$(O_{\lambda k}\ -\ M_{\lambda k}) / \epsilon_{\lambda k}$')
    ax.set_ylabel(r'Number of pixels')

    f.tight_layout()
    f.savefig(fileName)
    plt.close(f)

if __name__ == '__main__':
    for wei in ['w0', 'w1']:
        for rad in ['Galaxy', 'Nucleus', 'Bulge', 'Disc']:
        #for rad in ['Galaxy']:
            fileName = 'SpecResidStats4DR3_HistsRSU.%s.%s' % (wei, rad)
            binCenter, histoR, histoS, histoU = np.loadtxt(fileName, unpack = True)
            plotFitHisto(binCenter, histoR, '%s %s' % (rad, versionSuffix), 'histoR.%s.%s.%s.%s' % (versionSuffix, wei, rad, outputImgSuffix), fit = fit)
            plotFitHisto(binCenter, histoS, '%s %s' % (rad, versionSuffix), 'histoS.%s.%s.%s.%s' % (versionSuffix, wei, rad, outputImgSuffix), fit = fit)
            plotFitHisto(binCenter, histoU, '%s %s' % (rad, versionSuffix), 'histoU.%s.%s.%s.%s' % (versionSuffix, wei, rad, outputImgSuffix), fit = fit)