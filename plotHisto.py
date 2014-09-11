#!/usr/bin/python
#
# Lacerda@Saco - 19/Aug/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
#from astropy.modeling import models, fitting

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'
#versionSuffix = 'v20_q036.d13c'
#versionSuffix = 'v20_q043.d14a'
versionSuffix = 'v20_q046.d15a'
#versionSuffix = 'px1_q043.d14a'
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
    ax.bar(binCenter, histo, width = 0.05, edgecolor = 'black', color = 'lightgrey')
    #ax.set_title(title)
    
    print fileName, histo.sum()
    
    if fit:
        coeff, var_matrix = curve_fit(gauss, binCenter, histo, p0 = [1000., 0., 1.])
        y = gauss(binCenter, *coeff)
        ax.plot(binCenter, y, 'b-')
        txt = '$A=%.2f$\n$x_0=%.2f$\n$\sigma=%.2f$' % (coeff[0], coeff[1], np.abs(coeff[2]))
        print txt
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # ax.text(0.96, 0.96, txt,
        #         fontsize = 12, transform = ax.transAxes, va = 'top', ha = 'right',
        #         bbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.))
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    ax.set_xlabel(r'$(O_{\lambda k}\ -\ M_{\lambda k}) / \epsilon_{\lambda k}$')
    ax.set_ylabel(r'Number of pixels')

    f.tight_layout()
    f.savefig(fileName)
    plt.close(f)

if __name__ == '__main__':
    for wei in ['w0', 'w1']:
        for rad in ['Galaxy', 'Nucleus', 'Bulge', 'Disc']:
        #for rad in ['Galaxy']:
            fileName = 'SpecResidStats4DR2_HistsRSU.%s.%s' % (wei, rad)
            binCenter, histoR, histoS, histoU = np.loadtxt(fileName, unpack = True)
            #plotFitHisto(binCenter, histoR, '%s %s' % (rad, versionSuffix), 'histoR.%s.%s.%s.%s' % (versionSuffix, wei, rad, outputImgSuffix), fit = fit)
            #plotFitHisto(binCenter, histoS, '%s %s' % (rad, versionSuffix), 'histoS.%s.%s.%s.%s' % (versionSuffix, wei, rad, outputImgSuffix), fit = fit)
            plotFitHisto(binCenter, histoU, '%s %s' % (rad, versionSuffix), 'histoU.%s.%s.%s.%s' % (versionSuffix, wei, rad, outputImgSuffix), fit = fit)