#!/usr/bin/python
#
# Lacerda@Saco - 19/Aug/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
#from scipy.optimize import curve_fit
from astropy.modeling import models, fitting

mpl.rcParams['font.family'] = 'sans-serif'

fit = False
#fit = True

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# # Define model function to be used to fit to the data above:
# def gauss(x, *p):
#     A, mu, sigma = p
#     return A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

def plotFitHisto(binCenter, histo, title, fileName, fit = False):
    f = plt.figure()
    f.set_dpi(92)
    f.set_size_inches(2 * 3.45, 6)
    
    ax = plt.gca()
    ax.bar(binCenter, histo, width = 0.08, edgecolor = 'black', color = 'lightgrey')
    ax.set_title(title)
    ax.set_xlabel(r'$(O_{\lambda k}\ -\ M_{\lambda k}) / \epsilon_{\lambda k}$')
    ax.set_ylabel(r'Number of pixels')

    if fit:
        # Fit data
        g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g_init, binCenter, histo)
        ax.plot(binCenter, g(binCenter))
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # coeff, var_matrix = curve_fit(gauss, binCenter, histo, p0 = [1., 0., 1.])
        # ax.plot(binCenter, gauss(binCenter, *coeff), 'b-')
        # ax.text(0.96, 0.96, '$A=%.2f$\n$x_0=%.2f$\n$\sigma=%.2f$' % (coeff[0], coeff[1], coeff[2]),
        #         fontsize = 12, transform = ax.transAxes, va = 'top', ha = 'right',
        #         bbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.))
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    f.tight_layout()
    f.savefig(fileName)
    plt.close(f)

if __name__ == '__main__':
    for wei in ['w0', 'w1']:
        for rad in ['Galaxy', 'Nucleus', 'Bulge', 'Disc']:
            fileName = 'SpecResidStats4DR2_HistsRSU.%s.%s' % (wei, rad)
            binCenter, histoR, histoS, histoU = np.loadtxt(fileName, unpack = True)
            #fit = True
            plotFitHisto(binCenter, histoR, rad, 'histoR-%s-%s.png' % (wei, rad), fit = fit)
            plotFitHisto(binCenter, histoS, rad, 'histoS-%s-%s.png' % (wei, rad), fit = fit)
            plotFitHisto(binCenter, histoU, rad, 'histoU-%s-%s.png' % (wei, rad), fit = fit)