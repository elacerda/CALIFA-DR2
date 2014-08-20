#!/usr/bin/python
#
# Lacerda@Saco - 19/Aug/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

mpl.rcParams['font.family'] = 'sans-serif'

radCode = ['Galaxy', 'Nucleus', 'Bulge', 'Disc']

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A * np.exp(-(x - mu)**2 / (2. * sigma**2))

def plotFitHisto(binCenter, histo, title, fileName, fit = False):
    f = plt.figure()
    f.set_dpi(92)
    ratio = 3.45 / 4.
    f.set_size_inches(2 * 3.45, 6)
    
    ax = plt.gca()
    ax.bar(binCenter, histo, width = 0.08, edgecolor = 'black', color = 'lightgrey')
    ax.set_title(title)
    ax.set_xlabel(r'$(O_{\lambda k}\ -\ M_{\lambda k}) / \epsilon_{\lambda k}$')
    ax.set_ylabel(r'Number of pixels')

    if fit:
        # Fit data
        coeff, var_matrix = curve_fit(gauss, binCenter, histo, p0=[1., 0., 1.])
        ax.plot(binCenter, gauss(binCenter, *coeff), 'b-')
        ax.text(0.96, 0.96,'$A=%.2f$\n$x_0=%.2f$\n$\sigma=%.2f$' % (coeff[0], coeff[1], coeff[2]), 
                fontsize=12, transform = ax.transAxes, va = 'top', ha = 'right', 
                bbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.))

    f.tight_layout()
    f.savefig(fileName)

if __name__ == '__main__':
    for rad in radCode:
        fileName = 'SpecResidStats4DR2_HistsRSU.w1.%s.csv' % rad
        binCenter, histoR, histoS, histoU = np.loadtxt(fileName, unpack = True)
        
        plotFitHisto(binCenter, histoR, rad, 'histoR-%s.png' % rad)
        plotFitHisto(binCenter, histoS, rad, 'histoS-%s.png' % rad)
        plotFitHisto(binCenter, histoU, rad, 'histoU-%s.png' % rad)
