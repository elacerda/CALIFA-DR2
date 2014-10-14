import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

mpl.rcParams['font.size']       = 14
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

outputImgSuffix = 'pdf'

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

data_v14 = np.loadtxt('v20_q043.d14a_185gal/SpecResidStats4DR2_RestFrame.w0.Galaxy', unpack = True)
data_v15 = np.loadtxt('v20_q046.d15a/SpecResidStats4DR2_RestFrame.w0.Galaxy', unpack = True)

llow = 3700
lup = 6850

# X and Y
x = getAttribFromData(data_v14, attrib = 'l')
y_v14 = getAttribFromData(data_v14, attrib = 'aveR')
y_v15 = getAttribFromData(data_v15, attrib = 'aveR')

# masks
m = (x >= llow) & (x <= lup)
N_v14 = getAttribFromData(data_v14, attrib = 'N')
N_v15 = getAttribFromData(data_v15, attrib = 'N')
yv14 = np.ma.masked_array(y_v14, mask = (N_v14 == 0))
yv15 = np.ma.masked_array(y_v15, mask = (N_v15 == 0))

# plot
f = plt.figure()
f.set_size_inches(15, 5)
ax = f.gca()
ax.plot(x[m], yv14[m] * 100., lw = 1.5, label = 'v14')
ax.plot(x[m], yv15[m] * 100., lw = 1.5, label = 'v15')
ax.legend()
ax.set_xlim([llow, lup])
ax.set_ylim([-7, 7])
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.axhline(-3, color='k', ls = '--', lw = 0.1)
ax.axhline(0, color='k', ls = '--', lw = 0.1)
ax.axhline(3, color='k', ls = '--', lw = 0.1)
ax.set_xlabel(r'rest-frame wavelength $[\AA]$')
ax.set_ylabel(r'$\langle R_\lambda  / O_{\lambda5635} \rangle\ [\%]$')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
f.savefig('aveR_v14_v15.png')
