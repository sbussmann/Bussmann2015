"""

2014 July 16
Shane Bussmann

Plot the distribution of fluxdensities for the ALMA sample.  Compare total
observed flux (what a single-dish telescope with 20" FWHM resolution would see)
with the individual observed flux (accounting for blending) and with the
individual intrinsic flux (accounting for lensing).

"""

import matplotlib.pyplot as plt
import getfluxes
import numpy
from pylab import savefig


def getstat(fluxdist):
    fmean = fluxdist.mean()
    fstd = fluxdist.std()
    return fmean, fstd

# f1: total observed flux density
# f2: per-source observed flux density
# f3: per-source intrinsic flux density
f1, f2, f3 = getfluxes.all('uvfit25', mu_estimate=True)
fnu = f1['fnu']
fnu = numpy.append(fnu, f2['fnu'])
fnu = numpy.append(fnu, f3['fnu'])
m1 = 0
m2 = 1 + int(fnu.max())

binwidth = 2
bins = numpy.arange(m1, m2 + binwidth, binwidth)
fluxarr = f1['fnu']
plt.hist(fluxarr, bins = bins, histtype='stepfilled', edgecolor='black',
        hatch='///', facecolor='none', label='Observed, Unresolved')
fmean, fstd = getstat(fluxarr)
strmean = '{0:.1f}'.format(fmean)
strstd = '{0:.1f}'.format(fstd)
s1 = 'Observed, Unresolved <S_870> = ' + strmean + ' +/- ' + strstd
plt.text(65, 13, s1, ha='right', fontsize='large')

binwidth = 2
bins = numpy.arange(m1, m2 + binwidth, binwidth)
fluxarr = f2['fnu']
plt.hist(fluxarr, bins = bins, alpha=0.5, histtype='stepfilled', 
        color='red', label='Observed, Resolved')
fmean, fstd = getstat(fluxarr)
strmean = '{0:.1f}'.format(fmean)
strstd = '{0:.1f}'.format(fstd)
s1 = 'Observed, Resolved <S_870> = ' + strmean + ' +/- ' + strstd
plt.text(65, 11, s1, ha='right', fontsize='large')

binwidth = 2
fluxarr = f3['fnu']# / mu
bins = numpy.arange(m1, m2 + binwidth, binwidth)
plt.hist(fluxarr, bins = bins, alpha=0.2, histtype='stepfilled',
        color='blue', label='Intrinsic, Resolved')
fmean, fstd = getstat(fluxarr)
strmean = '{0:.1f}'.format(fmean)
strstd = '{0:.1f}'.format(fstd)
s1 = 'Intrinsic, Resolved <S_870> = ' + strmean + ' +/- ' + strstd
plt.text(65, 9, s1, ha='right', fontsize='large')

plt.xlabel(r'$S_{870} \, ({\rm mJy})$', fontsize='x-large')
plt.ylabel('N', fontsize='x-large')
plt.tight_layout()
plt.legend()
savefig('fluxdist.png')
