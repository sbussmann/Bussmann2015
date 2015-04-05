"""

Shane Bussmann

2014 August 22

Plot magnification factor vs. S500.  Compare with Schechter function and broken
power-law predictions for the same.

"""

from astropy.table import Table
import matplotlib.pyplot as plt
from pylab import savefig
import matplotlib
import numpy
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from sklearn.linear_model import LogisticRegression


def TrimMu(Stot, oldmu, thresh, norm=True):

    """ Remove mu measurements below the threshold and adjust mu measurements
    above the threshold to unity if we seek P(mu > thresh). """

    trim = oldmu > thresh
    if norm:
        oldmu[:] = 0
        oldmu[trim] = 1
    else:
        oldmu = oldmu[trim]
        Stot = Stot[trim]

    return Stot, oldmu

def model(x):
    return 1 / (1 + numpy.exp(-x))

def LogReg(trimset):

    """ Classify P(mu>mu_min) based on total apparent flux using Logistic
    Regression. """

    Stot = numpy.log10(trimset[0])
    #Stot = numpy.log(trimset[0])
    #npre = 100
    #prevector = numpy.linspace(-100, -1, npre)
    #prevector1 = numpy.zeros(npre)
    #Stot = trimset[0]#numpy.append(prevector, trimset[0])
    newmu = trimset[1]#numpy.append(prevector1, trimset[1])
    clf = LogisticRegression(C=1e5)
    Stotmatrix = Stot.reshape(len(Stot), 1)
    clf.fit(Stotmatrix, newmu)
    #X_test = numpy.linspace(-5, 5, 300)
    #X_test = numpy.linspace(0, 150, 300)
    X_test = numpy.linspace(-2, 3, 300)

    loss = model(X_test * clf.coef_ + clf.intercept_).ravel()
    #import pdb; pdb.set_trace()
    #plt.plot(X_test, loss, color='blue', linewidth=3)
    #plt.show()
    #return X_test, loss
    #return numpy.exp(X_test), loss
    return 10 ** (X_test), loss

def SmoothMu(trimset, smoothing):

    """ Smooth mu or P(mu>thresh). """

    Stot = trimset[0]
    newmu = trimset[1]

    maxflux = Stot.max()
    xfine = numpy.exp(numpy.linspace(0, numpy.log(maxflux), 300))
    #tck = interpolate.splrep(Stot, newmu, s=1, k=5)
    #yfine = interpolate.splev(xfine, tck, der=0)

    yfine = griddata(Stot, newmu, xfine, method='nearest')
    bad = yfine * 0 != 0
    good = yfine * 0 == 0
    replacement = yfine[good][0]
    yfine[bad] = replacement

    if newmu.max() == 1:
        buffmax = 1.
    else:
        buffmax = 10.
    buff0 = numpy.zeros(100) + newmu.min()
    buff1 = numpy.zeros(100) + buffmax
    yfinebig = numpy.append(buff0, yfine)
    yfinebig = numpy.append(yfinebig, buff1)
    yfinebig = numpy.log10(yfinebig)
    smoothedyfine = gaussian_filter(yfinebig, smoothing)
    smoothedyfine = 10 ** (smoothedyfine)

    buff0 = numpy.arange(100) - 100
    buff1 = numpy.arange(100) + xfine.max()
    xfinebig = numpy.append(buff0, xfine)
    xfinebig = numpy.append(xfinebig, buff1)
    #plt.clf()
    #plt.plot(Stot, newmu, 'o')
    #plt.plot(xfine, yfine, '.')
    #plt.plot(xfinebig, smoothedyfine, '.')
    #plt.loglog()

    return xfinebig, smoothedyfine

def MonteCarloMu(Stot, e_Stot, oldmu, e_oldmu, smoothing, norm=True, nsim=1000):

    """ Monte Carlo over the possible mu values from uncertainties. """

    nmu = len(oldmu)
    Stotcube = numpy.zeros([nmu, nsim])
    oldmucube = numpy.zeros([nmu, nsim])
    for imu in range(nmu):
        Stotval = Stot[imu]
        Stotstd = e_Stot[imu]
        muval = oldmu[imu]
        mustd = e_oldmu[imu]
        Stotcube[imu, :] = numpy.random.normal(loc=Stotval, scale=Stotstd,
                size=nsim)
        simmu = numpy.random.normal(loc=muval, scale=mustd, size=nsim)
        oldmucube[imu, :] = simmu

    #print(numpy.min(Stotcube, axis=0))
    Stot_0 = Stotcube[:, 0].copy()
    newmu_0 = oldmucube[:, 0].copy()
    trimset = TrimMu(Stot_0, newmu_0, thresh, norm)
    if norm:
        fineStot, finemu = LogReg(trimset)
    else:
        fineStot, finemu = SmoothMu(trimset, smoothing)
    nfine = len(fineStot)
    fineStotcube = numpy.zeros([nfine, nsim])
    finemucube = numpy.zeros([nfine, nsim])
    for isim in range(nsim):
        Stot_isim = Stotcube[:, isim]
        newmu_isim = oldmucube[:, isim]
        trimset = TrimMu(Stot_isim, newmu_isim, thresh, norm)
        if norm:
            fineStot, finemu = LogReg(trimset)
        else:
            fineStot, finemu = SmoothMu(trimset, smoothing)
        fineStotcube[:, isim] = fineStot
        finemucube[:, isim] = finemu

    bestmu = numpy.median(finemucube, axis=1)
    rmsmu = numpy.std(finemucube, axis=1)

    return fineStot, bestmu, rmsmu

def find_nearest(array, value):
    idx = (numpy.abs(array-value)).argmin()
    return idx

def Pmu_model(modloc, linestyle='solid', lf='Powerlaw', fill=False,
        truncate=0):

    # read in Fialkov data
    fractiondata = Table.read(modloc, format='ascii')

    xx = fractiondata['S870']
    min

    addlabel = True
    #if modloc[8:10] == 'Mu':
    #    addlabel = True
    #elif modloc[8:19] == 'Fractions_1':
    #    addlabel = True
    #else:
    #    addlabel = False

    # predicted P(mu) for nominal best-fit LF
    if lf == 'Schechter':
        key = 'steep'

        yy = fractiondata[key]
        color = 'purple'
        plt.plot(xx, yy, color=color, linestyle=linestyle, linewidth=1.5)
        print(xx[find_nearest(yy, 0.5)])
        if fill:
            pass
            #sigmadata = Table.read('../Data/Sig_1p1_30sh.txt',
            #        format='ascii')
            #y1 = fractiondata['steep'] - sigmadata['steep']
            #y2 = fractiondata['steep'] + sigmadata['steep']
            #plt.fill_between(xx, y1, y2, facecolor='none', zorder=-3,
            #        hatch='/', edgecolor=color)
        if addlabel:
            plt.plot(xx, yy, color=color, linestyle=linestyle, 
                    label='Schechter', linewidth=1.5)
    else:
        key = 'karim'
        ok = fractiondata['S870'] > truncate
        fractiondata = fractiondata[ok]
        xx = fractiondata['S870']
        yy = fractiondata[key]
        color = 'magenta'
        #plt.plot(xx, yy, color=color, linestyle=linestyle, linewidth=1.5)
        #print(key, xx[find_nearest(yy, 0.5)])
        if fill:
            pass
            #sigmadata = Table.read('../Data/Sig_1p1_30powerlaw.txt',
            #        format='ascii')
            #y1 = fractiondata['karim'] - sigmadata['karim']
            #y2 = fractiondata['karim'] + sigmadata['karim']
            #plt.fill_between(xx, y1, y2, facecolor=color, zorder=-4, alpha=0.3)
            #plt.fill_between(xx, y1, y2, facecolor='none', zorder=-4,
            #        edgecolor=color, hatch='\\\\')
        if addlabel:
            plt.plot(xx, yy, color=color, linestyle=linestyle, 
                    label='Karim+ 2013', linewidth=1.5)

        # predicted P(mu) for Steep LF
        key = 'break15'
        ok = fractiondata['S870'] > truncate
        fractiondata = fractiondata[ok]
        xx = fractiondata['S870']
        yy = fractiondata[key]
        color = 'blue'
        plt.plot(xx, yy, color=color, linestyle=linestyle, linewidth=1.5)
        print(key, xx[find_nearest(yy, 0.5)])
        if addlabel:
            plt.plot(xx, yy, color=color, linestyle=linestyle, 
                    label=r'$S_\star = 15\,$mJy', linewidth=1.5)

        if fill:
            pass
            #sigmadata = Table.read('../Data/Sig_1p1_30powerlaw.txt',
            #        format='ascii')
            #y1 = fractiondata['break15'] - sigmadata['break15']
            #y2 = fractiondata['break15'] + sigmadata['break15']
            #plt.fill_between(xx, y1, y2, facecolor=color, zorder=-5)
            #plt.fill_between(xx, y1, y2, facecolor=color, zorder=-5)

        # predicted P(mu) for Flat LF
        #yy = fractiondata['flat']
        #color = 'green'
        #plt.plot(xx, yy, color=color, linestyle=linestyle)

# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

fig = plt.figure(figsize=(5.0, 4.5))

# ALMA plotting parameters
acolor = 'green'
ams = 4
afmt = 's'

# SMA plotting parameters
bcolor = 'red'
bms = 4
bfmt = 'o'

# Fialkov plotting parameters
fcolor = 'purple'

mudat = Table.read('../Data/table_observed.dat', format='ascii')
smadat = Table.read('../Data/bussmann2013_muflux.dat', format='ascii')

# remove poorly measured magnifications
#goodmu = mudat['mu870'] > 2 * mudat['e_mu870']
#goodmu = mudat['e_mu870'] > 0
#mudat = mudat[goodmu]
nomu = mudat['e_remu'] == 0
mudat['e_remu'][nomu] = 0.001
nomu = smadat['e_mu'] == 0
smadat['e_mu'][nomu] = 0.001
#goodmu = smadat['mu'] > 2 * smadat['e_mu']
#goodmu = smadat['e_mu'] > 0
#smadat = smadat[goodmu]

# plot the ALMA and SMA data
mutot = mudat['remu']
mutot1 = mudat['remu']
mutot2 = smadat['mu']
e_mutot1 = mudat['e_remu']
e_mutot2 = smadat['e_mu']
mutot = numpy.append(mutot, smadat['mu'])
e_mutot = mudat['e_remu']
e_mutot = numpy.append(e_mutot, smadat['e_mu'])
Stot = mudat['f870']
Stot = numpy.append(Stot, smadat['fnu'])
e_Stot = mudat['e_f870']
e_Stot = numpy.append(e_Stot, smadat['e_fnu'])
Stot1 = mudat['f870']
e_Stot1 = mudat['e_f870']
Stot2 = smadat['fnu']
e_Stot2 = smadat['e_fnu']
indxsort = numpy.argsort(Stot)
Stot = Stot[indxsort]
e_Stot = e_Stot[indxsort]
mutot = mutot[indxsort]
e_mutot = e_mutot[indxsort]
notlensed = mutot == e_mutot
e_mutot[notlensed] = 0.001

ok = mutot > 2.0
Shigh = Stot[ok]
Smin = Shigh.min()

# tack on mu=1, Stot=0 at beginning
#Stot = numpy.append(-100+numpy.arange(100), Stot)
#Stot = numpy.append(Stot, numpy.arange(100) + Stot.max())
#mutot = numpy.append(numpy.zeros(100)+1.1, mutot)
#mutot = numpy.append(mutot, numpy.zeros(100) + 15)

# Fialkov lens predictions
#mudata = Table.read('../Data/Mu_1_30.txt', format='ascii')
#x = mudata['flux']
#y1 = mudata['bestfit-meanmu'] - mudata['bestfit-std']
#y2 = mudata['bestfit-meanmu'] + mudata['bestfit-std']
##plt.fill_between(x, y1, y2, edgecolor=fcolor, hatch='//', facecolor='none')
#plt.plot(x, mudata['bestfit-meanmu'], color=fcolor, label='Best-fit LF')

#y1 = mudata['flat-meanmu'] - mudata['flat-std']
#y2 = mudata['flat-meanmu'] + mudata['flat-std']
##plt.fill_between(x, y1, y2, facecolor='orange', alpha=0.1)
#plt.plot(x, mudata['flat-meanmu'], color='orange', label='Flat LF')

#y1 = mudata['steep-meanmu'] - mudata['steep-std']
#y2 = mudata['steep-meanmu'] + mudata['steep-std']
#plt.fill_between(x, y1, y2, facecolor='blue', alpha=0.1)
#plt.plot(x, mudata['steep-meanmu'], color='blue', label='Steep LF')

# read in Fialkov data

# Schechter luminosity function
#Pmu_model('../Data/Mu_1p1_30Sc.txt', linestyle='dotted', lf='Schechter')
#Pmu_model('../Data/Mu_1p2_30Sc.txt', linestyle='--', lf='Schechter')
#Pmu_model('../Data/Mu_SH_2_30_z1p5.txt', linestyle='solid', lf='Schechter')

# Broken power-law luminosity function
#Pmu_model('../Data/Mu_1p1_30powerlaw.txt', linestyle='solid', fill=True)
truncateval = Smin
Pmu_model('../Data/Mu_PL_2_30_z1p5.txt', linestyle='solid', fill=True,
        truncate=truncateval)

#Pmu_model('../Data/Mu_1p2_30powerlaw.txt', linestyle='--')
#Pmu_model('../Data/Mu_2_30powerlaw.txt', linestyle='solid')

# Schechter luminosity function
#Pmu_model('../Data/Mu_1p1_30sh.txt', linestyle='solid', lf='Schechter', fill=True)
#Pmu_model('../Data/Mu_SH_2_30_z1p5.txt', linestyle='solid', lf='Schechter', fill=True)
#Pmu_model('../Data/Mu_1p2_30sh.txt', linestyle='--', lf='Schechter')
#Pmu_model('../Data/Mu_2_30sh.txt', linestyle='solid', lf='Schechter')

#Stottmp = Stot.copy()
#mutottmp = mutot.copy()
#print('<mu>')
thresh = 2.0
smoothing = 30.0
Sfull, mubest, mustd = MonteCarloMu(Stot, e_Stot, mutot, e_mutot, smoothing, norm=False)
ok = Sfull > Smin
Sfull = Sfull[ok]
mubest = mubest[ok]
mustd = mustd[ok]
label = r'$\mu_{\rm min} = ' + str(thresh) + '$'
plt.plot(Sfull, mubest, label=label, color='black', linewidth=2.0)
y1 = mubest - mustd
y2 = mubest + mustd
plt.fill_between(Sfull, y1, y2=y2, color='gray', zorder=1)

yyy = mudat['remu']
xxx = mudat['f870']
yyyerr = mudat['e_remu']
xxxerr = mudat['e_f870']

plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=',', ecolor='0.5',
        capsize=0)
plt.plot(xxx, yyy, afmt, color=acolor, ms=ams, label='Herschel-ALMA',
    mec='black')
#plt.axhline(y=1.1)

yyy = smadat['mu']
xxx = smadat['fnu']
yyyerr = smadat['e_mu']
xxxerr = smadat['e_fnu']

plt.errorbar(xxx, yyy, fmt=',', yerr=yyyerr, xerr=xxxerr, ecolor='0.5',
        capsize=0)
plt.plot(xxx, yyy, bfmt, color=bcolor, ms=bms, label='Herschel-SMA', mec='black')

plt.loglog()

xmin = 0.8
xmax = 1.2e2
ymin = 0.9
ymax = 35
plt.axis([xmin, xmax, ymin, ymax])

plt.xlabel(r'$S_{\rm 870-observed}\,({\rm mJy})$', fontsize='x-large')
plt.ylabel(r'$\mu_{870}$', fontsize='x-large')
plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

plt.legend(loc='upper left', numpoints=1, handletextpad=0.35, borderpad=0.6,
        handlelength=1.5)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='medium')
plt.subplots_adjust(left=0.13, right=0.97, top=0.97, bottom=0.14, wspace=0.39)

savefig('../Figures/s870_mu_z1p5.pdf')

# ======================================
# P(mu) vs. S_870

plt.clf()

# read in Fialkov data

# Schechter luminosity function
#Pmu_model('../Data/Fractions_SH_1p1_30_z1p5.txt', linestyle='solid', lf='Schechter')
#Pmu_model('../Data/Fractions_1p2_30Sc.txt', linestyle='--', lf='Schechter')
#Pmu_model('../Data/Fractions_SH_2_30_z1p5.txt', linestyle='solid', lf='Schechter')

# Broken power-law luminosity function
#Pmu_model('../Data/Fractions_PL_1p1_30_z1p5.txt', linestyle='solid')
#Pmu_model('../Data/Fractions_1p2_30powerlaw.txt', linestyle='--')
Pmu_model('../Data/Fractions_PL_2_30_z1p5.txt', linestyle='solid')

#Pmu(Stot, mutot, 1.0, linestyle='-.')
#thresh = 1.1
#Sfull, mubest, mustd = MonteCarloMu(Stot, e_Stot, mutot, e_mutot, smoothing, norm=True)
#label = r'$\mu_{\rm min} = ' + str(thresh) + '$'
#plt.plot(Sfull, mubest, label=label, color='black', linewidth=2.0)
#y1 = mubest - mustd
#y2 = mubest + mustd
#plt.fill_between(Sfull, y1, y2=y2, color='gray', zorder=1)
thresh = 2.0
Sfull, mubest, mustd = MonteCarloMu(Stot, e_Stot, mutot, e_mutot, smoothing, norm=True)
label = 'Logistic Fit to Data'
plt.plot(Sfull, mubest, label=label, color='black', linewidth=2.0)
y1 = mubest - mustd
y2 = mubest + mustd
plt.fill_between(Sfull, y1, y2=y2, color='gray', zorder=1)
print(Sfull[find_nearest(mubest, 0.5)])
#thresh = 1.1
#Stottmp = Stot.copy()
#mutottmp = mutot.copy()
#smoother = 8.0
#Pmu(Stottmp, mutottmp, thresh, smoother, linestyle='solid', norm=True)
#thresh = 1.2
#Stottmp = Stot.copy()
#mutottmp = mutot.copy()
#Pmu(Stot, mutot, thresh, smoother, linestyle='dashed', norm=True)
#thresh = 2.0
#Stottmp = Stot.copy()
#mutottmp = mutot.copy()
#Pmu(Stottmp, mutottmp, thresh, smoother, linestyle='--', norm=True)

classified = TrimMu(Stot1, mutot1, thresh, norm = True)
xxx = classified[0]
yyy = classified[1]
xxxerr = e_Stot1
#plt.errorbar(xxx, yyy, xerr=xxxerr, fmt=',', ecolor='0.5',
#        capsize=0)
plt.plot(xxx, yyy, afmt, color=acolor, ms=ams, label='Herschel-ALMA',
    mec='black')

classified = TrimMu(Stot2, mutot2, thresh, norm = True)
xxx = classified[0]
yyy = classified[1]
#plt.errorbar(xxx, yyy, fmt=',', yerr=yyyerr, xerr=xxxerr, ecolor='0.5',
#        capsize=0)
plt.plot(xxx, yyy, bfmt, color=bcolor, ms=bms, label='Herschel-SMA', mec='black')

plt.semilogx()

plt.axis([1,150,-0.1,1.1])
plt.legend(loc=(0.05,0.5), numpoints=1, handletextpad=0.35, borderpad=0.4,
        handlelength=1.0, labelspacing=0.18)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='medium')

plt.xlabel(r'$S_{\rm 870-observed}\,({\rm mJy})$', fontsize='x-large')
plt.ylabel(r'$P\,(\mu_{870} > 2)$', fontsize='x-large')
plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

savefig('../Figures/Pmu_S870_z1p5.pdf')

import pdb; pdb.set_trace()
