"""
Created 2012 September 3
Shane Bussmann

Purpose: Plot SigmaFIR vs. LFIR
"""

# import modules and define plot appearance
from astropy.table import Table
import matplotlib
import numpy
import matplotlib.pyplot as plt
from pylab import savefig
#from scipy.interpolate import griddata
#from scipy import interpolate
from scipy.interpolate import Rbf
from astropy import cosmology


def interpolateLFIR(S870, z, S870_ref, z_ref, LFIR_ref):
    #points, antipoints = numpy.meshgrid(z_ref, S870_ref)
    #values = LFIR_ref
    #xi = (z, S870)
    #import pdb; pdb.set_trace()
    #LFIR = griddata(points, values, xi)
    #tck = interpolate.bisplrep(S870_ref, z_ref, LFIR_ref, s=0)
    #LFIR = interpolate.bisplev(S870, z, tck)
    rbf = Rbf(S870_ref, z_ref, LFIR_ref, epsilon=2, function='linear')
    LFIR = rbf(S870, z)
    return LFIR

def simLFIR(S870dist, zdist, nsim, S870_ref, z_ref, LFIR_ref):
    LFIRdist = numpy.zeros(nsim)
    index = numpy.random.uniform(size=nsim)
    intindex = (index * len(S870dist)).astype(int)
    S870rand = S870dist[intindex]
    intindex = (index * len(zdist)).astype(int)
    zrand = zdist[intindex]
    for isim in range(nsim):
        iz = zrand[isim]
        iS870 = S870rand[isim]
        LFIRdist[isim] = interpolateLFIR(iS870, iz, S870_ref, z_ref, LFIR_ref)
    return LFIRdist

def sampleLFIR(sample_S870, sample_e_S870, sample_z, sample_e_z, S870_ref,
        z_ref, LFIR_ref):
    nsample = sample_S870.size
    sample_LFIR = numpy.zeros(nsample)
    sample_e_LFIR = numpy.zeros(nsample)
    for isample in range(nsample):
        loc = sample_S870[isample]
        scale = sample_e_S870[isample]
        S870 = numpy.random.normal(loc=loc, scale=scale, size=1000)
        loc = sample_z[isample]
        scale = sample_e_z[isample]
        z = numpy.random.normal(loc=loc, scale=scale, size=1000)
        LFIRdist = simLFIR(S870, z, 1000, S870_ref, z_ref, LFIR_ref)
        sample_LFIR[isample] = LFIRdist.mean()
        sample_e_LFIR[isample] = LFIRdist.std()
    return sample_LFIR, sample_e_LFIR


def physicalDistance(reff, z):
    angdist = cosmology.angular_diameter_distance(z)
    Dist = reff / 206265 * angdist.value * 1e3
    return Dist

def simDist(rdist, zdist, nsim):
    distdist = numpy.zeros(nsim)
    index = numpy.random.uniform(size=nsim)
    intindex = (index * len(rdist)).astype(int)
    rrand = rdist[intindex]
    intindex = (index * len(zdist)).astype(int)
    zrand = zdist[intindex]
    for isim in range(nsim):
        ir = rrand[isim]
        iz = zrand[isim]
        distdist[isim] = physicalDistance(ir, iz)
    return distdist

def sampleDist(sample_r, sample_e_r, sample_z, sample_e_z):
    nsample = sample_r.size
    sample_dist = numpy.zeros(nsample)
    sample_e_dist = numpy.zeros(nsample)
    for isample in range(nsample):
        loc = sample_r[isample]
        scale = sample_e_r[isample]
        r = numpy.random.normal(loc=loc, scale=scale, size=1000)
        loc = sample_z[isample]
        scale = sample_e_z[isample]
        z = numpy.random.normal(loc=loc, scale=scale, size=1000)
        distdist = simDist(r, z, 1000)
        sample_dist[isample] = distdist.mean()
        sample_e_dist[isample] = distdist.std()
    return sample_dist, sample_e_dist

def sampleSFIR(sample_LFIR, sample_e_LFIR, sample_Dist, sample_e_Dist):
    nsample = sample_LFIR.size
    sample_SFIR = numpy.zeros(nsample)
    sample_e_SFIR = numpy.zeros(nsample)
    for isample in range(nsample):
        loc = sample_LFIR[isample]
        scale = sample_e_LFIR[isample]
        iLFIR = numpy.random.normal(loc=loc, scale=scale, size=1000)
        loc = sample_Dist[isample]
        scale = sample_e_Dist[isample]
        iDist = numpy.random.normal(loc=loc, scale=scale, size=1000)
        iSFIR = 0.5 * iLFIR / numpy.pi / iDist ** 2
        sample_SFIR[isample] = iSFIR.mean()
        sample_e_SFIR[isample] = iSFIR.std()

    return sample_SFIR, sample_e_SFIR

font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

# read in the plotting color scheme
colorschemeloc = 'colorscheme.txt'
colorscheme = Table.read(colorschemeloc, format='ascii')
ecolor = str(colorscheme['ecolor'][0])
color3 = colorscheme['color3'][0]
color4 = colorscheme['color4'][0]
gray = colorscheme['color6'][0]

# set ALMA plotting parameters
acolor = 'green'
ams = 6
afmt = 's'

# set SPT plotting parameters
hcolor = colorscheme['color2'][0]
hms = 7
hfmt = '*'

# set SMA plotting parameters
bcolor = colorscheme['color1'][0]
bms = 6
bfmt = 'o'

# get reference sample luminosities and S870 fluxes and redshifts
referencedat = Table.read('../Data/magnelli_fluxes.dat', format='ascii')

reflfirdat = Table.read('../../greybody03/greybody_params.dat', format='ascii')
nref = len(referencedat)
flagger = []
for iref in range(nref):
    shortname = referencedat['shortname'][iref]
    match = reflfirdat['shortname'] == shortname
    if len(reflfirdat[match]) == 0:
        flagger.append(1)
    else:
        flagger.append(0)

ok = numpy.array(flagger) == 0
referencedat = referencedat[ok]
S870_ref = referencedat['f850']
z_ref = referencedat['z_source']
LFIR_ref = 10 ** reflfirdat['bestlfir'] / 1e12

# get ALMA size and luminosity measurements
almatarget = Table.read('../Data/targetlist.dat', format='ascii')
almadat = Table.read('../Data/source_intrinsic.dat', format='ascii')
alma_S870 = almadat['f870']
alma_e_S870 = almadat['e_f870']
nsource = alma_S870.size
alma_z = numpy.zeros(nsource)
alma_e_z = numpy.zeros(nsource)
for isource in range(nsource):
    target = almadat['target'][isource]
    index = almatarget['dataname'] == target
    alma_z[isource] = almatarget['z_source'][index][0]
    alma_e_z[isource] = almatarget['e_z_source'][index][0]
alma_LFIR, alma_e_LFIR = sampleLFIR(alma_S870, alma_e_S870, alma_z, alma_e_z,
        S870_ref, z_ref, LFIR_ref)
alma_LFIR = alma_LFIR * 1e12

alma_r = almadat['reff']
alma_e_r = almadat['e_reff']
alma_Dist, alma_e_Dist = sampleDist(alma_r, alma_e_r, alma_z, alma_e_z)

alma_SFIR, alma_e_SFIR = sampleSFIR(alma_LFIR, alma_e_LFIR, alma_Dist,
        alma_e_Dist)

#plt.errorbar(alma_Dist, alma_LFIR, yerr=alma_e_LFIR, xerr=alma_e_Dist,
#        fmt=afmt, ecolor='gray', capsize=0)
#savefig('testsigfir.pdf')
#import pdb; pdb.set_trace()

# get redshift and flux data -- SMA
sma_greybody_loc = '../../greybody04'
sma_spirephotloc = sma_greybody_loc + '/bussmann_fluxes.dat'
sma_spirephot = Table.read(sma_spirephotloc, format='ascii')

# get mu data -- SMA
sma_lensproploc = sma_greybody_loc + '/bussmann13_mu.dat'
sma_lensprop = Table.read(sma_lensproploc, format='ascii')

# get IR luminosities -- SMA
sma_firdatloc = sma_greybody_loc + '/greybody_params.dat'
sma_firdat = Table.read(sma_firdatloc, format='ascii')

# get IR luminosity surface densities -- SMA
sma_sigfirloc = sma_greybody_loc + '/sigfir.dat'
sma_sigfir = Table.read(sma_sigfirloc, format='ascii')
sma_nsigfir = len(sma_sigfir)
sma_redshifts = numpy.zeros(sma_nsigfir)
sma_lfir = numpy.zeros(sma_nsigfir)
sma_e_lfir = numpy.zeros(sma_nsigfir)
sma_mu = numpy.zeros(sma_nsigfir)
e_sma_mu = numpy.zeros(sma_nsigfir)
for i in numpy.arange(sma_nsigfir):
    iauname = sma_sigfir['iauname'][i]
    indx = sma_spirephot['iauname'] == iauname
    shortname = sma_spirephot['shortname'][indx][0]
    sma_redshifts[i] = sma_spirephot['z_source'][indx][0]
    indx1 = sma_lensprop['iauname'] == iauname
    sma_mu[i] = numpy.mean(sma_lensprop['mu'][indx1])
    e_sma_mu[i] = numpy.mean(sma_lensprop['e_mu'][indx1])
    indx2 = sma_firdat['iauname'] == iauname
    sma_lfir[i] = sma_firdat['bestlfir'][indx2][0]
    sma_e_lfir[i] = sma_firdat['rmslfir'][indx2][0]

# get redshift and flux data -- SPT
spt_greybody_loc = '../../greybody05'
spt_spirephotloc = spt_greybody_loc + '/hezaveh_fluxes.dat'
spt_spirephot = Table.read(spt_spirephotloc, format='ascii')

# get mu data -- SPT
spt_muloc = spt_greybody_loc + '/hezaveh_mu.dat'
spt_mu = Table.read(spt_muloc, format='ascii')

# get IR luminosity surface densities -- SPT
spt_sigfirloc = spt_greybody_loc + '/sigfir.dat'
spt_sigfir = Table.read(spt_sigfirloc, format='ascii')
spt_nsigfir = len(spt_sigfir)
spt_redshifts = numpy.zeros(spt_nsigfir)
for i in numpy.arange(spt_nsigfir):
    iauname = spt_sigfir['iauname'][i]
    indx = spt_spirephot['iauname'] == iauname
    spt_redshifts[i] = spt_spirephot['z_source'][indx][0]

# get IR luminosities -- SPT
spt_lfirloc = spt_greybody_loc + '/greybody_params.dat'
spt_lfir = Table.read(spt_lfirloc, format='ascii')

nxplot = 1
nyplot = 1
xwidth = 05.
ywidth = 04.
fig = plt.figure(figsize=(xwidth, ywidth))
axB = fig.add_subplot(nyplot, nxplot, 1)
axL = fig.add_subplot(nyplot, nxplot, 1)
plt.subplots_adjust(left=0.19, right=0.84, top=0.95, bottom=0.14, hspace=0.2,
        wspace=0.23)

# convert LFIR to SFR assuming Salpeter IMF
lfirtosfr = 1.91 * 170 / 1e12

# plot median value and range from Tacconi et al. 2006
tacconi_rhalf = numpy.array([2 - 0.3, 2 + 0.3])
tacconi_lfir = numpy.array([900 - 400, 900 + 400]) / lfirtosfr
tacconi_siglfir = numpy.array([80 -20, 80 + 20]) / lfirtosfr
tacconi_z = numpy.array([2.4 - 0.25, 2.4 + 0.25])
y1 = numpy.array([tacconi_siglfir[0], tacconi_siglfir[0]])
y2 = numpy.array([tacconi_siglfir[1], tacconi_siglfir[1]])
plt.fill_between(tacconi_z, y1, y2, facecolor=color3, edgecolor='none', \
        alpha=0.3)

# plot typical values from Younger et al. studies
#sfr_younger = numpy.array([1300., 6000.])
#rhalf_younger = numpy.array([1.0, 4.])
#sigmasfr_younger = 0.5 * sfr_younger / math.pi / rhalf_younger**2
#y1 = [sigmasfr_younger[0], sigmasfr_younger[0]]
#y2 = [sigmasfr_younger[1], sigmasfr_younger[1]]
#plt.fill_between(sfr_younger, y1, y2=y2, color='red', alpha=0.5)
    
# plot SigmaSFR vs. SFR

# plot theoretical limit for Sigma_SFR
plt.axhline(y=1e3 / lfirtosfr, linestyle='--', color=color4, zorder=0)

# plot Sigma_LFIR vs. z -- SPT
yyy = 10 ** spt_sigfir['bestsigfir']
yyyerrplus = 10**(spt_sigfir['bestsigfir'] + spt_sigfir['rmssigfir'])
yyyerrminus = 10**(spt_sigfir['bestsigfir'] - spt_sigfir['rmssigfir'])
yyyerr = numpy.abs(yyyerrplus - yyyerrminus) / 2.
xxx = spt_redshifts
xxxerr = spt_redshifts / 1e4
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=hfmt, color=hcolor, \
        ecolor=ecolor, capsize=0, ms=8)
plt.plot(xxx, yyy, hfmt, ms=8, color=hcolor, label='SPT', zorder=5)

# plot Sigma_LFIR vs. z -- SMA
yyy = 10 ** sma_sigfir['bestsigfir']
yyyerrplus = 10**(sma_sigfir['bestsigfir'] + sma_sigfir['rmssigfir'])
yyyerrminus = 10**(sma_sigfir['bestsigfir'] - sma_sigfir['rmssigfir'])
yyyerr = numpy.abs(yyyerrplus - yyyerrminus) / 2.
xxx = sma_redshifts
xxxerr = sma_redshifts / 1e4
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=bfmt, color='red', \
        ecolor=ecolor, capsize=0, ms=6)
plt.plot(xxx, yyy, bfmt, ms=6, color='red', label='SMA', zorder=10)

# plot Sigma_LFIR vs. z -- ALMA
ok = alma_SFIR > alma_e_SFIR
yyy = alma_SFIR[ok]
yyyerr = alma_e_SFIR[ok]
xxx = alma_z[ok]
xxxerr = alma_e_z[ok]
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=afmt, \
        ecolor=ecolor, capsize=0, ms=ams, zorder=10)
plt.plot(xxx, yyy, afmt, ms=ams, color=acolor, label='ALMA', zorder=11)

xmin = 1
xmax = 6
ymin = 2e10
ymax = 5e13
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel(r'$z_{\rm source}$', fontdict={'fontsize':'x-large'})
plt.ylabel(r'$\Sigma_{\rm FIR} \, ({\rm L}_\odot \, {\rm kpc}^{-2})$', 
        fontdict={'fontsize':'x-large'})
plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
    #plt.text(xmin*1.1, ymin*1.5, spirephot['shortname'][i], fontsize=12)

#plt.xticks(range(0.01, 10, 10))

#plt.loglog()
plt.semilogy()

plt.legend(loc='lower right', numpoints=1, handletextpad=0.00, borderpad=0.3,
        labelspacing=0.2, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small') 

axR = axB.twinx()
axR.plot(xxx, yyy * lfirtosfr, '.', mfc='none', zorder=-1)
axR.set_ylabel(r'$\Sigma_{\rm SFR} \, ({\rm M}_\odot \, {\rm yr}^{-1} \, ' + \
        r'{\rm kpc}^{-2}$)', fontdict={'fontsize':'x-large'})
#axR.axhline(y=1e3, linestyle='--', color=color4, zorder=-10)

axR.minorticks_on()
axR.tick_params(width=1.5, which='both')
axR.tick_params(length=2, which='minor')
axR.tick_params(length=4, which='major')
axR.semilogy()
axR.axis(numpy.array([xmin, xmax, ymin * lfirtosfr, ymax * lfirtosfr]))

saveloc = '../Figures/redshift_sigfir.pdf'
savefig(saveloc)

import pdb; pdb.set_trace()

plt.clf()

plt.subplots_adjust(left=0.19, right=0.84, top=0.85, bottom=0.18, hspace=0.2,
        wspace=0.23)
axB = fig.add_subplot(nyplot, nxplot, 1)

y1 = numpy.array([tacconi_siglfir[0], tacconi_siglfir[0]])
y2 = numpy.array([tacconi_siglfir[1], tacconi_siglfir[1]])
plt.fill_between(tacconi_lfir, y1, y2, facecolor=color3, edgecolor='none', \
        alpha=0.3)

# plot theoretical limit for Sigma_SFR
plt.axhline(y=1e3 / lfirtosfr, linestyle='--', color=color4, zorder=0)

# plot Sigma_LFIR vs. z -- SPT
yyy = 10 ** spt_sigfir['bestsigfir']
yyyerrplus = 10**(spt_sigfir['bestsigfir'] + spt_sigfir['rmssigfir'])
yyyerrminus = 10**(spt_sigfir['bestsigfir'] - spt_sigfir['rmssigfir'])
yyyerr = numpy.abs(yyyerrplus - yyyerrminus) / 2.
xxx = 10 ** spt_lfir['bestlfir']
xxxerrplus = 10**(spt_lfir['bestlfir'] + spt_lfir['rmslfir'])
xxxerrminus = 10**(spt_lfir['bestlfir'] - spt_lfir['rmslfir'])
xxxerr = numpy.abs(xxxerrplus - xxxerrminus) / 2.
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=hfmt, color=hcolor, \
        ecolor=ecolor, capsize=0, ms=hms)
plt.plot(xxx, yyy, hfmt, ms=hms, color=hcolor, label='SPT', \
        zorder=5)

# plot Sigma_LFIR vs. z -- SMA
yyy = 10 ** sma_sigfir['bestsigfir']
yyyerrplus = 10**(sma_sigfir['bestsigfir'] + sma_sigfir['rmssigfir'])
yyyerrminus = 10**(sma_sigfir['bestsigfir'] - sma_sigfir['rmssigfir'])
yyyerr = numpy.abs(yyyerrplus - yyyerrminus) / 2.
xxx = 10 ** sma_lfir
xxxerrplus = 10**(sma_lfir + sma_e_lfir)
xxxerrminus = 10**(sma_lfir - sma_e_lfir)
xxxerr = numpy.abs(xxxerrplus - xxxerrminus) / 2.
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=bfmt, color=bcolor, \
        ecolor=ecolor, capsize=0, ms=bms)
plt.plot(xxx, yyy, bfmt, ms=bms, color=bcolor, label='Herschel', \
        zorder=10)

plt.loglog()
plt.axis([xmin, xmax, ymin, ymax])

xmin = 5e10
xmax = 1e14
ymin = 1e9
ymax = 5e13
plt.xlabel(r'$L_{\rm FIR} \, ({\rm L}_\odot)$', fontdict={'fontsize':'x-large'})
plt.ylabel(r'$\Sigma_{\rm FIR} \, ({\rm L}_\odot \, {\rm kpc}^{-2})$', 
        fontdict={'fontsize':'x-large'})

plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
plt.loglog()
plt.axis(numpy.array([xmin, xmax, ymin, ymax]))

#plt.legend(loc='lower right', numpoints=1, handletextpad=0.00, borderpad=0.3,
#        labelspacing=0.15, handlelength=1.0)
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='small') 

axT = axB.twiny()
axT.plot(xxx * lfirtosfr, yyy, bfmt, ms=bms, mfc='none', zorder=0)
axT.set_xlabel(r'SFR (${\rm M}_\odot \, {\rm yr}^{-1}$)' ,
        fontdict={'fontsize':'x-large'})

axT.loglog()

#plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
#plt.loglog()
plt.axis(numpy.array([xmin * lfirtosfr, xmax * lfirtosfr, ymin, ymax]))

axR = axB.twinx()
axR.plot(xxx, yyy * lfirtosfr, bfmt, ms=bms, mfc='none', zorder=0)
axR.set_ylabel(r'$\Sigma_{\rm SFR} \, ({\rm M}_\odot \, {\rm yr}^{-1} \, $' +\
        r'${\rm kpc}^{-2}$)', fontdict={'fontsize':'x-large'})
#axT.xaxis.tick_top()
#axT.xaxis.set_label_position("top")
#plt.plot(xxx * lfirtosfr, yyy * lfirtosfr, '.', mec='none', mfc='none')
#plt.xlabel(r'SFR (${\rm M}_\odot \, {\rm yr}^{-1}$)')

#axR = fig.add_subplot(1, 1, 1, sharex=axL, frameon=False)
#axR.yaxis.tick_right()
#axR.yaxis.set_label_position("right")
#plt.plot(xxx * lfirtosfr, yyy * lfirtosfr, '.', mec='none', mfc='none')
#plt.xlabel(r'SFR (${\rm M}_\odot \, {\rm yr}^{-1}$)')
    #plt.text(xmin*1.1, ymin*1.5, spirephot['shortname'][i], fontsize=12)

plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
plt.loglog()
plt.axis(numpy.array([xmin, xmax, ymin * lfirtosfr, ymax * lfirtosfr]))

saveloc = 'lfir_sigfir.pdf'
#saveloc = 'submit/f8b.eps'
savefig(saveloc)

plt.clf()
axB = fig.add_subplot(nyplot, nxplot, 1)

plt.subplots_adjust(left=0.16, right=0.84, top=0.85, bottom=0.18, hspace=0.2,
        wspace=0.23)

y1 = numpy.array([tacconi_rhalf[0], tacconi_rhalf[0]])
y2 = numpy.array([tacconi_rhalf[1], tacconi_rhalf[1]])
plt.fill_between(tacconi_lfir, y1, y2, facecolor=color3, edgecolor='none', \
        alpha=0.3)

# plot theoretical limit for Sigma_SFR
#plt.axhline(y=1e3 / lfirtosfr, linestyle='--', color=color4, zorder=0)

# plot Rhalf vs. LFIR -- SPT
yyy = spt_sigfir['bestrhalf']
yyyerr = spt_sigfir['rmsrhalf']
xxx = 10 ** spt_lfir['bestlfir']
xxxerrplus = 10**(spt_lfir['bestlfir'] + spt_lfir['rmslfir'])
xxxerrminus = 10**(spt_lfir['bestlfir'] - spt_lfir['rmslfir'])
xxxerr = numpy.abs(xxxerrplus - xxxerrminus) / 2.
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=hfmt, color=hcolor, \
        ecolor=ecolor, capsize=0, ms=hms)
plt.plot(xxx, yyy, hfmt, ms=hms, color=hcolor, label='SPT', \
        zorder=5)

# plot rhalf vs. LFIR -- SMA
yyy = sma_sigfir['bestrhalf']
yyyerr = sma_sigfir['rmsrhalf']
xxx = 10 ** sma_lfir
xxxerrplus = 10**(sma_lfir + sma_e_lfir)
xxxerrminus = 10**(sma_lfir - sma_e_lfir)
xxxerr = numpy.abs(xxxerrplus - xxxerrminus) / 2.
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=bfmt, color=bcolor, \
        ecolor=ecolor, capsize=0, ms=bms, zorder=0, mec='none')
plt.plot(xxx, yyy, bfmt, ms=bms, color=bcolor, label='Herschel', \
        zorder=10)

plt.legend(loc='lower right', numpoints=1, handletextpad=0.00, borderpad=0.3,
        labelspacing=0.15, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='small') 

plt.loglog()
plt.axis([xmin, xmax, ymin, ymax])

xmin = 5e10
xmax = 1e14
ymin = 0.01
ymax = 10
plt.xlabel(r'$L_{\rm FIR} \, ({\rm L}_\odot)$', fontdict={'fontsize':'x-large'})
plt.ylabel(r'$r_{\rm half} \,{\rm (kpc)}$', fontdict={'fontsize':'x-large'})

plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
plt.loglog()
plt.axis(numpy.array([xmin, xmax, ymin, ymax]))

#plt.legend(loc='lower right', numpoints=1, handletextpad=0.00, borderpad=0.3,
#        labelspacing=0.15, handlelength=1.0)
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='small') 

axT = axB.twiny()
axT.plot(xxx * lfirtosfr, yyy, bfmt, ms=bms, mec='none', zorder=0,
        color=bcolor)
axT.set_xlabel(r'SFR (${\rm M}_\odot \, {\rm yr}^{-1}$)' ,
        fontdict={'fontsize':'x-large'})

axT.semilogx()

#plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
#plt.loglog()
plt.axis(numpy.array([xmin * lfirtosfr, xmax * lfirtosfr, ymin, ymax]))

saveloc = 'lfir_rhalf.pdf'
#saveloc = 'submit/f8a.eps'
savefig(saveloc)

plt.clf()

plt.subplots_adjust(left=0.19, right=0.84, top=0.85, bottom=0.18, hspace=0.2,
        wspace=0.23)
axB = fig.add_subplot(nyplot, nxplot, 1)

# plot theoretical limit for Sigma_SFR
plt.axhline(y=1e3 / lfirtosfr, linestyle='--', color=color4, zorder=0)

# plot Sigma_LFIR vs. mu -- SPT
yyy = 10 ** spt_sigfir['bestsigfir']
yyyerrplus = 10**(spt_sigfir['bestsigfir'] + spt_sigfir['rmssigfir'])
yyyerrminus = 10**(spt_sigfir['bestsigfir'] - spt_sigfir['rmssigfir'])
yyyerr = numpy.abs(yyyerrplus - yyyerrminus) / 2.
xxx = spt_mu['mu']
xxxerr = spt_mu['e_mu']
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=hfmt, color=hcolor, \
        ecolor=ecolor, capsize=0, ms=hms)
plt.plot(xxx, yyy, hfmt, ms=hms, color=hcolor, label='SPT', \
        zorder=5)

# plot Sigma_LFIR vs. mu -- SMA
yyy = 10 ** sma_sigfir['bestsigfir']
yyyerrplus = 10**(sma_sigfir['bestsigfir'] + sma_sigfir['rmssigfir'])
yyyerrminus = 10**(sma_sigfir['bestsigfir'] - sma_sigfir['rmssigfir'])
yyyerr = numpy.abs(yyyerrplus - yyyerrminus) / 2.
xxx = sma_mu
xxxerr = e_sma_mu
plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=bfmt, color=bcolor, \
        ecolor=ecolor, capsize=0, ms=bms)
plt.plot(xxx, yyy, bfmt, ms=bms, color=bcolor, label='Herschel', \
        zorder=10)

plt.semilogy()
plt.axis([xmin, xmax, ymin, ymax])

xmin = 0.
xmax = 22.
ymin = 2e10
ymax = 1e13
plt.xlabel(r'$\mu_{\rm 880 \, \mu m}$', fontdict={'fontsize':'x-large'})
plt.ylabel(r'$\Sigma_{\rm FIR} \, ({\rm L}_\odot \, {\rm kpc}^{-2})$', 
        fontdict={'fontsize':'x-large'})

plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
plt.semilogy()
plt.axis(numpy.array([xmin, xmax, ymin, ymax]))

#plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
#plt.loglog()
plt.axis(numpy.array([xmin * lfirtosfr, xmax * lfirtosfr, ymin, ymax]))

axR = axB.twinx()
axR.plot(xxx, yyy * lfirtosfr, bfmt, ms=bms, mec='none', mfc='none')
axR.set_ylabel(r'$\Sigma_{\rm SFR} \, ({\rm M}_\odot \, {\rm yr}^{-1} \, $' +\
        r'${\rm kpc}^{-2}$)', fontdict={'fontsize':'x-large'})
#axT.xaxis.tick_top()
#axT.xaxis.set_label_position("top")
#plt.plot(xxx * lfirtosfr, yyy * lfirtosfr, '.', mec='none', mfc='none')
#plt.xlabel(r'SFR (${\rm M}_\odot \, {\rm yr}^{-1}$)')

#axR = fig.add_subplot(1, 1, 1, sharex=axL, frameon=False)
#axR.yaxis.tick_right()
#axR.yaxis.set_label_position("right")
#plt.plot(xxx * lfirtosfr, yyy * lfirtosfr, '.', mec='none', mfc='none')
#plt.xlabel(r'SFR (${\rm M}_\odot \, {\rm yr}^{-1}$)')
    #plt.text(xmin*1.1, ymin*1.5, spirephot['shortname'][i], fontsize=12)

plt.minorticks_on()
plt.tick_params(width=1.5, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')
plt.semilogy()
plt.axis(numpy.array([xmin, xmax, ymin * lfirtosfr, ymax * lfirtosfr]))

saveloc = 'mu_sigfir.pdf'
savefig(saveloc)

