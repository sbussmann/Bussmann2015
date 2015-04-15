"""

Shane Bussmann

2014 August 20

Plot dN/dA as a function of angular separation from the center of light.  dN =
number of objects between radius 1 and radius 2.  dA = area between radius 1
and radius 2.

"""

from astropy.table import Table
from astropy.io import ascii
import matplotlib
import matplotlib.pyplot as plt
from pylab import savefig
import numpy
import sep_util


# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.2

fig = plt.figure(figsize=(5.0, 4.5))

nbins = 10
binwidth = 1.0
bin_edges = numpy.arange(0, nbins + 1, binwidth)

# ***********
# ALMA sample with ALESS uv coverage and 1.4 mJy flux limit
# ***********

# NOTE: scaling factor by which to reduce flux of each ALMA source in my
# sample should be determined by comparing total 870um flux densities for each
# Herschel source and each LESS source.  So far, the factor of 3 was meant to
# reproduce the median S500 flux density ratio, but S870 will be more accurate.

simcolor = 'MediumSpringGreen'
simms = 4
simfmt = 's-'

fluxcomponent_file = '../Data/table_ALESSsim.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='target')

# filter out objects below a given flux threshold
fluxcomponent = sep_util.setThresh(fluxcomponent, 1.4, fluxstring='peakflux')
nmultiples = len(fluxcomponent)

simalma = sep_util.getSeparation(fluxcomponent, rastring='ra_alma', \
        decstring='dec_alma', fluxstring='peakflux')
avgsep_simalma, wmeansep_simalma, ra_simalma, dec_simalma = simalma

sep_util.histArea(avgsep_simalma, nbins, color=simcolor, fmt=simfmt, ms=simms,
        norm=nmultiples, label='Herschel-ALMA ALESS-Sim')

#sep_util.simArea(fluxcomponent, nsim, bin_edges, fluxstring='S_870_observed',
#        edgecolor=asimcolor, facecolor='none', hatch='\\', norm=nmultiples)

# *********
# Cowley simulation
# *********

c15color = 'Orange'
c15ms = 6
c15fmt = 'x-'
c15mew = 1.5

c15 = ascii.read('../Data/SPIRE_ALMA_Cat_v4.txt')
s500_c15 = c15['SourceS500']
zc = c15['z']
hithresh = (s500_c15 > 50) & (zc > 1)
c15 = c15[hithresh]
c15 = sep_util.rmSingles(c15, targetstring='SurveyID')
nmultiples = len(c15)

simc15 = sep_util.getSeparation(c15, degrees=True, rastring='GalaxyX', \
        decstring='GalaxyY', fluxstring='GalaxyS850', targetstring='SurveyID')
avgsep_c15, wmeansep_c15, ra_c15, dec_c15 = simc15
sep_util.histArea(avgsep_c15, nbins, color=c15color, fmt=c15fmt,
        ms=c15ms, norm=nmultiples, showerror=False, label='C15 Simulation',
        mew=c15mew)

# *********
# ALESS
# *********

# Plotting parameters
hodgecolor = 'LightPink'
hodgesimcolor = 'LightPink'
hodgems = 4
hodgefmt = 'D-'

# Load the data
fluxcomponent_file = '../Data/hodge2013.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='lessid')
nmultiples = len(fluxcomponent)

hodge = sep_util.getSeparation(fluxcomponent, rastring='ra_alma', \
        decstring = 'dec_alma', targetstring='lessid')
avgsep_hodge, wmeansep_hodge, ra_hodge, dec_hodge = hodge

deltasep = avgsep_hodge.max() - avgsep_hodge.min()
#nbins = deltasep / binwidth
sep_util.histArea(avgsep_hodge, nbins, color=hodgecolor, fmt=hodgefmt,
        ms=hodgems, norm=nmultiples, label='ALESS')

indexsort = numpy.argsort(avgsep_hodge)
avgsep_hodge = avgsep_hodge[indexsort]
flux_hodge = fluxcomponent['f880'][indexsort]
nflux = flux_hodge.size
sumflux_hodge = numpy.zeros(nflux)
for i in range(nflux):
    sumflux_hodge[i] = flux_hodge[0:i].sum()
#plt.plot(avgsep_hodge, sumflux_hodge)

# plot simulated positions
nsim = 1000
#sep_util.simArea(fluxcomponent, nsim, bin_edges, targetstring='lessid',
#        edgecolor=hodgesimcolor, facecolor='none', hatch='//', norm=nmultiples)

# ***********
# ALMA sample
# ***********

# plotting parameters
acolor = 'green'
asimcolor = '0.2'
ams = 5
afmt = 's-'

fluxcomponent_file = '../Data/table_intrinsic.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='target')
nmultiples = len(fluxcomponent)

alma = sep_util.getSeparation(fluxcomponent, fluxstring='f870')
avgsep_alma, wmeansep_alma, ra_alma, dec_alma = alma

indexsort = numpy.argsort(avgsep_alma)
avgsep_alma = avgsep_alma[indexsort]
flux_alma = fluxcomponent['f870'][indexsort]
nflux = flux_alma.size
sumflux_alma = numpy.zeros(nflux)
for i in range(nflux):
    sumflux_alma[i] = flux_alma[0:i].sum()
#plt.plot(avgsep_alma, sumflux_alma)
#plt.show()

#import pdb; pdb.set_trace()

sep_util.histArea(avgsep_alma, nbins, color=acolor, fmt=afmt, ms=ams,
        norm=nmultiples, label='Herschel-ALMA')

sep_util.simArea(fluxcomponent, nsim, bin_edges, fluxstring='f870',
        edgecolor=asimcolor, facecolor='none', hatch='\\', norm=nmultiples,
        label='Randomly Distributed')

hayward = Table.read('../Data/dNdA_40arcsec_bright.txt', format='ascii')
xxx = hayward['separation']
yyy = hayward['dNdA']

plt.plot(xxx, yyy, '+-', ms=8, mew=1.5, color='blue', linewidth=1.5, 
        label='HB13 Simulation')

#import pdb; pdb.set_trace()

#plt.semilogy()

xmin = 0
ymin = 0
xmax = 6
ymax = 0.15
plt.axis([xmin, xmax, ymin, ymax])

plt.xlabel(r'${\rm Radial\,Offset\,from\,Centroid\,(arcsec)}$', fontsize='large')
plt.ylabel(r'$dN/dA \, ({\rm arcsec}^{-2}$)', fontsize='large')
plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

fake = numpy.arange(2) + 1e5
#plt.plot(fake, color=hodgecolor, label='ALESS')
#plt.plot(fake, 'x-', ms=6, mew=2, linewidth=1.5, color=c15color, label='C15 Simulation')
##plt.plot(fake, color=hodgesimcolor, linestyle='--', 
##        label='Hodge+13 Random')
#plt.plot(fake, color=acolor, label='Herschel-ALMA')
#plt.plot(fake, color=simcolor, label='Herschel-ALMA ALESS-Sim')
#plt.plot(fake, color=asimcolor, linestyle='--', 
#        label='Randomly Distributed')
plt.legend(loc='upper right', numpoints=1, handletextpad=0.35, borderpad=0.4,
        labelspacing=0.18, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='medium')
plt.subplots_adjust(left=0.14, right=0.95, top=0.97, bottom=0.13, wspace=0.39)

savefig('../Figures/dNdA.pdf')
import pdb; pdb.set_trace()
