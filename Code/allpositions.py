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
from scipy.io import readsav


model = 'Hayward'

# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.2

fig = plt.figure(figsize=(4.5, 4.5))
plt.clf()

# Plotting parameters
hodgecolor = 'LightPink'
hodgesimcolor = 'LightPink'
hodgems = 4
hodgefmt = 'D'

# Load the data
fluxcomponent_file = '../Data/hodge2013.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='lessid')
nmultiples = len(fluxcomponent)

sep_hodge = sep_util.getSeparation(fluxcomponent, targetstring='lessid',
        rastring='ra_alma', decstring='dec_alma')
avgsep_hodge = sep_hodge[0]
wmeansep_hodge = sep_hodge[1]
ra_hodge = numpy.array(sep_hodge[2]) * 3600
dec_hodge = numpy.array(sep_hodge[3]) * 3600

plt.plot(ra_hodge, dec_hodge, hodgefmt, alpha=0.8, color=hodgecolor,
        ms=hodgems, label='ALESS')

# ***********
# ALMA sample
# ***********

# plotting parameters
acolor = 'green'
asimcolor = 'green'
ams = 5
afmt = 's'

fluxcomponent_file = '../Data/table_intrinsic.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

# filter out single source systems
fluxcomponent = sep_util.rmSingles(fluxcomponent, targetstring='target')
nmultiples = len(fluxcomponent)

sep_alma = sep_util.getSeparation(fluxcomponent, fluxstring='f870')
avgsep_alma = sep_alma[0]
wmeansep_alma = sep_alma[1]
ra_alma = numpy.array(sep_alma[2]) * 3600
dec_alma = numpy.array(sep_alma[3]) * 3600

plt.plot(ra_alma, dec_alma, afmt, alpha=0.8, color=acolor, ms=ams, 
        label='Herschel-ALMA')

#if model == 'Hayward':
# add Hayward et al. simulated galaxies
h13 = readsav('../Data/100_brightest_shane_data.sav')
flux_h13 = h13['s850']
#hithresh = flux_h13 > 2
flux_h13 = h13['bright_s850']#[hithresh]
ra_h13 = h13['bright_dra']#[hithresh]
dec_h13 = h13['bright_ddec']#[hithresh]
smod = 'HB13'
#else:
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

ra_c15 = numpy.array(ra_c15) * 3600
dec_c15 = numpy.array(dec_c15) * 3600

#flux_c15 = c15['GalaxyS850']
#ra_c15 = c15['GalaxyX']
#dec_c15 = c15['GalaxyY']
#smod = 'C15'

#if model == 'Hayward':
plt.plot(ra_h13, dec_h13, '+', color='blue', label='HB13')
#else:
plt.plot(ra_c15, dec_c15, 'x', color='orange', label='C15')
#plt.hexbin(ra_c15, dec_c15, cmap='rainbow', gridsize=15)
#cbar = plt.colorbar()
#if model == 'Hayward':
#    cbar.set_label(r'$N_{\rm HB13}$')
#else:
#cbar.set_label(r'$N_{\rm C15}$')

xmin = -6
ymin = -6
xmax = 6
ymax = 6
plt.axis([xmin, xmax, ymin, ymax])

plt.xlabel(r'${\rm RA\,Offset\,from\,Centroid\,(arcsec)}$', fontsize='large')
plt.ylabel(r'${\rm Dec\,Offset\,from\,Centroid\,(arcsec)}$', fontsize='large')
plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

fake = numpy.arange(2) + 1e5
plt.legend(loc='lower left', numpoints=1, handletextpad=0.35, borderpad=0.4,
        labelspacing=0.18, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='medium')
plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.12, wspace=0.39)

savefig('../Figures/AllPositions.pdf')

import pdb; pdb.set_trace()
