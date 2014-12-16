"""
Created 2013 March 12
Shane Bussmann

Purpose: Make color-flux plots of Herschel SPIRE photometry
for lensed sample and compare with H-ATLAS galaxies
From Wardlow et al. 2012: 
    - F350/F500 vs. F500
"""

# import modules and define plot appearance
from astropy.io import fits
import matplotlib
from astropy.table import Table
from astropy.io import ascii
#matplotlib.use('pdf')
import numpy
import matplotlib.pyplot as mpl
from pylab import savefig
import pdb


# set font properties
font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 14}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

# read in the plotting color scheme
colorschemeloc = 'colorscheme.txt'
colorscheme = Table.read(colorschemeloc, format='ascii')
ecolor = str(colorscheme['ecolor'][0])
color3 = colorscheme['color3'][0]
color4 = colorscheme['color4'][0]
gray = colorscheme['color6'][0]

# set SMA plotting parameters
acolor = colorscheme['color5'][0]
ams = 6
afmt = 's'

# set SMA plotting parameters
bcolor = colorscheme['color1'][0]
bms = 6
bfmt = 'o'

# set SPT plotting parameters
hcolor = colorscheme['color2'][0]
hms = 7
hfmt = '*'

# set the level at which a Herschel source is a lens candidate
candidate1 = 0.10

# read in the SPIRE photometry for the ALMA sample
almaspirephotloc = '../Data/targetlist.dat'
#print spirephotloc
almaspirephot = Table.read(almaspirephotloc, format='ascii')

alma_f250 = almaspirephot['f250']
alma_f350 = almaspirephot['f350']
alma_f500 = almaspirephot['f500']

alma_e250 = almaspirephot['e250']
alma_e350 = almaspirephot['e350']
alma_e500 = almaspirephot['e500']

# read in the SPIRE photometry for lensed SMGs
spirephotloc = '../Data/table_sma.dat'
#print spirephotloc
spirephot = Table.read(spirephotloc, format='ascii')
#spirephot = asciidata.open(spirephotloc)

f250 = spirephot['f250']
f350 = spirephot['f350']
f500 = spirephot['f500']

e250 = spirephot['e250']
e350 = spirephot['e350']
e500 = spirephot['e500']

# read in SPIRE photometry for all H-ATLAS galaxies
atlasloc = '../../HATLASCatalog/'
g09 = fits.getdata(atlasloc+'G09_Phase1_release_v2_allids.fits.gz', 1)
g12 = fits.getdata(atlasloc+'G12_Phase1_release_v2_allids.fits.gz', 1)
g15 = fits.getdata(atlasloc+'G15_Phase1_release_v2_allids.fits.gz', 1)
ngpaloc = atlasloc+'NGP-h1-v1-v2_madx5_19-Feb-2011-110650_id_ac_names.dat'
ngpa = Table.read(ngpaloc, format='ascii')
ngpbloc = atlasloc+'NGP-h5-h7-v5_madx5_19-Feb-2011-104234_id_ac_names.dat'
ngpb = Table.read(ngpbloc, format='ascii')
ngpcloc = atlasloc+'NGP-h6-h8-v7_madx5_19-Feb-2011-014156_id_ac_names.dat'
ngpc = Table.read(ngpcloc, format='ascii')

# read in lensed candidates in GAMA + NGP fields
mattialoc = '../Data/mattialenses_2013june28.dat'
mattia = Table.read(mattialoc, format='ascii')

fig = mpl.figure(figsize=(5.0, 3.5))
ax = fig.add_subplot(1, 1, 1)
#mpl.subplots_adjust(left=0.16, right=1.0, top=0.96, bottom=0.21, wspace=0.0)

# signal-to-noise requirement for H-ATLAS targets to be plotted
snr = 3.0
hatlas250 = numpy.append(g09['f250_best'], g12['f250_best'])
hatlas350 = numpy.append(g09['f350_best'], g12['f350_best'])
hatlas500 = numpy.append(g09['f500_best'], g12['f500_best'])
eatlas250 = numpy.append(g09['e250_best'], g12['e250_best'])
eatlas350 = numpy.append(g09['e350_best'], g12['e350_best'])
eatlas500 = numpy.append(g09['e500_best'], g12['e500_best'])
hatlas250 = numpy.append(hatlas250, g15['f250_best'])
hatlas350 = numpy.append(hatlas350, g15['f350_best'])
hatlas500 = numpy.append(hatlas500, g15['f500_best'])
eatlas250 = numpy.append(eatlas250, g15['e250_best'])
eatlas350 = numpy.append(eatlas350, g15['e350_best'])
eatlas500 = numpy.append(eatlas500, g15['e500_best'])
hatlasbig250 = numpy.append(hatlas250, ngpa['f250'])
hatlasbig350 = numpy.append(hatlas350, ngpa['f350'])
hatlasbig500 = numpy.append(hatlas500, ngpa['f500'])
hatlasbig250 = numpy.append(hatlasbig250, ngpb['f250'])
hatlasbig350 = numpy.append(hatlasbig350, ngpb['f350'])
hatlasbig500 = numpy.append(hatlasbig500, ngpb['f500'])
hatlasbig250 = numpy.append(hatlasbig250, ngpc['f250'])
hatlasbig350 = numpy.append(hatlasbig350, ngpc['f350'])
hatlasbig500 = numpy.append(hatlasbig500, ngpc['f500'])
#hatlas250 = numpy.append(hatlas250, ngpa['f250'])
#hatlas350 = numpy.append(hatlas350, ngpa['f350'])
#hatlas500 = numpy.append(hatlas500, ngpa['f500'])
#eatlas250 = numpy.append(eatlas250, ngpa['e250'])
#eatlas350 = numpy.append(eatlas350, ngpa['e350'])
#eatlas500 = numpy.append(eatlas500, ngpa['e500'])
#hatlas250 = numpy.append(hatlas250, ngpb['f250'])
#hatlas350 = numpy.append(hatlas350, ngpb['f350'])
#hatlas500 = numpy.append(hatlas500, ngpb['f500'])
#eatlas250 = numpy.append(eatlas250, ngpb['e250'])
#eatlas350 = numpy.append(eatlas350, ngpb['e350'])
#eatlas500 = numpy.append(eatlas500, ngpb['e500'])
#hatlas250 = numpy.append(hatlas250, ngpc['f250'])
#hatlas350 = numpy.append(hatlas350, ngpc['f350'])
#hatlas500 = numpy.append(hatlas500, ngpc['f500'])
#eatlas250 = numpy.append(eatlas250, ngpc['e250'])
#eatlas350 = numpy.append(eatlas350, ngpc['e350'])
#eatlas500 = numpy.append(eatlas500, ngpc['e500'])
hatlasra = numpy.append(g09['RA_J2000'], g12['RA_J2000'])
hatlasra = numpy.append(hatlasra, g15['RA_J2000'])
hatlasbigra = numpy.append(hatlasra, ngpa['ra_new'])
hatlasbigra = numpy.append(hatlasbigra, ngpb['ra_new'])
hatlasbigra = numpy.append(hatlasbigra, ngpc['ra_new'])
hatlasdec = numpy.append(g09['DEC_J2000'], g12['DEC_J2000'])
hatlasdec = numpy.append(hatlasdec, g15['DEC_J2000'])
hatlasbigdec = numpy.append(hatlasdec, ngpa['dec_new'])
hatlasbigdec = numpy.append(hatlasbigdec, ngpb['dec_new'])
hatlasbigdec = numpy.append(hatlasbigdec, ngpc['dec_new'])
#hatlasdec = numpy.append(hatlasdec, ngpa['dec_new'])
#hatlasdec = numpy.append(hatlasdec, ngpb['dec_new'])
#hatlasdec = numpy.append(hatlasdec, ngpc['dec_new'])
hatlasr = numpy.append(g09['sdss_rmodelmag'], g12['sdss_rmodelmag'])
hatlasr = numpy.append(hatlasr, g15['sdss_rmodelmag'])
hatlassep = numpy.append(g09['sdss_sep'], g12['sdss_sep'])
hatlassep = numpy.append(hatlassep, g15['sdss_sep'])

# replace fluxes in Mattia's lens catalog with those in the most recent H-ATLAS
# data release
#nmattia = mattia.size
#for i in numpy.arange(nmattia):
#    rai = mattia['ra'][i]
#    deci = mattia['dec'][i]
#    offra = rai - hatlasbigra
#    offdec = deci - hatlasbigdec
#    offset = numpy.sqrt(offra ** 2 + offdec ** 2) * 3600
#    match = offset < 4.
#    nmatch = offset[match].size
#    print mattia['iauname'][i], offset[match]
#    if nmatch > 0:
#        print mattia['f250'][i], hatlasbig250[match][0] * 1e3, \
#                mattia['f350'][i], hatlasbig350[match][0] * 1e3, \
#                mattia['f500'][i], hatlasbig500[match][0] * 1e3
#        mattia['f250'][i] = hatlasbig250[match][0] * 1e3
#        mattia['f350'][i] = hatlasbig350[match][0] * 1e3
#        mattia['f500'][i] = hatlasbig500[match][0] * 1e3

# plot the H-ATLAS galaxies
good = (hatlas250 > eatlas250) & \
    (hatlas350 > eatlas350) & \
    (hatlas500 > eatlas500) 
yyy = hatlas350[good] / hatlas500[good]
xxx = hatlas500[good] * 1e3
mpl.hexbin(xxx, yyy, cmap='gray_r', label='H-ATLAS Phase I Galaxies',
    xscale='log', bins='log')
#mpl.plot(xxx, yyy, '.', color='brown')

# plot the lensed candidates
hgood = (hatlas500 > candidate1) 
yyy = hatlas350[hgood] / hatlas500[hgood]
xxx = hatlas500[hgood] * 1e3
#mpl.plot(xxx, yyy, '.', color='brown')
#mpl.plot(xxx, yyy, ',', color='blue', label='Local Spirals')

# plot the vetted lens candidates
gmattia = (mattia['f500'] > candidate1*1e3)# & \
#        (mattia['f350']/mattia['f500'] > 0.5)
yyyup = numpy.round(mattia['f350'][gmattia])
yyydown = numpy.round(mattia['f500'][gmattia])
yyy = yyyup / yyydown
xxx = numpy.round(mattia['f500'][gmattia])
mpl.plot(xxx, yyy, 'o', mfc='yellow', label='Herschel',
        markersize=bms/2., zorder=6)
mhihi = mattia['f500'] > 170
mattiahi = mattia[mhihi]

# print the RA and Dec of the lensed candidates
#hra = hatlasra[good]
#hdec = hatlasdec[good]
#for i in numpy.arange(hra.size): print i,',', hra[i], ',', hdec[i]

# plot the lensed systems
good = (f250 > e250) & \
    (f350 > e350) & \
    (f500 > e500) 
yyy = f350[good]/f500[good]
xxx = f500[good]
mpl.plot(xxx, yyy, bfmt, color=bcolor, label='SMA Subsample', ms=bms, zorder=5)
shihi = f500 > 170
smahi = spirephot[shihi]

# plot the ALMA sample
yyy = alma_f350 / alma_f500
xxx = alma_f500
mpl.plot(xxx, yyy, afmt, color=acolor, label='ALMA Subsample', ms=ams, zorder=4)

# add representative error bars
err_rep = numpy.median(e500)
err_col = numpy.sqrt((e350/f350)**2 + (e500/f500)**2) * f350/f500
err_col_rep = numpy.median(err_col)
xspot = numpy.array([20])
yspot = numpy.array([0.3])
mpl.errorbar(xspot, yspot, yerr=err_col_rep, xerr=err_rep, ecolor='black',
        capsize=0)
nlens = f250[good].size
#for i in numpy.arange(nlens):
#    namei = spirephot['shortname'][good][i]
#    xi = xxx[i]
#    yi = yyy[i]
#    if namei[0] == 'N':
#        mpl.annotate(namei, xy=(xi, yi), fontsize='x-small')

# plot SPT lensed SMGs
h13loc = '../Data/hezaveh_fluxes.dat'
h13dat = Table.read(h13loc, format='ascii')
h13dat['f350'][-2:] = h13dat['f350'][-2:].sum()
h13dat['f500'][-2:] = h13dat['f500'][-2:].sum()
xxx = h13dat['f500']
yyy = h13dat['f350'] / h13dat['f500']
mpl.plot(xxx, yyy, hfmt, color=hcolor, label='SPT', ms=hms, zorder=7)

# plot the stacked ALESS SED
xxx = [18.5]
yyy = [1.1]
#mpl.plot(xxx, yyy, 'D', color='magenta', ms=bms, zorder=3)
#mpl.text(21, 0.9, 'ALESS', color='white', fontsize='x-small')

# plot the simulated SMGs from Cowley et al. 2014
cowleycat = ascii.read('../Data/SPIRE_ALMA_Cat_v4.txt')
xxx = cowleycat['SourceS500']
yyy = cowleycat['SourceS350'] / cowleycat['SourceS500']
mpl.plot(xxx, yyy, '.', color='teal', zorder=8, label='Cowley+2014')


mpl.semilogx()

xmin = 10
xmax = 1000.
ymin = 0.0
ymax = 7.0
mpl.axis([xmin, xmax, ymin, ymax])
cb = mpl.colorbar(use_gridspec=True)
cb.set_label('log10(N)')
mpl.xlabel(r'$S_{500}\,({\rm mJy})$', fontsize='large')
mpl.ylabel(r'$S_{350}/S_{500}$', fontsize='large')
mpl.minorticks_on()
mpl.tick_params(width=1.5, which='both')
mpl.tick_params(length=2, which='minor')
mpl.tick_params(length=4, which='major')
mpl.subplots_adjust(left=0.12, right=0.96, top=0.96, bottom=0.18, wspace=0.39)

mpl.legend(loc='upper right', numpoints=1, handletextpad=0.00, borderpad=0.3,
        labelspacing=0.15, handlelength=1.0)
leg = mpl.gca().get_legend()
ltext  = leg.get_texts()
mpl.setp(ltext, fontsize='x-small')

#mpl.axvline(x=170)

#mpl.tight_layout(pad=0.3)

# rename files to remove '.' from the name
saveloc = '../Figures/spirecolflux_noaless.pdf'
savefig(saveloc)
pdb.set_trace()
