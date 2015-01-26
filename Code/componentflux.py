"""

Shane Bussmann

2014 August 15

Plot total flux of an object vs. the individual component fluxes of that
object.  Like Figure 5 in Hodge+13.

"""

from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
from pylab import savefig
import numpy


# set font properties
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

# *****************
# Hodge+2013 sample
# *****************

# Plotting parameters
hodgecolor = 'LightPink'
hodgems = 3
hodgefmt = 'D'

# Load the data
fluxcomponent_file = '../Data/hodge2013.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

nindiv = len(fluxcomponent)

superfluxtotal = []
e_superfluxtotal = []
for icomp in range(nindiv):
    target = fluxcomponent['lessid'][icomp]

    match = fluxcomponent['lessid'] == target
    fluxsources = fluxcomponent['best880'][match]
    e_fluxsources = fluxcomponent['e_best880'][match]
    ftot = fluxsources.sum()
    e_ftot = e_fluxsources.sum()

    superfluxtotal.append(ftot)
    e_superfluxtotal.append(e_ftot)

fluxcomponent_hodge = fluxcomponent['best880']
fluxtotal_hodge = numpy.array(superfluxtotal)

xxxhodge = numpy.array(superfluxtotal)
xxxerrhodge = numpy.array(e_superfluxtotal)
yyyhodge = fluxcomponent['best880']
yyyerrhodge = fluxcomponent['e_best880']

# **********
# Cowley simulations
# **********

# Plotting parameters
bargcolor = 'LightBlue'
bargms = 4
bargfmt = '*'

# Load the data
fluxcomponent_file = '../Data/SPIRE_ALMA_Cat_v2.txt'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')
hiz = fluxcomponent['z'] > 1.5
fluxcomponent = fluxcomponent[hiz]
hiflux = fluxcomponent['S500'] > 50
fluxcomponent = fluxcomponent[hiflux]

nindiv = len(fluxcomponent)

superfluxtotal = []
e_superfluxtotal = []
for icomp in range(nindiv):
    target = fluxcomponent['TargetID'][icomp]

    match = fluxcomponent['TargetID'] == target
    fluxsources = fluxcomponent['S850'][match]
    e_fluxsources = fluxcomponent['S850'][match] / 10.
    ftot = fluxsources.sum()
    e_ftot = e_fluxsources.sum()

    superfluxtotal.append(ftot)
    e_superfluxtotal.append(e_ftot)

fluxcomponent_barger = fluxcomponent['S850']
fluxtotal_barger = numpy.array(superfluxtotal)

xxxbarger = numpy.array(superfluxtotal)
xxxerrbarger = numpy.array(e_superfluxtotal)
yyybarger = fluxcomponent['S850']
yyyerrbarger = fluxcomponent['S850'] / 10.


# **********
# Barger sample
# **********

# Plotting parameters
bargcolor = 'Black'
bargms = 4
bargfmt = '.'

# Load the data
#fluxcomponent_file = '../Data/barger2012.dat'
#fluxcomponent = Table.read(fluxcomponent_file, format='ascii')
#
#nindiv = len(fluxcomponent)
#
#superfluxtotal = []
#e_superfluxtotal = []
#for icomp in range(nindiv):
#    target = fluxcomponent['name'][icomp]
#
#    match = fluxcomponent['name'] == target
#    fluxsources = fluxcomponent['f860'][match]
#    e_fluxsources = fluxcomponent['e_f860'][match]
#    ftot = fluxsources.sum()
#    e_ftot = e_fluxsources.sum()
#
#    superfluxtotal.append(ftot)
#    e_superfluxtotal.append(e_ftot)
#
#fluxcomponent_barger = fluxcomponent['f860']
#fluxtotal_barger = numpy.array(superfluxtotal)
#
#xxxbarger = numpy.array(superfluxtotal)
#xxxerrbarger = numpy.array(e_superfluxtotal)
#yyybarger = fluxcomponent['f860']
#yyyerrbarger = fluxcomponent['e_f860']

# ***********
# ALMA sample
# ***********

# plotting parameters
acolor = colorscheme['color5'][0]
ams = 4
afmt = 's'

fluxcomponent_file = '../Data/table_intrinsic.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

nindiv = len(fluxcomponent)

superfluxtotal = []
e_superfluxtotal = []
for icomp in range(nindiv):
    target = fluxcomponent['target'][icomp]

    match = fluxcomponent['target'] == target
    fluxsources = fluxcomponent['f870'][match]
    e_fluxsources = fluxcomponent['e_f870'][match]
    ftot = fluxsources.sum()
    e_ftot = e_fluxsources.sum()

    superfluxtotal.append(ftot)
    e_superfluxtotal.append(e_ftot)

fluxcomponent_alma = fluxcomponent['f870']
fluxtotal_alma = numpy.array(superfluxtotal)

xxxalma = superfluxtotal
xxxerralma = e_superfluxtotal
yyyalma = fluxcomponent['f870']
yyyerralma = fluxcomponent['e_f870']

# **********
# SMA sample
# **********

# Plotting parameters
bcolor = colorscheme['color1'][0]
bms = 3
bfmt = 'o'

# Load the data
fluxcomponent_file = '../Data/bussmann2013_muflux.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

nindiv = len(fluxcomponent)

superfluxtotal = []
e_superfluxtotal = []
for icomp in range(nindiv):
    target = fluxcomponent['dataname'][icomp]

    match = fluxcomponent['dataname'] == target
    fluxsources = fluxcomponent['fnu'][match] / fluxcomponent['mu'][match]
    e_fluxsources = fluxcomponent['e_fnu'][match] / fluxcomponent['mu'][match]
    ftot = fluxsources.sum()
    e_ftot = e_fluxsources.sum()

    superfluxtotal.append(ftot)
    e_superfluxtotal.append(e_ftot)


goodmodel = fluxcomponent['mu'] > fluxcomponent['e_mu']

fluxcomponent_sma = fluxcomponent['fnu'][goodmodel] / fluxcomponent['mu'][goodmodel]
fluxtotal_sma = numpy.array(superfluxtotal)[goodmodel]

xxxsma = numpy.array(superfluxtotal)[goodmodel]
xxxerrsma = numpy.array(e_superfluxtotal)[goodmodel]
yyysma = fluxcomponent['fnu'][goodmodel] / fluxcomponent['mu'][goodmodel]
yyyerrsma = fluxcomponent['e_fnu'][goodmodel] / fluxcomponent['mu'][goodmodel]

# set the plotting window size
#fig = plt.figure(figsize=(5.0, 4.5))
fig, ax1 = plt.subplots(figsize=(5,4.5))
#ax2 = ax1.twinx()

# plot the average S_component/S_total vs. S_total
#plt.clf()

def gethist(xflux, componentflux, totalflux):
    nflux = xflux.size
    avgratio = numpy.zeros(nflux)
    for i in range(nflux - 1):
        iflux = xflux[i]
        iflux1 = xflux[i + 1]
        match = (totalflux > iflux) & (totalflux < iflux1)
        if totalflux[match].size > 0:
            ratioflux = componentflux[match] / totalflux[match]
            avgratio[i] = numpy.mean(ratioflux)

    return avgratio

maxflux = 45
minflux = 0

xflux = numpy.arange(minflux, maxflux, 1)

ratio_hodge = gethist(xflux, fluxcomponent_hodge, fluxtotal_hodge)
ratio_barger = gethist(xflux, fluxcomponent_barger, fluxtotal_barger)
ratio_sma = gethist(xflux, fluxcomponent_sma, fluxtotal_sma)
ratio_alma = gethist(xflux, fluxcomponent_alma, fluxtotal_alma)
ratio_alma[0] = 0

#ax2.plot(xflux, ratio_hodge, color=hodgecolor, linewidth=2, zorder=-1,
#        drawstyle='steps-post')
#plt.plot(xflux, ratio_barger, ls='steps', color=bargcolor)
#plt.plot(xflux, ratio_sma, ls='steps', color=bcolor)
#ax2.plot(xflux, ratio_alma, color=acolor, linewidth=2, zorder=0, 
#        drawstyle='steps-post')

#ax2.fill_between(xflux, ratio_hodge, hatch='//', edgecolor=hodgecolor,
#        facecolor='none', zorder=0)
#ax2.fill_between(xflux, ratio_alma, hatch='\\', edgecolor=acolor,
#        facecolor='none', zorder=1)
#ax2.fill_between(xflux, ratio_barger, facecolor='gray', alpha=0.5)

xmin = 0
xmax = 20
ymin = 0
ymax = 1.1
#ax2.axis([xmin, xmax, ymin, ymax])

# customize plot
#ax2.set_xlabel(r'$S_{\rm total} \,({\rm mJy}$)', fontsize='x-large')
#ax2.set_ylabel(r'$<S_{\rm component}/S_{\rm total}>$', fontsize='large')
#ax2.minorticks_on()
#ax2.tick_params(width=1.5, which='both')
#ax2.tick_params(length=2, which='minor')
#ax2.tick_params(length=4, which='major')
#plt.subplots_adjust(left=0.14, right=0.96, top=0.97, bottom=0.13, wspace=0.39)

#plt.legend(loc='upper right', numpoints=1, handletextpad=0.35, borderpad=0.4,
#        labelspacing=0.18, handlelength=1.0)
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='small')

# plot the 1:1 line
ax1.plot([0,45], [0,45], color='gray')
#plt.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=',', ecolor='gray',
#        capsize=0)
ax1.plot(xxxhodge, yyyhodge, hodgefmt, ms=hodgems, mfc=hodgecolor,
        mec='black', label='ALESS')
#ax1.plot(xxxhodge, yyyhodge, hodgefmt, ms=hodgems, mfc=hodgecolor,
#        mec='none', label='Hodge+13')

#ax1.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=',', ecolor='gray',
#        capsize=0)
#ax1.plot(xxxbarger, yyybarger, bargfmt, ms=bargms, color=bargcolor, label='Barger+12')

ax1.plot(xxxbarger, yyybarger, bargfmt, ms=bargms, color=bargcolor,
        label='Cowley+2014 SAM')
#ax1.plot(xxxbarger, yyybarger, bargfmt, ms=bargms, color=bargcolor,
#        mec='none', label='Cowley+2014')

#ax1.errorbar(xxx, yyy, yerr=yyyerr, xerr=xxxerr, fmt=',', ecolor='gray',
#        capsize=0)
#ax1.plot(xxxsma, yyysma, bfmt, ms=bms, color=bcolor, label='SMA Sample')

ax1.errorbar(xxxalma, yyyalma, yerr=yyyerralma, xerr=xxxerralma, fmt=',',
        ecolor='gray', capsize=0)
ax1.plot(xxxalma, yyyalma, afmt, ms=ams, color=acolor, label='Herschel-ALMA')

# customize plot
#plt.xlabel(r'${\rm Total} \, S_{870}\,({\rm mJy}$)', fontsize='x-large')
#plt.ylabel(r'${\rm Component} \, S_{870}\,({\rm mJy}$)', fontsize='x-large')
ax1.set_xlabel(r'$S_{\rm total}\,({\rm mJy}$)', fontsize='large')
ax1.set_ylabel(r'$S_{\rm component}\,({\rm mJy}$)', fontsize='large')
ax1.minorticks_on()
ax1.tick_params(width=1.5, which='both')
ax1.tick_params(length=2, which='minor')
ax1.tick_params(length=4, which='major')
plt.subplots_adjust(left=0.12, right=0.87, top=0.97, bottom=0.12, wspace=0.39)

xmin = 0
xmax = 25
ymin = 0
ymax = 25
ax1.axis([xmin, xmax, ymin, ymax])

ax1.legend(loc='upper left', numpoints=1, handletextpad=0.00, borderpad=0.4,
        labelspacing=0.35, handlelength=1.0)
#leg = ax1.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='small')
#ax2.set_zorder(0)
ax1.set_zorder(1)
#ax1.set_frame_on(False)
#ax2.set_frame_on(True)

outfile = '../Figures/fluxtotalcomponent.pdf'
savefig(outfile)

#savefig('avgratio')
import pdb; pdb.set_trace()
