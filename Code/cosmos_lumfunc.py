"""

Estimate luminosity function in COSMOS from interferometric follow-up of
Miettinen+ 2014, Younger+ 2007, and Younger+2009.

"""

import numpy
import matplotlib.pyplot as plt
from pylab import savefig
from astropy.table import Table
import matplotlib


def MonteCarloCounts(fluxes, errors):

    hist890, bin_edges = numpy.histogram(fluxes)
    nbins = bin_edges.size - 1

    nsource = fluxes.size
    nsim = 1000
    obsPDF = numpy.zeros([nsource, nsim])
    for i in range(nsource):
        imean = fluxes[i]
        irms = errors[i]
        obsPDF[i, :] = numpy.random.normal(loc=imean, scale=irms, size=nsim)

    histPDF = numpy.zeros([nbins, nsim])

    for isim in range(nsim):

        hist, bedge = numpy.histogram(obsPDF[:, isim], bins=bin_edges)
        histPDF[:, isim] = hist

    histmean = numpy.mean(histPDF, axis=1)
    histrms = numpy.std(histPDF, axis=1)

    return histmean, histrms, bin_edges

S_1100 = [10.7, 9.0, 7.6, 6.8, 7.6, 7.9, 8.3, 5.5, 5.8, 4.7, 4.7, 4.5, 4.4,
        4.3, 4.3, 4.2]
S_1100 = numpy.array(S_1100)

S_890 = [15.6, 12.4, 8.7, 14.4, 9.3, 8.6, 12.0, 19.7, 9.0, 5.3, 14.4,
        13.5, 8.2, 5.0, 3.9, 4.4]
e_S_890 = [1.1, 1.0, 1.5, 1.9, 1.3, 1.3, 1.5, 1.8, 2.2, 1.0, 2.9, 1.8, 1.8,
        1.0, 1.0, 1.0]
S_890 = numpy.array(S_890)
e_S_890 = numpy.array(e_S_890) * 1.2

plt.clf()
plt.plot(S_1100, S_890/S_1100, 'o')

# This plot illustrates that the typical correction factor from total 1.1mm
# flux density to total 890um flux density is ~1.5

Oskari1300 = [2.07, 2.15, 1.58, 1.53, 1.78, 1.04, 4.82, 5.72, 1.85, 3.37, 2.19,
        1.27, 1.82, 0.99, 1.41, 1.79, 1.72, 2.85, 0.98, 0.90, 3.36, 2.38, 2.45,
        9.01, 1.53]
e_Oskari1300 = [0.62, 0.63, 0.43, 0.46, 0.54, 0.36, 1.33, 1.85, 0.49, 1.03,
        0.83, 0.40, 0.59, 0.29, 0.42, 0.53, 0.53, 0.78, 0.36, 0.28, 0.97, 0.77,
        0.67, 2.39, 0.45]

Oskari890 = numpy.array(Oskari1300) * 2.5
e_Oskari890 = numpy.array(e_Oskari1300) * 2.5 * 1.5

cosmos890 = numpy.append(S_890, Oskari890)
e_cosmos890 = numpy.append(e_S_890, e_Oskari890)
#cosmos890 = S_890
#e_cosmos890 = e_S_890

completeness = [0.28, 0.5, 0.8, 0.9, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0]

# AzTEC coverage in COSMOS is 0.15 deg^2, centered on z=0.7 overdensity
MCresult = MonteCarloCounts(cosmos890, e_cosmos890)
hist890mean = MCresult[0]
hist890rms = MCresult[1]
bin_edges = MCresult[2]
bin_width = bin_edges[1] - bin_edges[0]
area_aztec = 0.15
norm890 = hist890mean / area_aztec
e_norm890 = hist890rms / area_aztec
nbins = norm890.size

cum890 = numpy.zeros(nbins)
bin_centers = numpy.zeros(nbins)
for ibin in range(nbins):
    bin_centers[ibin] = (bin_edges[ibin] + bin_edges[ibin + 1]) / 2
    cum890[ibin] = norm890[ibin:].sum()

diff890 = norm890 / bin_width / completeness #/ bin_centers
e_diff890 = e_norm890 / bin_width / completeness #/ bin_centers

# Barger catalog
bargerloc = '../Data/barger_catalog.txt'
bargercat = Table.read(bargerloc, format='ascii')
bargerfluxes = bargercat['S860']
e_bargerfluxes = bargercat['e_S860'] * 1.2

MCresult = MonteCarloCounts(bargerfluxes, e_bargerfluxes)
barger890mean = MCresult[0]
barger890rms = MCresult[1]
bin_edges = MCresult[2]
bin_width = bin_edges[1] - bin_edges[0]
area_barger = 0.09
barger890 = barger890mean / area_barger
e_barger890 = barger890rms / area_barger

diffbarger890 = barger890 / bin_width# / completeness #/ bin_centers
e_diffbarger890 = e_barger890 / bin_width# / completeness #/ bin_centers

nbins = barger890.size
cum890 = numpy.zeros(nbins)
barger_bin_centers = numpy.zeros(nbins)
for ibin in range(nbins):
    barger_bin_centers[ibin] = (bin_edges[ibin] + bin_edges[ibin + 1]) / 2

# Smolcic catalog
smolcicloc = '../Data/smolcic_catalog.txt'
smolciccat = Table.read(smolcicloc, format='ascii')
smolcicfluxes = smolciccat['S1300'] * 2.5
e_smolcicfluxes = smolciccat['e_S1300'] * 2.5 * 1.5

MCresult = MonteCarloCounts(smolcicfluxes, e_smolcicfluxes)
smolcic890mean = MCresult[0]
smolcic890rms = MCresult[1]
bin_edges = MCresult[2]
bin_width = bin_edges[1] - bin_edges[0]
area_smolcic = 0.7 / 3.5
smolcic890 = smolcic890mean / area_smolcic
e_smolcic890 = smolcic890rms / area_smolcic

diffsmolcic890 = smolcic890 / bin_width / completeness #/ bin_centers
e_diffsmolcic890 = e_smolcic890 / bin_width / completeness #/ bin_centers

nbins = smolcic890.size
cum890 = numpy.zeros(nbins)
smolcic_bin_centers = numpy.zeros(nbins)
for ibin in range(nbins):
    smolcic_bin_centers[ibin] = (bin_edges[ibin] + bin_edges[ibin + 1]) / 2

# ALESS number counts from Karim et al. 2013
alesscounts = [52.3, 32.3, 24.9, 15.6, 1.6]#, 0.0, 0.0]
e_alesscounts = [18.2, 13.6, 7.9, 12.2, 7.2]#, 0.0, 0.0]
alessfluxes = [4.8, 5.9, 7.5, 8.8, 9.7]#, 11.0, 14.0]

#shadescounts = [2506, 844, 362, 150, 68, 33, 15, 7.4, 3.9, 2.0]
shadescounts = [831, 240, 106, 41, 17, 8.8, 3.9, 1.8, 1.0, 0.6]
shadescounts = numpy.array(shadescounts)
shadesfluxes = numpy.array([2.77, 4.87, 6.90, 8.93, 10.94, 12.95, 14.96, 16.96,
    18.96, 20.97]) / 1.5

# Aretxaga luminosity function

true_centers = [1.41, 2.44, 3.44, 4.45, 5.45, 6.46, 7.46, 8.46, 9.46, 10.46,
        11.46]
true_centers = numpy.array(true_centers)
true_edges = true_centers - 0.5
true_edges = numpy.append(true_edges, true_centers[-1] + 0.5)
true_diffaretxaga = [394, 269, 176, 99.5, 49.9, 22.3, 10.3, 5.83, 4.07, 2.94,
        1.87]
true_diffaretxaga = numpy.array(true_diffaretxaga)

aretxaga = Table.read('../Data/aretxagacatalog.fits')
aretxaga_S1100 = aretxaga['S1_1mm']
hist_aretxaga, edge_aretxaga = numpy.histogram(aretxaga_S1100, bins=true_edges)
nbins = hist_aretxaga.size

#cum890 = numpy.zeros(nbins)
aretxaga_centers = numpy.zeros(nbins)
for ibin in range(nbins):
    aretxaga_centers[ibin] = (edge_aretxaga[ibin] + edge_aretxaga[ibin + 1]) / 2
    #cum890[ibin] = norm890[ibin:].sum()

area_aretxaga = 0.71
normaretxaga = hist_aretxaga / area_aretxaga
aretxaga_completeness = [0.5, 0.85, 0.92, 0.95, 0.97, 0.98, 0.99, 1.0, 1.0,
        1.0, 1.0]
aretxaga_completeness = numpy.array(aretxaga_completeness)
diffaretxaga = normaretxaga / aretxaga_completeness #/ aretxaga_centers

# set font properties
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

fig = plt.figure(figsize=(5.0, 4.5))

plt.clf()

# plot the intrinsic luminosity functions used to make predictions shown in
# mag-flux.py
Sstar = 7.
nstar = 424.
alpha = 1.9

Svector = numpy.arange(1e3)/10
dndS1 = nstar / Sstar
dndS2 = (Svector / Sstar) ** (-alpha)
dndS3 = numpy.exp(-Svector / Sstar)
dndS = dndS1 * dndS2 * dndS3

#line1, = plt.plot(Svector, dndS, color='blue', label='Schechter')

Sstar = 8.
Nstar = 20.
beta1 = 2.0
beta2 = 6.9
dndS1 = Nstar * (Svector / Sstar) ** (-beta1)
dndS2 = Nstar * (Svector / Sstar) ** (-beta2)
dndS = dndS1
high = Svector > Sstar
dndS[high] = dndS2[high]

#line2 = plt.plot(Svector, dndS, color='black', lw=1.5, label='Karim+ 2013')
line2 = plt.plot(Svector, dndS, color='magenta', lw=1.5, label='Karim+ 2013')

Sstar = 15.
Nstar = 5.
beta1 = 2.0
beta2 = 6.9
dndS1 = Nstar * (Svector / Sstar) ** (-beta1)
dndS2 = Nstar * (Svector / Sstar) ** (-beta2)
dndS = dndS1
high = Svector > Sstar
dndS[high] = dndS2[high]

line3, = plt.plot(Svector, dndS, color='blue', lw=1.5,
        label=r'PL, $S_\star = 15\,{\rm mJy}$')

data1, = plt.plot(bin_centers, diff890, 'o', label='COSMOS', color='#A9D0F5')
plt.errorbar(bin_centers, diff890, yerr=e_diff890, fmt='o',
        ecolor='gray', capsize=0, color='#A9D0F5')
data2, = plt.plot(alessfluxes, alesscounts, 'D', label='ALESS', color='pink')
plt.errorbar(alessfluxes, alesscounts, yerr=e_alesscounts, fmt='D',
        ecolor='gray', capsize=0, color='pink')
#data3, = plt.plot(barger_bin_centers, diffbarger890, 's', label='Barger',
#        color='orange')
#plt.errorbar(barger_bin_centers, diffbarger890, yerr=e_diffbarger890, 
#        fmt='s', ecolor='gray', capsize=0, color='orange')
#data4, = plt.plot(smolcic_bin_centers, diffsmolcic890, 's', label='Smolcic',
#        color='orange')
#plt.errorbar(smolcic_bin_centers, diffsmolcic890, yerr=e_diffsmolcic890, 
#        fmt='s', ecolor='gray', capsize=0, color='orange')
#plt.plot(shadesfluxes, shadescounts, 's', label='SHADES')
#plt.hist(cosmos890, cumulative=-1)
#plt.plot(aretxaga_centers, diffaretxaga, '+', label='Aretxaga+ 2011: Me')
#plt.plot(true_centers, true_diffaretxaga, 'x', label='Aretxaga+ 2011: True')
#plt.loglog()
plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')

plt.minorticks_on()
plt.tick_params(width=1.2, which='both')
plt.tick_params(length=2, which='minor')
plt.tick_params(length=4, which='major')

plt.axis([01., 120, .001, 300])
first_legend = plt.legend(loc='lower left', numpoints=1, handletextpad=0.35, 
        borderpad=0.4, labelspacing=0.18, handlelength=1.0)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()

#ax = plt.gca().add_artist(first_legend)

# Create another legend for the second line.
#plt.legend(handles=[line2], loc=4)

#plt.setp(ltext, fontsize='medium')
plt.subplots_adjust(left=0.15, right=0.98, top=0.97, bottom=0.13, wspace=0.39)

plt.ylabel(r'$dN/dS\;{\rm (mJy}^{-1} \, {\rm deg}^{-2})$', fontsize='large')
plt.xlabel(r'$S_{870}\;{\rm (mJy)}$', fontsize='large')
savefig('../Figures/DifferentialNumberCounts.pdf')
import pdb; pdb.set_trace()
