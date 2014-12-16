"""
2014 February 11

Shane Bussmann

Compare uvfit10 and uvfit11 model results.

uvfit10 = MLE + statwt + CASA-bin
uvfit11 = chi2 + statwt + CASA-bin
"""

from astropy.table import Table
import glob
import modifypdf
import matplotlib.pyplot as plt
from pylab import savefig
import matplotlib

modelloc = '../../ModelFits/'
uvfitdirs = ['uvfit08', 'uvfit09', 'uvfit10', 'uvfit11']
pdflocs = 'posteriorpdf_burnin.*.dat'

ncol = 2

# set font properties
font = {'family' : 'Arial Narrow',
        'weight' : 'bold',
        'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

targlistloc = 'targetlist.dat'
targlist = Table.read(targlistloc, format='ascii')
ntarg = len(targlist)

for itarg in range(0, ntarg):
    dataname = targlist['dataname'][itarg]
    print "Working on " + dataname

    for uvfitdir in uvfitdirs:

        fitdir = modelloc + dataname + '/' + uvfitdir + '/'
        fitfiles = glob.glob(fitdir + pdflocs)
        fitfile = fitfiles[-1]
        fitresults = Table.read(fitfile, format='ascii', data_start=-5000)

        fitresults = modifypdf.prune(fitresults)
        columnheaders = fitresults.keys()
        columnheaders.sort()
        nparam = len(columnheaders)
        nrow = nparam / ncol + 1
        if uvfitdir == uvfitdirs[0]:
            fig = plt.figure(figsize=(8.0, 2.0 * nrow))

        counter = 1
        for key in columnheaders:
            plt.subplot(nrow, ncol, counter)
            if uvfitdir == 'uvfit08':
                hatch = '//'
                fc = 'none'
                ec = 'red'
                alpha = 1.0
            if uvfitdir == 'uvfit09':
                hatch = '\\'
                fc = 'none'
                ec = 'blue'
                alpha = 1.0
            if uvfitdir == 'uvfit10':
                hatch = ''
                fc = 'red'
                ec = 'black'
                alpha = 0.5
            if uvfitdir == 'uvfit11':
                hatch = ''
                fc = 'blue'
                ec = 'black'
                alpha = 0.5
            plt.hist(fitresults[key], histtype='stepfilled', alpha=alpha,
                    hatch=hatch, ec=ec, fc=fc)
            plt.xlabel(key)
            plt.ylabel('N')
            counter += 1

    plt.subplots_adjust(left=0.1, right=0.97, top=0.99, bottom=0.15, hspace=0.4)
    savefig('modelcompare/' + dataname + '.modelcompare.pdf')
