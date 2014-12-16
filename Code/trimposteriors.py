"""
2013 July 10
Shane Bussmann

Trim finalposterior.dat down to the final N iterations
"""

import glob
from astropy.table import Table
import numpy
import os
import modifypdf


# set the uvfit ID to be trimmed and pruned
fitid = 'uvfit11'

# set the original iterations
burninfile = 'posteriorpdf_burnin.*.dat'

cwd = os.getcwd()
anloc = '../../ModelFits/'
dataphotloc = 'targetlist.dat'
dataphot = Table.read(dataphotloc, format='ascii')
dataname = dataphot['dataname']
ntarg = dataname.size
dataflag = numpy.zeros(ntarg)
for itarg in range(0, ntarg):

    dirloc = anloc + dataname[itarg] + '/' + fitid
    print dirloc
    #if dataname[itarg] == 'XMM01':
    #    continue
    os.chdir(dirloc)
    burninfiles = glob.glob(burninfile)
    burninlast = burninfiles[-1]
    fitresults = Table.read(burninlast, format='ascii', data_start=-5000)
    prunedpdf = modifypdf.prune(fitresults)
    import matplotlib.pyplot as plt
    plt.clf()
    maxprob = fitresults['lnprob'].max() + 0.1
    hist1 = numpy.log10(maxprob - fitresults['lnprob'])
    plt.hist(hist1, bins=15, log=True)
    hist1 = numpy.log10(maxprob - prunedpdf['lnprob'])
    plt.hist(hist1, bins=15, log=True)
    plt.title(dataname[itarg])
    from pylab import savefig
    savefig('prunecheck')
    os.chdir(cwd)
