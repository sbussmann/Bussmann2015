"""
 Plot posterior distributions of model parameters
 Shane Bussmann
 2014 July 14

"""

import numpy
import os
from astropy.table import Table
import modifypdf
import setuputil
import yaml
from astropy.coordinates import SkyCoord
from astropy import units as u


# get S500 and S850 values from Magnelli
magnelli = Table.read('magnelli_fluxes.dat', format='ascii')
ok = (magnelli['f850'] > 0) & (magnelli['f500'] > 0)
magnelli = magnelli[ok]
avg850to500 = numpy.median(magnelli['f850'] / magnelli['f500'])

# get S500 values for ALMA sample from Herschel
targetlist = Table.read('targetlist.dat', format='ascii')

fluxcomponent_file = 'source_observed.dat'
fluxcomponent = Table.read(fluxcomponent_file, format='ascii')

nindiv = len(fluxcomponent)

superfluxtotal = []
e_superfluxtotal = []
s500a = []
s500b = []
for icomp in range(nindiv):
    target = fluxcomponent['target'][icomp]

    match = fluxcomponent['target'] == target
    fluxsources = fluxcomponent['f870'][match]
    e_fluxsources = fluxcomponent['e_f870'][match]
    ftot = fluxsources.sum()
    e_ftot = e_fluxsources.sum()

    superfluxtotal.append(ftot)
    e_superfluxtotal.append(e_ftot)

    match2 = targetlist['dataname'] == target
    s500tot = targetlist['f500'][match2][0]
    iflux = fluxcomponent['f870'][icomp]
    s500a.append(iflux / ftot * s500tot)

    s500b.append(iflux / avg850to500)

fluxcomponent_alma = fluxcomponent['f870']
fluxtotal_alma = numpy.array(superfluxtotal)
s500a_alma = numpy.array(s500a)
s500b_alma = numpy.array(s500b)

goodfitloc = 'uvfitlist.dat'
goodfitdat = Table.read(goodfitloc, format='ascii')
mu_estimate = False

cwd = os.getcwd()

# read in the sample targets
anloc = '../../ModelFits/'
dataphotloc = 'targetlist.dat'
dataphot = Table.read(dataphotloc, format='ascii')
dataname = dataphot['dataname']
ntarg = dataname.size
dataflag = numpy.zeros(ntarg)
#for i in range(ntarg):
#    datanamei = dataname[i]
#    dataloc = anloc + datanamei + '/' + goodfit
#    if os.path.exists(dataloc):
#        dataflag[i] = 1

#yesmodel = dataflag == 1
#dataphot = dataphot[yesmodel]
dataname = dataphot['dataname']
ntarg = len(dataphot)

counter = 0
for itarg in range(0, ntarg):

    # extract object name from directory
    objname = dataname[itarg]
    iauname = dataphot['iauname'][itarg]
    goodfit = goodfitdat['intrinsic'][itarg]
    fitdir = anloc + dataname[itarg] + '/' + goodfit + '/'
    os.chdir(fitdir)

    # ------------------------------------------------------------------------
    # Read posterior probability distributions
    fitloc = 'posteriorpdf.fits'
    fitresults = Table.read(fitloc)
    fitresults = fitresults[-5000:]
    fitresultsgood = modifypdf.cleanColumns(fitresults)
    fitresultsgood = modifypdf.bigScoop(fitresultsgood)
    fitresultsgood = modifypdf.prune(fitresultsgood, quiet=True)

    configfile = open('config.yaml')
    config = yaml.load(configfile)
    paramData = setuputil.loadParams(config)

    nregions = paramData['nregions']

    flux_target = 0

    os.chdir(cwd)
    for iregion in range(nregions):
        sr = str(iregion)
        region = 'Region' + sr
        nsource = paramData['nsource_regions'][iregion]
        radeg_centroid = config[region]['RACentroid']
        decdeg_centroid = config[region]['DecCentroid']
        for isource in range(nsource):
            ss = str(isource)
            datanamesource = dataname[itarg] + '.Source' + ss
            tag1 = 'IntrinsicFlux_Source' + ss + '_Region' + sr
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Source' + ss  + ' IntrinsicFlux'
            fluxdist = fitresultsgood[tag1]
            tag2 = 'mu_tot.Source' + ss + '.Region' + sr

            if tag2 in fitresults.keys():
                if fitresultsgood[tag2].mean() == 100:
                    continue
                observedflux = fluxdist * fitresultsgood[tag2]
                mu880 = fitresultsgood[tag2].mean()
                e_mu880 = fitresultsgood[tag2].std()
            else:
                observedflux = fluxdist
                mu880 = 1.0
                e_mu880 = 0.0
            flux_target += observedflux
            flux_indiv = observedflux.mean()
            e_flux_indiv = observedflux.std()

            S500a = s500a[counter]
            S500b = s500b[counter]
            e_S500a = s500a[counter] / 5.
            e_S500b = s500b[counter] / 5.

            msg = '{0:15} {1:.3f} {2:.3f} {3:.1f} {4:.1f} {5:.1f} {6:.1f}'
            msgfmt = msg.format(datanamesource, mu880, e_mu880, S500a,
                e_S500a, S500b, e_S500b)
            print(msgfmt)
            counter += 1
        
    flux_measure = flux_target.mean()
    flux_error = flux_target.std()
    #import matplotlib.pyplot as plt
    #plt.hist(flux_target)
    #plt.show()
    #import pdb; pdb.set_trace()
    radeg = radeg_centroid
    decdeg = decdeg_centroid
    c = SkyCoord(ra=radeg*u.degree, dec=decdeg*u.degree)
    fullRADec = c.to_string('hmsdms', sep=':')
    #msg = '{0:20} {1} {2} {3:.2f} {4:.2f}'
    #msgfmt = msg.format(objname, fullRADec, flux_measure, flux_error)
    #print(msgfmt)

