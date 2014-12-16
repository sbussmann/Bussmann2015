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


def all(goodfit, update=False, mu_estimate=False):

    # Check to see whether we've already loaded the data
    if not update:
        ffile = 'flux_total_observed.dat'
        if os.path.exists(ffile):
            ffile = 'flux_total_observed.dat'
            flux_total_observed = Table.read(ffile, format='ascii')
            ffile = 'flux_indiv_observed.dat'
            flux_indiv_observed = Table.read(ffile, format='ascii')
            ffile = 'flux_indiv_intrinsic.dat'
            flux_indiv_intrinsic = Table.read(ffile, format='ascii')
        else:
            update = True
    if update:
        cwd = os.getcwd()

        # read in the sample targets
        anloc = '../../ModelFits/'
        dataphotloc = 'targetlist.dat'
        dataphot = Table.read(dataphotloc, format='ascii')
        dataname = dataphot['dataname']
        ntarg = dataname.size
        dataflag = numpy.zeros(ntarg)
        for i in range(ntarg):
            datanamei = dataname[i]
            dataloc = anloc + datanamei + '/' + goodfit
            if os.path.exists(dataloc):
                dataflag[i] = 1

        yesmodel = dataflag == 1
        dataphot = dataphot[yesmodel]
        dataname = dataphot['dataname']
        ntarg = len(dataphot)

        dataname_indiv = []
        flux_total = []
        e_flux_total = []
        flux_indiv = []
        e_flux_indiv = []
        flux_indiv_int = []
        e_flux_indiv_int = []
        for itarg in range(0, ntarg):

            # extract object name from directory
            objname = dataname[itarg]
            iauname = dataphot['iauname'][itarg]
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
                nsource = paramData['nsource_regions'][iregion]
                for isource in range(nsource):
                    ss = str(isource)
                    datanamesource = dataname[itarg] + '.Source' + ss
                    tag1 = 'IntrinsicFlux_Source' + ss + '_Region' + sr
                    if tag1 not in fitresults.keys():
                        tag1 = 'Region' + sr + ' Source' + ss  + ' IntrinsicFlux'
                    fluxdist = fitresultsgood[tag1]
                    tag2 = 'mu_tot.Source' + ss + '.Region' + sr

                    # get the by-eye magnification estimates
                    if mu_estimate:
                        datamu = Table.read('mu_estimate.dat', format='ascii')
                        muindx = datamu['target'] == datanamesource
                        muerr = datamu['muerr'][muindx][0]
                        if muerr == 0.25:
                            mu_by_eye = datamu['mu'][muindx][0]
                        else:
                            mu_by_eye = 1.
                    else:
                        mu_by_eye = 1.
                    flux_indiv_int.append(fluxdist.mean() / mu_by_eye)
                    e_flux_indiv_int.append(fluxdist.std() / mu_by_eye)
                    if tag2 in fitresults.keys():
                        if fitresultsgood[tag2].mean() == 100:
                            continue
                        observedflux = fluxdist * fitresultsgood[tag2]
                        print(datanamesource,
                                fitresultsgood[tag2].mean(),
                                fitresultsgood[tag2].std())
                    else:
                        observedflux = fluxdist
                        print(datanamesource, 1.0)
                    flux_target += observedflux
                    flux_indiv.append(observedflux.mean())
                    e_flux_indiv.append(observedflux.std())
                    dataname_indiv.append(datanamesource)
                
            flux_measure = flux_target.mean()
            flux_total.append(flux_measure)
            flux_error = flux_target.std()
            e_flux_total.append(flux_error)
            #import matplotlib.pyplot as plt
            #plt.hist(flux_target)
            #plt.show()
            #import pdb; pdb.set_trace()
            #print(objname, flux_measure, ' +/- ', flux_error)

        flux_total_observed = {'targets': dataname, 'fnu': flux_total, 
                'e_fnu': e_flux_total}
        t1 = Table(flux_total_observed)
        t1.write('flux_total_observed.dat', format='ascii')
        flux_indiv_observed = {'targets': dataname_indiv, 'fnu': flux_indiv, 
                'e_fnu': e_flux_indiv}
        t1 = Table(flux_indiv_observed)
        t1.write('flux_indiv_observed.dat', format='ascii')
        flux_indiv_intrinsic = {'targets': dataname_indiv, 'fnu': flux_indiv_int, 
                'e_fnu': e_flux_indiv_int}
        t1 = Table(flux_indiv_intrinsic)
        t1.write('flux_indiv_intrinsic.dat', format='ascii')

        #import matplotlib.pyplot as plt
        #m1 = 0
        #m2 = 1 + int(numpy.array(flux_total + flux_indiv + flux_indiv_int).max())
        #binwidth = 2
        #bins = numpy.arange(m1, m2 + binwidth, binwidth)
        #plt.hist(flux_total, bins = bins, histtype='stepfilled')
        #binwidth = 1.8
        #bins = numpy.arange(m1, m2 + binwidth, binwidth)
        #plt.hist(flux_indiv, bins = bins, alpha=0.5, histtype='stepfilled')
        #binwidth = 1.6
        #bins = numpy.arange(m1, m2 + binwidth, binwidth)
        #plt.hist(flux_indiv_int, bins = bins, alpha=0.2, histtype='stepfilled')
        #plt.show()
    return flux_total_observed, flux_indiv_observed, flux_indiv_intrinsic
