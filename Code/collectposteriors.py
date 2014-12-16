"""
 Plot posterior distributions of model parameters
 Shane Bussmann
 2014 July 14

"""

#def plot(objname, redo, deltachi2):

import math
import numpy
import os
from astropy.table import Table
#import aplpy
import matplotlib.pyplot as mpl
from pylab import savefig
#from matplotlib import rcParams.update({'font.size': 8})
from matplotlib import rc
import pdb

rc('font',**{'family':'sans-serif','sans-serif':['Arial'],'size':'6'})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['New Century Schoolbook']})
#rc('text', usetex=True)

#deltachi2 = 100

# read in the sample targets
anloc = '../../ModelFits/'
dataphotloc = 'targetlist.dat'
dataphot = Table.read(dataphotloc, format='ascii')
dataname = dataphot['dataname']
ntarg = dataname.size
dataflag = numpy.zeros(ntarg)
for i in numpy.arange(ntarg):
    datanamei = dataname[i]
    dataloc = anloc + datanamei + '/uvfit25'
    if os.path.exists(dataloc):
        dataflag[i] = 1

yesmodel = dataflag == 1
dataphot = dataphot[yesmodel]
dataname = dataphot['dataname']
nmodel = dataphot.size

# examine each inputs_lens.dat file to determine total number of lenses
nsources = 0
for i in range(nmodel):

    # Grab constraints on lens(es)
    dir = anloc + dataname[i] + '/uvmcbest/'
    bestfitloc = dir + 'inputs_source.dat'
    fitinputs_source = Table.read(bestfitloc)
    nsource = len(fitinputs_source['flux']) / 2
    nsources = nsources + nsource

# store each source separately
dt = {'names':['iauname', 'ra_s', 'e_ra_s', 'dec_s', 'e_dec_s', 
    'flux_s', 'e_flux_s', 'n_s', 'e_n_s', 'r_s', 'e_r_s', 'ell_s', 'e_ell_s',
    'phi_s', 'e_phi_s', 'mu', 'e_mu', 'minchi2', 'rmschi2', 'ndof'],
    'formats':['a30', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 
    'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 'f4', 'f4', 'f4']}
modfits_source = numpy.zeros((nsources,), dtype=(dt))

# examine each inputs_lens.dat file to determine total number of lenses
nlenses = 0
for i in range(nmodel):

    # Grab constraints on lens(es)
    dir = anloc + dataname[i] + '/uvmcbest/'
    bestfitloc = dir + 'inputs_lens.dat'
    fitinputs_lens = asciitable.read(bestfitloc)
    nlens = len(fitinputs_lens['theta_E']) / 2
    nlenses = nlenses + nlens

# store each lens separately
dt = {'names':['iauname', 'ra_l', 'e_ra_l', 'dec_l', 'e_dec_l', 'theta_E',
    'e_theta_E', 'ell_l', 'e_ell_l', 'phi_l', 'e_phi_l', 'minchi2', 'rmschi2', 'ndof'],
    'formats':['a30', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2', 'f2',
    'f2', 'f2', 'f2', 'f4', 'f4', 'f4']}
modfits_lens = numpy.zeros((nlenses,), dtype=(dt))

counter_l = 0
counter_s = 0
for iii in numpy.arange(nmodel):

    dir = anloc + dataname[iii] + '/uvmcbest/'

    # extract object name from directory
    objname = dataname[iii]
    iauname = dataphot['iauname'][iii]
    burninloc = dir + 'finalposterior.dat'
    muloc = dir + 'mu.dat'

    #----------------------------------------------------------------------------
    # Grab constraints on lens(es)
    bestfitloc = dir + 'inputs_lens.dat'
    fitinputs_lens = asciitable.read(bestfitloc)
    nlens = len(fitinputs_lens['theta_E']) / 2
    evenindx = numpy.arange(nlens) * 2
    oddindx = evenindx + 1
    p_u = fitinputs_lens['theta_E'][evenindx]
    p_l = fitinputs_lens['theta_E'][oddindx]
    p_u = numpy.append(p_u, fitinputs_lens['x_l'][0])
    p_l = numpy.append(p_l, fitinputs_lens['x_l'][1])
    p_u = numpy.append(p_u, fitinputs_lens['y_l'][0])
    p_l = numpy.append(p_l, fitinputs_lens['y_l'][1])
    p_u = numpy.append(p_u, fitinputs_lens['ell_l'][evenindx])
    p_l = numpy.append(p_l, fitinputs_lens['ell_l'][oddindx])
    p_u = numpy.append(p_u, fitinputs_lens['phi_l'][evenindx])
    p_l = numpy.append(p_l, fitinputs_lens['phi_l'][oddindx])

    # get offset of other lenses from primary lens
    x_l_off = fitinputs_lens['x_l'][evenindx]
    y_l_off = fitinputs_lens['y_l'][evenindx]
    x_l_off[0] = 0.0
    y_l_off[0] = 0.0

    # don't use shear in fits unless it's listed in constraints
    yes_shear = False
    if fitinputs_lens['shear'].sum() != 0.0: yes_shear = True

    print "Reading output from emcee"

    #------------------------------------------------------------------------------
    # Read posterior probability distributions
    muresults = asciitable.read(muloc)
    fitresults = asciitable.read(burninloc)
    print objname, 'prior to pruning: ', fitresults['chi2'].mean()

    # identify the good fits
    chi2 = fitresults['chi2']
    fitresults['chi2'] = chi2 / 2.
    minchi2 = chi2.max()
    dchi2 = numpy.abs(chi2 - minchi2)
    medchi2 = numpy.median(dchi2)
    avgchi2 = numpy.mean(dchi2)
    skewchi2 = numpy.abs(avgchi2 - medchi2)
    rmschi2 = numpy.std(dchi2)
    scaler = 1.5
    if objname == 'NA.v1.177':
        scaler = 0.5
    if objname == 'ID081':
        scaler = 0.5
    if objname == 'ID130':
        scaler = 0.5
    if objname == 'NA.v1.56':
        scaler = 1.0
    if objname == 'NC.v1.268':
        scaler = 0.5
    inliers = (dchi2 < scaler*rmschi2)
    fitresults = fitresults[inliers]
    muresults = muresults[inliers]
    while skewchi2 > 0.1*medchi2:
        chi2 = fitresults['chi2']
        minchi2 = chi2.max()
        dchi2 = numpy.abs(chi2 - minchi2)
        rmschi2 = numpy.std(dchi2)
        scaler = 1.5
        if objname == 'NA.v1.177':
            scaler = 0.5
        if objname == 'NC.v1.268':
            scaler = 0.5
        inliers = (dchi2 < scaler*rmschi2)
        fitresultstmp = fitresults[inliers]
        if fitresultstmp.size == fitresults.size:
            inliers = (dchi2 < scaler/2.*rmschi2)
        fitresults = fitresults[inliers]
        muresults = muresults[inliers]
        chi2 = fitresults['chi2']
        dchi2 = numpy.abs(chi2 - minchi2)
        medchi2 = numpy.median(dchi2)
        avgchi2 = numpy.mean(dchi2)
        skewchi2 = numpy.abs(avgchi2 - medchi2)
        print medchi2, avgchi2, skewchi2

    fitresultsgood = fitresults.copy()
    # identify the good fits
    goodfits = (fitresults['chi2'] <= minchi2)
    fitresultsgood = fitresults[goodfits].copy()

    print "Making histograms for the lens(es)"
    ndim = 5 + (nlens - 1) * 3 + 1
    nrow = numpy.ceil(math.sqrt(ndim))
    ncol = numpy.ceil(math.sqrt(ndim))
    j = 1
    ax = mpl.subplot(nrow, ncol, j)
    #good = muresults['chi2'] > minchi2 - 150

    mpl.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95,
        wspace=0.4, hspace=0.55)

    # Chi-squared: identify best-fit model and range of best-fit models
    minchi2 = fitresultsgood['chi2'].max()
    avgchi2 = fitresultsgood['chi2'].mean()
    rmschi2 = fitresultsgood['chi2'].std()
    indx = fitresultsgood['chi2'] == minchi2
    bestfit = fitresultsgood[indx]
    #print minchi2, avgchi2, rmschi2, deltachi2
    mpl.clf()

    # Cycle through the lens(es)
    for i in numpy.arange(nlens)+1:

        # Position of primary lens
        p_u = fitinputs_lens['x_l'][2 * (i - 1)]
        p_l = fitinputs_lens['x_l'][2 * (i - 1) + 1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['x_l'])
            avgval = -1 * numpy.mean(fitresultsgood['x_l'])
            modfits_lens[counter_l]['ra_l'] = avgval
            modfits_lens[counter_l]['e_ra_l'] = rmsval
            print 'x_l = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['x_l']
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('x_l')
        p_u = fitinputs_lens['y_l'][2 * (i - 1)]
        p_l = fitinputs_lens['y_l'][2 * (i - 1) + 1]
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['y_l'])
            avgval = numpy.mean(fitresultsgood['y_l'])
            modfits_lens[counter_l]['dec_l'] = avgval
            modfits_lens[counter_l]['e_dec_l'] = rmsval
            print 'y_l'+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['y_l']
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('y_l')

        p_u = fitinputs_lens['theta_E'][2*(i-1)]
        p_l = fitinputs_lens['theta_E'][2*(i-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['theta_E'+str(i)])
            avgval = numpy.mean(fitresultsgood['theta_E'+str(i)])
            modfits_lens[counter_l]['theta_E'] = avgval
            modfits_lens[counter_l]['e_theta_E'] = rmsval
            print 'theta_E'+str(i)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['theta_E'+str(i)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('theta_E'+str(i))
        p_u = fitinputs_lens['ell_l'][2*(i-1)]
        p_l = fitinputs_lens['ell_l'][2*(i-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            tobeplot = fitresultsgood['ell_l'+str(i)]
            rmsval = numpy.std(tobeplot)
            avgval = numpy.mean(tobeplot)
            modfits_lens[counter_l]['ell_l'] = avgval
            modfits_lens[counter_l]['e_ell_l'] = rmsval
            print 'ell_l'+str(i)+' = ', avgval, ' +/- ', rmsval
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('ell_l'+str(i))
        p_u = fitinputs_lens['phi_l'][2*(i-1)]
        p_l = fitinputs_lens['phi_l'][2*(i-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            tobeplot = fitresultsgood['phi_l'+str(i)]
            rmsval = numpy.std(tobeplot)
            avgval = numpy.mean(tobeplot)
            modfits_lens[counter_l]['phi_l'] = avgval
            modfits_lens[counter_l]['e_phi_l'] = rmsval
            # if phi=0 and phi=180 are the best-fit value, re-center on 0
            ind1 = tobeplot < 90
            ind2 = tobeplot > 90
            rms1 = numpy.std(tobeplot[ind1])
            rms2 = numpy.std(tobeplot[ind2])
            rmsfull = numpy.std(tobeplot)
            if rms1+rms2 < rmsfull:
                plus90 = tobeplot > 90
                tobeplot[plus90] = tobeplot[plus90] - 180.
                rmsval = numpy.std(tobeplot)
                avgval = numpy.mean(tobeplot)
            print 'phi_l'+str(i)+' = ', avgval, ' +/- ', rmsval
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('phi_l'+str(i))
        if yes_shear:
            p_u = fitinputs_lens['shear'][2*(i-1)]
            p_l = fitinputs_lens['shear'][2*(i-1)+1]
            totalwidth = p_u - p_l
            if totalwidth > 0:
                rmsval = numpy.std(fitresultsgood['shear'+str(i)])
                avgval = numpy.mean(fitresultsgood['shear'+str(i)])
                print 'shear'+str(i)+' = ', avgval, ' +/- ', rmsval
                tobeplot = fitresultsgood['shear'+str(i)]
                nbins = totalwidth / rmsval * 5
                ax = mpl.subplot(nrow, ncol, j)
                j += 1
                mpl.hist(tobeplot, nbins, edgecolor='blue')
                mpl.ylabel('n')
                mpl.xlabel('shear'+str(i))
            p_u = fitinputs_lens['phi_shear'][2*(i-1)]
            p_l = fitinputs_lens['phi_shear'][2*(i-1)+1]
            totalwidth = p_u - p_l
            if totalwidth > 0:
                rmsval = numpy.std(fitresultsgood['phi_shear'+str(i)])
                avgval = numpy.mean(fitresultsgood['phi_shear'+str(i)])
                print 'phi_shear'+str(i)+' = ', avgval, ' +/- ', rmsval
                tobeplot = fitresultsgood['phi_shear'+str(i)]
                nbins = totalwidth / rmsval * 5
                ax = mpl.subplot(nrow, ncol, j)
                j += 1
                mpl.hist(tobeplot, nbins, edgecolor='blue')
                mpl.ylabel('n')
                mpl.xlabel('phi_shear'+str(i))

        # store chi^2 data
        modfits_lens[counter_l]['minchi2'] = -1 * minchi2
        modfits_lens[counter_l]['rmschi2'] = rmschi2
        modfits_lens[counter_l]['iauname'] = iauname
        counter_l += 1
    mpl.suptitle(iauname)

    # plot the histogram of chi^2 values
    ax = mpl.subplot(nrow, ncol, j)
    mpl.hist(fitresultsgood['chi2'].max() - fitresultsgood['chi2'], bins=100)
    savefig(dir + objname + '_posterior_pdfs_lens.png')
    print objname, minchi2, avgchi2, rmschi2#, deltachi2

    mpl.clf()

    #----------------------------------------------------------------------------
    # Grab constraints on source(es)
    inputloc = dir + 'inputs_source.dat'
    table = asciitable.get_reader(Reader=asciitable.Basic)
    fitinputs_source = table.read(inputloc)
    nsource = len(fitinputs_source['flux']) / 2
    evenindx = numpy.arange(nsource) * 2
    oddindx = evenindx + 1
    p_u = numpy.append(p_u, fitinputs_source['flux'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['flux'][oddindx])
    p_u = numpy.append(p_u, fitinputs_source['x_s'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['x_s'][oddindx])
    p_u = numpy.append(p_u, fitinputs_source['y_s'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['y_s'][oddindx])
    p_u = numpy.append(p_u, fitinputs_source['n_s'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['n_s'][oddindx])
    p_u = numpy.append(p_u, fitinputs_source['r_s'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['r_s'][oddindx])
    p_u = numpy.append(p_u, fitinputs_source['ell_s'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['ell_s'][oddindx])
    p_u = numpy.append(p_u, fitinputs_source['phi_s'][evenindx])
    p_l = numpy.append(p_l, fitinputs_source['phi_s'][oddindx])

    print "Making histograms for the source(s)"
    ndim = nsource * 8 + 2
    nrow = numpy.ceil(math.sqrt(ndim))
    ncol = numpy.ceil(math.sqrt(ndim))
    j = 1
    ax = mpl.subplot(nrow, ncol, j)
    #good = muresults['chi2'] > minchi2 - 150
    mu = muresults['mu_aper']
    rmsval = numpy.std(mu[goodfits])
    tobeplot = mu[goodfits]
    mpl.hist(tobeplot, 100, edgecolor='blue')
    modfits_source[counter_s]['mu'] = tobeplot.mean()
    modfits_source[counter_s]['e_mu'] = rmsval
    j = j + 1
    print 'mu_aper ', tobeplot.mean(), ' +/- ', rmsval

    #for i in numpy.arange(nsource)+1:
    #    mu = muresults['mu'+str(i)]
    #    rmsval = numpy.std(mu[goodfits])
    #    tobeplot = mu[goodfits]
    #    mpl.hist(tobeplot, 100, edgecolor='green')
    #    print 'mu'+str(i) + ' = ', tobeplot.mean(), ' +/- ', rmsval

    mpl.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95,
        wspace=0.4, hspace=0.55)

    # Chi-squared: identify best-fit model and range of best-fit models
    minchi2 = fitresultsgood['chi2'].max()
    avgchi2 = fitresultsgood['chi2'].mean()
    rmschi2 = fitresultsgood['chi2'].std()
    indx = fitresultsgood['chi2'] == minchi2
    bestfit = fitresultsgood[indx]
    #print minchi2, avgchi2, rmschi2, deltachi2

    # Cycle through the source(s)
    for si in numpy.arange(nsource)+1:
        #si = 1
        mu = muresults['mu_aper'+str(si)]
        rmsval = numpy.std(mu[goodfits])
        avgval = numpy.mean(mu[goodfits])
        tobeplot = mu[goodfits]
        ax = mpl.subplot(nrow, ncol, j)
        j += 1
        mpl.hist(tobeplot, 100, edgecolor='green')
        mpl.ylabel('N')
        mpl.xlabel('mu_aper'+str(si))
        print 'mu_aper'+str(si) + ' = ', avgval, ' +/- ', rmsval
        modfits_source[counter_s]['mu'] = avgval
        modfits_source[counter_s]['e_mu'] = rmsval
        p_u = fitinputs_source['flux'][2*(si-1)]
        p_l = fitinputs_source['flux'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['flux'+str(si)])
            avgval = numpy.mean(fitresultsgood['flux'+str(si)])
            modfits_source[counter_s]['flux_s'] = avgval
            modfits_source[counter_s]['e_flux_s'] = rmsval
            print 'flux'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['flux'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('flux'+str(si))
        p_u = fitinputs_source['x_s'][2*(si-1)]
        p_l = fitinputs_source['x_s'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['x_s'+str(si)])
            avgval = -1 * numpy.mean(fitresultsgood['x_s'+str(si)])
            modfits_source[counter_s]['ra_s'] = avgval
            modfits_source[counter_s]['e_ra_s'] = rmsval
            print 'x_s'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['x_s'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('x_s'+str(si))
        p_u = fitinputs_source['y_s'][2*(si-1)]
        p_l = fitinputs_source['y_s'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['y_s'+str(si)])
            avgval = numpy.mean(fitresultsgood['y_s'+str(si)])
            modfits_source[counter_s]['dec_s'] = avgval
            modfits_source[counter_s]['e_dec_s'] = rmsval
            print 'y_s'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['y_s'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('y_s'+str(si))
        p_u = fitinputs_source['n_s'][2*(si-1)]
        p_l = fitinputs_source['n_s'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['n_s'+str(si)])
            avgval = numpy.mean(fitresultsgood['n_s'+str(si)])
            modfits_source[counter_s]['n_s'] = avgval
            modfits_source[counter_s]['e_n_s'] = rmsval
            print 'n_s'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['n_s'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('n_s'+str(si))
        p_u = fitinputs_source['r_s'][2*(si-1)]
        p_l = fitinputs_source['r_s'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['r_s'+str(si)])
            avgval = numpy.mean(fitresultsgood['r_s'+str(si)])
            modfits_source[counter_s]['r_s'] = avgval
            modfits_source[counter_s]['e_r_s'] = rmsval
            print 'r_s'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['r_s'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('r_s'+str(si))
        p_u = fitinputs_source['ell_s'][2*(si-1)]
        p_l = fitinputs_source['ell_s'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['ell_s'+str(si)])
            avgval = numpy.mean(fitresultsgood['ell_s'+str(si)])
            modfits_source[counter_s]['ell_s'] = avgval
            modfits_source[counter_s]['e_ell_s'] = rmsval
            print 'ell_s'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['ell_s'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('ell_s'+str(si))
        p_u = fitinputs_source['phi_s'][2*(si-1)]
        p_l = fitinputs_source['phi_s'][2*(si-1)+1]
        totalwidth = p_u - p_l
        if totalwidth > 0:
            rmsval = numpy.std(fitresultsgood['phi_s'+str(si)])
            avgval = numpy.mean(fitresultsgood['phi_s'+str(si)])
            modfits_source[counter_s]['phi_s'] = avgval
            modfits_source[counter_s]['e_phi_s'] = rmsval
            print 'phi_s'+str(si)+' = ', avgval, ' +/- ', rmsval
            tobeplot = fitresultsgood['phi_s'+str(si)]
            nbins = totalwidth / rmsval * 5
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(tobeplot, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel('phi_s'+str(si))

        # store chi^2 data
        modfits_source[counter_s]['minchi2'] = -1 * minchi2
        modfits_source[counter_s]['rmschi2'] = rmschi2
        modfits_source[counter_s]['iauname'] = iauname
        counter_s += 1
    mpl.suptitle(iauname)

    # plot the histogram of chi^2 values
    ax = mpl.subplot(nrow, ncol, j)
    mpl.hist(fitresultsgood['chi2'].max() - fitresultsgood['chi2'], bins=100)
    savefig(dir + objname + '_posterior_pdfs_source.png')
    print objname, minchi2, avgchi2, rmschi2#, deltachi2


    # identify best-fit model and range of best-fit models
    #minchi2 = fitresultsgood['col1'].max()
    #avgchi2 = fitresultsgood['col1'].mean()
    #rmschi2 = fitresultsgood['col1'].std()
    #indx = fitresultsgood['col1'] == minchi2
    #bestfit = fitresultsgood[indx]

    #print objname, minchi2, avgchi2, rmschi2#, deltachi2


asciitable.write(modfits_lens, output='table3_lenses.dat', formats={'ra_l': '%0.3f', 'e_ra_l': '%0.3f', 'dec_l': '%0.3f', 'e_dec_l': '%0.3f', 'theta_E': '%0.3f', 'e_theta_E': '%0.3f', 'ell_l': '%0.2f', 'e_ell_l': '%0.2f', 'phi_l': '%0.0f', 'e_phi_l': '%0.0f', 'minchi2': '%0.1f', 'rmschi2': '%0.1f', 'ndof': '%0.0f'})

asciitable.write(modfits_source, output='table3_sources.dat', formats={'ra_s': '%0.3f', 'e_ra_s': '%0.3f', 'dec_s': '%0.3f', 'e_dec_s': '%0.3f', 'flux_s': '%0.3f', 'e_flux_s': '%0.3f', 'n_s': '%0.1f', 'e_n_s': '%0.1f', 'r_s': '%0.2f', 'e_r_s': '%0.2f', 'ell_s': '%0.2f', 'e_ell_s': '%0.2f', 'phi_s': '%0.0f', 'e_phi_s': '%0.0f', 'mu': '%0.1f', 'e_mu': '%0.1f', 'minchi2': '%0.1f', 'rmschi2': '%0.1f', 'ndof': '%0.0f'})

pdb.set_trace()
