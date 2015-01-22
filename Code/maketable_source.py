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
from astropy.io import ascii

goodfitloc = '../Data/uvfitlist.dat'
goodfitdat = Table.read(goodfitloc, format='ascii')
goodfitdat.sort('iauname')
mu_estimate = False

cwd = os.getcwd()

observed = ascii.read('../Data/table_observed.dat')

# read in the sample targets
anloc = '../../../../ModelFits/'
dataphotloc = '../Data/targetlist.dat'
dataphot = Table.read(dataphotloc, format='ascii')
dataphot.sort('ra_250')
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
shortname = dataphot['shortname']
ntarg = len(dataphot)
objname_previous = ''

sourcecounter = 0

for itarg in range(0, ntarg):

    # extract object name from directory
    objname = dataname[itarg]
    iauname = dataphot['iauname'][itarg]
    goodfit = goodfitdat['intrinsic'][itarg]
    #goodfit = goodfitdat['observed'][itarg]
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
    nfits = len(fitresultsgood)

    configfile = open('config.yaml')
    config = yaml.load(configfile)
    paramData = setuputil.loadParams(config)

    if objname != objname_previous:
        radeg_centroid = config['Region0']['RACentroid']
        decdeg_centroid = config['Region0']['DecCentroid']
        c = SkyCoord(ra=radeg_centroid*u.degree, 
                dec=decdeg_centroid*u.degree)
        fullRA = c.ra.to_string(unit='hour', pad=True, sep=':', precision=3)
        fullDec = c.dec.to_string('deg', pad=True, alwayssign=True, sep=':', precision=2)
        msg = '${0:9}$ ${1:12}$ ${2:12}$'
        msgfmt = msg.format(objname, fullRA, fullDec)
        #print(msgfmt)
    objname_previous = objname

    nregions = paramData['nregions']

    # RA and Dec info
    pname = paramData['pname']
    poff = paramData['poff']

    # RA and Dec of primary lens (default is no lens)
    RA0 = 0.0
    Dec0 = 0.0

    flux_target = 0

    os.chdir(cwd)
    for iregion in range(nregions):
        sr = str(iregion)
        region = 'Region' + sr
        nsource = paramData['nsource_regions'][iregion]
        nlens = paramData['nlens_regions'][iregion]
        radeg_centroid = config[region]['RACentroid']
        decdeg_centroid = config[region]['DecCentroid']
        for ilens in range(nlens):

            # Lens RA
            sl = str(ilens)
            datanamelens = shortname[itarg] + '.Lens' + sl
            tag1 = 'DeltaRA_Lens' + sl + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Lens' + sl  + ' DeltaRA'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_lens = fitresultsgood[tag1]
            else:
                # this parameter was fixed
                parname = tag1.split(" ")
                if len(parname) == 1:
                    parname = tag1.split("_")
                    parname.reverse()
                crlp = config[parname[0]][parname[1]][parname[2]]
                dist_lens = numpy.zeros(nfits) + crlp['Limits'][1]
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Lens' + sl  + ' DeltaRA':
                    tag1 = 'DeltaRA_Lens' + sl + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Lens' + sl  + ' DeltaRA'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_lens = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_lens = fitresultsgood[fixtag]
                else:
                    # the reference parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_lens = crlp['Limits'][1]
            dist = dist_lens + off_lens
            dradeg_lens = dist.mean() / 3600
            eRA = dist.std() / 15
            if ilens == 0:
                RA0 = dist

            # Lens Dec
            tag1 = 'DeltaDec_Lens' + sl + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Lens' + sl  + ' DeltaDec'
                if tag1 not in fitresults.keys():
                    thisparam = 0
            if thisparam == 1:
                dist_lens = fitresultsgood[tag1]
            else:
                # this parameter was fixed
                parname = tag1.split(" ")
                if len(parname) == 1:
                    parname = tag1.split("_")
                    parname.reverse()
                crlp = config[parname[0]][parname[1]][parname[2]]
                dist_lens = numpy.zeros(nfits) + crlp['Limits'][1]
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Lens' + sl  + ' DeltaDec':
                    tag1 = 'DeltaDec_Lens' + sl + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Lens' + sl  + ' DeltaDec'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_lens = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_lens = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_lens = crlp['Limits'][1]
            dist = dist_lens + off_lens
            ddecdeg_lens = dist.mean() / 3600
            eDec = dist.std()
            if ilens == 0:
                Dec0 = dist

            radeg = radeg_centroid + dradeg_lens
            decdeg = decdeg_centroid + ddecdeg_lens
            c = SkyCoord(ra=radeg*u.degree, dec=decdeg*u.degree)
            fullRA = c.ra.to_string(unit='hour', pad=True, sep=':', precision=3)
            fullDec = c.dec.to_string('deg', pad=True, alwayssign=True, sep=':', precision=2)
            #fullRA = fullRA + '(' + eRA + ')'
            #fullDec = fullRA + '(' + eDec + ')'

            # Einstein radius
            tag1 = 'EinsteinRadius_Lens' + sl + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Lens' + sl  + ' EinsteinRadius'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_lens = fitresultsgood[tag1]
            else:
                # this parameter was fixed
                parname = tag1.split(" ")
                if len(parname) == 1:
                    parname = tag1.split("_")
                    parname.reverse()
                crlp = config[parname[0]][parname[1]][parname[2]]
                dist_lens = numpy.zeros(nfits) + crlp['Limits'][1]
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Lens' + sl  + ' EinsteinRadius':
                    tag1 = 'EinsteinRadius_Lens' + sl + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Lens' + sl  + ' EinsteinRadius'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_lens = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_lens = fitresultsgood[fixtag]
                else:
                    # the reference parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_lens = crlp['Limits'][1]
            dist = dist_lens + off_lens
            einstein_indiv = dist.mean()
            e_einstein_indiv = dist.std()

            # Axial Ratio
            tag1 = 'AxialRatio_Lens' + sl + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Lens' + sl  + ' AxialRatio'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_lens = fitresultsgood[tag1]
            else:
                # this parameter was fixed
                parname = tag1.split(" ")
                if len(parname) == 1:
                    parname = tag1.split("_")
                    parname.reverse()
                crlp = config[parname[0]][parname[1]][parname[2]]
                dist_lens = numpy.zeros(nfits) + crlp['Limits'][1]
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Lens' + sl  + ' AxialRatio':
                    tag1 = 'AxialRatio_Lens' + sl + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Lens' + sl  + ' AxialRatio'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_lens = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_lens = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_lens = crlp['Limits'][1]
            dist = dist_lens + off_lens
            q_indiv = dist.mean()
            e_q_indiv = dist.std()

            # Position Angle
            tag1 = 'PositionAngle_Lens' + sl + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Lens' + sl  + ' PositionAngle'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_lens = fitresultsgood[tag1]
            else:
                # this parameter was fixed
                parname = tag1.split(" ")
                if len(parname) == 1:
                    parname = tag1.split("_")
                    parname.reverse()
                crlp = config[parname[0]][parname[1]][parname[2]]
                dist_lens = numpy.zeros(nfits) + crlp['Limits'][1]
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Lens' + sl  + ' PositionAngle':
                    tag1 = 'PositionAngle_Lens' + sl + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Lens' + sl  + ' PositionAngle'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_lens = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_lens = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_lens = crlp['Limits'][1]
            dist = dist_lens + off_lens
            phi_indiv = dist.mean()
            e_phi_indiv = dist.std()

            msg = '{0:15} & ${1:12}\pm{2:5.3f}$ & ${3:12}\pm{4:4.2f}$ & ${5:.3f}\pm{6:.3f}$ & ${7:.3f}\pm{8:0.3f}$ & ${9:3.0f}\pm{10:3.0f}$'
            msgfmt = msg.format(datanamelens, fullRA, eRA, fullDec, eDec, 
                    einstein_indiv, e_einstein_indiv, q_indiv, e_q_indiv,
                    phi_indiv, e_phi_indiv)
            #msg = '{0:3.0f} {1:3.0f}'
            #msgfmt = msg.format(phi_indiv, e_phi_indiv)            
            #msg = '{0:.2f} {1:.2f}'
            #msgfmt = msg.format(mu_indiv, e_mu_indiv)            
            #print(msgfmt)

        for isource in range(nsource):

            pbcorr = observed['pbcorr'][sourcecounter]
            remu = observed['remu'][sourcecounter]
            e_remu = observed['e_remu'][sourcecounter]
            RA1 = observed['ra_alma'][sourcecounter]
            Dec1 = observed['dec_alma'][sourcecounter]
            fullRADec = RA1 + ' ' + Dec1
            sourcecounter += 1

            ss = str(isource)
            datanamesource = shortname[itarg] + '.' + ss
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
            if tag2 in fitresults.keys():
                if fitresultsgood[tag2].mean() == 100:
                    continue
                observedflux = fluxdist / pbcorr * fitresultsgood[tag2] / remu
                mu_indiv = remu#fitresultsgood[tag2].mean()
                e_mu_indiv = e_remu#fitresultsgood[tag2].std()
            else:
                observedflux = fluxdist / pbcorr
                mu_indiv = 1.0
                e_mu_indiv = 0.0
            flux_target += observedflux
            flux_indiv = observedflux.mean()
            e_flux_indiv = observedflux.std()

            # RA
            tag1 = 'DeltaRA_Source' + ss + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Source' + ss  + ' DeltaRA'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_source = fitresultsgood[tag1]
            else:
                dist_source = 0.0
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Source' + ss  + ' DeltaRA':
                    tag1 = 'DeltaRA_Source' + ss + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Source' + ss  + ' DeltaRA'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_source = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_source = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_source = crlp['Limits'][1]
            dist = dist_source + off_source - RA0
            dradeg_source = dist.mean()
            e_dradeg_source = dist.std()

            # Dec
            tag1 = 'DeltaDec_Source' + ss + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Source' + ss  + ' DeltaDec'
                if tag1 not in fitresults.keys():
                    thisparam = 0
            if thisparam == 1:
                dist_source = fitresultsgood[tag1]
            else:
                dist_source = 0.0
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Source' + ss  + ' DeltaDec':
                    tag1 = 'DeltaDec_Source' + ss + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Source' + ss  + ' DeltaDec'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_source = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_source = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_source = crlp['Limits'][1]
            dist = dist_source + off_source - Dec0
            ddecdeg_source = dist.mean()
            e_ddecdeg_source = dist.std()

            radeg = radeg_centroid + (dradeg_source) * \
                    numpy.cos(decdeg_centroid) / 3600
            decdeg = decdeg_centroid + (ddecdeg_source) / 3600
            c = SkyCoord(ra=radeg*u.degree, dec=decdeg*u.degree)
            #fullRADec = c.to_string('hmsdms', sep=':')

            # Effective Radius
            tag1 = 'Size_Source' + ss + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Source' + ss  + ' EffectiveRadius'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_source = fitresultsgood[tag1]
            else:
                dist_source = 0.0
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Source' + ss  + ' EffectiveRadius':
                    tag1 = 'Size_Source' + ss + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Source' + ss  + ' EffectiveRadius'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_source = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_source = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_source = crlp['Limits'][1]
            dist = dist_source + off_source
            size_indiv = dist.mean()
            e_size_indiv = dist.std()

            # Axial Ratio
            tag1 = 'AxialRatio_Source' + ss + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Source' + ss  + ' AxialRatio'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_source = fitresultsgood[tag1]
            else:
                dist_source = 0.0
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Source' + ss  + ' AxialRatio':
                    tag1 = 'AxialRatio_Source' + ss + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Source' + ss  + ' AxialRatio'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_source = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_source = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_source = crlp['Limits'][1]
            dist = dist_source + off_source
            q_indiv = dist.mean()
            e_q_indiv = dist.std()

            # Position Angle
            tag1 = 'PositionAngle_Source' + ss + '_Region' + sr
            thisparam = 1
            if tag1 not in fitresults.keys():
                tag1 = 'Region' + sr + ' Source' + ss  + ' PositionAngle'
                if tag1 not in fitresults.keys():
                    thisparam = 0
                
            if thisparam == 1:
                dist_source = fitresultsgood[tag1]
            else:
                dist_source = 0.0
            count = pname.count(tag1)
            if count == 0:
                if tag1 == 'Region' + sr + ' Source' + ss  + ' PositionAngle':
                    tag1 = 'AxialRatio_Source' + ss + '_Region' + sr
                else:
                    tag1 = 'Region' + sr + ' Source' + ss  + ' PositionAngle'
            indx = pname.index(tag1)
            fixtag = poff[indx]
            if fixtag == 'Free':
                off_source = 0.
            else:
                count = fitresultsgood.keys().count(fixtag)
                if count != 0:
                    off_source = fitresultsgood[fixtag]
                else:
                    # the foundational parameter was fixed
                    parname = fixtag.split(" ")
                    if len(parname) == 1:
                        parname = fixtag.split("_")
                        parname.reverse()
                    crlp = config[parname[0]][parname[1]][parname[2]]
                    off_source = crlp['Limits'][1]
            dist = dist_source + off_source
            phi_indiv = dist.mean()
            e_phi_indiv = dist.std()

            msg = '{0:15} ${1:6.3f}\pm{2:5.3f}$ & ${3:6.3f}\pm{4:5.3f}$ & ${5:5.2f}\pm{6:5.2f}$ & ${7:5.3f}\pm{8:5.3f}$ & ${9:5.2f}\pm{10:5.2f}$ & ${11:3.0f}\pm{12:3.0f}$ & ${13:5.2f}\pm{14:5.2f}$ \\\\'
            #msg = '{0:15} {1:6.3f} {2:5.3f} {3:6.3f} {4:5.3f}  {5:5.2f} {6:5.2f}  {7:5.3f} {8:5.3f}  {9:5.2f} {10:5.2f}  {11:3.0f} {12:3.0f}  {13:5.2f} {14:5.2f}'
            msgfmt = msg.format(datanamesource, dradeg_source, e_dradeg_source,
                ddecdeg_source, e_ddecdeg_source, 
                flux_indiv, e_flux_indiv, size_indiv, e_size_indiv, q_indiv, 
                e_q_indiv, phi_indiv, e_phi_indiv, mu_indiv, e_mu_indiv)
            #msg = '{0:3.0f} {1:3.0f}'
            #msgfmt = msg.format(phi_indiv, e_phi_indiv)            
            #msg = '{0:.2f} {1:.2f}'
            #msgfmt = msg.format(mu_indiv, e_mu_indiv)            
            #msg = '{0:15} {1:29} {2:.2f} {3:.2f} {4:.2f} {5:.2f} {6:.2f}'
            #msgfmt = msg.format(datanamesource, fullRADec, flux_indiv, e_flux_indiv, pbcorr, remu, e_remu)
            print(msgfmt)
        
    flux_measure = flux_target.mean()
    flux_error = flux_target.std()
    #import matplotlib.pyplot as plt
    #plt.hist(flux_target)
    #plt.show()
    #import pdb; pdb.set_trace()
    #radeg = radeg_centroid
    #decdeg = decdeg_centroid
    #c = SkyCoord(ra=radeg*u.degree, dec=decdeg*u.degree)
    #fullRADec = c.to_string('hmsdms', sep=':')
    msg = '{0:20} {1:.2f} {2:.2f}'
    msgfmt = msg.format(objname, flux_measure, flux_error)
    #print(msgfmt)

