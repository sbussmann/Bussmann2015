"""
Created by Shane Bussmann
2014 March 28

Purpose: compute aperture photometry and uncertainties within a
circular aperture

Inputs: 
    optical_band: instrument + filter combination of images
    apcentloc: ascii file containing aperture centers and radii in arcsec

Outputs: Aperture flux and 1-sigma uncertainty

"""

import numpy
from astropy.table import Table
from astropy.io import fits
import os


def rms(im, head, mask, apdiam, nrand=2000):

    # load HST image
    #im = hdulist[0].data
    #head = hdulist[0].header
    headkeys = head.keys()
    cd1_1 = headkeys.count('CD1_1')
    cd1_2 = headkeys.count('CD1_2')
    cd2_1 = headkeys.count('CD2_1')
    cd2_2 = headkeys.count('CD2_2')
    if cd1_1 == 0:
        cdelt1 = numpy.abs(head['CDELT1'] * 3600)
        cdelt2 = numpy.abs(head['CDELT2'] * 3600)
    else:
        cdelt1 = numpy.abs(head['CD1_1'] * 3600)
        cdelt2 = numpy.abs(head['CD2_2'] * 3600)
        cd11 = head['CD1_1']
        #cd12 = head['CD1_2']
        #cd21 = head['CD2_1']
        cd22 = head['CD2_2']
        cdelt1 = numpy.sqrt(cd11 ** 2) * 3600
        cdelt2 = numpy.sqrt(cd22 ** 2) * 3600
        #if cd12 == 0:
        #    cd12 = cd11 / 1e8
        #cdratio = numpy.abs(cd11 / cd12)
        #if cdratio < 1:
        #    cdratio = 1 / cdratio
        #if cdratio < 1e2:
            #print "This shit ain't rotated yo!"
            #import pdb; pdb.set_trace()
    
    # aperture radius in pixels
    pixelscale = numpy.sqrt(cdelt1 ** 2 + cdelt2 ** 2)
    radius = apdiam / 2. / pixelscale

    # Define the mask region for computing rms
    # We compute rms over the inner 1/2 region, excluding the source
    nx = im[:,0].size
    ny = im[0,:].size

    # construct meshgrid of x, y coordinates
    xarr = numpy.arange(nx)
    yarr = numpy.arange(ny)
    xmesh, ymesh = numpy.meshgrid(yarr, xarr)

    # compute fluxes within apertures randomly placed around the image
    x = numpy.random.uniform(size = nrand) * nx
    y = numpy.random.uniform(size = nrand) * ny
    #x = numpy.arange(nx / radius) * radius
    #y = numpy.arange(ny / radius) * radius

    # use only those random apertures that are far enough from the edge
    skyradius = 3. * radius
    maxdx = skyradius
    maxdy = skyradius
    goodrand = (x - maxdx > 0) & \
            (x + maxdx < nx) & \
            (y - maxdy > 0) & \
            (y + maxdy < ny)
    x = x[goodrand]
    y = y[goodrand]
    nrand = x.size

    # remove nearby neighbors
    nearby = numpy.zeros(nrand)
    neighbormask = numpy.zeros(im.shape)
    for irand in numpy.arange(nrand):
        xi = x[irand]
        yi = y[irand]
        maskspot = neighbormask[yi - radius:yi + radius, xi - radius:xi + radius]
        if maskspot.sum() > 0:
            nearby[irand] = 1
        else:
            neighbormask[yi - radius:yi + radius, xi - radius:xi + radius] = 1
    nonearby = nearby == 0
    x = x[nonearby]
    y = y[nonearby]
    nrand = x.size
    sourcemask = numpy.zeros(im.shape)

    # compute aperture fluxes
    noisefluxes = []#numpy.zeros(nrand)
    for irand in numpy.arange(nrand):
    #for xi in x:
    #    for yi in y:

        xi = x[irand]
        yi = y[irand]

        # create a temporary mask initialized to zero
        randommask = numpy.zeros(im.shape)

        # measure the distance from (xi,yi) to any point on the mesh
        dx2d = xmesh - xi
        dy2d = ymesh - yi
        dist2d = numpy.sqrt(dx2d**2 + dy2d**2)

        # make the temporary sky mask
        randomskyregion = (dist2d < skyradius) & \
                (dist2d > radius) & (mask > 0.5)

        # make the temporary aperture mask
        randomsourceregion = (dist2d < radius)
        randommask[randomsourceregion] = 1

        # multiply masks to determine where random aperture overlaps with
        # science aperture
        finalmask = randommask * mask
        fcheck = (finalmask > 0)

        # if there is no overlap, compute the aperture flux
        if finalmask[fcheck].size == 0:
            sourcesum = im[randomsourceregion].sum()
            Nsource = im[randomsourceregion].size
            skysum = numpy.mean(im[randomskyregion]) * Nsource
            noisefluxes.append(sourcesum - skysum)
            # track which random regions were used to compute rms
            sourcemask[randomsourceregion] = 2e3
    noisefluxes = numpy.array(noisefluxes)
    nonzero = (noisefluxes != 0) & (noisefluxes*0 == 0)
    noisefluxes = noisefluxes[nonzero]
    rmsnoisefluxes = numpy.std(noisefluxes)
    scaler = 2.5
    medsky = numpy.median(noisefluxes)
    inliers = (noisefluxes - medsky < scaler*rmsnoisefluxes)
    noisefluxes = noisefluxes[inliers]
    meansky = numpy.mean(noisefluxes)
    medsky = numpy.median(noisefluxes)
    noiseflux = numpy.std(noisefluxes)

    skewsky = numpy.abs(meansky - medsky)
    medsky_previous = 0.
    sky = noisefluxes.copy()
    while skewsky > 0.1 * noiseflux:
        medsky = numpy.median(sky)
        dsky = numpy.abs(sky - medsky)
        rmssky = numpy.std(dsky)
        inliers = (dsky < scaler*rmssky)
        noisefluxestmp = noisefluxes[inliers]
        if len(noisefluxestmp) == len(noisefluxes):
            inliers = (dsky < scaler/2.*rmssky)
        noisefluxes = noisefluxes[inliers]
        sky = noisefluxes.copy()
        medsky = numpy.median(sky)
        dsky = numpy.abs(sky - medsky)
        avgsky = numpy.mean(sky)
        skewsky = numpy.abs(avgsky - medsky)
        #print(medsky, avgsky, skewsky)
        if medsky == medsky_previous:
            scaler /= 1.5
        medsky_previous = medsky

    rms = numpy.std(sky)
    #import matplotlib.pyplot as plt
    #print rms
    #plt.clf()
    #plt.hist(sky)
    #plt.show()
    #plt.clf()
    #nim = im.size
    #imsort = numpy.sort(im.flatten())
    #vmax = imsort[0.8*nim]
    #vmin = imsort[0.2*nim]
    #plt.imshow(im, vmax=vmax, vmin=vmin, origin='lower', cmap='gray_r')
    #plt.contour(sourcemask)
    #plt.show()
    #import pdb; pdb.set_trace()
    return rms


gemzero = [24.91, 28.33, 28.33, 27.93, 26.84]
gemcolors = [0.38, 0.18, 0.1, 0.08, 0.05]
filterlist = ['u', 'g', 'r', 'i', 'z', 'Y', 'J', 'H', 'Ks']

nfilters = len(filterlist)

targetlistloc = 'targetlist.dat'
targetlist = Table.read(targetlistloc, format='ascii')
ntargets = len(targetlist)

fitsdir = 'fitscutouts'
os.chdir(fitsdir)

for itarg in range(ntargets):

    # get the target name from the filename
    dataname = targetlist['dataname'][itarg]
    iauname = targetlist['iauname'][itarg]

    for ifilt in filterlist:

        #dataname = 'ADFS01'
        #ifilt = 'i'
        imloc = dataname + '_' + ifilt + '_cut.fits'
        maskloc = dataname + '_' + ifilt + '_seg.fits'
        apdiam = 2.0
        if not os.path.exists(imloc):
            continue
        im = fits.getdata(imloc)
        head = fits.getheader(imloc)
        mask = fits.getdata(maskloc)
        imrms = rms(im, head, mask, apdiam)

        headerkeys = head.keys()
        cautionflag = True
        #if headerkeys.count('AIRMASS') > 0:
        #    gain = header['GAIN']
        #    exptime = header['EXPTIME']
        #    airmass = header['AIRMASS']
        #    zpt = gemzero[counter] - 2.5 * numpy.log10(gain / exptime) - \
        #        gemcolors[counter] * (airmass - 1.0)
        #    cautionflag = 0
        if headerkeys.count('PHOTZP') > 0:
            zpt = head['PHOTZP']
            cautionflag = False
        if headerkeys.count('MAGZPT') > 0:
            zpt = head['MAGZPT']
            cautionflag = False
        if cautionflag:
            # go back to Gemini reduction directory to get needed header info
            dataloc = '/Users/rbussman/Data/Gemini-South/'
            oldimloc = dataloc + dataname + '_' + ifilt + '.fits'
            oldheader = fits.getheader(oldimloc)
            gain = oldheader['GAIN']
            exptime = oldheader['EXPTIME']
            airmass = oldheader['AIRMASS']
            counter = filterlist.index(ifilt)
            zpt = gemzero[counter] - 2.5 * numpy.log10(gain / exptime) - \
                gemcolors[counter] * (airmass - 1.0)
            cautionflag = False

        depth5 = -2.5 * numpy.log10(5 * imrms) + zpt
        print dataname, ifilt, depth5
