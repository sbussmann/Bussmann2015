"""
Created by Shane Bussmann
2012 December 4

Purpose: compute aperture photometry and uncertainties

Inputs: Image, aperture center, aperture shape

Outputs: Aperture flux and 1-sigma uncertainty

"""


import math
import numpy
import os
from astropy.table import Table
from astropy import wcs
import glob
#import distutils.dir_util
from astropy.io import fits
import matplotlib.pylab as plt
from pylab import savefig

# specify how many targets to skip at start
skipnum = 14

#--------------------------------------------------------------------------
# Step 1: read in image and beam
datadir = 'linesearch/'
dataloc = datadir + '*linesearch.fits'
files = glob.glob(dataloc)

nfiles = len(files)
bmaj = numpy.zeros(nfiles)
bmin = numpy.zeros(nfiles)
extents = numpy.zeros(nfiles)

counter = 1
for ifile in files[skipnum:skipnum+2]:
    print ifile

    # get the dataname from the filename
    startname = 11
    endname = ifile.find('.linesearch')
    dataname = ifile[startname:endname]

    # aperture info
    imcentloc = datadir + dataname + '.apertures'
    imcentdat = Table.read(imcentloc, format='ascii')
    napers = len(imcentdat)

    # load datacube
    imfile = fits.open(ifile)
    limfile = len(imfile)
    smaimall = fits.getdata(ifile)
    nchan = smaimall[0, :, 0, 0].size

    # initialize the coadded spectrum
    spec_coadd = 0

    for iaper in range(napers):
        ra_centroid = imcentdat['ra_centroid'][iaper]
        dec_centroid = imcentdat['dec_centroid'][iaper]
        extent = imcentdat['extent'][iaper]

        spectrum = numpy.zeros(nchan)
        e_spectrum = numpy.zeros(nchan)

        for chan in range(nchan):

            # convert to mJy/beam
            smaim = smaimall[0,chan,:,:].copy() * 1e3
            if smaim.sum() == 0:
                continue
            headsmaim = fits.getheader(ifile)

            # get phase center of observations
            smawcs = wcs.WCS(headsmaim, naxis=2)
            nx = smaim[0,:].size
            ny = smaim[:,0].size

            # get x, y pixel location of emission centroid
            pixxy = smawcs.wcs_world2pix(ra_centroid, dec_centroid, 1)
            x_centroid = numpy.round(pixxy[0])
            y_centroid = numpy.round(pixxy[1])

            # grab image resolution from header data
            if limfile == 1:
                bmaj = headsmaim['BMAJ'] * 3600
                bmin = headsmaim['BMIN'] * 3600
            else:
                bmaj = imfile[1].data['BMAJ'].mean()
                bmin = imfile[1].data['BMIN'].mean()
            cdelt1 = headsmaim['CDELT1'] * 3600
            cdelt2 = headsmaim['CDELT2'] * 3600
            celldata = math.sqrt( abs(cdelt1) * abs(cdelt2) )
            cdelt = celldata / 3600 #* math.pi / 180
            strcdelt = str(cdelt)
            npix_sma = math.pi * bmaj/2 * bmin/2 / celldata**2 / math.log(2)

            # Define the mask region for computing rms
            # We compute rms over the inner 1/2 region, excluding the source
            mask = smaim.copy()
            mask[:] = 0
            nx = mask[:,0].size
            ny = mask[0,:].size
            innerhalf = nx / 4
            mask[innerhalf:ny-innerhalf,innerhalf:nx-innerhalf] = 1
            x_extent = extent / celldata
            y_extent = extent / celldata
            xlo = x_centroid - x_extent
            xhi = x_centroid + x_extent
            ylo = y_centroid - y_extent
            yhi = y_centroid + y_extent

            # set the source region in the mask file = 2
            mask[ylo:yhi, xlo:xhi] = 2

            rmsregion = mask == 1
            sourceregion = mask == 2
            finalmask = mask.copy()

            # subtract offset from SMA imaging to get a zero background level
            #bg = im1[goodregion].median()
            #smaim = smaim - bg
            #print, bg

            # compute sigma image from cutout of SMA flux image
            #Y = histogram( im1[goodregion], bin=0.0006, locations=X )
            #Result = GAUSSFIT( X, Y, A )
            #rms = A[2]
            rms = smaim[rmsregion].std()
            flux = smaim[sourceregion].sum() / npix_sma
            nbeams = smaim[sourceregion].size / npix_sma
            #print rms

            # compute primary beam correction factor
            freq = headsmaim['CRVAL3']
            c = 3e8
            lam = c / freq
            diam = 12.
            fwhm = 1.2 * lam / diam
            halfradius = fwhm / (2 * math.sqrt(2 * math.log(2))) * 206265

            # compute offset between phase center and emission centroid
            phasecent_ra = headsmaim['CRVAL1']
            phasecent_dec = headsmaim['CRVAL2']
            off_ra = ra_centroid - phasecent_ra
            off_dec = dec_centroid - phasecent_dec
            offtot = math.sqrt(off_ra**2 + off_dec**2) * 3600
            scalefac = math.exp(-0.5*(offtot/halfradius)**2)
            flux = flux / scalefac
            rms = rms / scalefac
            spectrum[chan] = flux
            e_spectrum[chan] = rms

            #mpl.clf()
            #mpl.imshow(smaim)
            #mpl.plot(x[nonzero], y[nonzero], '.')
            ##mpl.contour(mask)
            #plevs = numpy.arange(3)
            #mpl.contour(finalmask, colors='red', levels=plevs, linewidths=1.5)
            #print objname, flux, noiseflux, rms, nbeams, math.sqrt(nbeams)*rms, medsky
            #rename = objname.replace('.', '_')
            #savefig(saveloc + rename + '_photometry')

        spectrum = spectrum[spectrum.nonzero()]
        spec_coadd += spectrum
        plt.clf()
        fig = plt.figure(figsize=(12.0, 4.0))
        plt.subplots_adjust(left=0.05, right=0.97, top=0.98, bottom=0.09)
        #ax = fig.add_subplot(1, 2, 1)        
        ax = plt.subplot2grid((1,4), (0,0), colspan=3)
        plt.plot(spectrum, drawstyle='steps-mid')

        plt.text(0.05, 0.80, dataname, fontsize='xx-large', \
            color='black', transform=ax.transAxes)

        #ax2 = fig.add_subplot(1, 2, 2)        
        ax2 = plt.subplot2grid((1,4), (0,3), colspan=1)

        # load the corresponding continuum image
        ALMA_file = datadir + dataname + '_ALMA_870um.fits'
        ALMA_image = fits.getdata(ALMA_file)
        ALMA_image = ALMA_image[0, 0, :, :].copy() * 1e3
        ALMA_header = fits.getheader(ALMA_file)

        # compute the (x, y) center and pixel scale in the ALMA image
        wcs_ALMA = wcs.WCS(ALMA_header, naxis=2)
        pixxy = wcs_ALMA.wcs_world2pix(ra_centroid, dec_centroid, 1)
        x_ALMA = numpy.round(pixxy[0])
        y_ALMA = numpy.round(pixxy[1])

        # grab image resolution from header data
        bmaj = ALMA_header['BMAJ'] * 3600
        bmin = ALMA_header['BMIN'] * 3600
        cdelt1_ALMA = numpy.abs(ALMA_header['CDELT1'] * 3600)
        cdelt2_ALMA = numpy.abs(ALMA_header['CDELT2'] * 3600)
        pixels_per_beam = math.pi * bmaj / 2 * bmin / 2 / cdelt1_ALMA / \
            cdelt2_ALMA / math.log(2)
        bmajmin = math.sqrt(bmaj * bmin)

        x_extent = extent / cdelt1_ALMA
        y_extent = extent / cdelt2_ALMA
        xr0 = x_ALMA - x_extent
        yr0 = y_ALMA - y_extent
        xr1 = x_ALMA + x_extent
        yr1 = y_ALMA + y_extent
        ALMA_trimmed = ALMA_image[yr0:yr1, xr0:xr1]

        extent = [extent, -extent, -extent, extent]
        #pospix = (ALMA_trimmed < 1e5) & \
        #    (ALMA_trimmed > 0)        
        #imsort = numpy.sort(ALMA_trimmed[pospix]).flatten()
        #npix = ALMA_trimmed[pospix].size
        #vmax = imsort[0.99 * npix]        
        #vmin = imsort[0.5 * npix]        
        plt.imshow(ALMA_trimmed, interpolation='nearest', \
            extent=extent, origin='lower')        

        sri = str(iaper)
        plt.text(0.05, 0.80, dataname + '_' + sri, fontsize='xx-large', \
            color='black', transform=ax2.transAxes)
        
        savefig(datadir + dataname + '_' + sri + '_linesearch')

    plt.clf()
    fig = plt.figure(figsize=(9.0, 4.0))
    plt.subplots_adjust(left=0.05, right=0.97, top=0.98, bottom=0.09)
    ax = fig.add_subplot(1, 1, 1)        
    plt.plot(spec_coadd, drawstyle='steps-mid')
    plt.text(0.05, 0.80, dataname, fontsize='xx-large', \
        color='black', transform=ax.transAxes)
    
    savefig(datadir + dataname + '_linesearch')
    #import pdb; pdb.set_trace()

import pdb; pdb.set_trace()
