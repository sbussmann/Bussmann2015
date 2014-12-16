"""
Created by Shane Bussmann
2012 January 11

Purpose: compute aperture photometry and uncertainties within a
circular aperture; relies on PHOTFNU header keyword in HST images

Inputs: 
    optical_band: instrument + filter combination of images
    apcentloc: ascii file containing aperture centers and radii in arcsec

Outputs: Aperture flux and 1-sigma uncertainty

"""


import numpy
import asciitable
import pywcs
import glob
import pyfits
import matplotlib.pylab as plt

# specify HST and filter combination
optical_band = 'HST_F110W'

# specify how many targets to skip
skipnum = 0

# read in images
dataloc = 'images/'
optical_band_locs = dataloc + '*XMM05*' + optical_band + '.fits.gz'
optical_band_files = glob.glob(optical_band_locs)

# read in aperture centers and radial extents
centloc = 'apcent_' + optical_band + '.txt'
centers = asciitable.read(centloc)
target_names = centers['target_name']
ra_centers = centers['RA']
dec_centers = centers['Dec']
radial_extents = centers['radial_extent']

naper = len(centers)
fluxes = numpy.zeros(naper)
errors = numpy.zeros(naper)

# initiate a counter to track the aperture number
apnum = 0

for file in optical_band_files[skipnum:]:

    # load HST image
    im = pyfits.getdata(file)
    head = pyfits.getheader(file)

    # get WCS data
    wcs = pywcs.WCS(head)

    # get pixel scale in HST image
    CD1_1 = head['CD1_1'] * 3600
    CD1_2 = head['CD1_2'] * 3600
    CD2_1 = head['CD2_1'] * 3600
    CD2_2 = head['CD2_2'] * 3600
    cdelt1 = numpy.sqrt(CD1_1 ** 2 + CD1_2 ** 2)
    cdelt2 = numpy.sqrt(CD2_1 ** 2 + CD2_2 ** 2)
    pixelscale = numpy.sqrt(cdelt1 ** 2 + cdelt2 ** 2)

    # scale image by photometry conversion factor, Jy / (electron/sec)
    photfnu = head['PHOTFNU']
    im = im * photfnu

    # get the target name from the filename
    indx_start = file.find('images') + 7
    indx_end = file.find(optical_band) - 1
    target_name = file[indx_start:indx_end]

    # get the number of apertures associated with this target
    target_indx = target_names == target_name
    naper_target = target_names[target_indx].size

    # loop over each aperture for this target
    for targ in numpy.arange(naper_target):

        # get RA, Dec, and radial extent for this target/aperture
        ra = ra_centers[target_indx][targ]
        dec = dec_centers[target_indx][targ]
        radius = radial_extents[target_indx][targ]

        # get x, y pixel location of emission centroid
        pixxy = wcs.wcs_sky2pix(ra, dec, 1)
        if target_name == 'XMM05':
            pixxy = wcs.wcs_sky2pix(dec, ra, 1)
        x_center = numpy.round(pixxy[0][0])
        y_center = numpy.round(pixxy[1][0])

        # Define the mask region for computing rms
        # We compute rms over the inner 1/2 region, excluding the source
        mask = im.copy()
        mask[:] = 0
        nx = mask[:,0].size
        ny = mask[0,:].size

        # construct meshgrid of x, y coordinates
        xarr = numpy.arange(nx)
        yarr = numpy.arange(ny)
        xmesh, ymesh = numpy.meshgrid(yarr, xarr)
        dx2d = xmesh - x_center
        dy2d = ymesh - y_center
        dist2d = numpy.sqrt(dx2d**2 + dy2d**2)

        # make the sky aperture mask
        skyradius = 2.5 * radius
        skyradpix = skyradius / pixelscale
        skyregion = (dist2d < skyradpix)
        mask[skyregion] = 4

        # make the source aperture mask
        apradpix = radius / pixelscale
        sourceregion = (dist2d < apradpix) & \
                (im * 0 == 0)
        mask[sourceregion] = 2

        # assign the rms region to the sky aperture surrounding the source
        finalmask = mask.copy()
        rmsregion = (finalmask == 4)

        # compute the rms and background level in the sky aperture
        rms = im[rmsregion].std()
        bg = numpy.median(im[rmsregion])

        # measure photometry of target/aperture
        flux = im[sourceregion].sum()
        fluxes[apnum] = flux

        # compute fluxes within apertures randomly placed around the image
        nrand = 10
        x = numpy.random.uniform(size = nrand) * nx
        y = numpy.random.uniform(size = nrand) * ny

        # use only those random apertures that are far enough from the edge
        maxdx = skyradius / pixelscale
        maxdy = skyradius / pixelscale
        goodrand = (x - maxdx > 0) & \
                (x + maxdx < nx) & \
                (y - maxdy > 0) & \
                (y + maxdy < ny)
        x = x[goodrand]
        y = y[goodrand]
        nrand = x.size

        # compute aperture fluxes
        noisefluxes = numpy.zeros(nrand)
        for iter in numpy.arange(nrand):

            xi = x[iter]
            yi = y[iter]

            # initialize the temporary mask for the random pixel location
            tempmask = finalmask.copy()
            tempmask[:] = 0

            dx2d = xmesh - xi
            dy2d = ymesh - yi
            dist2d = numpy.sqrt(dx2d**2 + dy2d**2)

            # make the temporary sky mask
            tempregion = (dist2d < skyradpix)
            tempmask[tempregion] = 4

            # make the temporary aperture mask
            tempregion = (dist2d < apradpix)
            tempmask[tempregion] = 1

            # multiply masks to determine where random aperture overlaps with
            # science aperture
            productmask = tempmask * finalmask
            fcheck = (productmask > 4)

            # if there is no overlap, compute the aperture flux
            if productmask[fcheck].size == 0:
                noisefluxes[iter] = im[tempregion].sum()

        nonzero = (noisefluxes != 0)
        medsky = numpy.median(noisefluxes[nonzero])
        noiseflux = numpy.std(noisefluxes[nonzero])
        errors[apnum] = noiseflux

        # increment the aperture number by +1
        abmag = -2.5 * numpy.log10(flux / 1e23) - 48.6
        print target_name, abmag, flux, noiseflux, rms, medsky
        apnum += 1

    plt.clf()
    plt.imshow(im)
    plt.plot(x[nonzero], y[nonzero], '.')
    #plt.contour(mask)
    plevs = numpy.arange(4)
    plt.contour(finalmask, colors='red', levels=plevs, linewidths=1.5)
    #import pdb; pdb.set_trace()

import pdb; pdb.set_trace()
