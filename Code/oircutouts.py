"""
Created by Shane Bussmann
2012 January 9
Last modified: 2013 April 14

Purpose: plot cutouts of ALMA imaging + VIKING K + SPIRE 250um

"""

import math
import numpy
from astropy.table import Table
from astropy import wcs
from astropy.io import fits
#from astropy.coordinates import ICRS
from astropy.coordinates import Angle
#from astropy import units as u
from pylab import savefig
import matplotlib
#matplotlib.use('png')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
import yaml
#import phot


# Set font parameters
font = {'family': 'Arial',
        'weight': 'bold',
        'size': 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

filterlist = ['u', 'g', 'r', 'i', 'z', 
        'Y', 'F110W', 'J', 'H', 'Ks',
        'irac1', 'irac2', 'irac3', 'irac4', '24',
        'wise1', 'wise2', 'wise3', 'wise4']
#filterlist = ['i']
nfilters = len(filterlist)
#midirdata = Table.read('../Data/bestmidir.dat', format='ascii')

goodfitloc = '../Data/uvfitlist.dat'
goodfitdat = Table.read(goodfitloc, format='ascii')

targetlistloc = '../Data/targetlist.dat'
targetlist = Table.read(targetlistloc, format='ascii')
ntargets = len(targetlist)

# center and size of overlays
ra_centers = targetlist['ra_250']
dec_centers = targetlist['dec_250']

figwidth = 15.0#2.0 * nfilters
figlength = 12.0# * ntargets
fig = plt.figure(figsize=(figwidth, figlength))
plt.subplots_adjust(left=0.05, right=0.97, top=0.98, bottom=0.09)

# clear the plotting window
plt.clf()

for itarg in range(0, ntargets):

    counter = 1
    plt.clf()

    # get the target name from the filename
    dataname = targetlist['dataname'][itarg]
    iauname = targetlist['iauname'][itarg]

    # load the corresponding ALMA image
    submm_band = '870'
    ALMA_file = '../../fitscutouts/' + dataname + '_' + submm_band + '.fits'
    ALMA_image = fits.getdata(ALMA_file)
    ALMA_image = ALMA_image[0, 0, :, :].copy() * 1e3
    ALMA_header = fits.getheader(ALMA_file)

    # get the appropriate image center and radial extent
    goodfit = goodfitdat['intrinsic'][itarg]
    modfitloc = '../../../../ModelFits/' + dataname + '/' + goodfit
    configloc = modfitloc + '/config.yaml'
    configfile = open(configloc)
    config = yaml.load(configfile)
    
    ra_center = config['Region0']['RACentroid']#ra_centers[itarg]
    dec_center = config['Region0']['DecCentroid']#ra_centers[itarg]
    radial_extent = 9.0#config['Region0']['RadialExtent']
    print radial_extent
    #dec_center = dec_centers[itarg]
    #ra_center = Angle(ra_center + ' hours').degree
    #dec_center = Angle(dec_center + ' degrees').degree
    #ra_center = 34.826736
    #dec_center = -3.1809444
    #ra_center = 34.673403
    #dec_center = -3.8343472

    # set the radial extent to be equal to the primary beam FWHM
    ALMA_FWHM = 1.2 * 0.87e-3 / 12. * 206265
    #radial_extent = 0.6 * ALMA_FWHM

    oldradialextent = radial_extent
    for ifilt in filterlist:
        plt.subplot(4, 5, counter)

        # grab the axes instance
        ax = fig.add_subplot(4, 5, counter)
        #ax = fig.gca()

        if counter > 10:
            radial_extent = 9.

        # load the appropriate optical image
        optloc = '../../fitscutouts/' + dataname + '_' + ifilt + '.fits'
        if not os.path.exists(optloc):
            plt.plot(numpy.arange(100))
            plt.text(0.05, 0.05, ifilt, fontsize='x-large', \
                color='black', transform=ax.transAxes)
            if ifilt == 'u':
                plt.text(0.05, 0.80, dataname, fontsize='xx-large', \
                    color='black', transform=ax.transAxes)
            counter += 1
            continue

        optical_hdu = fits.open(optloc)
        #hoho = phot.rms(optical_hdu, 4.0)
        optical_image = fits.getdata(optloc)
        optical_header = fits.getheader(optloc)
        if dataname == 'XMM02':
            if ifilt == 'J':
                optical_header = optical_hdu[1].header
            if ifilt == 'K':
                optical_header = optical_hdu[1].header
        if dataname == 'XMM108':
            if ifilt == 'K':
                optical_header = optical_hdu[1].header

        # compute the (x, y) center and pixel scale in the optical image
        wcs_optical = wcs.WCS(optical_header, naxis=2)
        pixxy = wcs_optical.wcs_world2pix(ra_center, dec_center, 1)
        x_optical = numpy.round(pixxy[0])
        y_optical = numpy.round(pixxy[1])
        headerkeys = optical_header.keys()
        cd1_1 = headerkeys.count('CD1_1')
        if cd1_1 == 0:
            cdelt1_optical = numpy.abs(optical_header['CDELT1'] * 3600)
            cdelt2_optical = numpy.abs(optical_header['CDELT2'] * 3600)
        else:
            cdelt1_optical = numpy.abs(optical_header['CD1_1'] * 3600)
            cdelt2_optical = numpy.abs(optical_header['CD2_2'] * 3600)
            cd11 = optical_header['CD1_1']
            cd12 = optical_header['CD1_2']
            cd21 = optical_header['CD2_1']
            cd22 = optical_header['CD2_2']
            cdelt1_optical = numpy.sqrt(cd11 ** 2 + cd12 ** 2) * 3600
            cdelt2_optical = numpy.sqrt(cd21 ** 2 + cd22 ** 2) * 3600
            if cd12 == 0:
                cd12 = cd11 / 1e8
            cdratio = numpy.abs(cd11 / cd12)
            if cdratio < 1:
                cdratio = 1 / cdratio
            #if cdratio < 1e2:
                #print "This shit ain't rotated yo!"
                #import pdb; pdb.set_trace()
        x_extent = numpy.round(radial_extent / cdelt1_optical)
        y_extent = numpy.round(radial_extent / cdelt2_optical)
        nx = optical_image[:, 0].size
        ny = optical_image[0, :].size

        angularradius = 0.5 * 60
        pixscale = numpy.sqrt(cdelt1_optical * cdelt2_optical)
        pixradius = numpy.round(angularradius / pixscale).astype(int)
        x1 = str(x_optical.astype(int) - pixradius)
        x2 = str(x_optical.astype(int) + pixradius)
        y1 = str(y_optical.astype(int) - pixradius)
        y2 = str(y_optical.astype(int) + pixradius)
        region = '[' + x1 + ':' + x2 + ',' + y1 + ':' + y2 + ']'
        cutfile = 'fitscutouts/' + dataname + '_' + ifilt + '_cut.fits'
        cmd = 'imcopy ' + optloc + region + ' ' + cutfile
        #print cmd
        
        # make the trimmed optical image
        #if dataname == 'XMM109':
        #    optical_trimmed = numpy.ones([2 * y_extent, 2 * x_extent])
        #else:
        if (x_optical - x_extent < 0) | (x_optical + x_extent > nx) | \
                (y_optical - y_extent < 0) | (y_optical + y_extent > ny):
            # pad with zeros
            import pdb; pdb.set_trace()
            optical_trimmed = numpy.zeros([2 * y_extent, 2 * x_extent])
            nx_pad = optical_trimmed[0, :].size
            ny_pad = optical_trimmed[:, 0].size
            if x_optical - x_extent < 0:
                xr0 = 0
                xrr0 = x_extent - x_optical
            else:
                xr0 = x_optical - x_extent
                xrr0 = 0
            if y_optical - y_extent < 0:
                yr0 = 0
                yrr0 = y_extent - y_optical
            else:
                yr0 = y_optical - y_extent
                yrr0 = 0
            if x_optical + x_extent > nx:
                xr1 = nx
                xrr1 = nx_pad / 2 + (nx - x_optical)
            else:
                xr1 = x_optical + x_extent
                xrr1 = nx_pad
            if y_optical + y_extent > ny:
                yr1 = ny
                yrr1 = ny_pad / 2 + (ny - y_optical)
            else:
                yr1 = y_optical + y_extent
                yrr1 = ny_pad
            optical_trimmed[yrr0:yrr1, xrr0:xrr1] = optical_image[yr0:yr1, xr0:xr1]
        else:
            xr0 = x_optical - x_extent
            yr0 = y_optical - y_extent
            xr1 = x_optical + x_extent
            yr1 = y_optical + y_extent
            optical_trimmed = optical_image[yr0:yr1, xr0:xr1]

        # compute the (x, y) center and pixel scale in the ALMA image
        wcs_ALMA = wcs.WCS(ALMA_header, naxis=2)
        pixxy = wcs_ALMA.wcs_world2pix(ra_center, dec_center, 1)
        x_ALMA = numpy.round(pixxy[0])
        y_ALMA = numpy.round(pixxy[1])

        # grab image resolution from header data
        bmaj = ALMA_header['BMAJ'] * 3600
        bmin = ALMA_header['BMIN'] * 3600
        bpa = ALMA_header['BPA']
        cdelt1_ALMA = numpy.abs(ALMA_header['CDELT1'] * 3600)
        cdelt2_ALMA = numpy.abs(ALMA_header['CDELT2'] * 3600)
        pixels_per_beam = math.pi * bmaj / 2 * bmin / 2 / cdelt1_ALMA / \
            cdelt2_ALMA / math.log(2)
        bmajmin = math.sqrt(bmaj * bmin)

        # Define the mask region for computing rms
        # We compute rms over the entire ALMA region, excluding the source
        mask = ALMA_image.copy()
        mask[:] = 0
        nx = mask[:, 0].size
        ny = mask[0, :].size
        innerhalf = 0
        mask[innerhalf:ny - innerhalf, innerhalf:nx - innerhalf] = 1
        x_extent = radial_extent / 4. / cdelt1_ALMA
        y_extent = radial_extent / 4. / cdelt2_ALMA
        xr0 = x_ALMA - x_extent
        yr0 = y_ALMA - y_extent
        xr1 = x_ALMA + x_extent
        yr1 = y_ALMA + y_extent
        mask[yr0:yr1, xr0:xr1] = 0
        goodregion = mask == 1
        rms = ALMA_image[goodregion].std()
        print iauname, dataname, ifilt, rms

        # make trimmed ALMA image
        x_extent = radial_extent / cdelt1_ALMA
        y_extent = radial_extent / cdelt2_ALMA
        if (x_ALMA - x_extent < 0) | (x_ALMA + x_extent > nx) | \
                (y_ALMA - y_extent < 0) | (y_ALMA + y_extent > ny):
            # pad with zeros
            ALMA_trimmed = numpy.zeros([2 * y_extent, 2 * x_extent])
            nx_pad = ALMA_trimmed[0, :].size
            ny_pad = ALMA_trimmed[:, 0].size
            if x_ALMA - x_extent < 0:
                xr0 = 0
                xrr0 = x_extent - x_ALMA
            else:
                xr0 = x_ALMA - x_extent
                xrr0 = 0
            if y_ALMA - y_extent < 0:
                yr0 = 0
                yrr0 = y_extent - y_ALMA
            else:
                yr0 = y_ALMA - y_extent
                yrr0 = 0
            if x_ALMA + x_extent > nx:
                xr1 = nx
                xrr1 = nx_pad / 2 + (nx - x_ALMA)
            else:
                xr1 = x_ALMA + x_extent
                xrr1 = nx_pad
            if y_ALMA + y_extent > ny:
                yr1 = ny
                yrr1 = ny_pad / 2 + (ny - y_ALMA)
            else:
                yr1 = y_ALMA + y_extent
                yrr1 = ny_pad
            ALMA_trimmed[yrr0:yrr1, xrr0:xrr1] = ALMA_image[yr0:yr1, xr0:xr1]
        else:
            xr0 = x_ALMA - x_extent
            yr0 = y_ALMA - y_extent
            xr1 = x_ALMA + x_extent
            yr1 = y_ALMA + y_extent
            ALMA_trimmed = ALMA_image[yr0:yr1, xr0:xr1]
        #if x_extent > nx / 2:
        #    # pad with zeros
        #    ALMA_trimmed = numpy.zeros([2 * y_extent, 2 * x_extent])
        #    nx_pad = ALMA_trimmed[0, :].size
        #    ny_pad = ALMA_trimmed[:, 0].size
        #    xr0 = nx_pad / 2 - nx / 2
        #    yr0 = ny_pad / 2 - ny / 2
        #    xr1 = nx_pad / 2 + nx / 2
        #    yr1 = ny_pad / 2 + ny / 2
        #    ALMA_trimmed[xr0:xr1, yr0:yr1] = ALMA_image
        #else:
        #    xr0 = nx / 2 - x_extent
        #    yr0 = ny / 2 - y_extent
        #    xr1 = nx / 2 + x_extent
        #    yr1 = ny / 2 + y_extent
        #    ALMA_trimmed = ALMA_image[yr0:yr1, xr0:xr1]


        # print the target name in the upper left corner
        start = iauname.find('J')
        iauaddress = 'HerMES ' + iauname[start:]
        #plt.text(0.5, 1.05, iauaddress, fontsize='xx-large', \
        #    color='black', transform=ax.transAxes, ha='center')
        plt.text(0.05, 0.80, dataname, fontsize='xx-large', \
            color='black', transform=ax.transAxes)
        plt.text(0.05, 0.05, ifilt, fontsize='x-large', \
            color='black', transform=ax.transAxes)
        #plt.subplot(1, 1, 1)

        extent = [radial_extent, -radial_extent, -radial_extent, radial_extent]
        ok = optical_trimmed * 0 == 0
        minopt = optical_trimmed[ok].min()
        if counter < 10:
            optical_trimmed = optical_trimmed - minopt + 1
            optical_trimmed = numpy.log10(optical_trimmed)
        goodpix = (numpy.abs(optical_trimmed) < 1e5) & \
            (numpy.abs(optical_trimmed) > 0)
        pospix = goodpix.copy()#(optical_trimmed < 1e5) & \
            #(optical_trimmed > 0)
        #ALMA_trimmed = numpy.log10(ALMA_trimmed)
        imsort = numpy.sort(optical_trimmed[pospix]).flatten()
        npix = optical_trimmed[pospix].size
        vmax = imsort[0.999 * npix]
        if dataname == 'ADFS01':
            vmax = imsort[0.95 * npix]
        if dataname == 'XMM02':
            vmax = imsort[0.95 * npix]
        if dataname == 'XMM03':
            vmax = imsort[0.999 * npix]
        if dataname == 'XMM101':
            vmax = imsort[0.95 * npix]
        if dataname == 'XMM108':
            vmax = imsort[0.999 * npix]
        #if dataname == 'XMM109':
        #    vmax = imsort[0.9 * npix]
        vmin = imsort[0.5 * npix]
        if counter > 10:
            vmin = imsort.min()
            vmax = imsort.max()
        plt.imshow(optical_trimmed, cmap='gray_r', interpolation='nearest', \
            extent=extent, origin='lower', vmax=vmax, vmin=vmin)

        #plt.colorbar()

        # define the x- and y-vectors for the contour plot
        nx_ALMA = ALMA_trimmed[0, :].size
        ny_ALMA = ALMA_trimmed[:, 0].size
        x_vector = (nx_ALMA / 2 - numpy.arange(nx_ALMA)) * cdelt1_ALMA
        y_vector = (numpy.arange(ny_ALMA) - ny_ALMA / 2) * cdelt2_ALMA
        #cellplus = celldata*(2*xrad+1.1)/(2*xrad)
        #cmodx = ( numpy.arange(2*xrad) - xrad ) * (-cellplus) - celldata/2.
        #cmody = ( numpy.arange(2*yrad) - yrad ) * cellplus + celldata/2.

        # define contour level parameters
        #plevs = [4 * rms]
        #nlevs = [-4 * rms]
        plevs = 3 * rms * 2 ** (numpy.arange(10))
        nlevs = sorted(-3 * rms * 2 ** (numpy.arange(4)))
        #pcline = 'solid'
        #ncline = 'dashed'

        # draw the contours
        plt.contour(x_vector, y_vector, ALMA_trimmed, colors='red', levels=plevs, \
            linewidths=1.0)
        plt.contour(x_vector, y_vector, ALMA_trimmed, colors='blue', levels=nlevs, \
            linewidths=1.0)

        plt.minorticks_on()
        plt.tick_params(width=1.5, which='both')
        plt.tick_params(length=2, which='minor')
        plt.tick_params(length=4, which='major')
        #plt.xlabel(r'$\Delta$RA (arcsec)', fontsize='x-large')
        #plt.ylabel(r'$\Delta$Dec (arcsec)', fontsize='x-large')

        #xhi = xrad * cdelt1_optical
        #xlo = -xrad * celldata
        #yhi = yrad * celldata
        #ylo = -yrad * celldata
        #axisrange = [numpy.float(xhi), numpy.float(xlo), numpy.float(ylo),
        #    numpy.float(yhi)]
        plt.axis(extent)

        normx = extent[0] - extent[1]
        normy = extent[3] - extent[2]
        bparad = bpa / 180 * math.pi
        beamx = numpy.abs(numpy.sin(bparad) * bmaj) + numpy.abs(numpy.cos(bparad) *
            bmin)
        beamy = numpy.abs(numpy.cos(bparad) * bmaj) + numpy.abs(numpy.sin(bparad) *
            bmin)
        dx = normx  # xhi - xlo
        dy = normy  # yhi - ylo
        bufferx = 0.05
        buffery = 0.05
        xpos = 1 - beamx / dx / 2 - bufferx
        ypos = beamy / dy / 2 + buffery
        e = Ellipse((xpos, ypos), bmin / normx, bmaj / normy, angle=bpa,
            ec='black', lw=1.0, transform=ax.transAxes, fc='white',
            zorder=10, hatch='////')
        ax.add_artist(e)
        #plt.plot([-10,-8],[10,10], color='black')
        #ax.set_xticklabels([])
        #ax.set_yticklabels([])
        #plt.tight_layout()

        # save the figure to a file
        rename = dataname.replace('.', '_')
        counter += 1
        #savefig('oircutouts.pdf')
        #import pdb; pdb.set_trace()
    savefig('../Figures/cutouts_oirmir_alma_' + rename + '.pdf')
import pdb; pdb.set_trace()
