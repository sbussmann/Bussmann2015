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
import pdb
from scipy import interpolate
import aplpy
import sys
import imp
import yaml

# define the instrument and wavelength of the optical imaging
#optical_band = 'HST_F110W'
optical_band = 'SDSS_i'
#optical_band = 'KeckAO_Ks'
optical_filter = optical_band[-1]

submm_band = 'ALMA_870um'

# Set font parameters
font = {'family': 'Arial Narrow',
        'weight': 'bold',
        'size': 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

# set up the plot window
fig = plt.figure(figsize=(6.5, 5.5))
plt.subplots_adjust(left=0.05, right=1.03, top=0.98, bottom=0.05)

# identify ALMA targets
targetloc = '../Data/targetlist.dat'
targets = Table.read(targetloc, format='ascii')
ntargets = len(targets)

# get the optical image locations
opticallocs = '../Data/bestoptical.dat'
opticals = Table.read(opticallocs, format='ascii')

# center and size of overlays
ra_centers = targets['ra_250']
dec_centers = targets['dec_250']
#radial_extents = centers['radial_extent']

# define the location of the ALMA images
almaloc = '../../../../ModelFits/'

for targ in range(0, ntargets):

    plt.clf()
    # get the target name from the filename
    dataname = targets['dataname'][targ]
    iauname = targets['iauname'][targ]
    lensgrade = targets['lensgrade'][targ]

    herschel_band = str(targets['bestspire'][targ])
    herschel_filter = herschel_band[-1]

    # if there is ALMA data, proceed
    ALMA_file = '../../fitscutouts/' + dataname + '_870.1000.fits'
        
    # load the corresponding ALMA image
    ALMA_image = fits.getdata(ALMA_file)
    ALMA_image = ALMA_image[0, 0, :, :].copy()
    ALMA_header = fits.getheader(ALMA_file)
    RA_phscen = ALMA_header['OBSRA']
    Dec_phscen = ALMA_header['OBSDEC']

    herschelimloc = '../../fitscutouts/' + dataname + '_' + \
            herschel_band+ '.fits'

    # load the Herschel image
    herschel_file = fits.open(herschelimloc)
    herschel_image = herschel_file[0].data
    herschel_header = herschel_file[0].header

    # load the appropriate optical image
    #opt_indx = opticals['iauname'] == iauname
    #optloc = opticals['bestloc'][opt_indx][0]
    #opttel = opticals['telescope'][opt_indx][0]
    #optfilt = opticals['filter'][opt_indx][0]
    #opticalinfo = opttel + ' ' + optfilt
    #optical_image = fits.getdata(optloc)
    #optical_header = fits.getheader(optloc)

    # get the appropriate image center and radial extent
    cwd = almaloc + dataname + '/uvfit25/'
    #print cwd
    #sys.path.append(cwd)
    #if targ == 0:
    #    import config
    #config = imp.load_source("config", cwd + "/config.py")
    configloc = cwd + 'config.yaml'
    configfile = open(configloc, 'r')
    config = yaml.load(configfile)
    #reload(config)
    ralist = config['Region0']['RACentroid']
    declist = config['Region0']['DecCentroid']
    extents = config['Region0']['RadialExtent']
    ra_center = ralist#(ralist.min() + ralist.max()) / 2.
    dec_center = declist#(declist.min() + declist.max()) / 2.
    
    #ra_center = ra_centers[targ]
    #dec_center = dec_centers[targ]
    #ra_center = Angle(ra_center + ' hours').degree
    #dec_center = Angle(dec_center + ' degrees').degree
    #radial_extent = radial_extents[target_indx][0]

    # set the radial extent to be equal to the primary beam FWHM
    ALMA_FWHM = 1.2 * 0.87e-3 / 12. * 206265
    radial_extent = 0.75*ALMA_FWHM
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
    #print iauname, dataname, rms

    gc = aplpy.FITSFigure(ALMA_file, figure=fig, subplot=[0.07,0.10,0.95,0.88])
    gc.show_colorscale(cmap='jet', pmax=100., vmin=-1.0)
    plevs = sorted(4 * 0.007 * 2 ** (numpy.arange(8) / 2.))
    gc.show_contour(herschelimloc, colors='white', linewidths=2.5, levels=plevs,
            kernel='box', linestyles='solid', ignore_missing_end=True)
    gc.show_contour(herschelimloc, colors='black', linewidths=1,
            levels=plevs, kernel='box', linestyles='solid',
            ignore_missing_end=True)
    gc.recenter(RA_phscen, Dec_phscen, radius=radial_extent/3600)
    nregions = 1#ralist.size
    for iregion in range(nregions):
        ra_center = ralist
        dec_center = declist
        extent = 2 * extents / 3600
        gc.show_rectangles(ra_center, dec_center, extent, extent,
            edgecolor='white', linestyle='dashed', linewidth=2)
    gc.show_circles(RA_phscen, Dec_phscen, ALMA_FWHM / 7200, 
            edgecolor='white', linewidth=2)
    gc.add_beam()
    gc.beam.set_color('white')
    gc.beam.set_corner('bottom left')
    gc.add_colorbar()
    #gc.colorbar.set_axis_label_text(r'Intensity (mJy / beam)')
    gc.axis_labels.hide()
    gc.ticks.show()
    gc.ticks.set_color('black')
    gc.ticks.set_linewidth(1.5)
    gc.frame.set_linewidth(1.5)
    #gc.axis_labels.set_font(size='x-large', weight='bold', \
    #                     stretch='normal', #family='sans-serif', \
    #                     style='normal', variant='normal')
    gc.tick_labels.set_font(size='large', weight='bold', \
                         stretch='normal', #family='sans-serif', \
                         style='normal', variant='normal')
    gc.colorbar.set_font(size='large', weight='bold', \
                      stretch='normal', #family='sans-serif', \
                      style='normal', variant='normal')
    gc.colorbar.set_axis_label_font(size='x-large', weight='bold')
    rename = dataname.replace('.', '_')
    gc.add_label(0.06, 0.9, dataname, relative=True, size='24', 
            color='white', horizontalalignment='left')
    gc.add_label(0.94, 0.9, 'Lens Grade: '+lensgrade, relative=True, size='24', 
            color='white', horizontalalignment='right')
    gc.add_label(0.94, 0.1, 'SPIRE '+herschel_band, relative=True, size='24', 
            color='white', horizontalalignment='right')
    savefig('../Figures/overlays/' + rename + '_870_' + herschel_band + '.pdf')
    print(dataname)
    continue

    # compute the (x, y) center and pixel scale in the optical image
    #wcs_optical = wcs.WCS(optical_header, naxis=2)
    #pixxy = wcs_optical.wcs_world2pix(ra_center, dec_center, 1)
    #x_optical = numpy.round(pixxy[0])
    #y_optical = numpy.round(pixxy[1])
    #if (opttel == 'IRAC') | (opttel == 'CFHT') | (opttel == 'CTIO'):
    #    cdelt1_optical = numpy.abs(optical_header['CDELT1'] * 3600)
    #    cdelt2_optical = numpy.abs(optical_header['CDELT2'] * 3600)
    #else:
    #    cdelt1_optical = numpy.abs(optical_header['CD1_1'] * 3600)
    #    cdelt2_optical = numpy.abs(optical_header['CD2_2'] * 3600)
    #    cd11 = optical_header['CD1_1']
    #    cd12 = optical_header['CD1_2']
    #    cd21 = optical_header['CD2_1']
    #    cd22 = optical_header['CD2_2']
    #    if cd12 == 0:
    #        cd12 = cd11 / 1e8
    #    cdratio = numpy.abs(cd11 / cd12)
    #    if cdratio < 1:
    #        cdratio = 1 / cdratio
    #    if cdratio < 1e2:
    #        print "This shit ain't rotated yo!"
    #        import pdb; pdb.set_trace()
    #x_extent = numpy.round(radial_extent / cdelt1_optical)
    #y_extent = numpy.round(radial_extent / cdelt2_optical)
    #nx = optical_image[:, 0].size
    #ny = optical_image[0, :].size
    
    ## make the trimmed optical image
    ##if dataname == 'XMM109':
    ##    optical_trimmed = numpy.ones([2 * y_extent, 2 * x_extent])
    ##else:
    #if (x_optical - x_extent < 0) | (x_optical + x_extent > nx) | \
    #        (y_optical - y_extent < 0) | (y_optical + y_extent > ny):
    #    # pad with zeros
    #    optical_trimmed = numpy.zeros([2 * y_extent, 2 * x_extent])
    #    nx_pad = optical_trimmed[0, :].size
    #    ny_pad = optical_trimmed[:, 0].size
    #    if x_optical - x_extent < 0:
    #        xr0 = 0
    #        xrr0 = x_extent - x_optical
    #    else:
    #        xr0 = x_optical - x_extent
    #        xrr0 = 0
    #    if y_optical - y_extent < 0:
    #        yr0 = 0
    #        yrr0 = y_extent - y_optical
    #    else:
    #        yr0 = y_optical - y_extent
    #        yrr0 = 0
    #    if x_optical + x_extent > nx:
    #        xr1 = nx
    #        xrr1 = nx_pad / 2 + (nx - x_optical)
    #    else:
    #        xr1 = x_optical + x_extent
    #        xrr1 = nx_pad
    #    if y_optical + y_extent > ny:
    #        yr1 = ny
    #        yrr1 = ny_pad / 2 + (ny - y_optical)
    #    else:
    #        yr1 = y_optical + y_extent
    #        yrr1 = ny_pad
    #    optical_trimmed[yrr0:yrr1, xrr0:xrr1] = optical_image[yr0:yr1, xr0:xr1]
    #else:
    #    xr0 = x_optical - x_extent
    #    yr0 = y_optical - y_extent
    #    xr1 = x_optical + x_extent
    #    yr1 = y_optical + y_extent
    #    optical_trimmed = optical_image[yr0:yr1, xr0:xr1]

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
    print iauname, dataname, rms*1e3

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

    # grab image resolution from header data
    cdelt1_herschel = numpy.abs(herschel_header['CD1_1'] * 3600)
    cdelt2_herschel = numpy.abs(herschel_header['CD2_2'] * 3600)

    # compute the (x, y) center and pixel scale in the Herschel image
    wcs_herschel = wcs.WCS(herschel_header, naxis=2)
    pixxy = wcs_herschel.wcs_world2pix(ra_center, dec_center, 1)
    x_center = pixxy[0]
    y_center = pixxy[1]

    # re-sample the Herschel image to match the resolution of the ALMA image
    nx_herschel = herschel_image[0, :].size
    ny_herschel = herschel_image[:, 0].size
    x_extent_herschel = nx_herschel / 2 * cdelt1_herschel
    y_extent_herschel = ny_herschel / 2 * cdelt2_herschel
    x_vector = (nx_herschel / 2 - numpy.arange(nx_herschel)) * cdelt1_herschel
    y_vector = (numpy.arange(ny_herschel) - ny_herschel / 2) * cdelt2_herschel
    xex1 = -1. * x_extent_herschel
    xex2 = x_extent_herschel
    xexstep = -1j * nx_herschel
    yex1 = -1. * y_extent_herschel
    yex2 = y_extent_herschel
    yexstep = -1j * ny_herschel
    xold,yold = numpy.mgrid[xex1:xex2:xexstep, yex1:yex2:yexstep]
    res_ratio = cdelt1_herschel / cdelt1_ALMA
    xexstep = -1j * nx_herschel * res_ratio
    yexstep = -1j * ny_herschel * res_ratio
    xnew,ynew = numpy.mgrid[xex1:xex2:xexstep, yex1:yex2:yexstep]
    tck = interpolate.bisplrep(xold, yold, herschel_image, s=0)
    herschel_resample = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    nx_resample = herschel_resample[0, :].size
    ny_resample = herschel_resample[:, 0].size
    xr0 = numpy.round(x_center * res_ratio - radial_extent / cdelt1_ALMA)
    yr0 = numpy.round(y_center * res_ratio - radial_extent / cdelt2_ALMA)
    xr1 = numpy.round(x_center * res_ratio + radial_extent / cdelt1_ALMA)
    yr1 = numpy.round(y_center * res_ratio + radial_extent / cdelt2_ALMA)

    # make sure the trimmed image is square
    if xr1 - xr0 != yr1 - yr0:
        dx = yr1 - yr0 - (xr1 - xr0)
        xr1 += dx
        print dx

    # make the trimmed sub-mm image
    herschel_trimmed = herschel_resample[yr0:yr1, xr0:xr1]

    # make sure it is the same shape as the trimmed ALMA image
    if herschel_trimmed.shape != ALMA_trimmed.shape:
        herschel_trimmed = herschel_resample[yr0:yr1 - 1, xr0:xr1 - 1]
    if herschel_trimmed.shape != ALMA_trimmed.shape:
        herschel_trimmed = herschel_resample[yr0:yr1 + 1, xr0:xr1 + 1]
    if herschel_trimmed.shape != ALMA_trimmed.shape:
        pdb.set_trace()


    # clear the plotting window
    plt.clf()

    # grab the axes instance
    ax = fig.add_subplot(1, 1, 1)
    #ax = fig.gca()

    # print the target name in the upper left corner
    plt.text(0.05, 0.90, iauname, fontsize='xx-large', \
        color='black', transform=ax.transAxes)
    plt.text(0.05, 0.80, dataname, fontsize='xx-large', \
        color='black', transform=ax.transAxes)
    #plt.text(0.05, 0.05, opticalinfo, fontsize='x-large', \
    #    color='black', transform=ax.transAxes)
    #plt.subplot(1, 1, 1)

    extent = [radial_extent, -radial_extent, -radial_extent, radial_extent]
    #optical_trimmed = optical_trimmed + optical_trimmed.min() + 1
    #optical_trimmed = numpy.log10(optical_trimmed)
    #goodpix = (numpy.abs(optical_trimmed) < 1e5) & \
    #    (numpy.abs(optical_trimmed) > 0)
    #pospix = (optical_trimmed < 1e5) & \
    #    (optical_trimmed > 0)
    #ALMA_trimmed = numpy.log10(ALMA_trimmed)
    vmax = ALMA_trimmed.max()
    vmin = ALMA_trimmed.min()
    #imsort = numpy.sort(optical_trimmed[pospix]).flatten()
    #npix = optical_trimmed[pospix].size
    #vmax = imsort[0.99 * npix]
    #if dataname == 'XMM02':
    #    vmax = imsort[0.95 * npix]
    #if dataname == 'XMM03':
    #    vmax = imsort[0.999 * npix]
    #if dataname == 'XMM101':
    #    vmax = imsort[0.95 * npix]
    #if dataname == 'XMM108':
    #    vmax = imsort[0.999 * npix]
    #if dataname == 'XMM109':
    #    vmax = imsort[0.9 * npix]
    #vmin = imsort[0.5 * npix]
    plt.imshow(ALMA_trimmed, cmap='jet', interpolation='nearest', \
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
    plevs = [4 * rms]
    nlevs = [-4 * rms]
    #plevs = 2 * rms * 2 ** (numpy.arange(10) / 2.)
    #nlevs = sorted(-2 * rms * 2 ** (numpy.arange(4) / 2.))
    #pcline = 'solid'
    #ncline = 'dashed'

    # draw the contours
    #plt.contour(x_vector, y_vector, ALMA_trimmed, colors='red', levels=plevs, \
    #    linewidths=1.0)
    #plt.contour(x_vector, y_vector, ALMA_trimmed, colors='red', levels=nlevs, \
    #    linewidths=1.0)

    # define the x- and y-vectors for the Herschel contour plot
    nx_herschel = herschel_resample[0, :].size
    ny_herschel = herschel_resample[:, 0].size
    x_vector = (nx_ALMA / 2 - numpy.arange(nx_ALMA)) * cdelt1_ALMA
    y_vector = (numpy.arange(ny_ALMA) - ny_ALMA / 2) * cdelt2_ALMA
    #cellplus = celldata*(2*xrad+1.1)/(2*xrad)
    #cmodx = ( numpy.arange(2*xrad) - xrad ) * (-cellplus) - celldata/2.
    #cmody = ( numpy.arange(2*yrad) - yrad ) * cellplus + celldata/2.

    # define contour level parameters
    #plevs = numpy.array([0.6,0.85]) * herschel_trimmed.max()
    #plevs = [0.04, 0.080, .12]
    plevs = sorted(4 * 0.007 * 2 ** (numpy.arange(8) / 2.))
    #nlevs = sorted(-2 * rms * 2 ** (numpy.arange(4) / 2.))
    #pcline = 'solid'
    #ncline = 'dashed'

    # draw the contours
    plt.contour(x_vector, y_vector, herschel_trimmed, colors='black', levels=plevs, \
        linewidths=2.5, smooth=0, linestyles='solid')
    plt.contour(x_vector, y_vector, herschel_trimmed, colors='white', levels=plevs, \
        linewidths=1., smooth=0, linestyles='solid')
    #plt.contour(x_vector, y_vector, herschel_trimmed, colors='red', levels=nlevs, \
    #    linewidths=1.5)

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
        ec='white', lw=1.0, transform=ax.transAxes, fc='white',
        zorder=10)
    ax.add_artist(e)
    #plt.plot([-10,-8],[10,10], color='black')
    #ax.set_xticklabels([])
    #ax.set_yticklabels([])
    #plt.tight_layout()

    # save the figure to a file
    rename = dataname.replace('.', '_')
    savefig('../Figures/overlays/' + rename + '_870_' + herschel_band)
    pdb.set_trace()


# rename files to remove '.' from the name
pdb.set_trace()
