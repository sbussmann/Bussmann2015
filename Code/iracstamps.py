"""

2014 November 13

Shane Bussmann

Overlay ALMA contours on IRAC 3.6um, 4.5um, 8.0um 3-color image.

"""

from astropy.table import Table
#import aplpy
#from PIL import Image
import numpy
from astropy.io import fits
from astropy import wcs
from pylab import savefig
import img_scale
import yaml
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Ellipse
import matplotlib


def transform(imloc, ra_center, dec_center, radial_extent, vmax=1.5):

    """ Rescale and trim input image. """

    hdu = fits.open(imloc)
    im = hdu[0].data
    thisisalma = False
    if im.ndim == 4:
        thisisalma = True
        im = im[0, 0, :, :]
    optical_header = hdu[0].header
    bad = im * 0 != 0
    good = im * 0 == 0
    im[bad] = im[good].min()

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
    nx = im[:, 0].size
    ny = im[0, :].size

    #angularradius = 0.5 * 60
    #pixscale = numpy.sqrt(cdelt1_optical * cdelt2_optical)
    #pixradius = numpy.round(angularradius / pixscale).astype(int)
    #x1 = str(x_optical.astype(int) - pixradius)
    #x2 = str(x_optical.astype(int) + pixradius)
    #y1 = str(y_optical.astype(int) - pixradius)
    #y2 = str(y_optical.astype(int) + pixradius)
    #region = '[' + x1 + ':' + x2 + ',' + y1 + ':' + y2 + ']'
    #cutfile = 'fitscutouts/' + dataname + '_' + ifilt + '_cut.fits'
    #cmd = 'imcopy ' + optloc + region + ' ' + cutfile
    #print cmd
    
    # make the trimmed optical image
    #if dataname == 'XMM109':
    #    optical_trimmed = numpy.ones([2 * y_extent, 2 * x_extent])
    #else:
    if (x_optical - x_extent < 0) | (x_optical + x_extent > nx) | \
            (y_optical - y_extent < 0) | (y_optical + y_extent > ny):
        # pad with zeros
        import pdb; pdb.set_trace()
        trimmed = numpy.zeros([2 * y_extent, 2 * x_extent])
        nx_pad = trimmed[0, :].size
        ny_pad = trimmed[:, 0].size
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
        trimmed[yrr0:yrr1, xrr0:xrr1] = im[yr0:yr1, xr0:xr1]
    else:
        xr0 = x_optical - x_extent
        yr0 = y_optical - y_extent
        xr1 = x_optical + x_extent
        yr1 = y_optical + y_extent
        trimmed = im[yr0:yr1, xr0:xr1]
    #print(vmax)
    if not thisisalma:
        trimmed -= trimmed.min()# - 1
        trimmed = img_scale.sqrt(trimmed, scale_max=vmax)
        #trimmed *= 2.5
        #trimmed = numpy.sqrt(trimmed)
        #trimmed -= trimmed.min()
        #trimmedsort = numpy.sort(trimmed).flatten()
        #ntrimmed = trimmedsort.size
        #vmax = trimmedsort[0.99999 * ntrimmed]
        #vmin = trimmedsort[0.5 * ntrimmed]
        #print(vmin, vmax)
        #toohigh = trimmed > vmax
        #toolow = trimmed < vmin
        #trimmed[toolow] = vmin
        #trimmed[toohigh] = vmax
    return trimmed

# blue color is IRAC channel 1: 3.6um
blue = 'irac1'

# green is 4.5um
green = 'irac2'

# red is 8.0um
red = 'irac4'

# set font properties
font = {'family' : 'Arial Narrow',
        'weight' : 'bold',
        'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

iractargetloc = '../Data/iractargetlist.txt'
iractargetlist = Table.read(iractargetloc, format='ascii')

goodfitloc = '../Data/uvfitlist.dat'
goodfitdat = Table.read(goodfitloc, format='ascii')

ntarget = len(iractargetlist)

for itarget in range(ntarget):

    # get the appropriate image center and radial extent
    target = iractargetlist['target'][itarget]
    match = goodfitdat['dataname'] == target
    goodfit = goodfitdat['intrinsic'][match][0]
    modfitloc = '../../../../ModelFits/' + target + '/' + goodfit
    configloc = modfitloc + '/config.yaml'
    configfile = open(configloc)
    config = yaml.load(configfile)

    ra_center = config['Region0']['RACentroid']
    dec_center = config['Region0']['DecCentroid']
    radial_extent = 9.5#config['Region0']['RadialExtent']

    # Make the RGB image
    imdir = '../../fitscutouts/' + target
    blueimloc = imdir + '_' + blue + '.fits'
    blueim = transform(blueimloc, ra_center, dec_center, radial_extent,
            vmax=0.7)
    greenimloc = imdir + '_' + green + '.fits'
    greenim = transform(greenimloc, ra_center, dec_center, radial_extent,
            vmax=0.7)
    redimloc = imdir + '_' + red + '.fits'
    redim = transform(redimloc, ra_center, dec_center, radial_extent)
    optical_header = fits.getheader(blueimloc)
    naxis = redim[:, 0].size
    cubefits = '../Figures/' + target + '_rgb.fits'
    rgbArray = numpy.zeros((naxis, naxis, 3))
    rgbArray[:, :, 0] = redim
    rgbArray[:, :, 1] = greenim
    rgbArray[:, :, 2] = blueim
    
    plt.clf()
    fig = plt.figure(figsize=(3.0, 3.0))
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(left=0.08, right=0.97, top=0.97, 
            bottom=0.08, wspace=0.35)

    # load the corresponding ALMA image
    submm_band = '870'
    ALMA_file = '../../fitscutouts/' + target + '_' + submm_band + '.fits'
    ALMA_hdu = fits.open(ALMA_file)
    ALMA_fullimage = ALMA_hdu[0].data
    ALMA_trimmed = transform(ALMA_file, ra_center, dec_center, radial_extent)
    ALMA_header = ALMA_hdu[0].header
    bmaj = ALMA_header['BMAJ'] * 3600
    bmin = ALMA_header['BMIN'] * 3600
    bpa = ALMA_header['BPA']
    cdelt1_ALMA = numpy.abs(ALMA_header['CDELT1'] * 3600)
    cdelt2_ALMA = numpy.abs(ALMA_header['CDELT2'] * 3600)
    wcs_ALMA = wcs.WCS(ALMA_header, naxis=2)
    pixxy = wcs_ALMA.wcs_world2pix(ra_center, dec_center, 1)
    x_ALMA = numpy.round(pixxy[0])
    y_ALMA = numpy.round(pixxy[1])
    mask = ALMA_fullimage.copy()
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
    rms = ALMA_fullimage[goodregion].std()

    extent = [radial_extent, -radial_extent, -radial_extent, radial_extent]
    plt.imshow(rgbArray, interpolation='nearest', \
        extent=extent, origin='lower')

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
    plevs = 4 * rms * 2 ** (numpy.arange(10))
    nlevs = sorted(-4 * rms * 2 ** (numpy.arange(4)))
    #pcline = 'solid'
    #ncline = 'dashed'

    # draw the contours
    plt.contour(x_vector, y_vector, ALMA_trimmed, colors='white', levels=plevs, \
        linewidths=1.0)
    plt.contour(x_vector, y_vector, ALMA_trimmed, colors='white', levels=nlevs, \
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
        ec='black', lw=0.1, transform=ax.transAxes, fc='white',
        zorder=10)
    ax.add_artist(e)
    plt.text(0.95, 0.95, target, fontsize='xx-large', \
        color='white', va='top', ha='right', transform=ax.transAxes)
    

    #imshow(rgbArray, origin='lower')
    savefig('../Figures/' + target + '_rgb.pdf')
    #import pdb; pdb.set_trace()
    #aplpy.make_rgb_cube([redim, greenim, blueim], cubefits)
    #import pdb; pdb.set_trace()
    #cubeim = '../Figures/' + target + '_rgb.png'
    #aplpy.make_rgb_image([redimloc, greenimloc, blueimloc], cubeim, 
    #        stretch_r='sqrt', stretch_g='sqrt', stretch_b='sqrt', 
    #        pmax_r=99, pmax_g=99, pmax_b=99)
    #gc = aplpy.FITSFigure(redimloc)
    #gc.show_rgb(cubeim)
    #gc.save('../Figures/' + target + '_test.png')
    #import pdb; pdb.set_trace()
