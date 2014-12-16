"""
2014 March 27
Shane Bussmann

Generate input files for SExtractor and run it.
"""

from astropy.table import Table
from astropy.io import fits
import os
import numpy


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

    counter = 0
    for ifilt in filterlist:

        # load the appropriate optical image
        imloc = dataname + '_' + ifilt + '_cut.fits'
        if not os.path.exists(imloc):
            continue
        header = fits.getheader(imloc)

        # estimate the photometric zeropoint of the data
        gemzero = [24.91, 28.33, 28.33, 27.93, 26.84]
        gemcolors = [0.38, 0.18, 0.1, 0.08, 0.05]
        headerkeys = header.keys()
        cautionflag = 1
        #if headerkeys.count('AIRMASS') > 0:
        #    gain = header['GAIN']
        #    exptime = header['EXPTIME']
        #    airmass = header['AIRMASS']
        #    zpt = gemzero[counter] - 2.5 * numpy.log10(gain / exptime) - \
        #        gemcolors[counter] * (airmass - 1.0)
        #    cautionflag = 0
        if headerkeys.count('PHOTZP') > 0:
            zpt = header['PHOTZP']
            cautionflag = 0
        if headerkeys.count('MAGZPT') > 0:
            zpt = header['MAGZPT']
            cautionflag = 0
        if cautionflag:
            # go back to Gemini reduction directory to get needed header info
            dataloc = '/Users/rbussman/Data/Gemini-South/'
            oldimloc = dataloc + dataname + '_' + ifilt + '.fits'
            oldheader = fits.getheader(oldimloc)
            gain = oldheader['GAIN']
            exptime = oldheader['EXPTIME']
            airmass = oldheader['AIRMASS']
            zpt = gemzero[counter] - 2.5 * numpy.log10(gain / exptime) - \
                gemcolors[counter] * (airmass - 1.0)
            cautionflag = 0

        string0 = ' -MAG_ZEROPOINT ' + str(zpt)

        catalogname = dataname + '_' + ifilt + '_cut.cat'
        string1 = ' -CATALOG_NAME ' + catalogname
        segname = dataname + '_' + ifilt + '_seg.fits'
        string2 = ' -CHECKIMAGE_NAME ' + segname

        stringfinal = ' ' + imloc

        cmd = 'sex' + string0 + string1 + string2 + stringfinal
        print cmd
        os.system(cmd)

        counter += 1

