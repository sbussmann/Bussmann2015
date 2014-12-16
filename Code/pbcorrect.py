"""

2014 November 19

Shane Bussmann

Correct fluxes in source_observed.dat to account for primary beam fall-off.

"""

from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import yaml
from astropy.io import fits
import numpy


observedloc = '../Data/table_observed.dat'
observed = ascii.read(observedloc)
nobs = len(observed)

intrinsicloc = '../Data/table_intrinsic.dat'
intrinsic = ascii.read(intrinsicloc)

uvfitlistloc = '../Data/uvfitlist.dat'
uvfitlist = ascii.read(uvfitlistloc)

fitloc = '../../../../ModelFits/'

for iobs in range(nobs):
    ra = coord.Angle(observed['ra_alma'][iobs], unit=u.hour).degree
    dec = coord.Angle(observed['dec_alma'][iobs], unit=u.deg).degree

    target = observed['target'][iobs]

    match = uvfitlist['dataname'] == target
    fitdir = uvfitlist['observed'][match][0]

    fullfitloc = fitloc + target + '/' + fitdir + '/'

    configloc = fullfitloc + 'config.yaml'
    configfile = open(configloc)
    config = yaml.load(configfile)

    imloc = fullfitloc + config['ImageName']
    imhead = fits.getheader(imloc)

    ra_center = imhead['CRVAL1']
    dec_center = imhead['CRVAL2']

    offra = (ra - ra_center) * numpy.cos(dec_center * numpy.pi / 180)
    offdec = dec - dec_center
    offset = numpy.sqrt(offra ** 2 + offdec ** 2) * 3600


    diameter = 12.
    lambda_obs = 870e-6
    fwhm = 1.2 * lambda_obs / diameter * 206265
    halfpower = fwhm / 2.355

    pbcorr = numpy.exp(-0.5 * offset ** 2 / halfpower **2)
    #import pdb; pdb.set_trace()

    observedflux = observed['f870'][iobs]
    e_observedflux = observed['e_f870'][iobs]
    actualflux = observedflux / pbcorr
    e_actualflux = e_observedflux / pbcorr

    mu1 = 1.0
    mu2 = observed['mu870'][iobs]
    remu = (mu1 + mu2) / 2.
    e_remu = (mu2 - remu) / 2.

    if remu > 3:
        remu = mu2
        e_remu = observed['e_mu870'][iobs]

    intrinsicflux = intrinsic['fnu'][iobs]
    e_intrinsicflux = intrinsic['e_fnu'][iobs]

    msg = '{0:10} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f} {5:5.2f}'
    msgfmt = msg.format(target, pbcorr, actualflux, e_actualflux, remu, e_remu)
    print(msgfmt)
