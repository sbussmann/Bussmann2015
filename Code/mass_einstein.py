"""
Created 2013 March 13
Shane Bussmann

Purpose: Returns mass and velocity dispersion inside Einstein radius given an
Einstein radius (in arcsec), source redshift, and lens redshift

Inputs must be in numpy array format!

"""

from math import pi
from astropy.cosmology import WMAP9 as cosmo
import numpy

def me(theta_Es, e_theta_Es, zlenses, zsources):

    ntarg = theta_Es.size
    M_E = numpy.zeros(ntarg)
    e_M_E = numpy.zeros(ntarg)
    vdisp = numpy.zeros(ntarg)
    e_vdisp = numpy.zeros(ntarg)

    for targ in numpy.arange(ntarg):

        zsource = zsources[targ]
        zlens = zlenses[targ]
        theta_E = theta_Es[targ]
        e_theta_E = e_theta_Es[targ]

        # luminosity distances
        d_LS = cosmo.luminosity_distance(zsource).value * 3.08e24
        d_LL = cosmo.luminosity_distance(zlens).value * 3.08e24

        # comoving distances
        d_MS = d_LS / (1 + zsource)
        d_ML = d_LL / (1 + zlens)

        # angular diameter distances
        d_ALS = 1 / (1 + zsource) * ( d_MS - d_ML )
        d_AL = d_LL / (1 + zlens)**2
        d_AS = d_LS / (1 + zsource)**2

        # einstein radius in cm (7.1 kpc/" at z=0.7)
        theta_E_cm = theta_E / 206265. * d_AL
        e_theta_E_cm = e_theta_E / 206265. * d_AL

        # get a distribution of Einstein radii
        niters = 1e3
        theta_E_iters = numpy.random.normal(loc = theta_E_cm, \
                scale = e_theta_E_cm, size = niters)

        # compute the mass enclosed within the Einstein radius
        c = 3e10
        G = 6.67e-8
        sigma_crit = c**2 / 4 / pi / G * d_AS / d_AL / d_ALS
        M_E_iters = pi * sigma_crit * theta_E_iters**2 / 2e33
        M_E[targ] = numpy.mean(M_E_iters)
        e_M_E[targ] = numpy.std(M_E_iters)

        vdisp2 = theta_E_iters / d_AL / 4 / pi * c**2 * d_AS / d_ALS
        vdisp[targ] = numpy.mean(numpy.sqrt(vdisp2) / 1e5)
        e_vdisp[targ] = numpy.std(numpy.sqrt(vdisp2) / 1e5)

    return M_E, e_M_E, vdisp, e_vdisp
