"""

Shane Bussmann

2014 August 21

Various utilities for calculating angular separations.

"""


import matplotlib.pyplot as plt
import numpy
from astropy import units as u
from astropy.coordinates import SkyCoord


def rmSingles(fluxcomponent, targetstring='target'):

    """

    Filter out targets in fluxcomponent that have only one ALMA source.

    """

    nindiv = len(fluxcomponent)

    flagger = numpy.zeros(nindiv)
    for icomp in range(nindiv):
        target = fluxcomponent[targetstring][icomp]

        match = fluxcomponent[targetstring] == target
        nmatch = fluxcomponent[targetstring][match].size
        if nmatch == 1:
            flagger[icomp] = 1

    goodflag = flagger == 0
    fluxcomponent = fluxcomponent[goodflag]

    return fluxcomponent

def setThresh(fluxcomponent, thresh, fluxstring='f880'):

    """

    Filter out targets in fluxcomponent below a given `thresh` flux density.

    """

    if fluxstring == 'peakflux':
        thresh /= 1e3
    fluxhigh = fluxcomponent[fluxstring] > thresh
    fluxcomponent = fluxcomponent[fluxhigh]

    return fluxcomponent


def getSeparation(fluxcomponent, rastring='offra', decstring='offdec',
        targetstring='target', fluxstring='f880', degrees=False):

    nindiv = len(fluxcomponent)

    separation = []
    ra_vector = []
    dec_vector = []
    fluxweighted = []
    for icomp in range(nindiv):
        target = fluxcomponent[targetstring][icomp]

        match = fluxcomponent[targetstring] == target
        ras = fluxcomponent[rastring][match]
        decs = fluxcomponent[decstring][match]
        fluxs = fluxcomponent[fluxstring][match]
        nmatch = fluxcomponent[fluxstring][match].size
        try:
            if degrees:
                factor = 1.
            else:
                factor = 3600.
            avgra = ras.mean() / factor
            avgdec = decs.mean() / factor
            ra = fluxcomponent[rastring][icomp] / factor
            dec = fluxcomponent[decstring][icomp] / factor
            newra = ras / factor
            newdec = decs / factor
        except:
            newra = numpy.zeros(nmatch)
            newdec = numpy.zeros(nmatch)
            for imatch in range(nmatch):
                ira = ras[imatch]
                idec = decs[imatch]
                sc = SkyCoord(ira, idec, "icrs", unit=(u.hourangle, u.degree))
                newra[imatch] = sc.ra.deg
                newdec[imatch] = sc.dec.deg

            ra = fluxcomponent[rastring][icomp]
            dec = fluxcomponent[decstring][icomp]
            sc = SkyCoord(ra, dec, "icrs", unit=(u.hourangle, u.degree))
            ra = sc.ra.deg
            dec = sc.dec.deg
            avgra = newra.mean()
            avgdec = newdec.mean()
        
        raoffset = (ra - avgra) * numpy.cos(dec)
        decoffset = dec - avgdec
        offset = numpy.sqrt(raoffset ** 2 + decoffset ** 2) * 3600
        separation.append(offset)
        ra_vector.append(raoffset)
        dec_vector.append(decoffset)

        wmeanra_num = newra * fluxs
        wmeanra_den = fluxs
        wmeanra = wmeanra_num.sum() / wmeanra_den.sum()
        wmeandec_num = newdec * fluxs
        wmeandec_den = fluxs
        wmeandec = wmeandec_num.sum() / wmeandec_den.sum()
        raoffset = (ra - wmeanra) * numpy.cos(dec)
        decoffset = dec - wmeandec
        offset = numpy.sqrt(raoffset ** 2 + decoffset ** 2) * 3600
        fluxweighted.append(offset)

    return numpy.array(separation), numpy.array(fluxweighted), ra_vector, \
            dec_vector

def simPosition(fluxcomponent, distance=8.5, rastring='offra',
        decstring='offdec', targetstring='target'):

    """

    Replace the RA and Dec values for each target in fluxcomponent with random
    positions that are forced to be within "distance" arcseconds of the average
    position of the sources in that object.  Default distance is 8.5", aka one
    ALMA FOV.

    """

    from numpy.random import uniform

    nindiv = len(fluxcomponent)

    for icomp in range(nindiv):
        target = fluxcomponent[targetstring][icomp]
        match = fluxcomponent[targetstring] == target
        ras = fluxcomponent[rastring][match]
        decs = fluxcomponent[decstring][match]
        nmatch = fluxcomponent[rastring][match].size
        try:
            avgra = ras.mean() / 3600
            avgdec = decs.mean() / 3600
            ra = fluxcomponent[rastring][icomp] / 3600
            dec = fluxcomponent[decstring][icomp] / 3600
            newra = ras / 3600
            newdec = decs / 3600
            degflag = True
        except:
            degflag = False
            newra = numpy.zeros(nmatch)
            newdec = numpy.zeros(nmatch)
            for imatch in range(nmatch):
                ira = ras[imatch]
                idec = decs[imatch]
                sc = SkyCoord(ira, idec, "icrs", unit=(u.hourangle, u.degree))
                newra[imatch] = sc.ra.deg
                newdec[imatch] = sc.dec.deg

            ra = fluxcomponent[rastring][icomp]
            dec = fluxcomponent[decstring][icomp]
            sc = SkyCoord(ra, dec, "icrs", unit=(u.hourangle, u.degree))
            ra = sc.ra.deg
            dec = sc.dec.deg
            avgra = newra.mean()
            avgdec = newdec.mean()

        lowra = avgra - distance / 3600 / numpy.cos(avgdec)
        highra = avgra + distance / 3600 / numpy.cos(avgdec)
        lowdec = avgdec - distance / 3600
        highdec = avgdec + distance / 3600
        if lowdec < -0.025:
            lowdec = -0.025
        if highdec > 0.025:
            highdec = 0.025
        if degflag:
            randomra = uniform(low=lowra, high=highra) * 3600
            randomdec = uniform(low=lowdec, high=highdec) * 3600
        else:
            randomra = uniform(low=lowra, high=highra)
            randomdec = uniform(low=lowdec, high=highdec)
        #print(randomra, randomdec)
        #c = SkyCoord(ra=randomra * u.degree, dec=randomdec * u.degree)

        fluxcomponent[rastring][icomp] = randomra
        fluxcomponent[decstring][icomp] = randomdec
        #fluxcomponent[rastring][icomp] = c.ra.to_string('hourangle', sep=':')
        #fluxcomponent[decstring][icomp] = c.dec.to_string(sep=':')

    return fluxcomponent

def simArea(fluxcomponent, nsim, bin_edges, targetstring='target', 
        edgecolor='black', facecolor='none', hatch='', fluxstring='f880',
        norm=1.0, label=''):
    nbins = bin_edges.size - 1
    supersep = numpy.zeros([nsim, nbins])
    for isim in range(nsim):

        simP = simPosition(fluxcomponent, targetstring=targetstring)
        avgsep, wmeansep, ra1, dec1 = getSeparation(simP, 
                targetstring=targetstring, fluxstring=fluxstring)
        hist, bin_edges = numpy.histogram(avgsep, bins=bin_edges)
        area1 = numpy.pi * bin_edges ** 2
        area2 = numpy.pi * numpy.roll(bin_edges, -1) ** 2
        area = area2 - area1
        area = area[0:-1]
        xhist = (bin_edges + numpy.roll(bin_edges, -1)) / 2.
        xhist = xhist[0:-1]
        supersep[isim, :] = hist / area / norm

    mediansep = numpy.median(supersep, axis=0)
    stdsep = numpy.std(supersep, axis=0)
    uppersep = mediansep + stdsep
    lowersep = mediansep - stdsep

    plt.plot(xhist, mediansep, linestyle='--', color=edgecolor, label=label)
    plt.fill_between(xhist, lowersep, y2=uppersep, hatch=hatch,
            facecolor=facecolor, edgecolor=edgecolor)

    return mediansep, stdsep, uppersep, lowersep

def histArea(values, nbins, color='', norm=1.0, fmt='s', ms=6, linestyle='-',
        drawstyle='default', showerror=True, mew=0.5, label=''):
    hist, bin_edges = numpy.histogram(values, bins=nbins)
    area1 = numpy.pi * bin_edges ** 2
    area2 = numpy.pi * numpy.roll(bin_edges, -1) ** 2
    area = area2 - area1
    area = area[0:-1]
    xhist = (bin_edges + numpy.roll(bin_edges, -1)) / 2.
    xhist = xhist[0:-1]
    plt.plot(xhist, hist / area / norm, fmt, color=color, linestyle=linestyle,
            linewidth=2, drawstyle=drawstyle, ms=ms, mew=mew, label=label)
    #plt.fill_between(xhist, hist / area / norm, facecolor=color, alpha=0.2)
    #plt.plot(xhist, hist / area / norm, fmt, ms=ms, mew=mew, color=color)
    if showerror:
        plt.errorbar(xhist, hist / area / norm, yerr=numpy.sqrt(hist) / 
                area / norm, fmt=fmt, ms=ms, color=color, ecolor='gray', 
                capsize=0, mew=mew)
