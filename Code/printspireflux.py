"""

2015 March 31

Shane Bussmann

Print SPIRE flux densities from latest Monte Carlo simulations on Herschel maps
using ALMA and Spitzer positions as priors (fiducial), from lens candidate
catalog (StarFinder), and from original catalog used to select targets for ALMA
observations (SUSSEXtractor).

"""

from astropy.io import ascii


# read in fiducial flux densities
fiducial = ascii.read('../Data/targetlist.dat')

# read in StarFinder flux densities
starfinder = ascii.read('../Data/targetlist_starfinder.dat')

# read in SUSSEXtractor flux densities
sussex = ascii.read('../Data/targetlist_sussextractor.dat')

ntarg = len(fiducial)

for i in range(ntarg):
    ifid = fiducial[i]
    isf = starfinder[i]
    isussex = sussex[i]
    strfmt = '{0:12} {1:5.0f} {2:5.0f} {3:5.0f} {4:5.0f} ' + \
            '{5:5.0f} {6:5.0f} {7:5.0f} {8:5.0f} {9:5.0f}'
    print(strfmt.format(ifid['shortname'], 
        ifid['f250'], ifid['f350'], ifid['f500'],
        isf['f250'], isf['f350'], isf['f500'],
        isussex['f250'], isussex['f350'], isussex['f500']))
