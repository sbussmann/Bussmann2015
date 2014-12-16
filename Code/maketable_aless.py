"""

2014 September 24

Shane Bussmann

Take observed sources from ALMA sample and modify them to match ALESS sample
properties:

    1. The flux density of each source is reduced by a factor of 3
    2. Sources separated by less than 2" become a single source
    3. Sources with flux densities below 1.2 mJy are removed

"""

from astropy.table import Table

almasample = Table.read('table_positions.dat', format='ascii')

objname_previous = ''

nsources = len(almasample)

for isource in range(nsources):
    objname = almasample['target'][isource]

    if objname == objname_previous:
        match = almasample['target'] == objname
        ra_match = almasample['ra_alma'][match]
        dec_match = almasample['dec_alma'][match]

        # compute the separations between each source


    objname_previous = objname
