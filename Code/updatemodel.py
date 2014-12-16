"""

2014 July 1

Shane Bussmann

Copy best-fit model and residual images from uvfit* directory to modelfit/

"""

from astropy.table import Table
import glob
import os


# the uvfit directory containing the best-fit model
fitlist = Table.read('uvfitlist.dat', format='ascii')

# identify ALMA targets
targetloc = 'targetlist.dat'
targets = Table.read(targetloc, format='ascii')
ntargets = len(targets)

for itarg in range(ntargets):

    target = targets[itarg]['dataname']

    fitdir = fitlist['intrinsic'][itarg]
    fitloc = '../../ModelFits/' + target + '/' + fitdir + '/'

    pdfmap = 'LensedSBmap.optical.bestfit.pdf'
    origin = fitloc + pdfmap
    destination = 'modelfit/' + target + '.optical.bestfit.pdf'
    destination = destination.replace('.', '_', 2)
    cmd = 'cp ' + origin + ' ' + destination
    print(cmd)
    os.system(cmd)

    pdfmap = 'LensedSBmap.model.bestfit.pdf'
    origin = fitloc + pdfmap
    destination = 'modelfit/' + target + '.model.bestfit.pdf'
    destination = destination.replace('.', '_', 2)
    cmd = 'cp ' + origin + ' ' + destination
    print(cmd)
    os.system(cmd)

    pdfmap = 'LensedSBmap.residual.bestfit.pdf'
    origin = fitloc + pdfmap
    destination = 'modelfit/' + target + '.residual.bestfit.pdf'
    destination = destination.replace('.', '_', 2)
    cmd = 'cp ' + origin + ' ' + destination
    print(cmd)
    os.system(cmd)
