"""

2014 July 1

Shane Bussmann

Copy best-fit model and residual images from uvfit* directory to modelfit/

"""

from astropy.table import Table
import glob
import os


# store the current directory
cwd = os.getcwd()

# the uvfit directory containing the best-fit model
fitlist = Table.read('../Data/uvfitlist.dat', format='ascii')

# identify ALMA targets
targetloc = '../Data/targetlist.dat'
targets = Table.read(targetloc, format='ascii')
ntargets = len(targets)

i = 28

for itarg in range(i, i + 1):

    dataname = targets[itarg]['dataname']
    shortname = targets[itarg]['shortname']

    fitdir = fitlist['intrinsic'][itarg]
    fitloc = '../../../../ModelFits/' + dataname + '/' + fitdir + '/'

    #print(dataname + '/' + fitdir + '/config.yaml')

    # change into that directory and run visualize.bestFit()
    os.chdir(fitloc)
    import visualize
    visualize.bestFit(interactive=False, showOptical=True)
    os.chdir(cwd)

    pdfmap = 'LensedSBmap.optical.bestfit.pdf'
    origin = fitloc + pdfmap
    destination = '../Figures/modelfit/' + shortname + '.optical.bestfit.pdf'
    #destination = destination.replace('.', '_', 2)
    cmd = 'cp ' + origin + ' ' + destination
    print(cmd)
    os.system(cmd)

    pdfmap = 'LensedSBmap.model.bestfit.pdf'
    origin = fitloc + pdfmap
    destination = '../Figures/modelfit/' + shortname + '.model.bestfit.pdf'
    #destination = destination.replace('.', '_', 2)
    cmd = 'cp ' + origin + ' ' + destination
    print(cmd)
    os.system(cmd)

    pdfmap = 'LensedSBmap.residual.bestfit.pdf'
    origin = fitloc + pdfmap
    destination = '../Figures/modelfit/' + shortname + '.residual.bestfit.pdf'
    #destination = destination.replace('.', '_', 2)
    cmd = 'cp ' + origin + ' ' + destination
    print(cmd)
    os.system(cmd)
