"""

Shane Bussmann

2015 March 30

Plot SEDs for each Herschel source using new Herschel photometry from J.
Scudders.

"""

from astropy.table import Table
import matplotlib.pyplot as plt
from pylab import savefig


def plot(filename, color='black', label=''):

    data = Table.read(filename, format='ascii')
    
    ndata = len(data)

    for i in range(ndata):
        plt.figure(i + 1)
        flux = [data['250med'][i], data['350med'][i], data['500med'][i]]
        wave = [250, 350, 500]
        plt.plot(wave, flux, color=color, label=label)


# plot the total 870um fluxes without 24um data
plot('../Data/totalfluxes/total_870_fluxes.txt', color='red', 
        label='Ignoring 24um sources')

# plot the total 870um fluxes with 24um data where 24um source > 2.5" from
# 870um source
plot('../Data/totalfluxes/total_870_24_fluxes.txt', color='blue', 
        label='Using 24um sources if not near 870um')

# plot the total 870um fluxes with 24um data
plot('../Data/totalfluxes/total_870_24_nodel_fluxes.txt', color='green', 
        label='Using all 24um sources')

datadummy = Table.read('../Data/totalfluxes/total_870_fluxes.txt', 
        format='ascii')

ndata = len(datadummy)
fieldname = datadummy['Field']
for i in range(ndata):
    plt.figure(i + 1)
    plt.legend(loc='best')
    plt.title(fieldname[i])
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Flux density (mJy)')
    savefig('../Figures/SEDplot_' + fieldname[i] + '.png')

