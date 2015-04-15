"""

Add a source ID column to the Cowley mock catalog.

"""

from astropy.io import ascii


dataloc = '../Data/SPIRE_ALMA_Cat_v4.txt'
data = ascii.read(dataloc)

nrows = len(data)

sourcevector = [0]
sourceID = 0
sourceRAprevious = data['SourceX'][0]
for i in range(1, nrows):
    sourceRA = data['SourceX'][i]
    if sourceRA != sourceRAprevious:
        sourceID += 1
    sourcevector.append(sourceID)
    sourceRAprevious = sourceRA

data['sourceID'] = sourcevector

data.write('SPIRE_ALMA_Cat_v4_sourceID.txt', format='ascii')
