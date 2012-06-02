import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from csv2gslib import csv2gslib



indat = np.loadtxt('Mike_Sill_jump.dat')
csvfilename = 'Mike_Sill2.csv'
allrows = indat.shape[0]
allcols = indat.shape[1]

#resample the image
numc = 20
numr = 30

delr = np.round(allrows/numr)
delc = np.round(allcols/numc)

subrows = np.arange(0,allrows,delr)
subcols = np.arange(0,allcols,delc)
subrows = np.hstack((subrows,allrows-1))
subcols = np.hstack((subcols,allcols-1))

ofp = open(csvfilename,'w')
ofp.write('X,Y,Z\n')

for cr in subrows:
    for cc in subcols:
        ofp.write('%d,%d,%f\n' %(cc+1,cr+1,indat[cr,cc]))

ofp.close()

csv2gslib(csvfilename)


plt.figure()
plt.imshow(indat,cmap=cm.Greys)
plt.savefig('figure_orig.pdf')
plt.savefig('figure_orig.png')
