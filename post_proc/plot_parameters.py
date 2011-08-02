import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def read_data(infile,nrow,ncol):
    indat = np.genfromtxt(infile,names=True,dtype=None)
    pars = indat['ParamVal'].reshape(nrow,ncol)
    return pars
    


plt.close('all')
casename = "bp_test_pest"
nrow = 21
ncol = 21

parfiles = list()
allfiles = os.listdir(os.getcwd())
for cf in allfiles:
    if casename + ".bpp" in cf:
        parfiles.append(cf)
del allfiles

for cf in parfiles:
    cind = cf.strip().split('.')[-1]
    pars = read_data(cf,nrow,ncol)
    plt.figure()
    plt.imshow(pars,interpolation='nearest')
    plt.colorbar()
    plt.title('Parameters, iteration %s' %(cind))
    plt.savefig('parameters_' + casename + '_' + cind + '.png')
    
    

    