import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

def read_data_pars(infile,nrow,ncol):
    indat = np.genfromtxt(infile,names=True,dtype=None)
    pars = indat['ParamVal'].reshape(nrow,ncol)
    return pars
    

plt.close('all')
casename = "bp_test_pest"
nrow = 21
ncol = 21

parfiles = list()
residfiles = list()
allfiles = os.listdir(os.getcwd())
for cf in allfiles:
    if casename + ".bpp" in cf:
        parfiles.append(cf)
    elif casename + ".bre" in cf:
        residfiles.append(cf)
    elif '.png' in cf:
        os.remove(cf)
del allfiles

for cf in parfiles:
    cind = cf.strip().split('.')[-1]
    pars = read_data_pars(cf,nrow,ncol)
    plt.figure()
    plt.imshow(pars,interpolation='nearest')
    plt.colorbar()
    plt.title('Parameters, iteration %s' %(cind))
    plt.savefig('parameters_' + casename + '_' + cind + '.png')
    
for cf in residfiles:
    cind = cf.strip().split('.')[-1]
    resid_dat = np.genfromtxt(cf,names=True,dtype=None)
    inds = np.arange(len(resid_dat['ObsName']))
    width = 0.2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ms = ax.bar(inds,resid_dat['Measured'],width,color='b')
    ob = ax.bar(inds+width,resid_dat['Modeled'],width,color='r')
    ax.legend((ms[0],ob[0]),('Measured','Modeled'))
    plt.title('Observations, iteration %s' %(cind))
    plt.savefig('residuals_' + casename + '_' + cind + '.png')

    

    