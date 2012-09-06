import numpy as np
import subprocess as sub
import parallel_condor_Jacobian as pcj
import os, shutil





# initialize a Jacobian Master object to hold results
fulljack = pcj.Jacobian_Master()

fulljack.read_obs_names()

fulljack.read_pars()

fulljack.read_mio_ins()

fulljack.jacfolder = '#jacfilestmp#'

'''
# make a scratch directory for all the output files
fulljack.jacfolder = '#jacfilestmp#'
if os.path.exists(os.path.join(os.getcwd(),fulljack.jacfolder)):
    shutil.rmtree(os.path.join(os.getcwd(),fulljack.jacfolder))
os.mkdir(os.path.join(os.getcwd(),fulljack.jacfolder))

# run all the runs --> without CONDOR FOR NOW
for i in np.arange(fulljack.NPAR):
    print 'running run number ---> %d' %(i)
    sub.call('python condor_single_run.py %d' %(i))
    
# finally run the base case
sub.call('python condor_single_run.py -999')

# this is only in place until redirect is made (using python)
for cf in os.listdir(os.getcwd()):
    if '.obf.' in cf:
        shutil.move(cf,os.path.join(os.getcwd(),fulljack.jacfolder,cf))
'''

# read in the results and populate the Jacobian
fulljack.JAC = np.zeros((fulljack.NOBS,fulljack.NPAR))

fulljack.read_obs_files(-999)

for i in np.arange(fulljack.NPAR):
    fulljack.read_obs_files(i)
    

k=1