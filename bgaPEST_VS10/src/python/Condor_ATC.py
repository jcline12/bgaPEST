import numpy as np
import subprocess as sub
import parallel_condor_Jacobian as pcj
import os
import shutil





# initialize a Jacobian Master object to hold results
fulljack = pcj.Jacobian_Master()

# read in the metadat that will be necessary to perform a Jacobian run
fulljack.read_obs_names()

fulljack.read_pars()

fulljack.read_mio_ins()

fulljack.read_jacfle()

fulljack.update_condor_subfile()

# make a scratch directory for all the output files
fulljack.jacfolder = '#jacfilestmp#'

# if it exists, empty it -- else, create it
if os.path.exists(os.path.join(os.getcwd(),fulljack.jacfolder)):
    shutil.rmtree(os.path.join(os.getcwd(),fulljack.jacfolder))
os.mkdir(os.path.join(os.getcwd(),fulljack.jacfolder))

'''
# run all the runs --> without CONDOR FOR NOW
for i in np.arange(fulljack.NPAR):
    print '===============\n\nrunning run number ---> %d\n===============\n' %(i)

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

fulljack.read_obs_files(fulljack.NPAR)

for i in np.arange(fulljack.NPAR):
    fulljack.read_obs_files(i)
    

# get the derivative increments
fulljack.read_derinc()

# adjust from observation values to sensitivities in fulljack.JAC
fulljack.calc_JAC()

# finally, write out the Jacobian into a text file
fulljack.Jacobian2jac(fulljack.jacfle)
k=1