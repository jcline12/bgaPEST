import numpy as np
import os
'''
parallel_condor_Jacobian --> program for external bgaPEST derivatives using Condor.
a m!ke @ usgs joint

'''

class Jacobian_Master:
    # initialize the class
    def __init__(self):
        self.JAC = []
        self.tpl_pargp =  [] # dictionary with pargroups and tpl files 
        self.obs_names = []
        self.base_obs_vals = [] # base observations values
    
    def read_obs_names(self):
        # read in the base observations values
        indat = np.genfromtxt('bgaPEST.#obs',names=True,dtype=None)
        self.obs_names = indat['OBSNME']
        del indat

class Jacobian_one_run:
    # initialize the class
    def __init__(self):
        self.parvals = []
        self.pargroups = []
        self.pargpuniq = []
        self.parnames = []
        self.derinc = [] # dictionary with pargroups and derincs
        self.mio = [] # dictionary with tpl files/ins files and model/output files
        self.perturb = []
        self.base_obs = []
        self.JAC = []
        self.tpl_pargp = [] # dictionary with pargroups and tpl files 
        
    def make_model_input(self):
        # first make all the input files

        # determine which parameter must be perturbed and increment it
        # (if self.perturb == -999 then it's a base run --> no incrementing)
        if self.perturb > -999:
            pertgrp = self.pargroups[self.perturb]
            pertamt = self.derinc[pertgrp]
            self.parvals[self.perturb] += self.parvals[self.perturb]*pertamt

        for cg in self.pargpuniq:            
            # write a par file for the current group              
            cinds = np.where(self.pargroups == cg)[0]
            ofp = open(cg + '.#jacpars','w')
            ofp.write('single point\n')
            for cp in zip(self.parnames[cinds],self.parvals[cinds]):
                ofp.write('%s %f 1.0 0.0\n' %(cp[0],cp[1]))
            ofp.close()
            os.system('tempchek ' + 
                      self.tpl_pargp[cg] + ' ' +
                      self.mio[self.tpl_pargp[cg]] + ' ' +
                      cg + '.#jacpars')
                

            
    
    def read_parameters_and_meta_data(self):
        # read in the template and parameter group data and make a dictionary
        indat = np.genfromtxt('bgaPEST.#pgtpl', names=True,dtype=None)
        self.tpl_pargp = dict(zip(np.atleast_1d(indat['PARGROUP']),
                                  np.atleast_1d(indat['TPL_FILE'])))
        del indat
        
        # read in the parameter files and get group information
        indat = np.genfromtxt('bgaPEST.#pargp', names=True,dtype=None)
        self.derinc = dict(zip(np.atleast_1d(indat['PARGPNME']),
                                  np.atleast_1d(indat['DERINC'])))
        del indat
        
        # read in the parameter values file
        indat = np.genfromtxt('bgaPEST.#par', names=True,dtype=None)
        self.parnames = indat['PARNME']
        self.parvals = indat['PARVAL1']        
        self.pargroups  = indat['PARGP']
        self.pargpuniq  = np.unique(self.pargroups)
        del indat
        
        # read in the MIO information
        indat = np.genfromtxt('bgaPEST.#mio', names=True,dtype=None)
        self.mio = dict(zip(np.atleast_1d(indat['MIO_FILE']),
                                  np.atleast_1d(indat['MOD_FILE'])))
        del indat
        
    def write_par_files(self):
        # split out the par file into one for each group
        for cg in self.pargpuniq:
            ofp = open(cg + '.#par','w')
            ofp.write('%12s%12s%12s\n' %('PARGPNME','PARVAL1','PARGP'))

    def Jacobian2jac(Xtmp,outfile):
        # open the outfile
        ofp = open(outfile,'w')
        # write out the header
        ofp.write('%10d%10d%10d\n' %(nobs,npar,2))
        # write out the Jacobian
        k = 0
        cp = 0
        for i in Xtmp:
            k+=1
            cp+=1
            ofp.write('%16.8e ' %(i))
            if (k==8):
                k = 0
                ofp.write('\n')
            elif (cp==npar):
                cp = 0
                k = 0
                ofp.write('\n')
                
        # now write out the observation names
        ofp.write('* row names\n')
        for cobs in obsnames:
            ofp.write('%s\n' %(cobs))
        
        # now write out the parameter names
        ofp.write('* column names\n')
        for i in np.arange(npar):
            ofp.write('p%d\n' %(i+1))
        ofp.close()        
