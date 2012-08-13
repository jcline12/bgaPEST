import numpy as np
import sys

'''
parallel_condor_Jacobian --> program for external bgaPEST derivatives using Condor.
a m!ke @ usgs joint

'''

class Jacobian_calculator:
    # initialize the class
    def __init__(self):
        self.parameters = []
        self.pargroups = []
        self.pargpuniq = []
        self.derinc = [] # dictionary with pargroups and derincs
        self.mio = [] # dictionary with tpl files/ins files and model/output files
        self.perturb = []
        self.base_obs = []
        self.JAC = []
        self.tpl_pargp = [] # dictionary with pargroups and tpl files 
    
    def read_parameters_and_meta_data(self):
        # read in the template and parameter group data and make a dictionary
        indat = np.genfromtxt('bgaPEST.#pgtpl', names=True,dtype=None)
        self.tpl_pargp = dict(zip(np.atleast_1d(indat['PARGROUP']),
                                  np.atleast_1d(indat['TPL_FILE'])))
        
        # read in the parameter files and get group information
        indat = np.genfromtxt('bgaPEST.#pargp', names=True,dtype=None)
        self.derinc = dict(zip(np.atleast_1d(indat['PARGPNME']),
                                  np.atleast_1d(indat['DERINC'])))
        
        # read in the parameter values file
        indat = np.genfromtxt('bgaPEST.#par', names=True,dtype=None)
        self.parameters = indat['PARVAL1']        
        self.pargroups  = indat['PARGP']
        self.pargpuniq  = np.unique(self.pargroups)
        
        # read in the MIO information
        indat = np.genfromtxt('bgaPEST.#mio', names=True,dtype=None)
        self.mio = dict(zip(np.atleast_1d(indat['MIO_FILE']),
                                  np.atleast_1d(indat['MOD_FILE'])))
        
    def write_par_files(self):
        # split out the par file into one for each group
        for cg in self.pargpuniq:
            ofp = open(cg + '.#par','w')
            ofp.write('%12s%12s%12s\n' %('PARGPNME','))
# ####### #
 # M A I N #
  # ####### #
param_to_perturb = int(sys.argv[1])

jack = Jacobian_calculator()

jack.read_parameters_and_meta_data()

i=1