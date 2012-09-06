import numpy as np
import sys
import parallel_condor_Jacobian as pcj
'''
Single run configuration for external bgaPEST derivatives using Condor
'''

parind = int(sys.argv[1]) #index for which parameter to perturb

# ####### #
 # M A I N #
  # ####### #
# initialize a single model run object
jack_one_run = pcj.Jacobian_one_run()

# determine which parameter index to perturb
jack_one_run.perturb = int(sys.argv[1])

# read in the parameter values and meta data for the run
jack_one_run.read_parameters_and_meta_data()

# create the model input files using TEMPCHEK
jack_one_run.make_model_input()

# run the model by calling the command line argument passed from bgaPEST
jack_one_run.run_model()

# get the observation values and save to proper filename
jack_one_run.read_obs()


i=1
