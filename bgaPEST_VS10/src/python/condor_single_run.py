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
jack_one_run = pcj.Jacobian_one_run()

jack_one_run.perturb = int(sys.argv[1])

jack_one_run.read_parameters_and_meta_data()

jack_one_run.make_model_input()

i=1
