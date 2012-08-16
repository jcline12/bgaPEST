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
param_to_perturb = int(sys.argv[1])

jack = Jacobian_calculator()

jack.read_parameters_and_meta_data()

i=1
