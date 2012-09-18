import os
import subprocess as sub
import time


print 'checking for status file'
if os.path.exists(os.path.join(os.getcwd(),'jacdone.#stat')):
        print 'status file found and removed'
        os.remove(os.path.join(os.getcwd(),'jacdone.#stat'))
else:
        print 'no status file found'

if os.path.exists(os.path.join(os.getcwd(),'log')):
	for cf in os.listdir(os.path.join(os.getcwd(),'log')):
		print cf
		os.remove(os.path.join(os.getcwd(),'log',cf))
else:
	os.mkdir(os.path.join(os.getcwd(),'log'))


jac_in_proc = True
dag_sub = 'dag_jacobian.dag'
print 'Starting --> ' + dag_sub
for cf in os.listdir(os.getcwd()):
	if dag_sub + '.' in cf:
		os.remove(os.path.join(os.getcwd(),cf))
		print 'removing --> ' + os.path.join(os.getcwd(),cf)
sub.call('condor_submit_dag  -notification never ' + dag_sub, shell=True)
while jac_in_proc:
	if os.path.exists(os.path.join(os.getcwd(),'jacdone.#stat')):
		jac_in_proc = False
	time.sleep(10)


