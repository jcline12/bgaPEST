import os
import subprocess as sub

print 'checking for status file'
if os.path.exists(os.path.join(os.getcwd(),'jacdone.#stat')):
	print 'status file found and removed'
	os.remove(os.path.join(os.getcwd(),'jacdone.#stat'))
else:
	print 'no status file found'

print 'zipping up the data folder'
os.call('wine zip data.zip data')


