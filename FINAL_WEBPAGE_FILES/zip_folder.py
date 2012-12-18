import zipfile
import os
import re
import sys

zipdir = sys.argv[1] # get the name of the folder to zip up
currd = os.getcwd()
zf = zipfile.ZipFile(zipdir + '.zip','w')
for croot,cdir,cfiles in os.walk(os.path.join(os.getcwd(),zipdir)):
    for cf in cfiles:
        curroot = re.sub(currd,'',croot)
	zf.write(os.path.join(croot,cf),os.path.join(curroot,cf),zipfile.ZIP_DEFLATED)
zf.close()
