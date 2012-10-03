import os
import zipfile as zf

dir_name = 'bgaPEST_DAG'
currd = os.path.join(os.getcwd(),dir_name)
czip = zf.ZipFile('test.zip','w')
for root,dir,files in os.walk(currd):
    for f in files:
        fname = os.path.join(root,f)
        new_path = os.path.normpath(fname.replace(currd,''))
        czip.write(fname,new_path,zf.ZIP_DEFLATED)
czip.close()