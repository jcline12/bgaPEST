import numpy as np
import os

def csv2gslib(cf):
    
    
    print 'opening --> ' + cf
    indat = open(cf,'r').readlines()
    innames = indat.pop(0)
    innames = innames.strip().split(',')
    ofp = open(cf[:-4] + '.gslib','w')
    ofp.write('KirigingData\n%d\n' %(len(innames)))
    for cn in innames:
        ofp.write('%s\n' %(cn))
    for line in indat:
        for i in line.strip().split(','):
            ofp.write('%s ' %(i))
        ofp.write('\n')
    ofp.close()