import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 

def getHDF5HeaderInfo(pathname,filename):

    N = ['NPeriods','comments','countMode','countPeriodSeconds','counterAInput','counterTInput','daynight',
             'discAMode','filename','fixedDiscThresholdmV','freqband','pol','run','runTimeMinutes','timestamp','weather']

    df1 = pd.read_hdf(filename,key='header',names = N)

    runN = filename[3]
    
    newfilename = 'HeaderInfo_Run'+runN+'.txt'
    newf = open(dirc+newfilename,'w')

    for i in N:
        newf.write(i+': ')
        x = df1[i].values[0]
        if type(x) == str:
            newf.write(x)
        else:
            newf.write(str(x))
        newf.write('\n')
        
    newf.close()
    return

