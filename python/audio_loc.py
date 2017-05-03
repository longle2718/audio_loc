'''
Utility functions

Long Le <longle1@illinois.edu>
University of Illinois
'''
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from joblib import Parallel, delayed

import os,sys
os.system("taskset -p 0xff %d" % os.getpid())
sys.path.append(os.path.expanduser('~')+'/audio_class/python')
sys.path.append('../../audio_class/python')
import audio_class
import ridgeDTW

def ridgePar(k,data,fs,tBlk,tInc,tBkt):
    S,F,T,tBlk,tInc = audio_class.spectrographic(data[k,:],fs,tBlk,tInc)
    X = ridgeDTW.ridgeTracker(S,np.median(S.flatten()),tBkt,tInc,isMaxPool=False,supThresh=3.)
    return T,F,X

def segProc(data,fs,tSegBlk=2.0,tSegInc=1.0,tBlk=0.032,tInc=0.004,tBkt=0.008,debug=0):
    # segment-processing the multi-channel data at
    # frame scale (parameterized by tBlk and tInc in seconds) and
    # segment scale (parameterized by tSegBlk and tSegInc in seconds)
    #
    nCh,_ = np.shape(data)
    rvs = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(ridgePar)(k,data,fs,tBlk,tInc,tBkt) for k in range(nCh))

    # plot for debugging
    if debug >= 1:
        T,F,X = rvs[7]
        tInc = T[1]-T[0]
        n1 = int(35/tInc)#12000
        n2 = int(35/tInc)#12800
        print('n1 = %s, n2 = %s' % (n1,n2))
        plt.figure(figsize=(15,15))
        plt.pcolormesh(T[n1:n2],F,np.sqrt(X[:,n1:n2]))
        plt.show()
    if debug >= 2:
        plt.figure(figsize=(15,30))
        for k in range(nCh):
            T,F,X = rvs[k]
            plt.subplot(nCh,1,1+k)
            plt.pcolormesh(T[n1:n2],F,np.sqrt(X[:,n1:n2]))
        plt.show()

    # segmenting the frames

    return
