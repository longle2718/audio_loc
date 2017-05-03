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

def extrSpec(data,fs,tBlk,tInc):
    # routine for extracting spectrogram
    # tBlk: block time
    # tInc: increment time
    #
    S,F,T,tBlk,tInc = audio_class.spectrographic(data,fs,tBlk,tInc)
    return S,F,T,tBlk,tInc

def extrRidge(S,tInc,tBkt):
    # (parallel) routine for extracting ridges
    # tBkt: backtrack time
    #
    X = ridgeDTW.ridgeTracker(S,np.median(S.flatten()),tBkt,tInc,isMaxPool=False,supThresh=3.)
    return X,T,F

def pool(X,T,F,tSegBlk,tSegInc,tInc):
    # aggregate the input representation
    # to form the output representation
    #
    NF,NT = np.shape(X)
    fBlk = int(tSegBlk/tInc) # frame block in a segment
    fInc = int(tSegInc/tInc) # frame increment between consecutive segments
    NSeg = int(np.ceil(NT/fInc))

    FF = F
    # #-of-segs x #-of-frames-per-seg x #-of-seconds-per-frame
    TT = np.arange(NSeg)*fInc*tInc 
    XX = np.zeros((NF,nSeg))
    i = 0
    for k in range(nSeg):
        j = min(NT,i+fBlk)

        XX[:,k] = np.sum(np.pad(X[:,i:j],(NF,fBlk-(j-i)),'constant',constant_values=0),axis=1)

        i = min(NT,i+fInc)
    
    return XX,TT,FF

def segProc(data,fs,tSegBlk=2.0,tSegInc=1.0,tBlk=0.032,tInc=0.004,tBkt=0.008,debug=0):
    # segment-processing the multi-channel data at
    # frame scale (parameterized by tBlk and tInc in seconds) and
    # segment scale (parameterized by tSegBlk and tSegInc in seconds)
    #
    NCh,_ = np.shape(data)

    rvSpec = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrSpec)(data[k,:],fs,tBlk,tInc) for k in range(NCh))

    rvRidge = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrRidge)(rvSpec[k][0],tInc,tBkt) for k in range(NCh))

    # plot for debugging
    if debug >= 1:
        T,F,X = rvRidge[7]
        tInc = T[1]-T[0]
        n1 = int(35/tInc)#12000
        n2 = int(37/tInc)#12800
        print('n1 = %s, n2 = %s' % (n1,n2))
        plt.figure(figsize=(15,15))
        plt.pcolormesh(T[n1:n2],F,np.sqrt(X[:,n1:n2]))
        plt.show()
    if debug >= 2:
        plt.figure(figsize=(15,30))
        for k in range(NCh):
            T,F,X = rvRidge[k]
            plt.subplot(NCh,1,1+k)
            plt.pcolormesh(T[n1:n2],F,np.sqrt(X[:,n1:n2]))
        plt.show()

    # group frames into segments
    rvRidge1 = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(pool)(rvRidge[k][0],rvRidge[k][1],rvRidge[k][2],tSegBlk,tSegInc,tInc) for k in range(NCh))

    rvRidge2 = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrRidge)(rvRidge1[k][0],tInc,tBkt) for k in range(NCh))

    return
