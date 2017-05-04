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

def pool(S,F,T,tSegBlk,tSegInc,tInc):
    # aggregate the input representation
    # to form the output representation
    #
    NF,NT = np.shape(S)
    fBlk = int(tSegBlk/tInc) # # of frames in a segment
    tSegBlk = fBlk*tInc
    fInc = int(tSegInc/tInc) # # of frames incremented between consecutive segments
    tSegInc = fInc*tInc
    NSeg = int(np.ceil(NT/fInc))

    FF = F
    # #-of-segs x #-of-frames-per-seg x #-of-seconds-per-frame
    TT = np.arange(NSeg)*fInc*tInc 
    SS = np.zeros((NF,NSeg))
    i = 0
    for k in range(NSeg):
        j = min(NT,i+fBlk)

        SS[:,k] = np.mean(np.pad(S[:,i:j],((0,0),(0,fBlk-(j-i))),'constant',constant_values=0),axis=1)

        i = min(NT,i+fInc)
    
    return SS,FF,TT,tSegBlk,tSegInc

def extrRidge(S,tInc,bktRatio):
    # (parallel) routine for extracting ridges
    # bktRatio: backtrack ratio
    #
    tBkt = tInc*bktRatio # backtrack time
    X = ridgeDTW.ridgeTracker(S,np.median(S.flatten()),tBkt,tInc,isMaxPool=False,supThresh=5.5)
    return X

def hieProc(data,fs,tHieBlk=[0.032,2.0],tHieInc=[0.004,1.0]):
    # hierarchical-processing the multi-channel data at
    # frame scale (parameterized by tBlk and tInc in seconds) and
    # segment scale (parameterized by tSegBlk and tSegInc in seconds)
    #
    NCh,_ = np.shape(data)
    NHie = len(tHieBlk)
    hSpecs = [None]*NHie # hierarchical specs
    hRidges = [None]*NHie # hierarchical ridges
    bktRatio = 2

    hSpecs[0] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrSpec)(data[k,:],fs,tHieBlk[0],tHieInc[0]) for k in range(NCh))
    hRidges[0] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrRidge)(hSpecs[0][k][0],hSpecs[0][k][-1],bktRatio) for k in range(NCh))
    for l in range(1,NHie):
        hSpecs[l] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(pool)(hSpecs[l-1][k][0],hSpecs[l-1][k][1],hSpecs[l-1][k][2],tHieBlk[l],tHieInc[l],hSpecs[l-1][k][-1]) for k in range(NCh))
        hRidges[l] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrRidge)(hSpecs[l][k][0],hSpecs[l][k][-1],bktRatio) for k in range(NCh))

    return hRidges,hSpecs
