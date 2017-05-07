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
    X = ridgeDTW.ridgeTracker(S,np.median(S,axis=1),tBkt,tInc,isMaxPool=False,supThresh=5.5)
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

    hSpecs[0] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrSpec)\
            (data[k,:],fs,tHieBlk[0],tHieInc[0]) for k in range(NCh))
    hRidges[0] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrRidge)\
            (hSpecs[0][k][0],hSpecs[0][k][-1],bktRatio) for k in range(NCh))
    for l in range(1,NHie):
        hSpecs[l] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(pool)\
                (hSpecs[l-1][k][0],hSpecs[l-1][k][1],hSpecs[l-1][k][2],tHieBlk[l],tHieInc[l],hSpecs[l-1][k][-1]) for k in range(NCh))
        hRidges[l] = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(extrRidge)\
                (hSpecs[l][k][0],hSpecs[l][k][-1],bktRatio) for k in range(NCh))

    return hRidges,hSpecs

def segment(ts,nBlk,nInc):
    # segment a time series
    # ts: time series MxN
    # M: data dimension
    # N: length of the time series

    M,N = np.shape(ts)
    NSeg = int(np.ceil(N/nInc)) # total number of frames
    tsSeg = [None]*NSeg
    i = 0
    for k in range(NSeg):
        j = min(N,i+nBlk)

        tsSeg[k] = np.pad(ts[:,i:j],((0,0),(0,nBlk-(j-i))),'constant',constant_values=0)

        i = min(N,i+nInc)

    return tsSeg

def labelObjects(XX,mask=None):
    # count and label all TF objects in
    # a spectrographic image/matrix,
    # constrained by a mask
    #
    if mask == None:
        mask = np.ones(np.shape(XX))

    # make a copy of the input array
    X = np.array(XX)
    M,N = np.shape(X)
    
    def bfs(X,start,cnt):
        # breadth first search
        M,N = np.shape(X)
        TFObj = np.zeros((M,N))

        explored = set()
        frontierQ = []
        frontierQ.append(start)
        while len(frontierQ) > 0:
            node = frontierQ.pop(0)

            TFObj[node] = X[node]
            X[node] = 0.

            # visit neighbors
            for ngb in getNeighbor(node,X):
                if ngb not in explored:
                    explored.add(ngb)
                    frontierQ.append(ngb)

        return TFObj

    def getNeighbor(node,X):
        # define local constraints
        M,N = np.shape(X)
        ngb = []
        for d in [[0,1],[1,0],[1,1],[0,-1],[-1,0],[-1,-1],[1,-1],[-1,1]]:
            n = tuple(np.array(node)+d) 
            if n[0]>=0 and n[0]<M and n[1]>=0 and n[1]<N and X[n]>0 and mask[n]>0:
                ngb.append(n)

        return ngb

    TFObjs = []
    cnt = 0
    for k in range(M):
        for l in range(N):
            TFObj = bfs(X,(k,l),cnt)
            if np.any(TFObj):
                TFObjs.append(TFObj)
                cnt += 1

    '''
    plt.figure()
    plt.pcolormesh(X)
    plt.show()
    print(cnt)
    '''
    return TFObjs

def pruneObj(inObjs):
    N = len(inObjs)
    val = np.zeros(N)
    for k in range(N):
        val[k] = np.mean(inObjs[k])
    thresh = np.percentile(val,80)

    outObjs = []
    for k in range(N):
        if np.mean(inObjs[k]) > thresh:
            outObjs.append(inObjs[k])

    return outObjs
