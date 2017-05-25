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

        SS[:,k] = np.mean(np.pad(S[:,i:j],((0,0),(0,fBlk-(j-i))),'constant'),axis=1)

        i = min(NT,i+fInc)
    
    return SS,FF,TT,tSegBlk,tSegInc

def segment(X,tSegBlk,tSegInc,tInc):
    # segment a time series
    # X: time series MxN
    # M: data dimension
    # N: length of the time series

    M,N = np.shape(X)
    nBlk = int(tSegBlk/tInc)
    tSegBlk = nBlk*tInc
    nInc = int(tSegInc/tInc)
    tSegInc = nInc*tInc

    NSeg = int(np.ceil(N/nInc)) # total number of frames
    XSeg = [None]*NSeg
    i = 0
    for k in range(NSeg):
        j = min(N,i+nBlk)

        XSeg[k] = np.pad(X[:,i:j],((0,0),(0,nBlk-(j-i))),'constant')

        i = min(N,i+nInc)

    return XSeg,tSegBlk,tSegInc

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

def labelObjects(XX):
    # count and label all TF objects in
    # a spectrographic image/matrix,
    #
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
            if n[0]>=0 and n[0]<M and n[1]>=0 and n[1]<N and X[n]>0:
                ngb.append(n)

        return ngb

    # make a copy of the input array
    X = np.array(XX)
    M,N = np.shape(X)
    
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
    if N == 0:
        return inObjs

    val = np.zeros(N)
    for k in range(N):
        val[k] = np.mean(inObjs[k])
    thresh = np.percentile(val,80)

    outObjs = []
    for k in range(N):
        if np.mean(inObjs[k]) >= thresh:
            outObjs.append(inObjs[k])

    return outObjs

def seg2bounds(idx,tSegBlk,tSegInc,tInc):
    # return the sample bounds of a segment at an index
    n0 = int(idx*tSegInc/tInc)
    n1 = n0+int(tSegBlk/tInc)

    return n0,n1

def bestLink(grps):
    # grps is NCh x NGrp (length varies)
    # find the most overlapped set of groups
    NCh = len(grps)
    seq = []

    V = [None]*(NCh)
    backPtr = [None]*(NCh)
                    
    N = len(grps[0])
    if N == 0:
        return seq
    V[0] = np.zeros(N)
    backPtr[0] = -np.ones(N,dtype=int)
    for chIdx in range(NCh-1):
        N0 = N

        N = len(grps[chIdx+1]) # num of groups in a channel
        if N == 0:
            return seq
        V[chIdx+1] = np.zeros(N)
        backPtr[chIdx+1] = -np.ones(N,dtype=int)
        for k in range(N):
            link = V[chIdx]
            for l in range(N0):
                # compute the link weight as the geometric mean of two groups
                link[l] += max(gramCorr(grps[chIdx+1][k],grps[chIdx][l],NInc=50))
            # this looks like a max pooling
            V[chIdx+1][k] = np.max(link)
            backPtr[chIdx+1][k] = np.argmax(link)
            
    # backtracking
    seq.append(np.argmax(V[NCh-1]))
    for chIdx in range(NCh-1,0,-1):
        #print('chIdx = %s' % chIdx)
        #print('seq[-1] = %s' % seq[-1])
        seq.append(backPtr[chIdx][seq[-1]])

    return seq[::-1]

def gramCorrPar(n,spec1,spec2):
    _,NT1 = np.shape(spec1)
    _,NT2 = np.shape(spec2)
    return np.sum(np.sqrt(np.pad(spec1,((0,0),(max(0,NT2-n),n)),'constant')\
                         *np.pad(spec2,((0,0),(n,max(0,NT1-n))),'constant')))
def gramCorr(spec1,spec2,NInc=1):
    # AIgram/spectrogram correlation
    _,NT1 = np.shape(spec1)
    _,NT2 = np.shape(spec2)
    return Parallel(n_jobs=multiprocessing.cpu_count())(delayed(gramCorrPar)\
            (n,spec1,spec2) for n in range(0,NT1+NT2,NInc))

def rmGrp(grps,seq):
    for k in range(len(grps)):
        del grps[k][seq[k]]
    return
