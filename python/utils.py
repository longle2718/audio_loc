'''
Utility functions

Long Le <longle1@illinois.edu>
University of Illinois
'''
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from joblib import Parallel, delayed

cSound = 343 # speed of sound m/s

def delay2loc_grad(micsloc,meas_delayMat):
    # from delays to location using gradient descent
    #
    # micsloc: Mx3, 3-D locations of M mics
    # meas_delayMat: MxM, measured delay matrix
    #
    # loc: 1x3, output location of the sound source

    M = len(micsloc)
    nIter = 0

    locnow = np.mean(micsloc,axis=0)
    err0 = errEst(locnow,micsloc,meas_delayMat)
    while True
        grad = gradEst(locnow,miscloc,meas_delayMat,err0)
        # backtrack line search
        mu = 1. # gradient step size
        while True:
            # gradient descent update
            locNext = locnow-mu*grad
            err = errEst(locNext,micsloc,meas_delayMat)
            if err < .9*err0:
                locnow = locNext
                err0 = err
                break
            else:
                mu = mu/2

        nIter += 1
        # check terminal condition
        if np.norm(grad) < 1e-2 or nIter >= 1e3:
            break

    return locnow

def gradEst(loc,micsloc,meas_delayMat,err0):
    gradStep = .1 # meter
    N = len(loc)
    grad = [None]*N

    for k in range(N):
        locNext = loc
        locNext[k] += gradStep

        err = errEst(loc,micsloc,meas_delayMat)
        grad[k] = (err-err0)/gradStep
        
    return grad

def errEst(loc,micsloc,meas_delayMat):
    # estimate error at a location 
    # given measured delays and locations of array's mics

    delays,_ = estDelay(micsloc,locNext)
    delayMat = toMat(delays)
    err = np.mean((delayMat-meas_delayMat)**2)

    return err

def toMat(delays):
    # convert from delays vector to matrix
    delayMat = np.dot(np.ones((M,1)),np.expand_dims(delays,axis=0))-\
            np.dot(np.expand_dims(delays,axis=1),np.ones((1,M)))

    return delayMat

def estDelays(micsloc,srcloc):
    # estimate theoretical delays given a 
    # source location and an array
    #
    # micsloc: Mx3, 3-D locations of M mics in an array
    # srcloc: 1x3, 3-D location of source
    #
    # delays: Mx1, delays from srcloc to the array
    # gains: Mx1, gains from srcloc to the array
    
    M = len(micsloc)
    delays = [None]*M
    gains = [None]*M
    for k in range(M):
        dist = np.sqrt(np.sum((micsloc[k,:]-srcloc)**2))
        delays[k] = dist/cSound
        gains[k] = 1/dist
    
    # normalize
    delays = delays-np.mean(delays)
    gains = gains/np.norm(gains)

    return delays,gains
