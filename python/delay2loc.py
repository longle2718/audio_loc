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
    while True:
        grad = gradEst(locnow,micsloc,meas_delayMat,err0)
        # backtrack line search
        mu = 1. # gradient step size
        while True:
            # gradient descent update
            locNext = locnow-mu*grad
            err = errEst(locNext,micsloc,meas_delayMat)
            if err <= err0:
                locnow = locNext
                err0 = err
                break
            else:
                mu = mu/2
                if mu == 0.:
                    break

        nIter += 1

        # check terminal condition
        print('nIter = %s' % nIter)
        print('|grad| = %s' % np.linalg.norm(grad))
        print('mu = %s' % mu)
        if nIter >= 1e2:
            print('Done! Max # of iterations reached')
            break
        if np.linalg.norm(grad) < 1e-9:
            print('Done! Gradient is sufficiently small')
            break
        if mu < 1e-9:
            print('Done! Stepsize is sufficiently small')
            break

    return locnow,err0

def gradEst(loc,micsloc,meas_delayMat,err0):
    gradStep = .1 # meter
    N = len(loc)
    grad = np.zeros(N)

    for k in range(N):
        locNext = loc
        locNext[k] += gradStep

        err = errEst(loc,micsloc,meas_delayMat)
        grad[k] = (err-err0)/gradStep
        
    return grad

def errEst(loc,micsloc,meas_delayMat):
    # estimate error at a location 
    # given measured delays and locations of array's mics

    delays,_ = estDelay(micsloc,loc)
    delayMat = toMat(delays)
    err = np.mean((delayMat-meas_delayMat)**2)

    return err

def toMat(delays):
    # convert from delays vector to matrix
    M = len(delays)
    delayMat = np.dot(np.ones((M,1)),np.expand_dims(delays,axis=0))-\
            np.dot(np.expand_dims(delays,axis=1),np.ones((1,M)))

    return delayMat

def estDelay(micsloc,srcloc):
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
    gains = gains/np.linalg.norm(gains)

    return delays,gains
