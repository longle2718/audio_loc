'''
(Time) delays to location conversion functions

Long Le <longle1@illinois.edu>
University of Illinois
'''
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from joblib import Parallel, delayed

cSound = 340 # speed of sound m/s

def delay2loc_grad(micsloc,meas_delayMat,mu=1.,debug=False):
    # from delays to location using gradient descent
    #
    # micsloc: Mx3, 3-D locations of M mics
    # meas_delayMat: MxM, measured delay matrix
    # mu: gradient step size
    #
    # loc: 1x3, output location of the sound source

    M,N = np.shape(micsloc)
    nIter = 0

    locNow = np.mean(micsloc,axis=0)
    errNow = errEst(locNow,micsloc,meas_delayMat)
    grad = gradEst(locNow,errNow,micsloc,meas_delayMat)
    while True:
        # backtrack line search
        while True:
            # gradient descent update
            locNext = locNow-mu*grad
            #print('===== locNext = %s, locNow = %s, mu = %s, grad = %s' % (locNext,locNow,mu,grad))
            errNext = errEst(locNext,micsloc,meas_delayMat)
            if errNext <= errNow:
                locNow = locNext
                errNow = errNext
                grad = gradEst(locNow,errNow,micsloc,meas_delayMat)
                break
            else:
                mu = mu/2
                if mu == 0.:
                    break

        # check terminal condition
        nIter += 1
        if debug:
            print('nIter = %s, mu = %s, grad = %s, loc= %s, err = %s' % \
                    (nIter,mu,grad,locNow,errNow))
        if nIter >= 1e4:
            if debug:
                print('Done! Max # of iterations reached')
            break
        if np.linalg.norm(grad) < 1e-6:
            if debug:
                print('Done! Gradient is sufficiently small')
            break
        if mu < 1e-2:
            if debug:
                print('Done! Stepsize is sufficiently small')
            break

    return locNow,errNow,grad

def gradEst(loc,err,micsloc,meas_delayMat):
    # interior-point optimization
    # gradient of the ojective and the logarithmic barrier functions
    # of the constraint |loc-micsCen| <= 200,
    # i.e. the solution should be within the 200 m radius
    gradStep = .1 # meter
    N = len(loc)
    grad = np.zeros(N)
    micsCen = np.mean(micsloc,axis=0) # micArray centroid
    lambda = 1.

    for k in range(N):
        locNext = np.array(loc); locNext[k] += gradStep
        errNext = errEst(locNext,micsloc,meas_delayMat)

        d = np.linalg.norm(loc-micsCen)
        dNext = np.linalg.norm(locNext-micsCen)

        grad[k] = (errNext-err)/gradStep + \
                lambda/(200-d)*(dNext-d)/gradStep
        
    return grad

def errEst(loc,micsloc,meas_delayMat):
    # estimate error at a location 
    # given measured delays and locations of array's mics

    delays,_ = estDelay(micsloc,loc)
    delayMat = toMat(delays)
    err = 0.5*np.sum(np.abs(delayMat-meas_delayMat)) # total error in delay is the metric to minimize

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
