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

def delay2loc_grad(micsloc,meas_delays):
    # from delays to location using gradient descent
    # micsloc: Mx3, 3-D locations of M mics
    # meas_delays: MxM, measured delays

    mu = 1. # gradient step size

    startloc = np.mean(micsloc,axis=0)
    locnow = startloc

    for k in range(128):
        delays = estDelays(micsloc,locnow,fs)

    return

def estDelays(micsloc,srcloc,fs):
    # estimate theoretical delays given a 
    # source location srcloc and the array micsloc

    return
