#!/usr/bin/env python

import abcinfer_multi as abc
from SimpleRepressilator import HillRepressilator
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import integrate, array

def dx_dt(X, t, th):
    kA = 3.
    k2 = th[0]
    k3 = th[1]
    k4 = th[2]
    k5 = th[3]
    y = array([(kA- k4)*X[0] - k2*X[0]*X[1],
               -k3*X[1] + k5*X[2],
               k4*X[0] - k5*X[2]])
    return y
    
def main():
    theta = [1, 1, 1, 1]
    ds = abc.generate_dataset_full(dx_dt, theta)
    print ds
    plt.plot(ds)
    plt.show()
    print "Showing plots."

    
if __name__ == '__main__':
    main()
    print 'This I add now.'
