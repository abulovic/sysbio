#!/usr/bin/env python

from scipy import integrate
from scipy import array
from scipy import spatial
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import abcinfer_1 as abc
import math
import sys
from scipy.spatial.distance import euclidean
import stats_util as utils


#simple system that undergoes a hopf bifurcation
def dX_dt(X, t):
    ka = 2.5
    k2 = 1
    k3 = 1
    k4 = 1
    k5 = 1
    y = array([(ka-k4)*X[0] - k2*X[0]*X[1] + 1, -k3*X[1] + k3*X[2], k4*X[0] - k5*X[2]])
    return y

#simple harmonic oscillator
def dx_dt(X, t):
    freq = 1
    y = array([X[1] , -freq**2 * X[0]])
    return y

#positive simple harmonic oscillator with offsets from zero
def ds_dt(S, t):
    cs = 2
    cr = 1
    a = 1
    b = 5
    y = array([-cs + a**2 * S[1], cr - b**2 * S[0]])
    return y

def do_dt(O,t):
    y = array([5])
    return y

def lv(X,t, theta):
    a = theta[0]
    b = theta[1]
    y = array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1], b*X[0]*X[1] - X[1] ])
    return y

def eu_dist(v1, v2):
    return np.sqrt(math.fsum((v1-v2)**2))

#simple clock model with 2 genes A,B
def dA_dt(A, t, theta):
    k = theta[0]
    c = theta[1]
    b = theta[2]
    d = theta[3]
    b1 = theta[4]
    p = theta[5]
    k1 = theta[6]
    d1 = theta[7]
    y = array([(A[1]*A[0] / (k + A[0])) - c / (b + A[1]) - d*A[0],
               (b1*A[0]**p / (k1 + A[0]**p) - d1*A[1])])
    return y

#hes1 model as found in http://arxiv.org/abs/1106.6280
def hes1(M, t, theta):
    P0 = theta[0]
    v = theta[1]
    k1 = theta[2]
    h = theta[3]
    kdeg = 0.03
    p1 = 5.0
    p2 = 3.0
    y = array([-kdeg*M[0] + (1 / (1 + (M[2] / P0)**h)), -kdeg*M[1] + v*M[0] - k1*M[1],
                -kdeg*M[2] + k1*M[1]])
    return y

def switch(S, t, theta):
    m1 = theta[0]
    m2 = theta[1]
    k1 = theta[2]
    k2 = theta[3]
    n = theta[4]
    y = array([k1/(1 + S[1]**n) - m1*S[0], k2/(1 + S[0]**n) - m2*S[1], k2/(1 + S[0]**n) - m2*S[2]])
    return y

def main():
    theta = [2.4, 0.02, 0.2, 6.9 ]
    theta1 = [  1.37327745e+04   ,1.46981207e-02   ,1.99460716e-01  ,-8.26157999e+04]
    M0 = np.array([2.0, 5.0, 3.0])
    ds = abc.generate_dataset(hes1, theta)
    noisy_ds = abc.add_gaussian_noise(np.copy(ds))
    populations = (abc.smc(hes1, noisy_ds, [160.0, 120.0, 80.0,70.0,  60.0, 40.0, 30.0, 20.0, 15.0, 
                                           10.0, 9.5, 9.0, 8.5, 8.0, 7.5, 7.0,6.9,6.8,6.7,6.6, 6.5,
                                            6.0,5.9, 5.8, 5.7,5.6, 5.5,5.4, 5.3,5.2, 5.1,
                                            5.0,4.8, 4.5,4.3, 4.0,3.8, 3.5,3.2, 3.0, 2.8, 2.5, 2.3, 2.0]))
    sys.exit(0)
    t = np.arange(0, 15, 0.1)
    M = integrate.odeint(hes1, M0, t, args=(theta, ))
    M1 = integrate.odeint(hes1, M0, t, args=(theta1, ))
    m,p1,p2 = M.T
    m1, p11, p21 = M1.T
    plt.subplot(311)
    plt.plot(t, m, 'r-')
    plt.plot(t, m1, 'b-')
    plt.subplot(312)
    plt.plot(t, p1,'r-')
    plt.plot(t, p11,'b-')
    plt.subplot(313)
    plt.plot(t, p2,'r-')
    plt.plot(t, p21,'b-')
    plt.show()

def fourier_compare(x, x1):
    fx = np.fft.fft(x)
    fx1 = np.fft.fft(x1)
    main_ind = []
    main_vals = []
    for ind, x in enumerate(fx):
        if x > 10.0:
            main_ind.append(ind)
            main_vals.append(abs(x))

    main_vals1 = utils.select(fx1, main_ind)
    return euclidean(main_vals, main_vals1)
    
def fourier():
    X0 = np.array([1, 0.5])
    theta = np.array([1 ,1])
    theta1 = np.array([5, 8])
    t = np.linspace(0, 15, 100)
    X  = integrate.odeint(lv, X0, t, args=(theta,))
    X1 = integrate.odeint(lv, X0, t, args=(theta1,))
    x,y = X.T
    x1,y1 = X1.T
    plt.figure(1)
    plt.plot(t,x1)
    plt.plot(t,y1)
    plt.figure(2)
    plt.plot(t,x)
    plt.plot(t,y)
    print (fourier_compare(y, y1) + fourier_compare(x, x1)) / 2 #take the distance of the fourier spectra along all signals
    #sys.exit()
    #plt.figure(3)
    #plt.plot(mfx, 'ro')
    #plt.plot(mfx1, 'bo')
    #print "distance between transforms: ", euclidean(mfx, mfx1)
    plt.show()

def lv_test():
    theta = [1,1]
    ds = abc.generate_dataset(lv, theta)
    ds = abc.add_gaussian_noise(ds)
    population = abc.smc(lv, ds, [30.0, 16.0, 6.0, 5.0, 4.3])
    sys.exit(0)
    
if __name__ == '__main__':
    lv_test()
