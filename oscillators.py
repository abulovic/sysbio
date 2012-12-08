#!/usr/bin/env python

from scipy import integrate
from scipy import array
from scipy import spatial
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import abcinfer as abc
import math

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
    y = array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1]])
    return y

def eu_dist(v1, v2):
    return np.sqrt(math.fsum((v1-v2)**2))

if __name__ == '__main__':
    X0 = np.array([1, 0.5])
    theta = np.array([1 ,1])
    theta1 = np.array([5 ,1])
    t = np.linspace(0, 15, 100)
    X  = integrate.odeint(lv, X0, t, args=(theta,))
    X1 = integrate.odeint(lv, X0, t, args=(theta1,))
    x,y = X.T
    x1,y1 = X1.T
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #x.plot(x, y)
    plt.figure(1)
    plt.plot(t,x1)
    plt.plot(t,y1)
    plt.figure(2)
    plt.plot(t,x)
    plt.plot(t,y)
    fx = np.fft.fft(x)
    fx1 = np.fft.fft(x1)
    print fx
    #plt.plot(t,z)
    #plt.figure(3)
    #hx = np.fft.fft(x)
    #plt.plot(hx)
    plt.show()
