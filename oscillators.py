#!/usr/bin/env python

from scipy import integrate
from scipy import array
import numpy as np
import matplotlib.pyplot as plt


#simple system that undergoes a hopf bifurcation
def dX_dt(X, t):
    ka = 2.8
    k2 = 1
    k3 = 1
    k4 = 1
    k5 = 1
    y = array([(ka-k4)*X[0] - k2*X[0]*X[1], -k3*X[1] + k3*X[2], k4*X[0] - k5*X[2]])
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
    b = 1
    y = array([-cs + a**2 * S[1], cr - b**2 * S[0]])
    return y

def do_dt(O,t):
    y = array([5])
    return y

if __name__ == '__main__':
    X0 = np.array([1, 1, 1])
    t = np.linspace(0, 5 *np.pi, 100)
    X  = integrate.odeint(dX_dt, X0, t)
    x,y,z = X.T
    plt.plot(t, x, 'r-')
    plt.plot(t, y, 'b-')
    plt.plot(t, z, 'g-')
    plt.show()
