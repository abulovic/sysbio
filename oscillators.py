#!/usr/bin/env python

from scipy import integrate
from scipy import array
import numpy as np
import matplotlib.pyplot as plt


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

if __name__ == '__main__':
    X0 = np.array([1, 0.5])
    t = np.linspace(0, 5 *np.pi, 100)
    X  = integrate.odeint(ds_dt, X0, t)
    x,y = X.T
    plt.plot(t, x, 'r-')
    plt.plot(t, y, 'b-')
    plt.show()
