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

if __name__ == '__main__':
    X0 = np.array([1, 0.5])
    t = np.linspace(0, 5 *np.pi, 100)
    X  = integrate.odeint(dx_dt, X0, t)
    x,y = X.T
    plt.plot(t, x, 'r-')
    plt.plot(t, y, 'b-')
    plt.show()
