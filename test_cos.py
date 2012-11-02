#!/usr/bin/python

from scipy import integrate
from scipy import array
import matplotlib.pyplot as plt
import numpy as np
import math

def generate_dataset():
    pass

def dx_dt(X,t):
    a1 = 2
    a2 = 1
    b1 = 1
    b2 = 2
    b3 = 1
    c1 = 1
    d1 = 2
    d2 = 2
    y = array([((a1 * X[0]**2) / (b1 + X[0]**2 + c1 * X[1])) - d1 * X[0],
               ((a2 * X[0]**2) / (b2 + b3 * X[0]**2)) / d2 * X[1]])
    return y

if __name__ == "__main__":
    X0 = np.array([1, 0.5])
    t = np.linspace(0, 15, 1000)
    X, info = integrate.odeint(dx_dt, X0, t, full_output=True)
    x,y = X.T
    plt.plot(t, x, 'r-', label='x(t)')
    plt.plot(t, y, 'b-', label='y(t)')
    plt.plot(t, np.cos(t), 'g-', label='cos(t)')
    plt.xlabel('time')
    plt.show()
    
    
    
    
