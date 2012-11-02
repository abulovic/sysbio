#!/usr/bin/python

from scipy import integrate
from scipy import array
import matplotlib.pyplot as plt
import numpy as np
import math

def generate_dataset():
    pass

def dx_dt(X,t,th0, th1, th2, th3, th4, th5, th6, th7):
    y = array([((th0 * X[0]**2) / (th2 + X[0]**2 + th5 * X[1])) - th6 * X[0],
               ((th1 * X[0]**2) / (th3 + th4 * X[0]**2)) / th7 * X[1]])
    return y

if __name__ == "__main__":
    X0 = np.array([1, 0.5])
    t = np.linspace(0, 15, 1000)
    theta = (1,1,1,1,1,1,1,1)
    X, info = integrate.odeint(dx_dt, X0, t, args=theta, full_output=True)
    x,y = X.T
    plt.plot(t, x, 'r-', label='x(t)')
    plt.plot(t, y, 'b-', label='y(t)')
    plt.plot(t, np.cos(t), 'g-', label='cos(t)')
    plt.xlabel('time')
    plt.show()
    
    
    
    
