#!/usr/bin/env python

from scipy import array
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

#positive autoregulation
# bla bla bla
def dx_dt(X, t):
    tau = 10
    beta = 1
    beta1 = 0.1
    alpha = np.log(2) / tau
    y = array([beta + beta1*X[0] - alpha*X[0]])
    return y

#negative autoregulation
def nar(X,t):
    tau = 10
    beta = 2
    k = 1
    alpha = np.log(2) / tau
    y = array([(beta/(alpha+X[0]/k)) - alpha*X[0]])
    return y

def sr(X,t):
    beta = 2
    tau = 10
    alpha = np.log(2) / tau
    y = array([beta - alpha*X[0]])
    return y

def sr1(X,t):
    beta = 0.11
    tau = 10
    alpha = np.log(2) / tau
    y = array([beta - alpha*X[0]])
    return y

if __name__ == "__main__":
    t = np.arange(0, 40, 0.1)
    X0 = array([0.0])
    X = integrate.odeint(sr, X0, t)
    Y = integrate.odeint(dx_dt,X0, t)
    plt.plot(t, X)
    plt.plot(t,Y)
    plt.show()
    
    
