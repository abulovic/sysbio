#!/usr/bin/python

from scipy import integrate
from scipy import array
import matplotlib.pyplot as plt
import numpy as np
import math
import abcinfer as abc
import numpy as np
times = (1.1, 2.4, 3.9, 5.6, 7.5, 9.6, 11.9, 14.4)
datapoints = 8

def generate_data():
    th = 1
    dataset = np.zeros([datapoints, 2])
    for i in range(datapoints):
        dataset[i] = np.array([np.cos(times[i]), np.cos(times[i] + th)])
    return dataset
    
def dx_dt(X,t,th):
    a1 = th[0]
    a2 = th[1]
    b1 = th[2]
    b2 = th[3]
    b3 = th[4]
    c1 = th[5]
    d1 = th[6]
    d2 = th[7]
    y = array([((a1 * X[0]**2) / (b1 + X[0]**2 + c1 * X[1])) - d1 * X[1],
               ((a2 * X[0]**2) / (b2 + b3 * X[0]**2)) - d2 * X[0]])
    return y

#ode system for Lotka-Voltera model
def dy_dt(X,t,a,b):
    y = array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1]])
    return y

def plot_solution(population=[]):
    for i in range(50):
        theta =  np.random.uniform(-15,15,8)
    for p in population:
        theta.append(p[0])
        tpl_theta = tuple(theta)
        X0 = np.array([1, 0.5])
        t = np.arange(0, 15, 0.1)#np.linspace(0, 15, 1000)
        X= integrate.odeint(dx_dt, X0, t, args=(theta,))
        x,y = X.T
        plt.plot(t, x, 'r-', label='x(t)')
        plt.plot(t, y, 'b-', label='y(t)')
        plt.plot(t, np.cos(t+1), 'g-', label='cos(t+1)')
        plt.plot(t,np.cos(t), 'g-',label='cos(t)')
        plt.xlabel('time')
        plt.show()
    
if __name__ == "__main__":
    plot_solution()
    ds = generate_data()
    population = abc.rejector_algorithm(dx_dt, ds)
    plot_solution(population)
    
    
    
    
    
    
