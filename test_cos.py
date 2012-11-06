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

def dn_dt(N, t, theta):
    k = theta[0]
    a = theta[1]
    n = array([k * N[1] * N[0], -a * k * N[1] * N[0]])
    return n

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
def dy_dt(X,t,theta):
    a = theta[0]
    b = theta[1]
    y = array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1]])
    return y

def plot_solution(population=[]):
    theta = []
    for p in population:
        theta.append(math.fsum(p) / len(p))
    X0 = np.array([1.0, 0.5])
    t = np.arange(0, 15, 0.1)
    X= integrate.odeint(dx_dt, X0, t, args=(theta,))
    x,y = X.T
    plt.figure(1)
    plt.plot(t, x, 'r-', label='x(t)')
    plt.plot(t, y, 'b-', label='y(t)')
    plt.plot(t, np.cos(t), 'g-')
    plt.xlabel('time')
    plt.show()
    
if __name__ == "__main__":
#    print theta
#    n0 = [0.1,10]
#    t = np.arange(0, 480, 5)
#    theta1 = [0.005, 1]
#    r = integrate.odeint(dn_dt, n0, t, args=(theta,))
#    y = integrate.odeint(dn_dt, n0, t, args=(theta1,))
#    plt.plot(t, r, '-o')
#    plt.plot(t, y, 'r-')
#    plt.legend(['Bacteria', 'Nutrients'], loc='lower right')
#    plt.xlabel('Time')
#    plt.ylabel('Concentration')
#    plt.title('Simple ODE model for bacterial growth')
#    plt.show()

    ds = generate_data()
    population = abc.smc(dx_dt, ds, [20.0, 10.0, 5.0])
    plot_solution(population[len(population)-1])
    
    
    
    
    
    
    
