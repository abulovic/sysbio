#!/usr/bin/python

#simple abc rejector. samples parameter vector from uniform prior, simulates dataset and then compares simulated
#dataset with actual dataset and either rejects or accepts based on the euclidian distance between the two.
#we'll use Lotka-Voltera model -> get dataset with params a=b=1 and add gaussian noise

import scipy
import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import random
import warnings

#global vars used throughout
n = 1000
naccepted = 0
theta1 = []
theta2 = []
epsilon = 4.3
times = (11,24,39,56,75,96,119,144)
steps = 10000

def generate_dataset(theta):
    dataset = np.zeros([8,2])
    t = np.arange(0, 15, 0.1)
    X0 = scipy.array([1,0.5])
    X, info = scipy.integrate.odeint(dx_dt,X0, t,args=(theta[0],theta[1]), full_output=True)
    for i in range(8):
        el = np.array([X[times[i]][0], X[times[i]][1]])
        dataset[i] = el
    return dataset

def add_gaussian_noise(dataset):
    x_noise = np.random.normal(0,0.5,8)
    dataset[:,0] = dataset[:,0] + x_noise
    dataset[:,1] = dataset[:,1] + x_noise
    return dataset
    
#ode system for Lotka-Voltera model
def dx_dt(X,t,a,b):
    y = scipy.array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1]])
    return y

def euclidian_distance(dataset, sim_dataset):
    x = sim_dataset[:,0] - dataset[:,0]
    x **= 2
    x_error = x.sum()
    y = sim_dataset[:,1] - dataset[:,1]
    y **= 2
    y_error = y.sum()
    sq_error = y_error + x_error
    return sq_error

def rejector_algorithm(ds):
    #draw sample from uniform prior in the interval [-10,10]
    for i in range(steps):
        theta = np.random.uniform(-3,3,2)
        print i, theta
        sim_dataset = generate_dataset(theta)
        if euclidian_distance(ds, sim_dataset) <= epsilon:#accept
            theta1.append(theta[0])
            theta2.append(theta[1])
            #naccepted += 1

def calc_a():
    rv = uniform(0,1)
    
#def mcmc(ds):
#    theta = np.zeros(2)
#    init_theta = np.random.uniform(-3,3,2)
#    for i in range(steps):
#        sim_theta = np.random.normal((theta[0]+theta[1])/2,1,2)
#        sim_dataset = generate_dataset(theta)
#        if euclidean_distance(ds, sim_dataset) <= epsilon:
#            theta_prior = uniform(0,1)
#            prior_sim_theta = theta_prior(
        

if __name__ == "__main__":
    theta = np.array([1,1])
    ds = generate_dataset(theta)
    ds = add_gaussian_noise(ds)
    rejector_algorithm(ds)
    plt.hist(theta1)
    plt.hist(theta2)
    print "theta1 " , theta1
    print "theta2 ", theta2
    plt.show()

#    plt.plot(t, X[:,0],'ro',  t, X[:,1], 'bx')     # x = X[:,0] and y = X[:,1]
#    plt.legend(['x(t)','y(t)'])
#    plt.show()

