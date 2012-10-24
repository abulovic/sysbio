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
epsilon = 5
times = (11,24,39,56,75,96,119,144)
steps = 50000

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

def calc_a(sim_th, theta):
    from scipy.stats import uniform
    prior_sim = uniform(-3,3)
    pth_sim = prior_sim.pdf(sim_th)
    pth = prior_sim.pdf(theta)
    from scipy.stats import norm
    jumping_dist_sim = norm(sim_th,1)
    prop_th = jumping_dist_sim.pdf(theta)
    jumping_dist = norm(theta,1)
    prop_simth = jumping_dist_sim.pdf(sim_th)
    likelihood = (0.4 * prop_th) / (0.4 * prop_simth)
    print likelihood
    return min(1, likelihood)
    
def mcmc(ds):
    sigma = 3
    rej_streak = 0
    th1 = np.random.uniform(-3,3)
    th2 = np.random.uniform(-3,3)
    for i in range(5000):
        sim_th1 = np.random.normal(th1,sigma,1)
        sim_th2 = np.random.normal(th2,sigma,1)
        print i, sim_th1, sim_th2, sigma
        el = np.array([sim_th1[0], sim_th2[0]])
        sim_dataset = generate_dataset(el)
        e = euclidian_distance(ds, sim_dataset)
        if e <= epsilon:
            sigma = 0.1
            r = random.randint(0,1)
            a1 = calc_a(sim_th1[0],th1)
            a2 = calc_a(sim_th2[0],th2)
            if r > (1 - a1) or r > (1 - a2):
                th1 = sim_th1[0]
                theta1.append(th1)
                th2 = sim_th2[0]
                theta2.append(th2)
        else:
            rej_streak += 1
            if rej_streak > 50:
                rej_streak = 0
                sigma = 3

if __name__ == "__main__":
    theta = np.array([1,1])
    ds = generate_dataset(theta)
    ds = add_gaussian_noise(ds)
    mcmc(ds)
    print "theta1", theta1
    print "theta2 ", theta2

#    plt.plot(t, X[:,0],'ro',  t, X[:,1], 'bx')     # x = X[:,0] and y = X[:,1]
#    plt.legend(['x(t)','y(t)'])
#    plt.show()

