#!/usr/bin/python

#simple abc rejector. samples parameter vector from uniform prior, simulates dataset and then compares simulated
#dataset with actual dataset and either rejects or accepts based on the euclidian distance between the two.
#we'll use Lotka-Voltera model(http://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equation) -> get dataset
#with params a=b=1 and add gaussian noise
#also added a simple mcmc based on metropolis algorithm

from scipy import integrate
from scipy import array
from scipy.stats import uniform
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import random

#global vars used throughout
theta1 = []
theta2 = []
epsilon = 4.3
data_points = 8
times = (11,24,39,56,75,96,119,144)
steps = 60000
param_number = 2

def summary(theta):
    from scipy.stats.mstats import gmean
    from scipy.stats.mstats import mode
    return gmean(theta), mode(theta)

#ode system for Lotka-Voltera model
def dx_dt(X,t,a,b):
    y = array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1]])
    return y

def generate_dataset(theta):
    dataset = np.zeros([data_points, np.size(theta)])
    t = np.arange(0, 15, 0.1)
    X0 = array([1,0.5])
    X, info = integrate.odeint(dx_dt,X0, t,args=(theta[0],theta[1]), full_output=True)
    for i in range(data_points):
        dataset[i] = np.array([X[times[i]][0], X[times[i]][1]])
    return dataset

def add_gaussian_noise(dataset):
    x_noise = np.random.normal(0,0.5,data_points)
    for i in range(dataset.shape[1]):
        dataset[:,i] = dataset[:,i] + x_noise
    return dataset

def euclidian_distance(dataset, sim_dataset):
    sq_error = 0
    from scipy.spatial.distance import sqeuclidean
    for i in range(dataset.shape[1]):
        sq_error += sqeuclidean(sim_dataset[:,i],dataset[:,i])
    return sq_error

def rejector_algorithm(ds):
    naccepted = 0
    #draw sample from uniform prior in the interval [-10,10]
    for i in range(steps):
        theta = np.random.uniform(-5,5,2)
        print i, theta, naccepted
        sim_dataset = generate_dataset(theta)
        if euclidian_distance(ds, sim_dataset) <= epsilon:#accept
            naccepted += 1
            theta1.append(theta[0])
            theta2.append(theta[1])

def calc_a(sim_th, theta, sigma):
    prior_sim = uniform(-5,5)
    pth_sim = prior_sim.pdf(sim_th)
    pth = prior_sim.pdf(theta)
    jumping_dist_sim = norm(sim_th,sigma)
    prop_th = jumping_dist_sim.pdf(theta)
    jumping_dist = norm(theta,sigma)
    prop_simth = jumping_dist_sim.pdf(sim_th)
    likelihood = (0.1 * prop_th) / (0.1 * prop_simth)
    return min(1, likelihood)

#simple mcmc algorithm creating a separate chain for each parameter
def mcmc(ds):
    naccepted = 0
    sigma = 3
    rej_streak = 0
    i = 0
    #start of with random values for params taken from uniform prior
    th1 = np.random.uniform(-5,5)
    th2 = np.random.uniform(-5,5)
    for i in range(steps):
        sim_th1 = np.random.normal(th1,sigma,1)[0]
        sim_th2 = np.random.normal(th2,sigma,1)[0]
        print i,sim_th1, sim_th2, sigma, naccepted
        sim_dataset = generate_dataset(np.array([sim_th1, sim_th2]))
        if euclidian_distance(ds, sim_dataset) <= epsilon:
            rej_streak = 0
            sigma = 0.1
            r = random.randint(0,1)
            a1 = calc_a(sim_th1,th1,sigma)
            a2 = calc_a(sim_th2,th2,sigma)
            if r <= a1:
                naccepted += 1
                th1 = sim_th1
                theta1.append(th1)
            if r <= a2:
                th2 = sim_th2
                theta2.append(th2)
        else:
            rej_streak += 1
            if rej_streak > 10:
                rej_streak = 0
                sigma = 3
                th1 = np.random.uniform(-5,5)
                th2 = np.random.uniform(-5,5)

#returns a distribution from population and associated weights
def calc_weighted_distribution(population, weights):
    weighted_population = []
    for i in range(len(population)):
        for k in range(weights[i]):
            weighted_population.append(population[i])

    from scipy.stats.mstats import gmean
    mu = gmean(weighted_population)
    return norm(mu, 1)
    
def smc(ds, eps_seq=[30.0]):#, 16.0, 6.0, 5.0, 4.3]):
    current_theta = np.array(params_number)
    previous_population = []
    previous_weights = []
    current_population = []
    population_weights = []
    prior_dist = uniform(-1,1)
    for epsilon in eps_seq:
        if eps_seq.index(epsilon) == 0: #if first population draw from prior
            #get the simulated points
            for i in range(params_number):
                current_theta[i] = np.random.uniform(-1,1)
            sim_dataset = generate_dataset(current_theta)
            if euclidian_distance(sim_dataset, ds) >= epsilon:
                continue
            else:
              current_population.append(current_theta)
              population_weights.append(1)
        else: #draw from previous population
            pass

def write_to_file(filename,theta):
    f = open(filename, 'w')
    f.write("theta\n")
    for th in theta:
        f.write(str(th) + ",")
        
if __name__ == "__main__":
    theta = np.array([1,1])
    ds = generate_dataset(theta)
    ds = add_gaussian_noise(ds)
    mcmc(ds)
    write_to_file("theta1.txt", theta1)
    write_to_file("theta2.txt", theta2)
    print theta1
    print theta2



