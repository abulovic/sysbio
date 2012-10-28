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

#returns a weighted distribution from population and associated weights
def calc_weighted_distribution(population, weights):
    weighted_population = []
    for i in range(len(population)):
        for k in range(weights[i].astype(int)):
            weighted_population.append(population[i])
    #from scipy.stats.mstats import gmean
    #mu = gmean(weighted_population)
    return weighted_population#norm(mu, 1)

#initialises a list with a number of sublists
def init_list():
    lst = []
    for i in range(param_number):
        lst.append([])
    return lst

#returns an np.array with values drawn from uniform(start, end)
def draw_uniform(start, end):
    theta = np.array([])
    for i in range(param_number):
        theta = np.append(theta, np.random.uniform(start, end))
    return theta

#adds a particle (parameter vector) to corresponding sublists of current population
#th1 goes to sublist for th1, th2 goes to second sublist for th2 and so on
def add_particle_to_list(c_population, theta):
    for i in range(param_number):
        c_population[i].append(theta[i])
    return c_population

#adds a set of weights for current parameter vector to corresponding sublists of weights
#for example w1 of th1 goes to first sublist for weights associated with th1 and so on
def add_weights_to_list(c_weights, wei):
    for i in range(param_number):
        c_weights[i].append(wei[i])
    return c_weights

#adds a population to the general list of populations
#first it gets the weighted population  
def add_population_to_list(population, weights, current_population):
    for i in range(param_number):
        current_population[i] = calc_weighted_distribution(current_population, weights[i])
    population.append(current_population)

#sequential monte carlo
def smc(ds, eps_seq=[30.0]):#, 16.0, 6.0, 5.0, 4.3]):
    populations = []
    weights = init_list()
    current_population = init_list()
    for epsilon in eps_seq:
        if eps_seq.index(epsilon) == 0: #if first population draw from prior
            for i in range(5000):
                sim_theta = draw_uniform(-5,5)
                print sim_theta
                sim_dataset = generate_dataset(sim_theta)
                if euclidian_distance(sim_dataset, ds) < epsilon:
                    current_population = add_particle_to_list(current_population, sim_theta)
                    weights = add_weights_to_list(weights, np.ones(param_number))
        else: #draw from previous population
            pass
        add_population_to_list(populations, weights, current_population)
        current_population = []
        weights = []
    return populations
    
def write_to_file(filename,theta):
    f = open(filename, 'w')
    f.write("theta\n")
    for th in theta:
        f.write(str(th) + ",")
        
if __name__ == "__main__":
    theta = np.array([1,1])
    ds = generate_dataset(theta)
    ds = add_gaussian_noise(ds)
    populations = smc(ds)
    for pop in populations:
        print pop



