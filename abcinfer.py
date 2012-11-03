#!/usr/bin/python

#simple abc rejector. samples parameter vector from uniform prior, simulates dataset and then compares simulated
#dataset with actual dataset and either rejects or accepts based on the euclidian distance between the two.
#we'll use Lotka-Voltera model(http://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equation) -> get dataset
#with params a=b=1 and add gaussian noise
#also added a simple mcmc based on metropolis algorithm
#and a sequential monte carlo

from scipy import integrate
from scipy import array
from scipy.stats import uniform
from scipy.stats import norm
from itertools import izip
import numpy as np
np.seterr(all='ignore')
import matplotlib.pyplot as plt
import random

#global vars used throughout
theta1 = []
theta2 = []
epsilon = 4.3
data_points = 8
times = (11,24,39,56,75,96,119,144)
steps = 100000
param_number = 2

def summary(theta):
    from scipy.stats.mstats import gmean
    from scipy.stats.mstats import mode
    return gmean(theta), mode(theta)

#ode system for Lotka-Voltera model
def dx_dt(X,t,theta):
    a = theta[0]
    b = theta[1]
    y = array([a*X[0] - X[0]*X[1], b*X[0]*X[1] - X[1]])
    return y

def generate_dataset(dx_dt, theta):
    dataset = np.zeros([data_points, 2])
    t = np.arange(0, 15, 0.1)
    X0 = array([1,0.5])
    X= integrate.odeint(dx_dt,X0, t, args=(theta,),mxhnil=0,hmin=1e-20)
    for i in range(data_points):
        dataset[i] = create_datapoint(X[times[i]])
    return dataset

#create a datapoint from 
def create_datapoint(data):
    datapoint = np.array([])
    for x in data:
        datapoint = np.append(datapoint, x)
    return datapoint

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

def rejector_algorithm(dx_dt, ds):
    population = init_list()
    #draw sample from uniform prior in the interval [-10,10]
    for i in range(100000):
        theta = np.random.uniform(-10,10,param_number)
        sim_dataset = generate_dataset(dx_dt, array_to_tuple(theta))
        error = euclidian_distance(ds, sim_dataset)
        print i, theta, error
        if error <= epsilon:#accept
            naccepted += 1
            population = add_particle_to_list(population, theta)
    return population

def array_to_tuple(nparray):
    lst = [x for x in nparray]
    return lst

def calc_a(sim_th, theta, sigma):
    prior_sim = uniform(-10,10)
    pth_sim = prior_sim.pdf(sim_th)
    pth = prior_sim.pdf(theta)
    jumping_dist_sim = norm(sim_th,sigma)
    prop_th = jumping_dist_sim.pdf(theta)
    jumping_dist = norm(theta,sigma)
    prop_simth = jumping_dist_sim.pdf(sim_th)
    likelihood = (0.1 * prop_th) / (0.1 * prop_simth)
    return min(1, likelihood)

def add_particle(population, sim_theta, theta, sigma):
    for p, sth, th in izip(population, sim_theta, theta):
        a = calc_a(sth, th, sigma)
        r = random.randint(0,1)
        if r <= a:
            th = sth
            p.append(th)

def draw_from_jumping(theta, sigma):
    sim_theta = []
    for pth in theta:
        sim_theta.append(np.random.normal(pth, sigma))
    return sim_theta

#simple mcmc algorithm creating a separate chain for each parameter
def mcmc(ds):
    population = init_list()
    sigma = 3
    rej_streak = 0
    counter = 0
    #start of with random values for params taken from uniform prior
    theta = np.random.uniform(-5, 5, param_number)
    while counter < steps:
        counter += 1
        sim_theta = draw_from_jumping(theta, sigma)
        sim_dataset = generate_dataset(dx_dt, sim_theta)
        error = euclidian_distance(ds, sim_dataset)
        if error <= epsilon:
            print sim_theta
            rej_streak = 0
            sigma = 0.1
            add_particle(population, sim_theta, theta, sigma)
        else:
            rej_streak += 1
            if rej_streak > 10:
                theta = np.random.uniform(-5, 5, param_number)
                rej_streak = 0
                sigma = 3
    print "steps taken ", counter
    return population

def ftoi(w):
    if w < 1: return 1
    else : return w.astype(int)
    
#returns a weighted distribution from population and associated weights
def calc_weighted_mean(population, weights):
    wsum = 0.0
    sum_weights = 0.0
    wmu = 0.0
    for p, w in izip(population, weights):
        wsum += p*w
        sum_weights += w
    wmu = wsum / sum_weights
    return wmu

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
        current_population[i] = calc_weighted_distribution(current_population[i], weights[i])
    population.append(current_population)

def sample_from_previous(prev_population, weights):
    from scipy.stats import tstd
    theta = np.array([])
    for i in range(param_number):
        weighted_mu = calculate_weighted_mean(prev_population[i], weights[i])
        sigma = tstd(prev_population[i])
        particle = np.random.normal(weighted_mu, sigma)
        pert_particle = np.random.normal(particle, sigma)
        theta = np.append(theta, pert_particle)
    return theta

def calculate_weights(prev_population, prev_weights, sim_theta):
    from scipy.stats import tstd
    weights = np.array([])
    for i in range(param_number):
        rv = uniform(-5, 5)
        prior = rv.pdf(sim_theta[i])
        prod = []
        for w,th in izip(prev_weights[i], prev_population[i]):
            prod.append(w * norm(sim_theta[i], tstd(prev_population[i])).pdf(th))
            weights = np.append(weights, (0.1 / math.fsum(prod)))
    return weights

#sequential monte carlo
def smc(ds, eps_seq=[30.0, 16.0]):
    t = 0
    populations = []
    weights = []
    current_weights = init_list()
    current_population = init_list()
    for epsilon in eps_seq:
        print "population", t
        if eps_seq.index(epsilon) == 0: #if first population draw from prior
            for i in range(100):
                sim_theta = draw_uniform(-5,5)
                print i, sim_theta
                sim_dataset = generate_dataset(dx_dt, sim_theta)
                if euclidian_distance(sim_dataset, ds) < epsilon:
                    current_population = add_particle_to_list(current_population, sim_theta)
                    current_weights = add_weights_to_list(current_weights, np.ones(param_number))
        else: #draw from previous population
            for i in range(100):
                sim_theta = sample_from_previous(populations[t-1], weights[t-1])
                sim_dataset = generate_dataset(sim_theta)
                error = euclidian_distance(sim_dataset, ds)
                print i, sim_theta, error
                if error <= epsilon:
                    current_population = add_particle_to_list(current_population, sim_theta)
                    wei = calculate_weights(populations[t-1], weights[t-1], sim_theta)
                    current_weights = add_weights_to_list(current_weights, wei)
        print "current_weights ", t, " ", current_weights
        populations.append(current_population)
        weights.append(current_weights)
        current_population = init_list()
        current_weights = init_list()
        t += 1
    return populations
    
def write_to_file(filename,theta):
    f = open(filename, 'w')
    f.write("theta\n")
    for th in theta:
        f.write(str(th) + ",")

def plot_solution(population):
    theta1 = np.array([1,1])
    theta = []
    for p in population:
        theta.append(math.fsum(p) / len(p))
    X0 = np.array([1, 0.5])
    t = np.arange(0, 15, 0.1)
    X= integrate.odeint(dx_dt, X0, t, args=(theta,))
    Y= integrate.odeint(dx_dt, X0, t, args=(theta1,))
    x,y = X.T
    x1,y1 = Y.T
    plt.figure(1)
    plt.subplot(211)
    plt.plot(t, x, 'r-', label='x(t)')
    plt.plot(t, x1,'g-',label='x(t))')
    plt.subplot(212)
    plt.plot(t, y, 'b-', label='y(t)')
    plt.plot(t, y1, 'g-', label='y(t)')
    plt.xlabel('time')
    plt.show()
        
if __name__ == "__main__":
    theta1 = np.array([1,1])
    ds = generate_dataset(dx_dt,  theta1)
    ds = add_gaussian_noise(ds)
    population = mcmc(ds)
    plot_solution(population)
    


