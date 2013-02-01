#!/usr/bin/env python

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
import math
import stats_util as utils

#global vars used throughout
epsilon = 5.0
data_points = 8
#times = (0, 10, 20, 30, 40, 50, 60, 70)
times = (11, 24, 39, 56, 75, 96, 119, 144)
steps = 100000
param_number = 8
#weighted_mu = np.zeros(param_number)
#sigma = np.zeros((param_number, param_number))

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
    X0 = array([1.0, 0.5])
    X= integrate.odeint(dx_dt, X0, t, args=(theta,),mxhnil=0,hmin=1e-20)
    #plt.plot(t, X)
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
    x_noise = np.random.normal(0, 0.5, data_points)
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
    naccepted = 0
    i = 0
    population = init_list()
    #draw sample from uniform prior in the interval [-10,10]
    while naccepted < 200:
        i += 1
        theta = np.random.uniform(-10,10,param_number)
        sim_dataset = generate_dataset(dx_dt,theta)
        error = euclidian_distance(ds, sim_dataset)
        print i, theta, error, naccepted
        if error <= epsilon:#accept
            naccepted += 1
            population = add_particle_to_list(population, theta)
    return population

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
def mcmc(dx_dt, ds):
    naccepted = 0
    population = init_list()
    sigma = 5
    rej_streak = 0
    counter = 0
    #start of with random values for params taken from uniform prior
    theta = np.random.uniform(-5, 5, param_number)
    while counter < 50000:#steps:
        if len(population[1]) > 500: break
        counter += 1
        sim_theta = draw_from_jumping(theta, sigma)
        sim_dataset = generate_dataset(dx_dt, sim_theta)
        error = euclidian_distance(ds, sim_dataset)
        print counter, sim_theta, sigma, error, naccepted
        if error <= epsilon:
            rej_streak = 0
            sigma = 0.1
            add_particle(population, sim_theta, theta, sigma)
            naccepted += 1
        else:
            rej_streak += 1
            if rej_streak > 10:
                theta = np.random.uniform(-5, 5, param_number)
                rej_streak = 0
                sigma = 1
    print "steps taken ", counter
    return population

#returns a weighted distribution from population and associated weights
def calc_weighted_mean(population, weights):
    wsum = np.zeros(param_number)
    sum_weights = 0.0
    for i in range(len(population)):
        wsum += population[i]*weights[i]
        sum_weights += weights[i]
    return wsum / sum_weights

#returns an np.array with values drawn from uniform(start, end)
def draw_uniform(start, end):
    theta = np.array([])
    for i in range(param_number):
        theta = np.append(theta, np.random.uniform(start, end))
    return theta

def get_pert_sigma(prev_population, sim_theta):
    M = 20
    from scipy.spatial.distance import sqeuclidean
    distances = []
    for p in prev_population:
        distances.append(sqeuclidean(p, sim_theta))

    nearest = []
    indices = [i[0] for i in sorted(enumerate(distances), key=lambda x:x[1])]
    for index in indices[:M]:
        nearest.append(prev_population[index])
    return np.cov(np.vstack(nearest).T)

#to do: perturb particle before returning
def sample_from_previous(prev_population, weights):
    weighted_mu = calc_weighted_mean(prev_population, weights)
    sigma = np.cov(np.vstack(prev_population).T)
    particle = np.random.multivariate_normal(weighted_mu, sigma)
    pert_sigma = get_pert_sigma(prev_population, particle)
    pert_particle = np.random.multivariate_normal(particle, pert_sigma)
    return pert_particle

def calculate_weight(prev_population, prev_weights, sim_theta):
    wsum = 0.0
    sigma = np.cov(np.vstack(prev_population).T)
    mean = utils.colMeans(np.vstack(prev_population))
    for particle, weight in izip(prev_population, prev_weights):
        wsum += weight * utils.dmvnorm(sim_theta, mean, sigma)
    return 0.1 / wsum

def show_histogram(population):
    plt.figure(1)
    for pop in population:
        plt.hist(pop)
        plt.show()
        plt.figure(2)

def norm_weights(weights):
    sum_weights = math.fsum(weights)
    n_weights = [weight/sum_weights for weight in weights]
    return n_weights

def calc_pert_params(prev_population, weights):
    weighted_mu = calc_weighted_mean(prev_population, weights)
    sigma = np.cov(np.vstack(prev_population).T)
    return weighted_mu, sigma

#sequential monte carlo
def smc(dx_dt, ds, eps_seq):
    i = 0
    naccepted = 0
    t = 0
    populations = []
    weights = []
    current_weights = []
    current_population = []
    for epsilon in eps_seq:
        print "population", t
        if eps_seq.index(epsilon) == 0: #if first population draw from prior
            while naccepted < 100:
                i += 1
                sim_theta = draw_uniform(-10,10)
                print i, sim_theta, naccepted
                sim_dataset = generate_dataset(dx_dt, sim_theta)
                if euclidian_distance(sim_dataset, ds) < epsilon:
                    naccepted += 1
                    current_population.append(sim_theta)
                    current_weights.append(1)
        else: #draw from previous population
            #weighted_mu, sigma = calc_pert_params(populations[t-1], weights[t-1])
            while naccepted < 100:
                i += 1
                sim_theta = sample_from_previous(populations[t-1], weights[t-1])
                sim_dataset = generate_dataset(dx_dt, sim_theta)
                error = euclidian_distance(sim_dataset, ds)
                print i, sim_theta, error, naccepted
                if error <= epsilon:
                    naccepted += 1
                    current_population.append(sim_theta)
                    wei = calculate_weight(populations[t-1], weights[t-1], sim_theta)
                    current_weights.append(wei)
        populations.append(current_population)
        weights.append(norm_weights(current_weights))
        current_population = []
        current_weights = []
        t += 1
        naccepted = 0
    return populations
    
def write_to_file(filename,theta):
    f = open(filename, 'w')
    f.write("theta\n")
    for th in theta:
        f.write(str(th) + ",")
        
def plot_solution(population, ds):
    #ds = generate_dataset(dx_dt, theta)
    ti = [t/10 for t in times]
    theta1 = np.array([1,1])
    plt.figure(1)
    theta = utils.colMeans(np.vstack(population))
    X0 = np.array([1, 0.5])
    t = np.arange(0, 15, 0.1)
    X= integrate.odeint(dx_dt, X0, t, args=(theta,))
    Y= integrate.odeint(dx_dt, X0, t, args=(theta1,))
    x,y = X.T
    x1,y1 = Y.T
    plt.figure(3)
    plt.subplot(211)
    plt.plot(t, x, 'r-', label='x(t)')
    plt.plot(t, x1,'g-',label='x(t))')
    plt.plot(ti, ds[:, 0], marker='s', linestyle='', color='g')
    plt.subplot(212)
    plt.plot(t, y, 'b-', label='y(t)')
    plt.plot(t, y1, 'g-', label='y(t)')
    plt.plot(ti, ds[:, 1], marker='^', linestyle='', color='g')
    plt.xlabel('time')
    plt.show()

if __name__ == "__main__":
    theta = [1,1]
    ds = generate_dataset(dx_dt, theta)
    ds = add_gaussian_noise(ds)
    #population = smc(dx_dt, ds, [30.0, 16.0, 6.0, 5.0])
    population = mcmc(dx_dt,ds)
    last_population = population[len(population)-1]
    plot_solution(population, ds)
    

	   
