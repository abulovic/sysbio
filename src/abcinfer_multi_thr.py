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
import itertools
from scipy.stats.mstats import mquantiles
from multiprocessing import Process, Queue
import sys
import oscillators

#global vars used throughout
epsilon = 5.0
data_points = 8
#times = (0, 10, 20, 30, 40, 50, 60, 70)
times = (11, 24, 39, 56, 75, 96, 119, 144)
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
    #init = np.array([2.0, 5.0, 3.0])
    init = np.array([1, 0.5])
    t = np.arange(0, 15, 0.1)
    X= integrate.odeint(dx_dt, init, t, args=(theta,),mxhnil=0,hmin=1e-20)
    for i in xrange(data_points):
        dataset[i] = create_datapoint(X[times[i]])
    return dataset

#create a datapoint from 
def create_datapoint(data):
    data_n = len(data)
    datapoint = np.zeros(data_n)
    for i in xrange(data_n):
        datapoint[i] = data[i]
    return datapoint

def add_gaussian_noise(dataset):
    x_noise = np.random.normal(0, 0.5, data_points)
    for i in xrange(dataset.shape[1]):
        dataset[:,i] = dataset[:,i] + x_noise
    
    return dataset

def generate_dataset_full(dx_dt, theta):
    init = np.array([1, 0.5])
    #init = np.array([2.0, 5.0, 3.0])
    t = np.arange(0, 15, 0.1)
    X= integrate.odeint(dx_dt, init, t, args=(theta,),mxhnil=0,hmin=1e-20)
    return X

def add_gaussian_noise_full(dataset):
    x_noise = np.random.normal(0, 0.5, np.shape(dataset))
    return dataset + x_noise

def euclidian_distance(dataset, sim_dataset):
    sq_error = 0
    from scipy.spatial.distance import sqeuclidean
    for i in xrange(dataset.shape[1]):
        sq_error += sqeuclidean(sim_dataset[:,i],dataset[:,i])
    return sq_error

#returns the average distance between the signals in the 2 datasets
def fourier_distance(dataset, sim_dataset):
    sum_ferr = 0.
    signals = np.shape(dataset)[1]
    for i in xrange(signals):
        sum_ferr += oscillators.fourier_compare(dataset[:, i], sim_dataset[:, i])
    return sum_ferr / signals

def fitness(dataset, sim_dataset):
    global eta
    fitness = (eta*fourier_distance(dataset, sim_dataset) +
               (1-eta)*euclidean(dataset, sim_dataset))
    return fitness / 2.

#adds a particle (parameter vector) to corresponding sublists of current population
#th1 goes to sublist for th1, th2 goes to second sublist for th2 and so on
def add_particle_to_list(c_population, theta):
    for i in range(param_number):
        c_population[i].append(theta[i])
    return c_population

#initialises a list with a number of sublists
def init_list():
    lst = []
    for i in range(param_number):
        lst.append([])
    return lst

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
    for i in xrange(len(population)):
        wsum += population[i]*weights[i]
        sum_weights += weights[i]
    return wsum / sum_weights

#returns an np.array with values drawn from uniform(start, end)
def draw_uniform(start, end):
    theta = np.array([])
    for i in xrange(param_number):
        theta = np.append(theta, np.random.uniform(start, end))
    return theta

def get_pert_sigma(prev_population, sim_theta):
    M = (int) (len(prev_population) * 0.2)
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
    pert_sigma = 2 * sigma#get_pert_sigma(prev_population, particle)
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

def get_population_prior(num_particles,dx_dt, a, b, qd, qw, qp, epsilon, ds):
    i = 0
    naccepted = 0
    while naccepted < num_particles:
        i += 1
        sim_theta = draw_uniform(a,b)
        print i, sim_theta, naccepted
        sim_dataset = generate_dataset_full(dx_dt, sim_theta)
        error = fourier_distance(ds, sim_dataset)
        if error < epsilon:
            qd.put(error)
            naccepted += 1
            qp.put(sim_theta)
            qw.put(1)

def dump_queue(queue):
    result = []
    for i in range(queue.qsize()):
        result.append(queue.get())

    return result

def dump_list_queues(queue_list):
    current_populations = []
    current_population = []
    for pop in queue_list:
        current_populations.append(dump_queue(pop))
    all_pop = itertools.chain(*current_populations)
    for p in all_pop: current_population.append(p)
    return current_population

def get_population_previous(num_particles, prev_population, prev_weights, dx_dt, epsilon, qd, qp, qw, ds):
    i, naccepted = 0, 0
    while naccepted < num_particles:
        i += 1
        sim_theta = sample_from_previous(prev_population, prev_weights)
        sim_dataset = generate_dataset_full(dx_dt, sim_theta)
        error = fourier_distance(ds, sim_dataset)
        print i, sim_theta, error, naccepted, epsilon
        if error <= epsilon:
            qd.put(error)
            naccepted += 1
            qp.put(sim_theta)
            wei = calculate_weight(prev_population, prev_weights, sim_theta)
            qw.put(wei)
    
#sequential monte carlo
def smc(dx_dt, ds, eps_seq):
    i = 0
    naccepted = 0
    t = 0
    populations = []
    weights = []
    current_weights = []
    current_population = []
    distances_prev = []
    epsilon = eps_seq[0]
    prev_epsilon = eps_seq[0]
    num_particles = 100
    num_threads = 1
    while True:
        print "population", t
        distances_lst = []
        populations_lst = []
        weights_lst = []
        num_particles_thread = (int) (num_particles / num_threads) #divide the particles equally among threads
        threads = []
        for i in xrange(num_threads):
            qd = Queue()
            qp = Queue()
            qw = Queue()            
            distances_lst.append(qd)
            populations_lst.append(qp)
            weights_lst.append(qw)
            if t == 0:
                p = (Process(target=get_population_prior, args=(num_particles_thread,
                                                                dx_dt, -10, 10, qd, qw,
                                                                qp, epsilon, ds,)))
            else:
                p = (Process(target=get_population_previous, args=(num_particles_thread,
                                                                   populations[t-1],
                                                                   weights[t-1], dx_dt,
                                                                   epsilon, qd, qp, qw, ds,)))
            p.start()
            threads.append(p)
        for thread in threads: thread.join()
        distances_prev = dump_list_queues(distances_lst)
        current_weights = dump_list_queues(weights_lst)
        current_population = dump_list_queues(populations_lst)
        
        populations.append(current_population)
        weights.append(norm_weights(current_weights))
        epsilon = mquantiles(distances_prev, prob=[0.1, 0.25, 0.5, 0.75])[1]
        if prev_epsilon - epsilon < 0.05: break
        else: prev_epsilon = epsilon
        current_population = []
        current_weights = []
        distances_prev = []
        t += 1
        naccepted = 0
    return populations
    
def write_to_file(filename,theta):
    f = open(filename, 'w')
    f.write("theta\n")
    for th in theta:
        f.write(str(th) + ",")
        
def plot_solution(population, ds):
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

def main():
    theta = [1,1]
    ds = generate_dataset(dx_dt, theta)
    ds = add_gaussian_noise(ds)
    population = smc(dx_dt, ds, [300.0, 16.0, 6.0, 5.0, 4.3])
    last_population = population[-1:]
    plot_solution(last_population, ds)
    
if __name__ == "__main__":
    main()
    

	   
