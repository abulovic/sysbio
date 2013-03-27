#!/usr/bin/env python

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import abcinfer as abc
from scipy.spatial.distance import euclidean
from scipy import integrate, array
#from mpl_toolkits.mplot3d import Axes3D
import sys
import oscillators
import math

def log_likelihood(sim_ds, orig_ds):
    sum_dist = 0.
    for i in range(len(sim_ds)):
        dist = euclidean(sim_ds[i], orig_ds[i])
        sum_dist += dist**2
    return -0.5 * sum_dist

def log_likelihood_fourier(sim_ds, orig_ds):
    return -0.5*abc.fourier_distance(sim_ds, orig_ds)

def log_likelihood_distance(sim_ds, orig_ds):
    pass

    
def main_hopf():
    dx_dt = abc.dx_dt # simplest system with hopf bifurcation
    orig_theta = [3., 1.]
    orig_ds = abc.generate_dataset_full(dx_dt, orig_theta)
    #orig_ds = abc.add_gaussian_noise_full(orig_ds)
    param_range = np.arange(1., 4., 0.1)
    param_range_1 = np.arange(0.5, 3, 0.1)
    likelihood_vals = np.zeros((len(param_range), len(param_range_1)))
    for ind1, th in enumerate(param_range):
        for ind2, th1 in enumerate(param_range_1):
            
            sim_ds = abc.generate_dataset_full(dx_dt, [th, th1])
            l = log_likelihood(sim_ds, orig_ds)
            likelihood_vals[ind1, ind2] = l 
            print th, th1, l

    print "========================"
    i,j = np.unravel_index(likelihood_vals.argmax(), likelihood_vals.shape)
    print "max likelihood ", likelihood_vals[i, j], " at: ()",  i, " ,", j, ")"
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    m1, m2 = np.meshgrid(param_range, param_range_1)
    print np.shape(m1), np.shape(m2), np.shape(likelihood_vals)
    surf = ax.plot_surface(m1, m2, likelihood_vals.T, rstride=1,  cstride=1, cmap=cm.jet, linewidth=0.1, antialiased=True)


    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    ax.set_xlabel('kA')
    ax.set_ylabel('k2')
    ax.set_zlabel('likelihood')
    plt.figure()
    plt.imshow(likelihood_vals)
    plt.show()

def main_hopf_1():
    dx_dt = abc.dx_dt # simplest system with hopf bifurcation
    orig_theta = [2.2]
    orig_ds = abc.generate_dataset(dx_dt, orig_theta)
    #orig_ds = abc.add_gaussian_noise_full(orig_ds)
    param_range = np.arange(1, 4, 0.05)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset(dx_dt, [th])
        likelihood_vals.append(log_likelihood(sim_ds, orig_ds))

    print "max: ", max(likelihood_vals)
    print "max val: ", param_range[likelihood_vals.index(max(likelihood_vals))]
    plt.plot(param_range, likelihood_vals)
    plt.show()

def main_rep():
    orig_theta = [216., 205.]
    orig_ds = abc.generate_dataset_rep(orig_theta)
    param_range = np.arange(0, 250, 5.)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset_rep([th])
        likelihood_vals.append(log_likelihood_fourier(sim_ds, orig_ds))
    plt.plot(param_range, likelihood_vals)
    plt.show()

def main_rep_1():
    orig_theta = [5., 5.]
    orig_ds = abc.generate_dataset_rep(orig_theta)
    #orig_ds = abc.add_gaussian_noise_full(orig_ds)
    param_range = np.arange(3., 7., 0.1)
    param_range_1 = np.arange(3, 7, 0.1)
    likelihood_vals = np.zeros((len(param_range), len(param_range_1)))
    for ind1, th in enumerate(param_range):
        for ind2, th1 in enumerate(param_range_1):
            sim_ds = abc.generate_dataset_rep([th, th1])
            l = log_likelihood(sim_ds, orig_ds)
            likelihood_vals[ind1, ind2] = l 
            print th, th1, l

    print "========================"
    i,j = np.unravel_index(likelihood_vals.argmax(), likelihood_vals.shape)
    print "max likelihood ", likelihood_vals[i, j], " at: ()",  i, " ,", j, ")"
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    m1, m2 = np.meshgrid(param_range, param_range_1)
    print np.shape(m1), np.shape(m2), np.shape(likelihood_vals)
    surf = ax.plot_surface(m1, m2, likelihood_vals.T, rstride=1,  cstride=1, cmap=cm.jet, linewidth=0.1, antialiased=True)


    fig.colorbar(surf, shrink=0.5, aspect=5)
    
    ax.set_xlabel('alpha')
    ax.set_ylabel('beta')
    ax.set_zlabel('likelihood')
    plt.figure()
    plt.imshow(likelihood_vals)
    plt.show()
    
def fisher_information(kA):
    dx_dt = abc.dx_dt # simplest system with hopf bifurcation
    ka_range = np.linspace(2, 3.8, 19)
    datasets = []
    N = len(ka_range)
    h = 0.01
    #produce all observable datasets into D
    for ka in ka_range:
        datasets.append(abc.generate_dataset_full(dx_dt, [ka]))

    sum_partial = 0.
    for dataset in datasets:
        partial_deriv = (log_likelihood(dataset, abc.generate_dataset_full(dx_dt, [kA + h])) -
                         log_likelihood(dataset, abc.generate_dataset_full(dx_dt, [kA])))
        
        sum_partial += partial_deriv**2

    fi  = sum_partial/(4*N*h**2)
    return fi

def fisher_information_function():
    fi = []
    ka_range = np.linspace(2, 3.8, 19)
    for th in ka_range:
        print th
        fi.append(fisher_information(th))

    plt.plot(ka_range, fi)
    plt.show()


def switch(X, t, theta):
    k1 = theta[0]
    k2 = theta[1]
    k3 = theta[2]
    k4 = theta[3]
    y = array([2*k1*X[1]-k2*X[0]**2-k3*X[0]*X[1]-k4*X[0],
               k2*X[0]**2-k1*X[1]])
    return y

def main_lv():
    lv = abc.lv
    orig_theta = [1.]
    t = np.linspace(0, 90, 1000)
    orig_ds = abc.generate_dataset(lv, orig_theta)
    param_range = np.linspace(0, 4, 100)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset(lv, [th])
        likelihood_vals.append(log_likelihood(sim_ds, orig_ds))

    plt.plot(param_range, likelihood_vals)
    #plt.figure()
    #plt.plot(t, orig_ds)
    plt.show()
        

def get_freq_ind(tsignal):
    abs_freq = abs(np.fft.fft(tsignal))
    dom_freqs = []
    for ind, val in enumerate(abs_freq):
        if val > 100:
            dom_freqs.append(ind)

    return dom_freqs
    
def create_comparables(dom_freqs, dom_freqs_sim):
    dom_freqs_comp = np.zeros(len(dom_freqs))
    for val in dom_freqs_sim:
        ind = min(range(len(dom_freqs)), key=lambda i: abs(dom_freqs[i]-val))
        dom_freqs_comp[ind] = val

    return dom_freqs, dom_freqs_comp
    
def shape_sim_score(dom_freqs, dom_freqs_sim):
    gap_score = 5.
    sim_score = 0.
    if len(dom_freqs) > len(dom_freqs_sim):
        seq1, seq_2 = create_comparables(dom_freqs, dom_freqs_sim)
    else:
        seq1, seq_2 = create_comparables(dom_freqs_sim, dom_freqs)

    for val1, val2 in zip(seq1, seq_2):
        if (val1 == 0) ^ (val2 == 0):
            sim_score += gap_score
        else:
            sim_score += abs(val1 - val2)

    return (sim_score / len(seq1)) * 20

def get_shape_sim(orig_ds, sim_ds):
    sim_score = 0.
    #get similarity scores between all signal components of datasets
    num_comps = np.shape(orig_ds)[1]
    for i in range(1):
        dom_freqs = get_freq_ind(orig_ds[:,i])
        dom_freqs_sim = get_freq_ind(sim_ds[:,i])
        sign_sim = shape_sim_score(dom_freqs, dom_freqs_sim)
        sim_score += sign_sim
        
    return sim_score / num_comps
    
if __name__ == "__main__":
    dx_dt = abc.dx_dt
    orig_theta = [3.]
    orig_ds = abc.generate_dataset_full(dx_dt, orig_theta)
    param_range = np.arange(1., 4., .1)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset_full(dx_dt, [th])
        likelihood_vals.append(-get_shape_sim(orig_ds, sim_ds))

    print likelihood_vals
    plt.plot(param_range, likelihood_vals)
    plt.show()
    
