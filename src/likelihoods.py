#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import abcinfer as abc
from scipy.spatial.distance import euclidean
from mpl_toolkits.mplot3d import Axes3D
import sys

def log_likelihood(sim_ds, orig_ds):
    sum_dist = 0.
    for i in range(len(sim_ds)):
        dist = euclidean(sim_ds[i], orig_ds[i])
        sum_dist += dist**2
    return -0.5 * sum_dist

def main_hopf():
    dx_dt = abc.dx_dt # simplest system with hopf bifurcation
    orig_theta = [2.]
    orig_ds = abc.generate_dataset(dx_dt, orig_theta)
    param_range = np.arange(0, 5, 0.1)
    param_range_1 = np.arange(0, 5, 0.1)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset(dx_dt, [th])
        likelihood_vals.append(log_likelihood(sim_ds, orig_ds))
    #X,Y = np.meshgrid(param_range, param_range_1)
    #Z = np.array(likelihood_vals)
    #fig = plt.figure()
    #ax = Axes3D(fig)
    plt.plot(param_range, likelihood_vals)
    plt.show()

def main_rep():
    orig_theta = [216.]
    orig_ds = abc.generate_dataset_hp(orig_theta)
    param_range = np.arange(0, 250, 5.)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset_hp([th])
        likelihood_vals.append(log_likelihood(sim_ds, orig_ds))
    plt.plot(param_range, likelihood_vals)
    plt.show()

def fisher_information(kA):
    dx_dt = abc.dx_dt # simplest system with hopf bifurcation
    ka_range = np.arange(0, 4.5, 0.1)
    datasets = []
    N = len(ka_range)
    h = 0.01
    #produce all observable datasets into D
    for ka in ka_range:
        datasets.append(abc.generate_dataset(dx_dt, [ka]))

    sum_partial = 0.
    for dataset in datasets:
        partial_deriv = (log_likelihood(dataset, abc.generate_dataset(dx_dt, [kA + h])) -
                         log_likelihood(dataset, abc.generate_dataset(dx_dt, [kA])))
        
        sum_partial += partial_deriv**2

    fi  = sum_partial/(4*N*h**2)
    return fi
        
if __name__ == "__main__":
    print fisher_information(5.)