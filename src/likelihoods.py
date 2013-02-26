#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import abcinfer as abc
from scipy.spatial.distance import euclidean

def log_likelihood(sim_ds, orig_ds):
    sum_dist = 0.
    for i in range(len(sim_ds)):
        dist = euclidean(sim_ds[i], orig_ds[i])
        sum_dist += dist**2
    return -0.5 * sum_dist

def main():
    dx_dt = abc.dx_dt # simplest system with hopf bifurcation
    orig_theta = [3.]
    orig_ds = abc.generate_dataset(dx_dt, orig_theta)
    param_range = np.arange(0, 5, 0.1)
    likelihood_vals = []
    for th in param_range:
        sim_ds = abc.generate_dataset(dx_dt, [th])
        likelihood_vals.append(log_likelihood(sim_ds, orig_ds))
    plt.plot(param_range, likelihood_vals)
    plt.show()
    
if __name__ == "__main__":
    main()