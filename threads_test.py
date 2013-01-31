#!/usr/bin/env python

from multiprocessing import Process, Array
from itertools import izip
import numpy as np
import time
import math
from scipy.stats import norm
from scipy.stats import tstd
import sys

def calculate_var(index, weights, prev_pop, prev_wei):
    prod = []
    for w,th in izip(prev_wei, prev_pop):
        prod.append(w * norm(1.2, tstd(prev_pop)).pdf(th))
    weights[index] = 0.1 / math.fsum(prod)
        

def main():
    weights = [0, 0]
    prev_pop1 = np.random.randn(1000)
    prev_pop2 = np.random.randn(1000)
    prev_wei1 = np.random.randn(1000)
    prev_wei2 = np.random.randn(1000)
    prev_population = [prev_pop1, prev_pop2]
    prev_weights = [prev_wei1, prev_wei2]
    for i in range(len(weights)):
        calculate_var(i, weights, prev_population[i], prev_weights[i])
    for weight in weights:
        print weight

def main_threaded():
    prev_pop1 = np.random.randn(1000)
    prev_pop2 = np.random.randn(1000)
    prev_wei1 = np.random.randn(1000)
    prev_wei2 = np.random.randn(1000)
    weights = Array('d', [0.0, 0.0])
    p1 = Process(target=calculate_var, args=(0, weights,prev_pop1, prev_wei1))
    p2 = Process(target=calculate_var, args=(1, weights,prev_pop2, prev_wei2))
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    for weight in weights:
        print weight
    
if __name__ == "__main__":
    start = time.time()
    main_threaded()
    print "Elapsed Time: %s" % (time.time() - start)
