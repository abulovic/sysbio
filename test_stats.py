#!/usr/bin/python

import statistics

if __name__ == "__main__":
    a = [1, 2, 3]
    b = [3, 4, 5]
    m = statistics.multi_mu([a,b])
    print m
