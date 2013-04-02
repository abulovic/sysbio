import numpy as np
from scipy import integrate, array
import matplotlib.pyplot as plt

#genetic oscillator driven by a SNIC(Saddle Node Bifurcation on Invariant cycle)
class Oscillator:


    def __init__(self, deg_rate=12., basal_ratio_x=1.58, epsilon=0.05, prot_production=50., \
                 rep_strength = 1., x0=[0., 0.]):
        self.basal_ratio_y = epsilon * basal_ratio_x
        self.basal_ratio_x = basal_ratio_x
        self.epsilon = epsilon
        self.prot_production = prot_production
        self.deg_rate = deg_rate
        self.rep_strength = rep_strength
        self.times = None
        self.ds = None
        self.x0 = x0


    @staticmethod
    def dx_dt(X, t, deg_rate, basal_ratio_x, basal_ratio_y, prot_production, \
                  rep_strength):
        y = array([deg_rate*(basal_ratio_x*((1 + prot_production*X[0]**2)/(1 + X[0]**2 + rep_strength*X[1]**2)) - X[0]),
                   deg_rate*basal_ratio_y*((1 + prot_production*X[0]**2)/(1+X[0]**2)) - X[1]])
        return y

    def run(self, times=np.linspace(0, 15, 100)):
        self.times = times
        ds = integrate.odeint(self.dx_dt, self.x0, self.times, args=(self.deg_rate, self.basal_ratio_x,
                                                                     self.basal_ratio_y, self.prot_production,
                                                                     self.rep_strength))
        self.ds = ds

    def plot(self):
        plt.plot(self.times, self.ds)
        plt.show()

#genetic oscillator driven by a subcritical Hopf Bifurcation
class Oscillator_1:

    def __init__(self, deg_rate=100., basal_ratio_x=1.58, epsilon=0.05, prot_production=50., \
                 rep_strength = 5., x0=[1., 1.]):
        self.basal_ratio_y = epsilon * basal_ratio_x
        self.basal_ratio_x = basal_ratio_x
        self.epsilon = epsilon
        self.prot_production = prot_production
        self.deg_rate = deg_rate
        self.rep_strength = rep_strength
        self.times = None
        self.ds = None
        self.x0 = x0


    @staticmethod
    def dx_dt(X, t, deg_rate, basal_ratio_x, basal_ratio_y, prot_production, \
                  rep_strength):
        y = array([deg_rate*(basal_ratio_x*((1+prot_production*X[0]**2)/(1+X[0]**2))) - X[0]-rep_strength*X[0]*X[1],
                   deg_rate*basal_ratio_y*((1+prot_production*X[0]**2)/(1+X[0]**2))-X[1]])
        return y

    def run(self, times=np.linspace(0, 60, 1000)):
        self.times = times
        ds = integrate.odeint(self.dx_dt, self.x0, self.times, args=(self.deg_rate, self.basal_ratio_x,
                                                                     self.basal_ratio_y, self.prot_production,
                                                                     self.rep_strength))
        self.ds = ds

    def plot(self):
        plt.plot(self.times, self.ds)
        plt.show()

            
        
