import numpy as np
from scipy import integrate, array
import pylab as plt


class VanderpolOscillator:

    def __init__(self, B=10., d=2.):
        self.B = B
        self.d = d
        self.times = None
        self.ds = None
        self.x0 = [.01, .01]


    @staticmethod
    def dx_dt(X, t, B, d):
        y = array([X[1], -(B*X[0]**2 - d)*X[1] - X[0]])# + 2*np.cos(1.7*t)])
        return y

    def run(self, times=np.linspace(0, 50, 50)):
        self.times = times
        ds = integrate.odeint(self.dx_dt, self.x0, self.times, args=(self.B, self.d))
        self.ds = ds

    def plot(self):
        if self.ds == None or self.times == None:
            plt.plot(self.times, self.ds)
        plt.show()

    def add_noise(self, sigma=0.5):
        if not self.ds == None:
            x_noise = np.random.normal(0, sigma, np.shape(self.ds))
            self.ds += x_noise