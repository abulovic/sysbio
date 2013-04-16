import numpy as np
from scipy import integrate, array
import pylab as plt

class MitoticOscillator:

    def __init__(self,vi=0.025, Kd=0.01, kd=0.02, vd=0.25, k=0.005,
    v2=1.5, v4=0.5, vm1=3., vm3=1., kc=.5):
        self.vi = vi
        self.Kd = Kd
        self.vd = vd
        self.kd = kd
        self.k = k
        self.v2 = v2
        self.v4 = v4
        self.vm1 = vm1
        self.vm3 = vm3
        self.kc = kc
        self.x0 = [0.1, 0.1, 0.1]
        self.t = 0.
        self.ds = None
        self.times = None

    @staticmethod
    def dx_dt(X, t, vi, Kd, kd, vd, k, v2, v4, vm1, vm3, kc):
        V1 = X[0]/(kc + X[0])*vm1
        V3 = X[1]*vm3
        y = array([vi - vd*X[2]*(X[0]/(kd + X[0])) - kd*X[0], 
                   V1 * ((1-X[1])/(k+(1-X[1]))) - v2*(X[1]/(k+X[1])),
                   V3 * ((1-X[2])/(k+(1-X[2]))) - v4*(X[2]/(k+X[2]))])
        return y
            
            
    def run(self, times=np.linspace(0, 100, 1000)):
        self.times = times
        ds = integrate.odeint(self.dx_dt, self.x0, times, \
                              args=(self.vi, self.Kd, self.kd,\
                                    self.vd, self.k, self.v2, \
                                    self.v4, self.vm1, self.vm3, self.kc))
        self.ds = ds
        
    
    def plot(self):
        if self.ds == None or self.times == None:
            plt.plot(self.times, self.ds)
        plt.show()

    def add_noise(self, sigma=0.5):
        if not self.ds==None:
            x_noise = np.random.normal(0, sigma, np.shape(self.ds))
            self.ds += x_noise
            
            
        
    
