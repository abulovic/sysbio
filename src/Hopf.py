import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate, array


def dx_dt(X, t, theta):
    m = theta[0]
    omega = theta[1]
    b = theta[2]
    y = array([m*X[0] + X[0]**3 - X[0]**5, omega + b*X[0]**2])
    return y


def dy_dt(Y, t, theta):
    a = 0.1
    b = theta[0]
    y = array([-Y[0] + a*Y[1] + Y[0]**2*Y[1], b - a*Y[1] - Y[0]**2*Y[1]])
    return y

def dh_dt(X, t, theta):
    m = theta[0]
    y = array([m*X[0] - X[1] + X[0]*X[1]**2, X[0] + m*X[1] + X[1]**3])
    return y

def main_hopf_test():
    X0 = [.1, .1]
    t = np.linspace(0, 500, 1000)
    param_range = np.linspace(-2, 2., 50)
    for th in param_range:
        print th
        ds = integrate.odeint(dy_dt, X0, t,  args=([th],))
        plt.plot(ds[:, 0], ds[:, 1])
        plt.figure()
        plt.plot(t, ds)
        plt.show()
        plt.close()

def hopf(X, t, theta):
    mi = theta[0]
    omega = theta[1]
    b = theta[2]
    y = array([mi*X[0] - omega*X[1] - X[0]*X[1]**2 - b*X[1]*X[0]**2 - X[0]**3 - b*X[1]**3,
               mi*X[1] + omega*X[0] - X[1]*X[0]**2 + b*X[0]*X[1]**2 - X[1]**3 + b*X[0]**3])
    return y

def snic(X, t, theta):
    mi = theta[0]
    o = theta[1]
    b = theta[2]
    y = array([mi*X[0] - o*X[1] + X[0]*X[1]**2 - b*X[0]**2*X[1] + X[0]**3 - b*X[1]**3 - 2*X[0]**3*X[1]**2 - X[0]*X[1]**4 - X[0]**5,
               mi*X[1] + o*X[0] + X[1]*X[0]**2 + b*X[1]**2*X[0] - X[1]*X[0]**2 + b*X[0]**3 + X[1]**3 - 2*X[0]**2*X[1]**3 - X[1]**5])
    return y

def main_snic():
    X0 = [.5, .5]
    t = np.linspace(0, 100, 1000)
    theta = [-0.15, 1, 1]
    ds = integrate.odeint(snic, X0, t,  args=(theta,))
    plt.plot(ds[:, 0], ds[:, 1])
    plt.figure()
    plt.plot(t, ds)
    plt.show()
    
def main_hopf():
    X0 = [.1, .1]
    t = np.linspace(0, 100, 1000)
    theta = [-0.2, 1, 1]
    ds = integrate.odeint(hopf, X0, t,  args=(theta,))
    plt.plot(ds[:, 0], ds[:, 1])
    plt.figure()
    plt.plot(t, ds)
    plt.show()
    
def main_polar():
    X0 = [1, 1]
    t = np.arange(0, 15, 0.1)
    theta = [3, 1.5, 1.5]
    ds = integrate.odeint(dx_dt, X0, t,  args=(theta,))
    fig1 = plt.figure()
    ax1 = fig1.add_axes([0.1,0.1,0.8,0.8],polar=True)
    ax1.plot(ds[:, 0], ds[:, 1])
    plt.show()


if __name__ == "__main__":
    main_snic()
