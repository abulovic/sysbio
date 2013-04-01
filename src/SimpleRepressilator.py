import scipy, scipy.integrate, scipy.optimize, pylab, copy

class HillRepressilator:
    def __init__(self, alpha=216., alpha0=0.216, beta=5., n=2.,
                 K_m=scipy.log(2.)/120., K_p=scipy.log(2.)/600.,
                 T=50., K_b=1600.,
                 y1=[0.1,0.2,0.3,0.,0.,0.], y0=[0, 0, 0, 2, 1, 3]):
        """Initiate a HillRepressilator (i.e., a Repressilator where the
        protein-DNA interaction is subsumed in a simple Hill function).
        Initialize instance parameters based on the input parameters
        (e.g., self.alpha = alpha).  Also store the bare parameters
        (K_m, K_p, T, K_b) so that you can convert from scaled times
        and concentrations back to bare times and concentrations if
        desired.
        """
        self.alpha = alpha
        self.alpha0 = alpha0
        self.beta = beta
        self.n = n
        #
        self.K_m = K_m
        self.K_p = K_p
        self.T = T
        self.K_b = K_b
        #
        self.t = 0.
        self.traj = None
        self.y = y0
        self.times = None

    @staticmethod
    def dydt(y,t,alpha,n,alpha0,beta):
        """Define the right-hand-side of the Repressilator equations
        for use with scipy.integrate.odeint.  Since this needs to be called
        as if it were a standalone function, and not an instance method
        (which would implicitly assume the first argument is the instance
        self), we declare this as a staticmethod so that self is not passed
        as an implicit argument.  Alternatively, one could define dydt
        outside the HillRepressilator class, but this is cleaner."""
        m = y[:3]
        p = y[3:]
        dm = [(-m[i] + alpha/(1.+(p[(i-1)%3])**n) + alpha0) for i in range(3)]
        dp = [(-beta * (p[i]-m[i])) for i in range(3)]
        return scipy.array(dm+dp)

    def run(self, T, dT=None, nT=100, times=None):
        """Run the Repressilator for the specified amount of time T, returning
        output either for fixed time step dT, or over a specified number
        of timesteps nT, or for a specified array of times.  Store the
        trajectory returned by odeint in the instance variable self.traj,
        concatenating the result to the existing self.traj if a previous
        trajectory had been created."""
        if times is None:
            if dT is None:
                #times = scipy.linspace(self.t, self.t+T, nT)
                times = scipy.arange(0., 50., 0.1)
            else:
                times = scipy.arange(self.t, self.t+T, dT)
        traj = scipy.integrate.odeint(self.dydt, self.y, times, \
                                      args=(self.alpha, self.n, self.alpha0,
                                            self.beta), mxstep=1000)
        if self.traj is None:
            self.traj = traj
            self.y = self.traj[-1]
            self.times = times
            self.t = self.times[-1]
        else:
            self.traj = scipy.concatenate((self.traj, traj))
            self.y = self.traj[-1]
            self.times = scipy.concatenate((self.times, times))
            self.t = self.times[-1]
        return traj[:, :3]#assume only mRna conc known
        
    def reset(self):
        self.t = 0.
        self.traj = None
        self.y = [0.,0.,0.001,0.,0.,0.]
        self.times = None
        self.scale_parameters()
        
    def find_steady_state(self,alpha=None, n=None, alpha0=None):
        """Return the steady-state concentration of mRNAs and proteins
        (i.e., an array of length 6) for the specified parameters
        alpha, n, alpha0.  Use scipy.optimize.fsolve."""
        if alpha is None:
            alpha = self.alpha
        if n is None:
            n = self.n
        if alpha0 is None:
            alpha0 = self.alpha0
        #
        f = lambda m_, alpha_, n_, alpha0_: alpha_/(1.+m_**n_) + alpha0_ - m_
        cstar = scipy.optimize.fsolve(f, 1., args=(alpha, n, alpha0))
        return cstar * scipy.ones(6)
    

    def in_steady_state(self):
        """Return True or False indicating whether the current state of
        the system is in the steady-state.  Do this by evaluating the
        instantaneous time derivative of the equations of motion and
        checking whether the vector norm of the instantaneous velocity
        is sufficiently small (e.g., smaller than 1.0e-3)."""
        eps = 1.0e-03
        dy = self.dydt(self.y, self.t, self.alpha, self.n, self.alpha0,
                       self.beta)
        return scipy.sum(dy*dy) < eps

    def rescale_trajectory(self):
        """Return a scaled trajectory from the data contained in the
        self.traj variable.  This undoes the scaling transformation into
        natural units, producing a 6-component trajectory
        containing raw molecular concentrations."""
        straj = scipy.zeros_like(self.traj)
        for i in range(0,3):
            straj[:,i] = self.traj[:,i] * (self.K_p/(self.T*(self.K_b**(1./self.n))))
        for i in range(3,6):
            straj[:,i] = self.traj[:,i] * (self.K_b**(1./self.n))
        return straj

    def plot(self, show_proteins=True, show_mRNAs=False, rescaled=True):
        """Plot the trajectory of the Repressilator, optionally showing
        either the protein concentrations or the mRNA concentrations,
        in either dimensionless or rescaled units."""
        if rescaled:
            traj = self.rescale_trajectory()
        else:
            traj = self.traj
        if show_proteins:
            pylab.plot(self.times, traj[:,3])
            pylab.plot(self.times, traj[:,4])
            pylab.plot(self.times, traj[:,5])
        if show_mRNAs: 
            pylab.plot(self.times, traj[:,0])
            pylab.plot(self.times, traj[:,1])
            pylab.plot(self.times, traj[:,2])
        pylab.show()
        
    @staticmethod
    def f(beta, X):
        print (((beta+1.)**2)/beta)-((3.*X*X)/(4.+2.*X))
        return (((beta+1.)**2)/beta)-((3.*X*X)/(4.+2.*X))

    def find_unstable_beta(self, alpha):
        pstar = self.find_steady_state(alpha)[-1]
        X = alpha*self.n*(pstar**(self.n-1.))/((1.+pstar**self.n)**2.)
        print pstar, X
        #f = lambda beta, X=X: (((beta+1.)**2)/beta)-((3.*X*X)/(4.+2.*X))
        betastar = scipy.optimize.fsolve(self.f, self.beta, args=(X,))
        return betastar

    def plot_phase_boundary(self):
        alphas = 10.**scipy.linspace(0., 5, 20)
        betas = []
        for alpha in alphas:
            betastar = self.find_unstable_beta(alpha)
            betas.append(betastar)
        pylab.plot(alphas,betas)

    def compute_phase_diagram(self):
        results = {}
        for alpha in 10.**scipy.linspace(0., 5, 20):
            for beta in 10.**scipy.linspace(0.,4.,20):
                print scipy.log10(alpha), scipy.log10(beta)
                self.reset()
                self.alpha = alpha
                self.beta = beta
                self.run(T=200.,nT=10)
                in_ss = self.in_steady_state()
                results[(alpha, beta)] = in_ss
        return results

    def plot_phase_diagram(self, ph_diag):
        x = []
        y = []
        for k,v in ph_diag.items():
            if v:
                x.append(k[0])
                y.append(k[1])
        pylab.plot(x,y,'r.')
        pylab.loglog()
    
    def divide(self):
        new_hr = copy.deepcopy(self)
        mRNA_div = scipy.random.normal(0.5, 0.1, (3,))
        prot_div = scipy.random.normal(0.5, 0.1, (3,))
        self.y[0:3] *= mRNA_div
        self.y[3:6] *= prot_div
        new_hr.y[0:3] *= (1.-mRNA_div)
        new_hr.y[3:6] *= (1.-prot_div)
        return new_hr


def WriteStrajAsDataFile(straj, datafile, t0=1000.):
    output = open(datafile, 'w')
    pass



def L2Cost(p, repressilator, data):
    times = []
    for idx, dataset in enumerate(data):
        times.extend(dataset[:,0])
    model_traj = repressilator.run(max(times), times=times)
    for idx, dataset in enumerate(data):
        times.extend(dataset[:,0])
    
        

# Copyright (C) Cornell University
# All rights reserved.
# Apache License, Version 2.0


