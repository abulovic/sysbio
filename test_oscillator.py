#!/usr/bin/python

from scipy import array
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

def dx_dt(X, t):
    compartment_cell = 1.0
    global_par_VM1 = 3.0
    global_par_VM3 = 1.0
    global_par_Kc = 0.5
    global_par_V1=X[1]*global_par_VM1*(X[1]+global_par_Kc)**(-1)
    global_par_V3=x(2)*global_par_VM3;
    reaction_reaction1_vi=0.025;
    reaction_reaction2_kd=0.01;
    reaction_reaction2=x(1)*compartment_cell*reaction_reaction2_kd;
    reaction_reaction3_vd=0.25;
    reaction_reaction3_Kd=0.02;
    reaction_reaction3=x(1)*compartment_cell*reaction_reaction3_vd*x(3)*(x(1)+reaction_reaction3_Kd)**(-1);
    reaction_reaction4_K1=0.0050;
    reaction_reaction4=compartment_cell*(1- 1*x(2))*global_par_V1*(reaction_reaction4_K1- 1*x(2)+1)**(-1);
    reaction_reaction5_V2=1.5;
    reaction_reaction5_K2=0.0050;
    reaction_reaction5=compartment_cell*x(2)*reaction_reaction5_V2*(reaction_reaction5_K2+x(2))**(-1);
    reaction_reaction6_K3=0.0050;
    reaction_reaction6=compartment_cell*global_par_V3*(1- 1*x(3))*(reaction_reaction6_K3- 1*x(3)+1)**(-1);
    reaction_reaction7_K4=0.0050;
    reaction_reaction7_V4=0.5;
    reaction_reaction7=compartment_cell*reaction_reaction7_V4*x(3)*(reaction_reaction7_K4+x(3))**(-1);

    y = array([(1/(compartment_cell))*(( 1.0 * reaction_reaction1) + (-1.0 * reaction_reaction2) + (-1.0 * reaction_reaction3)), (1/(compartment_cell))*(( 1.0 * reaction_reaction4) + (-1.0 * reaction_reaction5)), (1/(compartment_cell))*(( 1.0 * reaction_reaction6) + (-1.0 * reaction_reaction7))])
    return y

if __name__ == "__main__":
    x0 = np.array([0.01, 0.01, 0.01])
    t = np.linspace(0, 90, 100)
    r = integrate.odeint(dx_dt, t, x0)
    plt.plot(r, 'r-')
    plt.show()
#    plt.plot(t, r, 'r-')
#    plt.plot(t, y, 'g-')
#    plt.plot(t, z, 'b-')
#    plt.show()
    
