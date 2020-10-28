# C1_average_charge.py : calculate C1 average charge using Kaneko model
# Ref: T. Kaneko et al,NIMB315(2013)76

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt, exp
from scipy import integrate

plt.rcParams["font.family"] = "Arial"
plt.rcParams["xtick.direction"] = "out"             
plt.rcParams["ytick.direction"] = "in"              
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True          
plt.rcParams["xtick.major.width"] = 1.5             
plt.rcParams["ytick.major.width"] = 1.5             
plt.rcParams["xtick.minor.width"] = 1.0             
plt.rcParams["ytick.minor.width"] = 1.0             
plt.rcParams["xtick.major.size"] = 6                
plt.rcParams["ytick.major.size"] = 6                
plt.rcParams["xtick.minor.size"] = 3
plt.rcParams["ytick.minor.size"] = 3
plt.rcParams["font.size"] = 18                      
plt.rcParams["axes.linewidth"] = 1.5                

#carbon cluster
Z = 6
E_C = 300
E_0 = 900
v_0 = sqrt(E_0/E_C)

def gauss(t):
    return exp(-t**2)

def main():
    #carbon ions
    Z = 6
    v_bind = 1.045*Z**(2/3)

    # calc average charge
    vs = np.linspace(0, 10, 100)
    ys = [Z*2/sqrt(pi)*integrate.quad(gauss, 0, sqrt(3/8)*v/v_bind)[0] for v in vs]
    #plot
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(vs, ys)    
    #x axis
    ax.set_xlim(0, 10)
    ax.set_xlabel('$v/v_0$')
    #y axis
    ax.set_ylim(0, 6)
    ax.set_ylabel('$q_{ave}$')
    #text
    ax.text(0.05, 0.9, 'average charge of C$_1$', transform=ax.transAxes)
    q_ave = Z*2/sqrt(pi)*integrate.quad(gauss, 0, sqrt(3/8)*v_0/v_bind)[0] 
    
    print(q_ave)
    
if __name__ == '__main__':
    main()
