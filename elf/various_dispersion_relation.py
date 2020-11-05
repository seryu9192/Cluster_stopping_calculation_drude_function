# various_dispersion_relation.py
# plot various dispersion relation between w and k (in atomic unit)

import os
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FuncFormatter
from scipy import integrate

#general parameters of gragh
plt.rcParams["font.family"] = "Arial"                
plt.rcParams["xtick.direction"] = "out"               
plt.rcParams["ytick.direction"] = "out"
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
plt.rcParams["font.size"] = 16
plt.rcParams["axes.linewidth"] = 1.5


def main():
    fig = plt.figure(figsize=(6, 5), dpi=300)
    ax = fig.add_subplot(1, 1, 1)

    #Fermi velocity (Glycine)
    v_f = 1.06979565420349
    #beam velocity
    v = 3**(1/2)
    
    #omega
    w = 1.5

    # x:k, y:w_k
    x = np.linspace(0, 10, 500)
    
    k_1 = v *(1-(1-2*w/v**2)**(1/2))
    k_2 = v *(1+(1-2*w/v**2)**(1/2))
    
    print(k_1)
    print(k_2)
    
    #calc various dispersion relation
    y_0 = np.array([w]*500) 
    y_1 = w + x**2/2
    y_2 = (w**2 + (x**2/2)**2)**(1/2)
    y_3 = (w**2 + 3/5*v_f**2 * x**2 + (x**2/2)**2)**(1/2)
    
    ax.plot(x, y_0, label="$\omega_k = \omega$", linewidth=1)
    ax.plot(x, y_1, label="$\omega_k = \omega + k^2/2$", linewidth=1)
    ax.plot(x, y_2, label="$\omega_k^2 = \omega^2 + (k^2/2)^2$", linewidth=1)
    ax.plot(x, y_3, label="$\omega_k^2 = \omega^2 +3/5v_F^2k^2 + (k^2/2)^2$", linewidth=1)
    
    ax.text(0.7, 0.9, f"$\omega$ = {w} au", transform=ax.transAxes)
    ax.legend(fontsize=13, frameon=False)
    
    #setup for xaxis
    x_min, x_max, x_major, x_minor = 0, 4, 1, 0.5
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.set_xlabel('momentum transfer $k$ (au)', fontsize=20)

    #setup for yaxis
    y_min, y_max, y_major, y_minor = 0, 10, 5, 1
    y_scale = 1
    y_max *= y_scale
    y_min *= y_scale
    y_major *= y_scale
    y_minor *= y_scale
    
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(MultipleLocator(y_major))
    ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
    # ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: '{:.2f}'.format(y/y_scale)))
    ax.set_ylabel('$\omega_k$(au)', fontsize=20)
    

if __name__ == '__main__':
    main()