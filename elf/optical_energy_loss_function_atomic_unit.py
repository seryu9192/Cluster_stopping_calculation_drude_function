# Drude-type_function_atomic_unit.py

import os
import numpy as np
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
plt.rcParams["font.size"] = 20
plt.rcParams["axes.linewidth"] = 1.5

#CONSTANT
Z_ave = 4 #for Glycine
E_0 = 27.2 # Hartree energy

omega_p = 19.6487041570916 / E_0

def drude_function(x):
    #parameter(from Z. Tan el al, Rad.Env.Biophys. 45 2(2006) 135-143)
    a = 0.2882925381282788
    b = (19.927 + 0.9807 * Z_ave) / E_0
    c = (13.741 + 0.3215 * Z_ave) / E_0
    return a*x/((x**2-b**2)**2 + (c*x)**2)
    

def main():
    fig = plt.figure(figsize=(6, 5), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    
    x = np.linspace(0, 10, 500)
    y = drude_function(x)
    ax.plot(x, y, color='black')
    
    #setup for xaxis
    x_min, x_max, x_major, x_minor = 0, 6, 1, 0.5
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.set_xlabel('energy transfer $\hbar\omega_p$ (au)', fontsize=20)

    #setup for yaxis
    y_min, y_max, y_major, y_minor = 0, 1.5, 0.5, 0.1
    y_scale = 1
    y_max *= y_scale
    y_min *= y_scale
    y_major *= y_scale
    y_minor *= y_scale
    
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(MultipleLocator(y_major))
    ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: '{:.1f}'.format(y/y_scale)))
    ax.set_ylabel('OELF (arb.unit)', fontsize=22)

    integ = integrate.quad(lambda x: drude_function(x)*x, 0, np.Inf)[0]
    print("integ = ", integ)
    print("pi/2*w**2 = ", np.pi/2*omega_p**2)
    print("pi/2*w**2 / integ = ", (np.pi/2*omega_p**2)/integ)
        
if __name__ == '__main__':
    main()