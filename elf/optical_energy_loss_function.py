# energy_loss_function.py
#calculate energy loss function of glycine

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
plt.rcParams["font.size"] = 20
plt.rcParams["axes.linewidth"] = 1.5

# glycine parameters
dens = 1.607 #g cm^-3
NA = 6.02E23 #Avogadro
M = 75 
Z = 30 # number of valence electrons per molecule
Z_ave = 4 #average atomic number

#parameter(from Z. Tan et al, Rad.Env.Biophys. 45 2(2006) 135-143 / Z. Tan et al, Rad. Env. Biophys 43 (2004) 173-182)
#f1 param
a = 8018.762592073095
b = 19.927 + 0.9807 * Z_ave
c = 13.741 + 0.3215 * Z_ave

#constants
hbar = 1.0546E-34/1.602E-19 #eV*s

#plasmon energy
w_p = 56410 * (NA*dens/M)**(1/2)
E_p = hbar*w_p

def elf(x):
    return a*x/((x**2-b**2)**2 + (c*x)**2)

def elf_x(x):
    return x*a*x/((x**2-b**2)**2 + (c*x)**2)

#ufunc version of elf
elf_ndarray = np.frompyfunc(elf, 1, 1)

def main():
    print(w_p)
    print(E_p)
    print(hbar)
    print(pi/2*E_p**2*Z)
    
    
    fig = plt.figure(figsize=(6, 5), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    
    x = np.linspace(0, 100, 500)
    y = elf_ndarray(x)
    ax.plot(x, y, color='black')
    
    #setup for xaxis
    x_min, x_max, x_major, x_minor = 0, 100, 20, 5
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.set_xlabel('energy transfer $\hbar\omega_p$ (eV)', fontsize=20)

    #setup for yaxis
    y_min, y_max, y_major, y_minor = 0, 2, 1, 0.2
    y_scale = 1
    y_max *= y_scale
    y_min *= y_scale
    y_major *= y_scale
    y_minor *= y_scale
    
    ax.set_ylim(y_min, y_max)
    # ax.yaxis.set_major_locator(MultipleLocator(y_major))
    # ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: '{:.2f}'.format(y/y_scale)))
    ax.set_ylabel('OELF Im[-1/$\epsilon$($\omega)$] (arb.unit)', fontsize=22)
    
    integ = integrate.quad(elf_x, 0, np.Inf)[0]
    print(integ)

if __name__ == '__main__':
    main()