# Drude-type_ELF_atomic_unit.py

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FuncFormatter
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm


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
plt.rcParams["font.size"] = 14
plt.rcParams["axes.linewidth"] = 1.5

#CONSTANT
Z_ave = 4 #for Glycine
E_0 = 27.2 # Hartree energy

omega_p = 19.6487041570916 / E_0

#ELF
def drude_function(k, w):
    #parameter(from Z. Tan et al, Rad.Env.Biophys. 45 2(2006) 135-143)
    a = 0.2882925381282788
    b = (19.927 + 0.9807 * Z_ave) / E_0
    c = (13.741 + 0.3215 * Z_ave) / E_0
    w_k = w - k**2/2
    if w_k >= 0:
        return w_k/w * a*w_k/((w_k**2-b**2)**2 + (c*w_k)**2)
    else:
        return 0

drude_ufunc = np.frompyfunc(drude_function, 2, 1)

def main():
    #calculation
    x_ma, y_ma, z_ma = 4, 6, 1.2
    x_ti, y_ti = 0.05, 0.05
    x = np.arange(0.01, x_ma, x_ti)
    y = np.arange(0.01, y_ma, y_ti)
    X, Y = np.meshgrid(x, y)

    Z = drude_ufunc(X, Y)
    Z = np.array(Z, np.float64) # cast object type to float64
            
    #plot 3D graph
    fig = plt.figure(figsize=(6, 5), dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm)
    ax.view_init(elev=20, azim=210)
    
    #setup for xaxis
    x_min, x_max, x_major, x_minor = 0, x_ma, 1, 0.5
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.set_xlabel('$k$ (au)', fontsize=16)

    #setup for yaxis
    y_min, y_max, y_major, y_minor = 0, y_ma, 1, 0.5
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(MultipleLocator(y_major))
    ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
    ax.set_ylabel('$\omega$ (au)', fontsize=16)

    #setup for zaxis
    z_min, z_max, z_major, z_minor = 0, z_ma, 1, 0.5
    ax.set_zlim(z_min, z_max)
    ax.zaxis.set_major_locator(MultipleLocator(z_major))
    ax.zaxis.set_minor_locator(MultipleLocator(z_minor))
    ax.set_zlabel('Im[$-1/\epsilon(k,\omega)$]', fontsize=16, )
        
if __name__ == '__main__':
    main()