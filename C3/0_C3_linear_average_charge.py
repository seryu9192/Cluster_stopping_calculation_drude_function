# C3_linear_average_charge.py : calculate linear C3 average charge using Kaneko model(Self-consistent calculation)
# Ref: T. Kaneko et al,NIMB315(2013)76, T. Kaneko, PRA 66 (2002) 052901

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from numpy import pi, sqrt, exp
from scipy import integrate
import os
import sys
import json
sys.path.append('../')
from common_library import *

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

#filepath
working_dir = './'
param_filename = 'param_C3.json'
param_path = os.path.join(working_dir, param_filename)
output_dir = 'results'

#parameters
E = 0
v = 0

#interatomic dist in A
r = 0 #in atomic unit
q_init = 2
qs = [q_init, q_init, q_init]

#NUMBER OF ITERATION
NUM_ITER = 100

def read_parameters(path):
    global E, v, r
    with open(path, 'r') as f:
        params = json.loads(f.read())
    E = params["E0"]
    v = sqrt(E/E_CARBON)
    r = params["r"]["linear"]
    return

def calc_q(i):
    global qs, v, r   
    v_b = 1.092 * Z_CARBON**(4/3)
    for j in range(len(qs)):
        if j == i:
            continue
        v_b += 2*qs[j]/(r[str_ind(i, j)]/a_0) #r:atomic unit!
    v_b = sqrt(v_b)
    y = sqrt(3/8) * v / v_b
    q_new = Z_CARBON * 2/sqrt(pi) * integrate.quad(lambda x: exp(-x**2), 0, y)[0]
    qs[i] = q_new

def main():
    read_parameters(param_path)

    #results for showing iteration
    ys_1 = [q_init]
    ys_2 = [q_init]
    ys_3 = [q_init]

    #self-consistent calculation
    for it in range(NUM_ITER):
        for i in range(len(qs)):
            calc_q(i)
        ys_1.append(qs[0])
        ys_2.append(qs[1])
        ys_3.append(qs[2])

    #show results
    print('E = {}'.format(E))
    print('v = {}'.format(v))
    print(r)
    print(qs)
    
    #save as a json file
    output_filename = 'C3_average_charge.json'
    output_path = os.path.join(output_dir, output_filename)
    #read json file if exists
    try:
        with open(output_path, 'r') as f:
            dat = json.loads(f.read())
    except Exception:
        dat = {}
        os.makedirs(output_dir, exist_ok=True)
    #update json
    #setdefault
    dat.setdefault(f"{E}",{})
    dat[f"{E}"].setdefault("linear", {})

    for i, q in enumerate(qs):
        dat[f"{E}"]["linear"][f"{i}"] = q
    with open(output_path, 'w') as f:
        f.write(json.dumps(dat, indent=4))
    
    # plot
    xs = np.arange(0, len(ys_1))
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xs, ys_1, label='ion1')    
    ax.plot(xs, ys_2, label='ion2')    
    ax.plot(xs, ys_3, label='ion3')
    #x axis
    x_min, x_max, x_major, x_minor = 0, 10, 2, 1
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.set_xlabel('iteration')
    #y axis
    ax.set_ylabel('$q$')
    ax.text(0.02, 0.9, 'linear C$_3$($E$ = {} keV/atom)'.format(E), transform=ax.transAxes)
    ax.legend(fontsize=15)

if __name__ == '__main__':
    main()
