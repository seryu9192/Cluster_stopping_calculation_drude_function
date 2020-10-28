# C1_average_charge.py : calculate C1 average charge using Kaneko model
# Ref: T. Kaneko et al, NIMB315(2013)76

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt, exp
from scipy import integrate
import os
import json

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
param_filename = 'param_C1.json'
param_path = os.path.join(working_dir, param_filename)
output_dir = 'results'

#CONSTANTS
E_C = 300
a_B = 0.529

#carbon cluster
Z = 6
E_0 = 900
v_0 = sqrt(E_0/E_C)

#interatomic dist in A
r = a_B #in au
q_init = 2
qs = [q_init, q_init]

#NUMBER OF ITERATION
NUM_ITER = 100

def gauss(t):
    return exp(-t**2)

def read_parameters(path):
    global E_0, v_0
    with open(path, 'r') as f:
        params = json.loads(f.read())
    E_0 = params["E0"]
    v_0 = sqrt(E_0/E_C)
    return

def main():
    read_parameters(param_path)
    v_bind = 1.045*Z**(2/3)

    # calc average charge
    q = Z*2/sqrt(pi)*integrate.quad(gauss, 0, sqrt(3/8)*v_0/v_bind)[0]
    print("E_0 = {} keV".format(E_0))    
    print("v_0 = {} au".format(v_0))    
    print("q_ave = ", q)

    #save as a text file
    output_filename = '0_average_charge.txt'
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as f:
        tmp = str(q) + '\n'
        f.write(tmp)

if __name__ == '__main__':
    main()
