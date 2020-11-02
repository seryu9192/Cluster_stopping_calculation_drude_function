# C1_average_charge.py : calculate C1 average charge using Kaneko model
# Ref: T. Kaneko et al, NIMB315(2013)76

import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sqrt, exp
from scipy import integrate
import os
import sys
import json
sys.path.append('../')
from common_library import *

#filepath
working_dir = './'
param_filename = 'param_C1.json'
param_path = os.path.join(working_dir, param_filename)
output_dir = 'results'

#parameters
E = 0
v = 0

#interatomic dist in A
r = 0 #in au
q_init = 2
qs = [q_init]

#NUMBER OF ITERATION
NUM_ITER = 100

def read_parameters(path):
    global E, v
    with open(path, 'r') as f:
        params = json.loads(f.read())
    E = params["E0"]
    v = sqrt(E/E_CARBON)
    return

def main():
    read_parameters(param_path)

    #average speed of bound electrons
    v_bind = 1.045*Z_CARBON**(2/3)

    # calc average charge
    q = Z_CARBON * 2/sqrt(pi) * integrate.quad(lambda x: exp(-x**2), 0, sqrt(3/8)*v/v_bind)[0]
    print("E = {} keV".format(E))
    print("v = {} au".format(v))    
    print("q_ave =", q)

    #save as a text file
    output_filename = '0_C1_average_charge.txt'
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as f:
        tmp = str(q) + '\n'
        f.write(tmp)

if __name__ == '__main__':
    main()
