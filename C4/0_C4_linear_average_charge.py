# C4_linear_average_charge.py : calculate linear C4 average charge using Kaneko model(Self-consistent calculation)
# Ref: T. Kaneko et al,NIMB315(2013)76, T. Kaneko, PRA 66 (2002) 052901

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
param_filename = 'param_C4_linear.json'
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
qs = [q_init, q_init, q_init, q_init]

#NUMBER OF ITERATION
NUM_ITER = 100

def str_ind(i, j):
    res = str(i) + str(j)
    if res > res[::-1]:res = res[::-1]
    return res

def gauss(t):
    return exp(-t**2)

def read_parameters(path):
    global E_0, v_0, r, target
    with open(path, 'r') as f:
        params = json.loads(f.read())
    E_0 = params["E0"]
    v_0 = sqrt(E_0/E_C)
    r = params["r"]
    return

def calc_q(i):
    global Z, qs, v_0, r   
    v_b = 1.092 * Z**(4/3)
    for j in range(len(qs)):
        if j == i:
            continue
        v_b += 2*qs[j]/(r[str_ind(i, j)]/a_B) #r:atomic unit!
    v_b = sqrt(v_b)
    y = sqrt(3/8) * v_0 / v_b
    q_new = Z*2/sqrt(pi)*integrate.quad(gauss, 0, y)[0]
    qs[i] = q_new

def main():
    global Z, qs, v_0, r   
    read_parameters(param_path)
    ys_1 = [q_init]
    ys_2 = [q_init]
    ys_3 = [q_init]
    ys_4 = [q_init]
    for it in range(NUM_ITER):
        for i in range(len(qs)):
            calc_q(i)
        ys_1.append(qs[0])
        ys_2.append(qs[1])
        ys_3.append(qs[2])
        ys_4.append(qs[3])
    print('E = {}'.format(E_0))
    print('v = {}'.format(v_0))
    print(r)
    print(qs)

    output_filename = '0_C4_linear_average_charge.txt'
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(output_dir, exist_ok=True)
    with open(output_path, 'w') as f:
        tmp = ''
        for q in qs:
            tmp += str(q) + '\n'
        f.write(tmp)
    
    # # plot
    # xs = np.arange(1, len(ys_1)+1)
    # fig = plt.figure(dpi=300)
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(xs, ys_1, label='C1')    
    # ax.plot(xs, ys_2, label='C2')    
    # ax.plot(xs, ys_3, label='C3')
    # ax.plot(xs, ys_4, label='C4')
    # ax.legend()
    # #x axis
    # ax.set_xlabel('$iteration$')
    # #y axis
    # ax.set_ylabel('$q$')
    # ax.text(0.1, 0.9, 'linear C$_4$', transform=ax.transAxes)
    # ax.legend()

if __name__ == '__main__':
    main()
