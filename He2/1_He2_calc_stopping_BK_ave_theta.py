#0_He2_calc_stopping_BK.py: Brandt-Kitagawaモデルを2原子クラスターの阻止能公式に適用し、イオンに対する阻止能を計算する(積分はscipy.integrateを利用)
# stopping averaged over theta

import numpy as np
from numpy import sqrt, sin, cos, pi, log, radians
from scipy.special import j0
from scipy import integrate
import os
import sys
import json
sys.path.append('../')
from common_library import *

#directory path
working_dir =  r'./'
input_dir = os.path.join(working_dir, 'results')
q_ave_path = os.path.join(input_dir, 'He2_average_charge.json')
param_path = r'./param_He2.json'

#experimental condition
# Energy(keV) and velocity(au) of projectile
E = 0
#target name
target = ''
#plasmon energy(atomic unit)
E_p = 0

#parameters for integration
r_close = 0
r_dist = 0
k_min = 0 #not kappa !!
k_max = 0
r = {"00":0}
N = [0]
q = [0]
lamb = [0]
w_max = 0

#beam velocity v (atomic unit)
v = 0

#被積分関数
def integrand(k, w):
    #Brandt-Kitagawaモデルの遮蔽関数（のFourier変換）
    zeta = [0] * len(q)
    
    #calc each zeta
    for i in range(len(q)):
        zeta[i] = (q[i] + (k*lamb[i])**2)/(1 + (k*lamb[i])**2)

    #calc diagonal and cross term
    diagonal_terms = 0
    cross_terms = 0
    #diag
    for i in range(len(q)):
       diagonal_terms += zeta[i]**2 
    #cross
    for i in range(len(q)):
        for j in range(len(q)):
            if i >= j:
                continue
            cross_terms += 2*zeta[i]*zeta[j]*sin(k * r[str_ind(i, j)])/(k * r[str_ind(i, j)]) #averaged over theta
    elf = drude_function(w)
    return w * (diagonal_terms + cross_terms)*elf/k

#calc stopping
def calc_stopping_BK():
    #integration result
    global w_max # w_max must be in [0, 1]
    w_max = 1
    res = integrate.dblquad(integrand, 0, w_max * v**2/2, lambda x: v*(1-(1-2*x/v**2)**(1/2)), lambda x: v*(1+(1-2*x/v**2)**(1/2)))[0]
    #calc stopping
    res *= Z_HELIUM**2
    res *= 2/pi/v**2
    res *= e2/a_0**2 # to eV/A
    return res

def set_parameters(path):
    global E, E_p, target
    global r_close, r_dist 
    global k_min, k_max
    global N
    global q
    global lamb
    global r
    global v

    #import parameters
    with open(param_path, 'r') as f:
        params = json.loads(f.read())
    #set projectile parameters
    E = params["E0"]
    v = sqrt(E/E_HELIUM)
    r = {k:v/a_0 for k, v in params["r"].items()} # to atomic unit

    #set target parameters
    target = params["target"]
    E_p = params["Ep"][target]/E_0 # to atomic unit
    
    #import average charge as np.array
    with open(q_ave_path, 'r') as f:
        dat = json.loads(f.read())
        q = np.array([v for k, v in dat[f"{E}"].items()])

    # N : 束縛電子の数 (Z - q)
    N = Z_HELIUM - q
    # q : イオンの価数を電離度(0 - 1)に変換 -> BKモデルのqに対応
    q /= Z_HELIUM

    #Brandt-Kitagawaモデルの遮蔽定数 lambda (atomic unit)
    lamb = [0] * len(q)
    for i in range(len(q)):
        lamb[i] = 2 * A * (N[i]/Z_HELIUM)**(2/3)/(Z_HELIUM**(1/3)*(1-N[i]/Z_HELIUM/7))
    
    #calc r_close, r_dist(atomic unit)
    r_close = 1/2/v
    r_dist = v/E_p

    #calc k range(atomic unit)
    k_min = 1/r_dist
    k_max = 1/r_close
    return

def main():
    set_parameters(param_path)
    #check parameters
    print('E = {} keV/atom'.format(E))
    print('v = {} au'.format(v))
    print('target: {}'.format(target))
    print('E_p = {} au'.format(E_p))
    print('N : ', N)
    print('q : ', q)
    print('r : ', r)
    print('lambda : ', lamb)
    print('r_dist = {} au'.format(r_dist))
    print('r_close = {} au'.format(r_close))
    print('k_min = {} au'.format(k_min))
    print('k_max = {} au'.format(k_max))
    
    stopping = calc_stopping_BK()
    print(f"S = {stopping} eV/A")
    #ファイルに書き込み
    output_filename = 'E={}keV_atom_He2_linear_{}_ave_theta.txt'.format(E, target)
    with open(os.path.join(input_dir, output_filename), 'w') as f:
        f.write(f"{stopping}")
    print('successfully written to {}'.format(os.path.join(input_dir, output_filename)))
    print('finished!')

if __name__ == '__main__':
    main()

