#0_C4_ring_calc_stopping_BK.py: Brandt-Kitagawaモデルを3原子クラスターの阻止能公式に適用し、イオンに対する阻止能を計算する(積分はscipy.integrateを利用)
# ring structure (n=4)

import numpy as np
from numpy import sqrt, sin, cos, pi, log, radians
from scipy.special import j0
from scipy import integrate
import os
import sys
import json

#directory path
working_dir =  r'./'
input_dir = os.path.join(working_dir, 'results')
q_ave_path = os.path.join(input_dir, '0_C4_ring_average_charge.txt')
param_path = r'./param_C4_ring.json'

#CONSTANTS
a_0 = 0.529 #Bohr radius(A)
E_Ryd = 13.6 #Rydberg energy(eV)
e2 = 14.4 #(eV A)
E_CARBON = 300 #(keV at 1 au velocity)
Z_CARBON = 6 #atomic number of carbon
A = 0.240 # BKモデルの遮蔽定数の計算のための定数

#experimental condition
# Energy(keV) and velocity(au) of projectile
E_0 = 900

#target name
target = ''

#plasmon energy(eV)
E_p = 20

#parameters for integration
r_close = 0
r_dist = 0
k_max = 0
N = [0, 0]
q = [0, 0]
lamb = [0]
r = {"00":0}
t = 0
v = sqrt(E_0/E_CARBON)


def str_ind(i, j):
    res = str(i) + str(j)
    if res > res[::-1]:res = res[::-1]
    return res

#被積分関数
def func(k):
    #Brandt-Kitagawaモデルの遮蔽関数（のFourier変換）
    zeta = [0] * len(q)
    
    #calc each zeta
    for i in range(len(q)):
        zeta[i] = (q[i] + ((1/r_dist)**2 + k**2)*(lamb[i]**2))/(1 + ((1/r_dist)**2 + k**2)*(lamb[i]**2))

    #calc diagonal and cross_terms term
    diagonal_terms = 0
    cross_terms = 0
    for i in range(len(q)):
       diagonal_terms += zeta[i]**2
    for i in range(len(q)):
        for j in range(len(q)):
            if i >= j:
                continue
            cross_terms += 2* zeta[i]*zeta[j]*j0(k * r[str_ind(i, j)])
    return k * (diagonal_terms + cross_terms)/(k**2 + 1/r_dist**2)

#calc stopping
def calc_stopping_BK():
    #integration result
    res, err = integrate.quad(func, 0, k_max)
    #calc stopping
    res *= Z_CARBON**2
    res *= e2
    res /= r_dist**2
    return res

def set_parameters(path):
    global E_0, E_p, target
    global r_close
    global r_dist 
    global k_max
    global N
    global q
    global lamb
    global r
    global t
    global v

    #import parameters
    with open(param_path, 'r') as f:
        params = json.loads(f.read())
    #set projectile parameters
    E_0 = params["E0"]
    v = sqrt(E_0/E_CARBON)
    r = params["r"]
    #set target parameters
    target = params["target"]
    E_p = params["Ep"][target]
    
    #import average charge as np.array
    with open(q_ave_path, 'r') as f:
        q = np.genfromtxt(f)

    # N : 束縛電子の数 (Z - q)
    N = Z_CARBON - q
    # q : イオンの価数を電離度(0 - 1)に変換 -> BKモデルのqに対応
    q /= Z_CARBON

    #Brandt-Kitagawaモデルの遮蔽定数 lambda
    lamb = [0] * len(q)
    for i in range(len(q)):
        lamb[i] = 2 * A * (N[i]/Z_CARBON)**(2/3)/(Z_CARBON**(1/3)*(1-N[i]/Z_CARBON/7)) * a_0
    
    #calc r_close, r_dist
    r_close = a_0/2/v
    r_dist = 2*v*E_Ryd/E_p*a_0
    k_max = sqrt((1/r_close)**2 - (1/r_dist)**2)

    return

def main():
    set_parameters(param_path)
    #check parameters
    print('target: {}'.format(target))
    print('E0 = {} keV/atom'.format(E_0))
    print('v0 = {} au'.format(v))
    print('N : ', N)
    print('q : ', q)
    print('r : ', r)
    print('lambda : ', lamb)
    print('r_dist = {} A'.format(r_dist))
    print('r_close = {} A'.format(r_close))
    print('k_max = {} A'.format(k_max))
    print('E_p = {} eV'.format(E_p))
    
    #calc stopping
    stopping = calc_stopping_BK()
            
    #ファイルに書き込み
    output_filename = 'E={}keV_atom_C4_ring_{}.txt'.format(E_0, target)
    with open(os.path.join(input_dir, output_filename), 'w') as f:
        f.write(str(stopping))
    print('successfully written to {}'.format(os.path.join(input_dir, output_filename)))
    print('finished!')

if __name__ == '__main__':
    main()

