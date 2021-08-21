# 1_H_calc_stopping_BK_low_velocity.py: Brandt-Kitagawaモデルを阻止能公式に適用し、イオンに対する阻止能を計算する(積分はscipy.integrateを利用)
# Fermi速度よりも低い入射エネルギーに対する阻止能の計算

import numpy as np
from numpy import sqrt, sin, cos, pi, log, radians, Inf
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
q_ave_path = os.path.join(input_dir, 'H_average_charge.json')

param_path = r'./param_H.json'

#experimental condition
# Energy(keV) and velocity(au) of projectile
E = 0
#target name
target = ''
#plasmon energy(atomic unit)
E_p = 0
#Fermi velocity(atomic unit)
v_F = 0

#parameters for integration
N = 0     #the number of bound electrons
q = 0     #average charge
qs = {}   #average charges
lamb = 0  #screening parameters of Brandt-Kitagawa model

#beam velocity v (atomic unit)
v = 0

#被積分関数
def integrand(w, k):
    #Brandt-Kitagawaモデルの遮蔽関数（のFourier変換）
    zeta = 0
    
    #calc each zeta
    zeta = (q + (k*lamb)**2)/(1 + (k*lamb)**2)

    #calc diagonal term
    # diagonal_terms = zeta**2
    diagonal_terms = 1

    elf = elf_low_v(k, w, target)
    return w * diagonal_terms * elf/k

#calc stopping
def calc_stopping_BK():
    #integration result
    v_F = AMINO_PROP[target]["vF"]
    k_F = v_F # k_F = m * v_F/h_bar = v_F (in atomic unit)

    res = integrate.dblquad(integrand, 0, 2*k_F, lambda k: 0, lambda k: k*v)[0]

    #calc stopping(eV/A)
    res *= Z_PROTON**2
    res *= 2/pi/v**2
    res *= e2/a_0**2 # to eV/A
    return res

def read_parmaters(path):
    global qs
    global E_p, target
    
    #import parameters
    with open(param_path, 'r') as f:
        params = json.loads(f.read())

    #set target parameters
    target = params["target"]
    E_p = AMINO_PROP[target]["Ep"]/E_0

    #import average charge as np.array
    with open(q_ave_path, 'r') as f:
        qs = json.loads(f.read())
    
    return

def set_parameters(ene):
    global E, v
    global N
    global q
    global lamb

    #set projectile parameters
    E = float(ene)
    v = sqrt(E/E_PROTON)

    # N : 束縛電子の数 (Z - q)
    q = qs[ene]
    N = Z_PROTON - q
    # q : イオンの価数を電離度(0 - 1)に変換 -> BKモデルのqに対応
    q /= Z_PROTON

    #Brandt-Kitagawaモデルの遮蔽定数 lambda (atomic unit)
    lamb = 2 * A * (N/Z_PROTON)**(2/3)/(Z_PROTON**(1/3)*(1-N/Z_PROTON/7))
    
    return

def main():
    read_parmaters(param_path)
    
    total_output = ""
    for ene in sorted(list(map(int,qs.keys()))):
        set_parameters(str(ene))

        #check parameters
        print('E = {} keV/atom'.format(ene))
        print('v = {} au'.format(v))
        print('target: {}'.format(target))
        print('E_p = {} au'.format(E_p))
        print('N : ', N)
        print('q : ', q)
        print('lambda(au) :', lamb)

        stopping = calc_stopping_BK()
        print(f"S = {stopping} eV/A")
        total_output += f"{ene}\t{str(v)}\t{stopping}\n"
        # #ファイルに書き込み
        # output_filename = 'E={}keV_atom_H_{}.txt'.format(ene, target)
        # with open(os.path.join(input_dir, output_filename), 'w') as f:
        #     f.write(f"{stopping}")
        # print('successfully written to {}'.format(os.path.join(input_dir, output_filename)))

    #まとめデータをファイルに書き込み
    output_filename = 'H_{}_low_velocity.txt'.format(target)
    with open(os.path.join(input_dir, output_filename), 'w') as f:
        f.write(total_output)
    print('successfully written to {}'.format(os.path.join(input_dir, output_filename)))
    print('finished!')

if __name__ == '__main__':
    main()

