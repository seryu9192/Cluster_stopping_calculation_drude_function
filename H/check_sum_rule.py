# check_sum_rule: 

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

def read_parmaters(path):
    global qs
    global E_p, target
    
    #import parameters
    with open(param_path, 'r') as f:
        params = json.loads(f.read())

    #set target parameters
    target = params["target"]
    print(target)
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

#被積分関数
def integrand(w):
    elf = elf_low_v(1, w, target)
    return elf


def main():
    read_parmaters(param_path)
    print(integrand(1))
    print(integrand(2))
    print(integrand(3))
    print(integrand(4))

if __name__ == '__main__':
    main()
