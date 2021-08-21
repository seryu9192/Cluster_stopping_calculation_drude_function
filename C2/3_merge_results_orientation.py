# 3_merge_results_orientation.py
# merge results for each energy to one file

import os
import sys
import json

working_dir =  r'results'
param_path = r'./param_C2.json'

target = ''
Es = []
input_dirname = os.path.join(working_dir, 'averaged')

def set_parameters(path):
    global Es, target

    #import parameters
    with open(param_path, 'r') as f:
        params = json.loads(f.read())
    #set projectile parameters
    Es = list(params["r"].keys())

    #set target parameters
    target = params["target"]

def main():
    set_parameters(param_path)

    # output data as string
    outputdat = ''
    # loop for incident energies
    for E in Es:
        input_filename = f'E={E}keV_atom_C2_linear_{target}_ave_para=0deg-.txt'
        input_path = os.path.join(input_dirname, input_filename)
        with open(input_path, 'r') as f:
            dat = [list(map(float, v.split())) for v in f.read().split('\n') if v != '']
            outputdat +=  '\t'.join([str(2), str(E), '-1', '0', str(dat[1][0])]) + '\n'
            outputdat +=  '\t'.join([str(2), str(E), '-1', '1', str(dat[1][2])]) + '\n'
    
    # write output file
    output_filename = f'C2_linear_{target}_stopping_orientation.txt'
    output_path = os.path.join(input_dirname, output_filename)
    with open(output_path, 'w') as f:
        f.write(outputdat)
    print("finished!")

if __name__ == '__main__':
    main()