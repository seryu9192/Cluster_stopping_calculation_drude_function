#_plot_theta_dependence.py: 

import os
import sys
library_path = r'C:\Users\RyuMurase\OneDrive\programs\Research\Experimental_data_analysis'
if library_path not in sys.path:sys.path(library_path)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#グラフ全体のセットアップパラメータ
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

target = 'Gly'
E = 900

#filepath
working_dir =  r'results'
fig_dir = r'../fig'
input_dir = working_dir
filename = 'E={}keV_atom_C2_linear_{}_wk=w.txt'.format(E, target)
inputfile_path = os.path.join(input_dir, filename)

def main():
    #import source data
    with open(os.path.join(inputfile_path), 'r') as f:
        dat = [v.split() for v in f.read().split('\n') if len(v) != 0]
    xs = [int(v) for v in dat[0]]
    dat = dat[1:]
    
    fig = plt.figure(figsize=(9, 6), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    for i in range(len(dat)):
        ys = [float(v) for v in dat[i]]
        
        #x-axis
        ax.set_xlabel('Orientation angle (degree)', fontsize=22)
        ax.set_xlim(0, 90)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))

        #y-axis
        y_min, y_max, y_major, y_minor = 0, 300, 50, 10
        ax.set_ylabel('', fontsize=22)
        ax.set_ylabel('SP (eV/$\AA$)', fontsize=22)
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(MultipleLocator(y_major))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
        
        #plot
        ax.plot(xs, ys, linewidth=2, color='black')
        
    # ax.text(0.02, 0.9, '{:.1f} MeV/atom C$_2^+$ in {}'.format(E/1000,target), fontsize=24, transform=ax.transAxes)
    ax.text(0.01, 0.9, filename, fontsize=24, transform=ax.transAxes)

    save_fig_on = False
    if save_fig_on:
        os.makedirs(fig_dir, exist_ok=True)
        fig_path = os.path.join(fig_dir, filename[:-4]+'.png')
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')        
        print(f"successfully saved to {fig_path}")

if __name__ == '__main__':
    main()