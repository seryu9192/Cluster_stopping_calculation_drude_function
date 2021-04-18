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

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

E = 1800
target = 'Gly'

#filepath
working_dir =  r'results'
fig_dir = r'.\figs'
input_dir = working_dir
filename = f'E={E}keV_atom_C2_linear_Gly_interfere.txt'
inputfile_path = os.path.join(input_dir, filename)

def main():
    #import source data
    with open(os.path.join(inputfile_path), 'r') as f:
        dat = [v.split() for v in f.read().split('\n') if len(v) != 0]
    xs = [int(v) for v in dat[0]]
    dat = dat[1:]
    
    fig = plt.figure(figsize=(9, 6), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    labels = ['$S_{\mathrm{A}}$', '$S_{\mathrm{B}}$', '$S_{\mathrm{A}}+S_{\mathrm{B}}$']
    for i in range(len(dat)):
        ys = [float(v) for v in dat[i]]
        
        #x-axis
        ax.set_xlabel('Orientation angle (degrees)', fontsize=22)
        x_min, x_max, x_major, x_minor = 0, 90, 10, 5
        ax.set_xlim(x_min, x_max)
        ax.xaxis.set_major_locator(MultipleLocator(x_major))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor))

        #y-axis
        y_min, y_max, y_major, y_minor = -20, 200, 50, 10
        ax.set_ylabel('', fontsize=22)
        ax.set_ylabel('SP (eV/$\AA$)', fontsize=22)
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(MultipleLocator(y_major))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
        
        #plot
        ax.plot(xs, ys, linewidth=2, color=colors[i], label=labels[i])
        
    ax.text(0.02, 0.9, f'{E/1000*2:.1f} MeV C$_2^+$ in {target}', fontsize=24, transform=ax.transAxes)
    ax.legend()
    # ax.text(0.01, 0.9, filename, fontsize=16, transform=ax.transAxes)

    save_fig_on = False
    if save_fig_on:
        os.makedirs(fig_dir, exist_ok=True)
        fig_path = os.path.join(fig_dir, filename[:-4]+'.png')
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')        
        print(f"successfully saved to {fig_path}")

if __name__ == '__main__':
    main()