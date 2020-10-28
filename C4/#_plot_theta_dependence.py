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

#filepath
working_dir =  r'results'
input_dir = working_dir
filename = 'E=900keV_atom_C4_linear_{}.txt'.format(target)
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
        ax.set_ylabel('', fontsize=22)
        ax.set_ylabel('Stopping(eV/$\AA$)', fontsize=22)
        ax.set_ylim(0, 700)
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(50))
        #plot
        ax.plot(xs, ys, linewidth=1)
        
    ax.text(0.02, 0.9, '0.9 MeV/atom C$_4^+$ in {}'.format(target), fontsize=24, transform=ax.transAxes)

    save_fig_on = False
    if save_fig_on:
        outputfig_dir = os.path.join(input_dir, 'fig')
        os.makedirs(outputfig_dir, exist_ok=True)
        fig_path = os.path.join(outputfig_dir, filename[:-4]+'.png')
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')        

if __name__ == '__main__':
    main()