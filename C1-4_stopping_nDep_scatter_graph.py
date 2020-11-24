#C1-4_stopping_nDep_scatter_graph.py

import os
import sys
import json
library_path = r'C:\Users\RyuMurase\OneDrive\programs\Research\Experimental_data_analysis'
if library_path not in sys.path:sys.path(library_path)
import numpy as np
import eda_library as lib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

#directory path
root_path = r'./'
inputdir_path = os.path.join(root_path, 'summary_of_results')
inputfile_name = r'E=900keV_w=0.0-1.0.txt'
inputdata_path = os.path.join(inputdir_path, inputfile_name)

#general parameters of gragh
plt.rcParams["font.family"] = "Arial"
plt.rcParams["xtick.direction"] = "out"        
plt.rcParams["ytick.direction"] = "out"
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
plt.rcParams["font.size"] = 20
plt.rcParams["axes.linewidth"] = 1.5

markers = ['o', 's', '^', 'v']
labels = ['single', 'linear-para.', 'linear-perp.','ring-perp.']

def main():
    #import source data as int array
    with open(inputdata_path, 'r') as f:
        dat = np.genfromtxt(f)

    fig = plt.figure(figsize=(7, 5), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
        
    for i in range(4):
        #plot
        if i == 0:#single ion
            xs = dat[:1, 0]
            ys = dat[:1, 1]
        else:     #para/perp/ring
            xs = dat[1:, 0]
            ys = dat[1:, i]
        ax.plot(xs, ys, markers[i], markersize=10, label=labels[i])

    #setup for xaxis
    x_min, x_max, x_major, x_minor = 0.5, 4.5, 1, 1
    ax.set_xlim(x_min, x_max)
    ax.xaxis.set_major_locator(MultipleLocator(x_major))
    ax.xaxis.set_minor_locator(MultipleLocator(x_minor))
    ax.set_xlabel('$n$ of C$_n^+$', fontsize=24)

    #setup for yaxis
    y_min, y_max, y_major, y_minor = 0, 300, 100, 20
    # y_min, y_max, y_major, y_minor = 0, 2, 0.5, 0.5
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(MultipleLocator(y_major))
    ax.yaxis.set_minor_locator(MultipleLocator(y_minor))
    ax.set_ylabel('calculated SP (eV/$\AA$)', fontsize=24)
    # ax.set_ylabel('$S_n/nS_1 (eV/\AA$)', fontsize=24)

    #labels
    ax.legend(loc='upper left', fontsize=16, frameon=False)
    
    #title  
    # ax.text(0.01, 1.05, "calculated stopping power of glycine to C$_n^+$", fontsize=20 ,transform=ax.transAxes)
            
if __name__ == '__main__':
    main()    