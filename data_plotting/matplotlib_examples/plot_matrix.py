#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Philipp Orekhov'
__copyright__ = 'Copyright (C) 2012 Philipp Orekhov'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os,optparse

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-m','--mtx', dest='mtx', help = 'file with matrix', metavar='matrix.dat')
    parser.add_option('-s','--step', dest='n_ticks', help = 'number of ticks',
    metavar='50')
            
    (options,args) = parser.parse_args()
    
    if options.mtx == None:
        print "Error: Marix file is not defined!\n"
        parser.print_help()
        exit(0)
    if options.n_ticks == None:
        n_ticks = 5
    else:
        n_ticks = int(options.n_ticks)

    absmatrix = os.path.abspath(options.mtx)
    
    
    with open(absmatrix, 'r') as enemat:
        header_items = enemat.readline().split() 
        header_list = header_items[:(len(header_items))] # to remove nonsense elements in the end
        enemat.close()
    
    mat = np.loadtxt(absmatrix, skiprows = 1)
    
    header_length = len(header_list)
    mat_max = mat.max()
    mat = mat/mat_max
    plt.imshow(mat, cmap = cm.jet, origin = 'lower', interpolation='nearest', aspect='equal')
    plt.colorbar()
    plt.title('Correlation map')
    plt.xlabel('Res. Indices as in PDB')
    plt.ylabel('Res. Indices as in PDB')
    plt.tick_params(direction = 'out')
    step = header_length/n_ticks
    range = np.arange(header_length, step=step)
    header_short_list = [] # a list with the only ticks to show
    for i in range: header_short_list.append(header_list[i])
    plt.yticks(range, header_short_list)
    plt.xticks(range, header_short_list, rotation=90)
    
    plt.show()
