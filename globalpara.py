#!/bin/usr/env python

'''
*******************************************************************************
Global parameters
*******************************************************************************
'''
import numpy as np
import scipy as sp

# atomic chains in y and x directions
# Honeycomb lattice: y is zigzag, while x is armchair
# Kagome lattice: Height = 3*ny where ny is number of 
#                 unit cell along [1/2, sqrt(3)/2] direction 
#                 while nx is that along [1, 0] direction
Height, Length = 3*2, 2 #3, 3

nNNcell = 9 # number of neighbor unit cells
pbcx, pbcy = 1, 1 # boundary condition: 1: true. 0: false
N_site = Height * Length
n_electron = N_site/3

pi = 3.1415926535897932

Lattice_type = 'kagome'
if (Lattice_type == 'honeycomb' or Lattice_type == 'kagome') :
  distNN, distNNN, distN3 = 1.0, np.sqrt(3.0), 2 # distance between NN/NNN sites


'''
*******************************************************************************
Here is the main program
*******************************************************************************
'''
if __name__ == '__main__':
  print('Run global parameters') 
else:
  print('\n Import global parameters \n ')
