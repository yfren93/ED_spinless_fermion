#!/usr/bin/env python
'''
-------------------------------------------------------------
This function block-diagonalize the basis functions by translation symmetry
The input parameters are N, n: construct full basis
N1, N2, and tr: define the translation operations
-------------------------------------------------------------
'''
def trans_symm_basis():
  import scipy as sp
  import numpy as np
  import scipy.special as sps
  import itertools
  from lattice import *
  from globalpara import *
  from ED_basis_fun import *

  "Initialize basis functions and their orders: return dob{}, num_stat1{}"
  N, n = N_site, n_electron
  
  num_basis = int(sps.comb(N,n)) # total number of basis functions
  order_basis = list( itertools.combinations(range(N),n) ) # basis functions presents in order
  
  dob = {} # dictionary of ordered basis function. key: value --> number: state tuple
  num_stat1 = {} # same as above dict with key <--> value
  
  for ii in range(0,num_basis):
    dob[ii] = order_basis[ii] # the value is tuple. Define dob
  num_stat1 = {v:k for k, v in dob.items()} # reverse the key and value
  
  del order_basis
 
  "Initialize the translation operations: return tr[]" 
  N1, N2, n_unit = Height/3, Length, 3

  tr = np.zeros((N1,N2,N_site),dtype = int) # define translation operator with torus geometry
  for ii1 in range(0,N1):
    for ii2 in range(0,N2):
      for jj in range(0,n_unit):
        for jj1 in range(0,N1):
          for jj2 in range(0,N2):
            tr[ii1,ii2,jj+jj1*n_unit+jj2*Height] = np.mod(jj+(jj1+ii1)*n_unit,Height) + np.mod((jj2+ii2),Length)*Height
  
  "Re-organize basis function according to the translation operations: return map_dob{}, dob1{}, renormalize_factor"
  map_dob = {ii:[] for ii in range(0,num_basis)} # new dict that map previous order_basis to new defined basis order
  dob1 = {}   # dict of block basis function
  renormalize_factor = [] # renormalization factor of each block basis.
  
  num_nb = -1 # number of basis function in reduced block
  bas1 = np.zeros((n,),dtype = int) # basis obtained after operation
  while dob != {}:
    bas = dob[next(iter(dob))]
    
    num_nb += 1
    #print '---------', num_nb
    dob1[num_nb] = []
    del_ele = []
    for ii1 in range(0,N1):
      for ii2 in range(0,N2):

        bas1[:] = tr[ii1,ii2,bas]
        t_bas = tuple(bas1)
  
        dob1[num_nb] += [t_bas]
        
        nbas1 = num_stat1[tuple(sorted(bas1))]
        map_dob[nbas1] += [num_nb]  # map old nbas1-th state to num_nb-th for block basis
        map_dob[nbas1] += [ii2+ii1*N2]     # denote the translations along x and y directions
        #map_dob[nbas1] += [jj1]

        del_ele += [nbas1]          # delete the used basis in dob
    for ii in set(del_ele):
      del dob[ii]
  
    renormalize_factor += [N1*N2/len(set(del_ele))] # define the appearance time of each basis
  
  return renormalize_factor, map_dob, dob1, num_stat1

'''
--------------------------------------------------------
               MAIN PROGRAM
--------------------------------------------------------
'''
if __name__ == '__main__':
# import numpy as np
  import scipy as sp
  import matplotlib.pyplot as plt

  from scipy.sparse.linalg import LinearOperator
  from scipy.sparse.linalg import eigsh
  import time
  from sparse_Ham_fun import *

  renormalize_factor, map_dob, dob1 = trans_symm_basis()

  for ii in range(0,len(dob1)):
    print 'dob1', ii, dob1[ii]
  print 'renormalize_factor', renormalize_factor
  print map_dob
  #for ii in map_dob:
  #  print map_dob[ii]
else:
  print ' Import rebasis function \n'
  
