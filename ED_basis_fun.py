#!/bin/usr/env python
from globalpara import *

'''
--------------------------------------------------------
 define permutation time /* from Fermion exchange */ for two states in different order
--------------------------------------------------------
'''
def permute_time3(a1,a12) : # a1 tuple, a12 list 
  #import numpy as np
  n = len(a1)

  da1 = {a1[ii]:ii for ii in range(0,n)} # define a map from site number to order
  ds1 = range(0,n)
  ds2 = range(0,n)
  for ii in range(0,n):    # redefine the order 
    ds2[ii] = da1[a12[ii]]

  permut = 0  # calculate permute time
  for ii in range(0,n-1) :  # move all right site to middle
    if ds2[ii] > ds1[ii] :
      permut += ds2[ii]-ds1[ii]
      for jj in range(ii+1,n):
        if (ds2[jj] < ds2[ii] and ds2[jj] >= ds1[ii]):
          ds2[jj] += 1      # when one site is moved to middle, the middle ones increase by one

      ds2[ii] = ds1[ii]

  return permut


'''
*******************************************************************************
Given a row of numbers "stat" and total sites "N",
Return the order of this state in the basis functions
Suitable for spinless particle in any lattice
There are n numbers in stat, range from 0 to N-1
*******************************************************************************
'''
def num_stat(stat,N):
  import scipy as sp

  n = len(stat) # number of particles

  # Initialize nstate, start from 1
  if (n > 1): # more than one particle
    nstat = stat[n-1] - stat[n-2]
  else:       # single particle case
    nstat = stat[n-1]+1

  # Define order
  for ii in range(1,n): # summation over each site
    if ii == 1:
      lb = 1
    else:
      lb = stat[ii-2]+1
      lb += 1 # number in stat is start from 0, while that in notes from 1

    ub = stat[ii-1]+1 # number in stat is start from 0, while that in notes from 1
    for jj in range(lb,ub):
      nstat += sp.special.factorial(N-jj)/(sp.special.factorial(n-ii)*sp.special.factorial(N+ii-n-jj))

  # Return the order, start from 0
  return int(nstat-1)


'''
*******************************************************************************
Based on the basis function | N1, N2, ..., Nn > with N1 < N2 < ... < Nn
of a spinless Fermionic system with N sites and n particles, define the Hamiltonian
matrix elements and diagonalize
*******************************************************************************
'''
def Ham_Operator(which = 'sp', Initialize_stat = 'Normal', vNN = 2.0, vNNN = 0.0, vN3 = 0.0) :
  # Initialize_stat: 'Normal': calculate all terms 
  #                  'Initialize': only hop, 'N_Initialize': only diagonal term 

# import numpy as np
  import scipy as sp
  import itertools
  import time
  # from scipy.special import factorial

  global Height, Length, nNNcell, pbcx, pbcy, N_site, n_electron, distNN, distNNN, distN3, Lattice_type


  ## ---------------------------------  Lattice Dependent ----------------------------------- ##
  "Initialize positions and Hamiltonian parameters"
  if Lattice_type == 'honeycomb' :
    distNc = pos_honeycomb() # Distance between sites in neighbor cells 
    tot_coordinate = 3       # number of coordinate sites with << hopping >>
  elif Lattice_type == 'kagome' :
    distNc,pos_x,pos_y = pos_kagome()
    tot_coordinate = 4


  ## ---------------------------------  Lattice InDependent ----------------------------------- ##
  N, n = N_site, n_electron # system size and electron number
  dist = np.zeros((N, N))
  # modification PBCs on atomic distances 
  if ((pbcx==1) and (pbcy==1)) :
    for ii0 in range(0,N):
      for jj0 in range(0,N):
        dist[ii0,jj0] = min(distNc[ii0,jj0,:])

  ts = time.time()
  "Initialize basis functions and their orders"
  num_basis = int(sp.special.comb(N,n)) # total number of basis functions
  order_basis = list( itertools.combinations(range(N),n) ) # basis functions presents in order
  te = time.time()
  print 'Bas num = ', num_basis
  print 'Bas time = ', te - ts

  # Interaction & hopping 
  # vNN, vNNN = 0.0, 0.3 # interactions
  tNN = -1.0          # hopping

  # Define on-site energy that depends on the boundary condition
  onsite = np.zeros((N,))
  if (Lattice_type == 'honeycomb') :
    uA1 = -3.0
    uA2, uB = 2.0, 0.
    site_A1 = [1, 10, 13, 22]
    for ii in range(0,Height):
      for jj in range(0,Length):
        if (np.mod((ii+jj),2) == 0):
          onsite[ii+jj*Height] = uB
        else:
          if ((ii+jj*Height) in site_A1):
            onsite[ii+jj*Height] = uA1
          else:
            onsite[ii+jj*Height] = uA2

  #for ii in range(0,N):
  #  print 'site ',ii+1, ', E =', onsite[ii]

  if (which == 'full'):
    Ham = np.zeros((num_basis,num_basis))

  '''
  Define sparsed Hamiltonian in triple format: row, column and data
  '''
  row_s = np.zeros(( num_basis*(tot_coordinate*n*2+1), ),dtype=np.int)
  col_s = np.zeros(( num_basis*(tot_coordinate*n*2+1), ),dtype=np.int)
  data_s = np.zeros(( num_basis*(tot_coordinate*n*2+1), ))
  if (Initialize_stat == 'Initialize') :
    num_elements = int(num_basis-1) # number of nonzero matrix elements
  else :
    num_elements = int(-1) # number of nonzero matrix elements

  for ii in range(0,num_basis):
    bas = list(order_basis[ii])  # call basis function
    n_bas = np.delete(np.arange(0,N),bas,None).tolist()

    if (np.mod(ii,10000) == 0):
      print '%-8.5e'%(float(ii)/num_basis)

    if (Initialize_stat != 'Initialize') :
      "# ---------------------- Diagonal terms ---------------------- #"
      # numbers of NN and NNN pairs of particles in each state 
      numNN, numNNN, numN3 = 0, 0, 0
      for kk in range(0,n):
        site1 = bas[kk]
        for kk1 in range(kk+1,n):
          site2 = bas[kk1]
          for kk2 in range(0,1):
            dist12 = dist[site1,site2]
            if (abs(dist12-distNN) < 0.01) :
              numNN += 1
            elif (abs(dist12-distNNN) < 0.01) :
              numNNN += 1
            elif (abs(dist12-distN3) < 0.01) :
              numN3 += 1
#x          for kk2 in range(0,nNNcell):
#x            dist12 = distNc[site1,site2,kk2]
#x            if (abs(dist12-distNN) < 0.01) :
#x              numNN += 1
#x            elif (abs(dist12-distNNN) < 0.01) :
#x              numNNN += 1
#x            elif (abs(dist12-distN3) < 0.01) :
#x              numN3 += 1

      num_elements += 1
      row_s[num_elements] = ii
      col_s[num_elements] = ii
      # on-site energy & Interaction energy
      data_s[num_elements] += sum(onsite[bas]) + numNN*vNN + numNNN*vNNN + numN3*vN3

      if (which == 'full'):
        Ham[ii,ii] += data_s[num_elements]
      #print 'bas[',ii,'] =', bas,'Eint =', Ham[ii,ii]

    if (Initialize_stat != 'N_Initialize') :
      "# ---------------------- Off-diagonal terms ------------------ #"
      # hopping energy
      for kk in range(0,n): # kk is index
        site1 = bas[kk]     # this the site considered
        # find neighbors
        for kk1 in n_bas:
          if (kk1 > site1): # define upper triangle elements
            if ( abs(dist[site1,kk1]-distNN) < 0.01 ):
              bas1 = np.array(bas)
              bas1[kk] = kk1 # change bas1 to new state by substitute

              # define permutation time /* from Fermion exchange */
              sign_permute = np.sign(bas1-sorted(bas1))
              permute_time = np.mod( np.sum(sign_permute), 2.0) - 1.0
              permute_time = permute_time * np.sign(np.sum(np.abs(sign_permute)))

              # upper triangle
              num_elements += 1
              row_s[num_elements] = ii
              col_s[num_elements] = num_stat(sorted(bas1),N)
              data_s[num_elements] = (-1.)**permute_time * tNN

              if (which == 'full'):
                Ham[row_s[num_elements], col_s[num_elements]] += data_s[num_elements]

              # lower triangle
              num_elements += 1
              row_s[num_elements] = col_s[num_elements-1]
              col_s[num_elements] = row_s[num_elements-1]
              data_s[num_elements] = np.conjugate(data_s[num_elements-1])

              if (which == 'full'):
                Ham[row_s[num_elements], col_s[num_elements]] += data_s[num_elements]

  print 'Finish calculation of Ham_sparse'

  if (which == 'full'):
    return Ham
  else:
    return num_basis, row_s, col_s, data_s
#  # Define the sparse Hamiltonian in triple format
#  if ((which == 'sp') and (Initialize_stat == 'Normal')):
#    Ham_triple = sp.sparse.coo_matrix((data_s,(row_s,col_s)), shape = (num_basis,num_basis))
#    return num_elements, Ham_triple
#  else: 
#    return num_elements, row_s, col_s, data_s


