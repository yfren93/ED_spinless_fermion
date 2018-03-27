#!/usr/bin/env python

from globalpara import *

'''
*******************************************************************************
Define the position of each site in a cell depending on the lattice type and
boundary conditions                         
The position of sites in nearest 8 cells are also defined for pbc case
*******************************************************************************
'''
def pos_honeycomb(shapeCell = 'square'):

# import numpy as np
  import scipy as sp
  import itertools
  import time

  global Height, Length, nNNcell, pbcx, pbcy, N_site, n_electron, distNN, distNNN

  '''
  Intialize the lattice and boundary conditions
  '''
  # define positions of each site
  N = N_site # Length*Height #Length*Height/2 # system size and electron number
  pos_x, pos_y = np.zeros((N,1)), np.zeros((N,1))

  # 'diamond' or 'square' shaped unit cell
  if (shapeCell == 'diamond'):
    for ii in range(0,Height):
      for jj in range(0,Length):
        numij = jj * Height + ii
        pos_x[numij,0] = jj*distNN*3.0/2.0 + np.mod(ii,2)*distNN/2.0
        pos_y[numij,0] = jj*distNNN/2.0    + ii*distNNN/2.0
    # define unit vectors for pbc
    a1 = np.array([0.0, distNNN])*Height/2; a2 = np.array([distNNN*np.sqrt(3.0)/2.0, distNNN/2.0])*Length
    
  elif (shapeCell == 'square'):
    print 'square unit cell'
    for ii in range(0,Height):
      for jj in range(0,Length):
        numij = jj * Height + ii
        pos_x[numij,0] = jj*distNN*3.0/2.0 - np.mod(ii+jj,2)*distNN/2.0
        pos_y[numij,0] = ii*distNNN/2.0
    # define unit vectors for pbc
    a1 = np.array([0.0, distNNN])*Height/2; a2 = np.array([distNN*3.0/2.0, 0])*Length

  "define distance between every two sites for both intra-cell and inter-cell cases"
  dist = np.sqrt((pos_x-pos_x.T)**2 + (pos_y-pos_y.T)**2) # Inside of unit cell
  distNc = np.zeros((N,N,nNNcell)) # Distance between sites in neighbor cells 
  for kk1 in range(0,3):
    for kk2 in range(0,3):
      pos_xb1 = pos_x + (kk1-1)*a1[0] + (kk2-1)*a2[0]; pos_yb1 = pos_y + (kk1-1)*a1[1] + (kk2-1)*a2[1]
      distNc[:,:,kk1*3+kk2] = np.sqrt( (pos_x-pos_xb1.T)**2 + (pos_y-pos_yb1.T)**2 )

#  # modification PBCs on atomic distances 
#  if ((pbcx==1) and (pbcy==1)) :    
#    for ii0 in range(0,N):
#      for jj0 in range(0,N):
#        dist[ii0,jj0] = min(distNc[ii0,jj0,:])

  return distNc


'''
*******************************************************************************
Define the position of each site in a cell depending on the lattice type and
boundary conditions                         
The position of sites in nearest 8 cells are also defined for pbc case
*******************************************************************************
'''
def pos_kagome(N_site0=N_site, Height0 = Height, Length0=Length):

# import numpy as np
  import scipy as sp
  import itertools
  import time

  global Height, Length, nNNcell, pbcx, pbcy, N_site, n_electron, distNN, distNNN, distN3

  '''
  Intialize the lattice and boundary conditions
  '''
  # define positions of each site
  N = N_site0 # Length*Height #Length*Height/2 # system size and electron number
  pos_x, pos_y = np.zeros((N,1)), np.zeros((N,1))

  pos_x[0:3,0] = np.array([0,distNN/2,distNN]) # define the first three sites
  pos_y[0:3,0] = np.array([0,distNNN/2,0])

  for ii in range(1,Height0/3): # define the first column
    pos_x[ii*3:ii*3+3,0] = pos_x[0:3,0]+ii*distN3/2.0;
    pos_y[ii*3:ii*3+3,0] = pos_y[0:3,0]+ii*distN3*np.sqrt(3.0)/2.0

  for ii in range(1,Length0): # define the others
    pos_x[ii*Height0:(ii+1)*Height0,0] = pos_x[0:Height0,0] + distN3*ii
    pos_y[ii*Height0:(ii+1)*Height0,0] = pos_y[0:Height0,0]

  # define unit vectors for pbc
  ay = np.array([distNN, distNNN])*Height0/3; ax = np.array([distN3, 0])*Length0

  "define distance between every two sites for both intra-cell and inter-cell cases"
  dist = np.sqrt((pos_x-pos_x.T)**2 + (pos_y-pos_y.T)**2) # Inside of unit cell
  distNc = np.zeros((N,N,nNNcell)) # Distance between sites in neighbor cells 
  for kky in range(0,3):    # for kky, kkx, ay, ax, xy indicate the two unit vectors
    for kkx in range(0,3):  # for posx, posy, xy indicate the x-y axis
      pos_xb1 = pos_x + (kky-1)*ay[0] + (kkx-1)*ax[0]; pos_yb1 = pos_y + (kky-1)*ay[1] + (kkx-1)*ax[1]
      distNc[:,:,kky*3+kkx] = np.sqrt( (pos_x-pos_xb1.T)**2 + (pos_y-pos_yb1.T)**2 )

  return distNc, pos_x, pos_y


'''
*******************************************************************************
Define the Hamiltonian of single particle kagome lattice for give momentum  
Momentum is defined by the periodic boundary condition of torus used
The eigenvalues and eigenstates of each state is given in UnitaryTransform
*******************************************************************************
'''
def GetKagomeUnitaryTransform(phi0 = 0.0):
  import numpy as np

  global Height, Length, nNNcell, pbcx, pbcy, N_site, n_electron, distNN, distNNN, distN3, Lattice_type
  # define constants
  pi = 3.1415926535897932

  N=N_site

  N0, Height0, Length0 = 3, 3, 1
  # define unit vectors
  ay = np.array([distNN, distNNN])*Height0/3; ax = np.array([distN3, 0])*Length0
  onsite=np.zeros((N0,))

  # hopping
  tNN = 1.0 #*np.sqrt(2);
  t2 = -0.1
  phi = phi0*pi;
  distNc, pos_x, pos_y = pos_kagome(N_site0=3, Height0=3, Length0=1)
  hopE = np.zeros([N0,N0,9])*0j

  for kky in range(0,3):
    for kkx in range(0,3):
      for ii1 in range(0,N0):
        for ii2 in range(0,N0):
          if (abs(distNc[ii1,ii2,kky*3+kkx]-distNN) < 0.01) :
            if np.mod(ii1-ii2,3) == 2:
              hopE[ii1,ii2,kky*3+kkx] += tNN*np.exp(phi*1j)
            else:
              hopE[ii1,ii2,kky*3+kkx] += tNN*np.exp(-phi*1j)
          elif (abs(distNc[ii1,ii2,kky*3+kkx]-distNNN) < 0.01) :
            hopE[ii1,ii2,kky*3+kkx] += t2

  # define the momentum 
  by = 2.0*pi/(distNNN)*np.array([0, 1.0])
  bx = 2.0*pi/(distNNN)*np.array([np.sqrt(3.0)/2.0, -1.0/2.0])

  UnitaryTransform = {}
  for iix in range(0, Length) :
    for iiy in range(0, Height/3) :
      UnitaryTransform[iiy*Length+iix]={} # count momentum from left to right & from bottom to top
      kx = iix*bx[0]/Length + iiy*by[0]/(Height/3)
      ky = iix*bx[1]/Length + iiy*by[1]/(Height/3)

      Ham = np.diag(onsite+0j)
      for kky in range(0,3) :
        for kkx in range(0,3) :
          posdx = (kky-1)*ay[0] + (kkx-1)*ax[0]; posdy = (kky-1)*ay[1] + (kkx-1)*ax[1]
          Ham += hopE[:,:,kky*3+kkx]*np.exp(1j*kx*posdx+1j*ky*posdy)
      Ham = (Ham + np.conjugate(Ham.T))/2
#      if (iiy*Length+iix) == 1 :
#        print 'Ham 1 \n', Ham
      eige,eigf = np.linalg.eigh(Ham)

      UnitaryTransform[iiy*Length+iix]['eige'] = eige
      UnitaryTransform[iiy*Length+iix]['eigf'] = eigf
#      print 'xx', iix, iiy, iiy*Length+iix 
  return UnitaryTransform

'''
*************************************************************
Define the hopping between unit cells of kagome lattice
*************************************************************
'''
def unit_kagome():
  import numpy as np
  
  nunit = 3
  hopu[]

  for ii in range(0,3):

    for jj in range(0,3):
    
      
