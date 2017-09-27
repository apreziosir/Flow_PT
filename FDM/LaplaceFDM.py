#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solver de Laplace con diferencias finitas - Se hace con matrices eficientes
Incluye el calculo de la red de flujo con los valores de presion (gradiente del 
campo escalar de presiones)
Created on Tue Sep 12 10:11:10 2017
@author: apreziosir
"""

import numpy as np
import scipy.sparse as scsp
import scipy.io 
import matplotlib.pyplot as plt
from Exp_cond import alt_media, inf_vel
from Bound_cond import fill_tbc, fill_bbc, fill_bbc_N, fill_rbc, fill_lbc
from FDM_Auxiliar import positions, nzero, comp  
from Builder import RHS_build, LHS_build, LHS_build_N 
from Velocity_prof import gw_vel

# =============================================================================
# Input variables - flow conditions
# =============================================================================

U = 0.35                    # Mean flow velocity in m/s
H = 0.015                   # Dune height
d = 0.10                    # Mean depth of flow
phi = 0.33                  # Porosity of material
q = -0.00015                     # Inflow or downflow velocity (+ up / - down)
K = 0.1195                  # Hydraulic conductivity
Neum = True                # Neumann condition at the bottom?

# =============================================================================
# Input variables - domain and bed characteristics
# =============================================================================

Lx = 6.40           # Length of the flume (m) (considered for numerical model) 
Ly = -0.25           # Depth of bed (m)
Lambda = 0.15       # Wavelength of bedform (m)
Dif = 1.0           # Diffusion coefficient (just for fun)

# =============================================================================
# Numerical model input parameters
# =============================================================================

Nx = 400              # Nodes in x direction (number)
Ny = 100              # Nodes in y direction  (number)

# Set up mesh - function that calculates everything
dx = np.abs(Lx / (Nx - 1))
dy = np.abs(Ly / (Ny - 1))

# =============================================================================
# Calculate hm value for the problem assigned, inflow Darcy velocity and set 
# a vector for the boundary condition values - outside functions
# =============================================================================

hm = alt_media(U, H, d)

if Neum == True : v = inf_vel(phi, q, K)
Tbc = fill_tbc(Lx, Nx, hm, Lambda)
if Neum == False : Bbc = fill_bbc(Tbc, Nx, Lx, Ly)
else : Bbc = fill_bbc_N(Nx, Dif, q, dy, K)
Lbc = fill_lbc(Ly, Ny, Tbc[0])
Rbc = fill_rbc(Lx, Ly, Ny, Tbc[Nx - 1])

# =============================================================================
# Node positions - mesh generation
# =============================================================================

# Matrix with nodes ID's and positions
xn = positions(Lx, Ly, Nx, Ny)

# =============================================================================
# Calculating the number of nonzero elements in matrix
# =============================================================================

nzero(Nx, Ny, Neum)

# =============================================================================
# Building RHS of system to solve Laplace Equation 
# (Results vector) - This vector is built the same way if there is a Dirichlet
# or Neumann BC since it is just a vector
# =============================================================================

RHS = RHS_build(Tbc, Bbc, Lbc, Rbc)

# =============================================================================
# Building LHS matrix to solve the system
# Coordinate system storage - Later transformed to CSR (by Python script)
# =============================================================================

if Neum == False : LHS = LHS_build(Nx, Ny, dx, dy, Dif)
else : LHS = LHS_build_N(Nx, Ny, dx, dy, Dif)
LHS = LHS.tocsr()
scipy.io.mmwrite('matrix_test', LHS)

# Checking matrix construction
plt.spy(LHS, markersize = 2)
plt.show()

# =============================================================================
# Solving linear system for the pressure field (Laplace's equation)
# =============================================================================

# Gradiente con jugado - no esta funcionando
#P = scsp.linalg.cg(LHS, RHS)
P = scsp.linalg.spsolve(LHS, RHS)

# =============================================================================
# Reshaping and plotting solution to check values and different Bc's
# =============================================================================

RTA = np.reshape(P, (Ny, Nx), order='C')

np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

print('*---------------------------*')
print('El valor máximo de RTA es:')
print(np.amax(RTA))

print('El valor mínimo de RTA es:')
print(np.amin(RTA))

print('El valor de RTA es... ')
print(RTA)
print('*---------------------------*')

# =============================================================================
# Plotting the solution for visual check
# =============================================================================

x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y)
CS4 = plt.contourf(X, Y, RTA)
cbar = plt.colorbar(CS4)
#CS4.set_clim(vmin=-10000, vmax=1000)
#plt.clim(-np.amax(RTA),np.amax(RTA))
plt.gca().set_aspect(9, adjustable='box')
plt.ylim((Ly, 0))
#cbar.Normalize(clip=False)
#plt.clabel(CS4, fmt='%2.1f', colors='w', fontsize=14)
plt.show()

# =============================================================================
# Comparing with test functions values (just for test cases, not real 
# conditions)
# =============================================================================

#err = comp(RTA, Nx, Ny, Lx, Ly)
##CS5 = plt.contourf(X, Y, err)
##cbar = plt.colorbar(CS5)
##plt.gca().set_aspect(9, adjustable='box')
##plt.ylim((Ly, 0))
##plt.show()
#
## Guardando archivos para comparar en csv
#np.savetxt('sol_numerica.csv', RTA, delimiter=' ', newline='\n')
#np.savetxt('sol_analit.csv', err, delimiter=' ', newline='\n')

# =============================================================================
# Calculating the velocity field with Darcy's law q = K grad(h)
# =============================================================================

