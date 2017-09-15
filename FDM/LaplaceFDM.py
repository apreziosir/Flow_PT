#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solver de Laplace con diferencias finitas - Se hace con matrices eficientes
Created on Tue Sep 12 10:11:10 2017
@author: apreziosir
"""

import numpy as np
import scipy.sparse as scsp
import matplotlib.pyplot as plt
from FDM_Auxiliar import alt_media, fill_tbc, fill_bbc, fill_rbc, fill_lbc
from FDM_Auxiliar import positions, nzero, RHS_build, LHS_build

# =============================================================================
# Input variables - flow conditions
# =============================================================================

U = 0.35            # Mean flow velocity in m/s
H = 0.05            # Dune height
d = 0.20            # Mean depth of flow

# =============================================================================
# Input variables - domain and bed characteristics
# =============================================================================

Lx = 1.00           # Length of the flume (m) (considered for numerical model) 
Ly = 1.00           # Depth of bed (m)
Lambda = 0.60        # Wavelength of bedform (m)
Dif = 1.0           # Diffusion coefficient (just for fun)

# =============================================================================
# Numerical model input parameters
# =============================================================================

Nx = 8           # Nodes in x direction (number)
Ny = 8           # Nodes in y direction  (number)

# =============================================================================
# Calculate hm value for the problem assigned and set a vector for the 
# boundary condition values - outside functions
# =============================================================================

hm = alt_media(U, H, d)
Tbc = fill_tbc(Lx, Nx, hm, Lambda)
Bbc = fill_bbc(Tbc, Nx, Lx, Ly)
Lbc = fill_lbc(Ly, Ny, Tbc[0])
Rbc = fill_rbc(Lx, Ly, Ny, Tbc[Nx - 1])

# =============================================================================
# Node positions - mesh generation
# =============================================================================

# Set up mesh - function that calculates everything
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)

# Matrix with nodes ID's and positions
xn = positions(Lx, Ly, Nx, Ny)

# =============================================================================
# Calculating the number of nonzero elements in matrix
# =============================================================================

nzero(Nx, Ny)

# =============================================================================
# Building RHS of system to solve Laplace Equation 
# (Results vector)
# =============================================================================

RHS = RHS_build(Tbc, Bbc, Lbc, Rbc)

# =============================================================================
# Building LHS matrix to solve the system
# Coordinate system storage - Later transformed to CSR (by Python script)
# =============================================================================

LHS = LHS_build(Nx, Ny, dx, dy, Dif)
LHS = LHS.tocsr()

# Checking matrix construction
# plt.spy(LHS, markersize = 4)

# =============================================================================
# Solving linear system for the pressure field (Laplace's equation)
# =============================================================================

P = scsp.linalg.spsolve(LHS, RHS)

# =============================================================================
# Reshaping and plotting solution to check values and different Bc's
# =============================================================================

RTA = np.reshape(P, (Ny, Nx), order='C')

print('El valor máximo de RTA es:')
print(np.amax(RTA))

print('El valor mínimo de RTA es:')
print(np.amin(RTA))

# ==============================================================================
# Plotting the solution for visual check
# ==============================================================================

x = np.arange(Nx)
y = np.arange(Ny)
X, Y = np.meshgrid(x, y)
CS4 = plt.contourf(X, Y, RTA)
cbar = plt.colorbar(CS4)
#plt.clabel(CS4, fmt='%2.1f', colors='w', fontsize=14)
plt.show()