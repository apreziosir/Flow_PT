# -*- coding: utf-8 -*-
"""
Script to calculate Laplace equation and then calculate flow for the particle 
tracking model
Antonio Preziosi-Ribero / Universidad Nacional de Colombia
August 2017
"""

import numpy as np
import scipy.sparse as scsp
import matplotlib.pyplot as plt
from var_funs import alt_media, fill_tbc, fill_bbc, fill_lbc, fill_rbc
from var_funs import positions 

# =============================================================================
# Input variables - flow conditions
# =============================================================================

U = 0.35            # Mean flow velocity in m/s
H = 0.05            # Dune height
d = 0.20            # Mean depth of flow

# =============================================================================
# Input variables - domain and bed characteristics
# =============================================================================

Lx = 6.00           # Length of the flume (m) (considered for numerical model) 
Ly = 0.30           # Depth of bed (m)
Lambda = 0.60       # Wavelength of bedform (m)

# =============================================================================
# Numerical model input parameters
# =============================================================================

Nx = 10            # Elements in x direction (number)
Ny = 10             # Elements in y direction  (number)
bbc = d + Ly

# =============================================================================
# Calculate hm value for the problem assigned and set a vector for the 
# top boundary condition values - outside functions
# =============================================================================

hm = alt_media(U, H, d)
Tbc = fill_tbc(Lx, Nx, hm, Lambda)
Bbc = fill_bbc(Lx, Nx, hm, Lambda, bbc)
Lbc = fill_lbc(Ly, Ny, Tbc[0])
Rbc = fill_rbc(Ly, Ny, Tbc[Nx - 1])

# =============================================================================
# Node positions - mesh generation
# =============================================================================

# Set up mesh - function that calculates everything
dx = Lx / Nx
dy = Ly / Ny

# Matrix with nodes ID's and positions
xn = positions(Lx, Ly, Nx, Ny)

# =============================================================================
# Calculation of nonzero entries in matrix
# =============================================================================

# Just checking
nzero = 12 + (Nx - 2) * 16 + (Nx - 2) * (Ny - 2) * 5
print('Total entries in matrix (size of matrix):')
print(Nx** 2 * Ny ** 2)
print('Nonzero entries:')
print(nzero)
print('Ratio of nonzero/total: ')
print(nzero / (Nx * Ny) ** 2)

# =============================================================================
# Creating vectors for storing matrix in COORD system
# THis parto of the script has been constructed due to the size of the problem
# and the neccesity of faster computation
# =============================================================================        

# =============================================================================
# Constructing vectors for top elements
# =============================================================================

# Top left corner data - 3 datos en matriz (es esquina)
# Matrix data 
TLcorn_data = np.array([-3 * (dy/dx + dx/dy), dy/dx, dx/dy])

# Matrix indexes
TLcorn_i = np.zeros(3)
TLcorn_j = np.array([0, 1, Nx])

# Top row of elements data - 4 datos en matriz por elemento (borde superior)
Trow_data = np.tile([-(2 * dy/dx + 3 * dx/dy), dy/dx, dy/dx, dx/dy], Nx - 2)
Trow_i = np.zeros((Nx - 2) * 4)
Trow_j = np.zeros((Nx - 2) * 4)

# Top row indexes for row and column  
for i in range(1, Nx - 1):    
    Trow_i[(i - 1) * 4:i * 4] = np.tile(i, 4)
    Trow_j[(i - 1) * 4:i * 4] = np.array([i, i + 1, i - 1, i + Nx])

# Top right corner data 
# Matrix data 
TRcorn_data = np.array([-3 * (dy/dx + dx/dy), dy/dx, dx/dy])
# Matrix indexes
TRcorn_i = np.array([Nx - 1, Nx - 1, Nx - 1])
TRcorn_j = np.array([Nx - 1, Nx - 2, 2 * Nx - 1])

Tel_data = np.concatenate((TLcorn_data, Trow_data, TRcorn_data), axis = 0) 
Tel_i = np.concatenate((TLcorn_i, Trow_i, TRcorn_i), axis = 0)
Tel_j = np.concatenate((TLcorn_j, Trow_j, TRcorn_j), axis = 0)

del(TLcorn_data, TLcorn_i, TLcorn_j, Trow_data, Trow_i, Trow_j, TRcorn_data,
    TRcorn_i, TRcorn_j)

# =============================================================================
# Constructing vectors of elements in the middle of the domain. Includes the 
# vertical boundary elements
# =============================================================================

Iel_data = np.zeros((Ny - 2) * (8 + (Nx - 2) * 5)) 
Iel_i = np.zeros((Ny - 2) * (8 + (Nx - 2) * 5))
Iel_j = np.zeros((Ny - 2) * (8 + (Nx - 2) * 5))



# =============================================================================
# Spying the matrix to see its construction
# =============================================================================

LHS = scsp.coo_matrix((Tel_data, (Tel_i, Tel_j)), shape = ((Nx * Ny), 
                       (Nx * Ny)))

plt.spy(LHS, markersize = 2)