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
from matplotlib import colors, ticker, cm
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

Nx = 300            # Elements in x direction (number)
Ny = 100            # Elements in y direction  (number)
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

LHS_data = np.array
LHS_i = np.array
LHS_j = np.array
RHS = np.zeros(Nx * Ny)

# =============================================================================
# Constructing vectors for top elements - data, rows and columns. Includes also 
# the construction of the RHS vector for each position
# =============================================================================

for i in range(0, Nx):
    
    if i % Nx == 0:        
        # Top left corner data - 3 datos en matriz (es esquina)
        # Matrix data 
        LHS_data = np.array([-3 * (dy/dx + dx/dy), dy/dx, dx/dy])        
        # Matrix indexes
        LHS_i = np.zeros(3)
        LHS_j = np.array([0, 1, Nx])
        RHS[i] = -2 * (Lbc[0] * dy/dx + Tbc[0] * dx/dy)
    
    elif i % Nx == Nx - 1:        
        # Top right corner data 
        # Matrix data 
        LHS_data = np.concatenate((LHS_data, [-3 * (dy/dx + dx/dy), dy/dx, 
                                              dx/dy]), axis = 0)
        # Matrix indexes
        LHS_i = np.concatenate((LHS_i, [i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [Nx - 1, Nx - 2, 2 * Nx - 1]), axis = 0)        
        RHS[i] = -2 * (Rbc[0] * dy/dx + Tbc[Nx - 1] * dx/dy)
    
    else:        
        # Top row of elements data - 4 datos en matriz por elemento (borde 
        # superior)
        LHS_data = np.concatenate((LHS_data, [-(2 * dy/dx + 3 * dx/dy), dy/dx, 
                                              dy/dx, dx/dy]), axis = 0)
        LHS_i = np.concatenate((LHS_i, [i, i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [i, i + 1, i - 1, i + Nx]), axis = 0)
        
        RHS[i] = -2 * (dx/dy) * Tbc[i]        
# =============================================================================
# Constructing vectors of elements in the middle of the domain. Includes the 
# vertical boundary elements. Includes also the construction of the RHS vector
# with values for the solution of the system
# =============================================================================

for i in range(Nx, (Ny - 1) * Nx):
    
    if i % Nx == 0:        
        # Left boundary matrix coefficients, rows and columns
        LHS_data = np.concatenate((LHS_data, [-(3 * dy/dx + 2 * dx/dy), dy/dx,
                                                dx/dy, dx/dy]), axis = 0)
        LHS_i = np.concatenate((LHS_i, [i, i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [i, i + 1, i - Nx, i + Nx]), axis = 0)
        
        # Right hand side values for Dirichlet boundary condition
        RHS[i] = - 2 * (dy/dx) * Lbc[int(i/Nx)]
        
    elif i % Nx == Nx - 1:
        # Right bounday matrix coefficients, rows and columns
        LHS_data = np.concatenate((LHS_data, [-(3 * dy/dx + 2 * dx/dy), dx/dy,
                                                dx/dy, dx/dy]), axis = 0)
        LHS_i = np.concatenate((LHS_i, [i, i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [i, i - 1, i - Nx, i + Nx]), axis = 0)
        # Right hand side for Dirichlet boundary condition
        RHS[i] = - 2 * (dy/dx) * Rbc[int(((i + 1)/Nx)-1)]
        
    else:
        # Most inner cells filling with rows and columns
        LHS_data = np.concatenate((LHS_data, [-2 * (dy/dx + dx/dy), dy/dx, 
                                              dy/dx, dx/dy, dx/dy]), axis = 0)
        LHS_i = np.concatenate((LHS_i, [i, i, i, i, i]))
        LHS_j = np.concatenate((LHS_j, [i, i + 1, i - 1, i + Nx, i - Nx]), 
                              axis = 0)
        # This part does not include RHS values since the equation equals zero
        # for the most inner cells
    
# =============================================================================
# Constructing vectors of elements in the bottom of the domain. Includes the 
# horizontal boundary elements. Includes also the RHS vector construction for 
# the solution of the linear system
# =============================================================================

count = 0
for i in range((Ny - 1) * Nx, Ny * Nx):
    
    if i % Nx == 0:        
        # Bottom left corner data - 3 datos en matriz (es esquina)
        # Matrix data 
        LHS_data = np.concatenate((LHS_data, [-3 * (dy/dx + dx/dy), dy/dx, 
                                              dx/dy]), axis = 0)        
        # Matrix indexes
        LHS_i = np.concatenate((LHS_i, [i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [i, i + 1, i - Nx]), axis = 0)
        
        RHS[i] = -2 * (Lbc[Ny - 1] * dy/dx + Bbc[0] * dx/dy)
    
    elif i % Nx == Nx - 1:        
        # Bottom right corner data 
        # Matrix data 
        LHS_data = np.concatenate((LHS_data, [-3 * (dy/dx + dx/dy), dy/dx, 
                                              dx/dy]), axis = 0)
        # Matrix indexes
        LHS_i = np.concatenate((LHS_i, [i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [i, i - 1, i - Nx]), axis = 0) 
        
        RHS[i] = -2 * (Rbc[Ny - 1] * dy/dx + Bbc[Nx - 1] * dx/dy)
    
    else:        
        # Bottom row of elements data - 4 datos en matriz por elemento (borde 
        # inferior)
        LHS_data = np.concatenate((LHS_data, [-(2 * dy/dx + 3 * dx/dy), dy/dx, 
                                              dy/dx, dx/dy]), axis = 0)
        LHS_i = np.concatenate((LHS_i, [i, i, i, i]), axis = 0)
        LHS_j = np.concatenate((LHS_j, [i, i + 1, i - 1, i - Nx]), axis = 0)
        
        RHS[i] = -2 * (dx/dy) * Bbc[count] 
        
    count +=1
# =============================================================================
# Spying the matrix to see its construction.
# =============================================================================

LHS = scsp.coo_matrix((LHS_data, (LHS_i, LHS_j)), shape = ((Nx * Ny), 
                       (Nx * Ny)))
#
## Transforming matrix to CSR format to perform operations quickly
LHS = LHS.tocsr()
#                
#plt.spy(LHS, markersize = 0.5)

# =============================================================================
# Solving linear system for the pressure field (Laplace's equation)
# =============================================================================

P = scsp.linalg.spsolve(LHS, RHS)

# =============================================================================
# Reshaping and plotting solution to check values and different Bc's
# =============================================================================

RTA = P.reshape((Nx, Ny))

X, Y = np.meshgrid(xn[:,1], xn[:,2])

im = plt.matshow(RTA, cmap=plt.cm.hot, aspect='auto')

#plt.figure()
#CS = plt.surf(RTA)
#plt.clabel(CS, inline=1, fontsize=10)
#plt.title('Simplest default with labels')
#plt.show()



