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

Nx = 700            # Elements in x direction (number)
Ny = 100             # Elements in y direction  (number)
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
# Initialize matrices with zeros (LHS and RHS)
# =============================================================================

LHS = scsp.csr_matrix((Nx * Ny, Nx * Ny))
RHS = np.zeros(Nx * Ny)

# =============================================================================
# Compute LHS and RHS elements 
# =============================================================================

# Inner cells with sides
count = 1
for i in range(Nx, (Ny - 1) * Nx):
    
    if i % Nx == 0:
        # Left Hand Side Matrix Coefficients (LEFT BOUNDARY)
        LHS[i, i] = -(3 * dy / dx + 2 * dx / dy)
        LHS[i, i + 1] = dy / dx
        LHS[i, i - Nx] = dx / dy
        LHS[i, i + Nx] = dx / dy     
        # Right Hand Side Vector Values
        RHS[i] = -2 * Lbc[count] * dy / dx # Meter LBC
        
    elif i % Nx == (Nx - 1):
        # Left Hand Side Matrix Coefficients (RIGHT BOUNDARY)
        LHS[i, i] = -(3 * dy / dx + 2 * dx / dy)
        LHS[i, i - 1] = dy / dx
        LHS[i, i - Nx] = dx / dy
        LHS[i, i + Nx] = dx / dy     
        # Right Hand Side Vector Values
        RHS[i] = -2 * Rbc[count] * dy / dx # Meter RBC
        count += 1
        
    else:
        # Left Hand Side Matrix Coeffcients
        LHS[i, i] = -2 * (dy / dx + dx / dy)
        LHS[i, i - 1] = dy / dx
        LHS[i, i + 1] = dy / dx
        LHS[i, i - Nx] = dx / dy
        LHS[i, i + Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = 0.
        
# Upper part of the domain
count = 1
for i in range(0, Nx):
    
    if i % Nx == 0:
        # Left Hand Side Matrix Coeffcients (TOP LEFT CORNER OF THE DOMAIN)
        LHS[i, i] = -3 * (dy / dx + dx / dy)
        LHS[i, i + 1] = dy / dx
        LHS[i, i + Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = -2 * (Lbc[0] * dy / dx + Tbc[0] * dx / dy) 
        
    elif i % Nx == (Nx - 1):
        # Left Hand Side Matrix Coeffcients (TOP RIGHT CORNER OF THE DOMAIN)
        LHS[i, i] = -3 * (dy / dx + dx / dy)
        LHS[i, i - 1] = dy / dx
        LHS[i, i + Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = -2 * (Rbc[0] * dy / dx + Tbc[Nx - 1] * dx / dy) 
    
    else: 
        # Left Hand Side Matrix Coeffcients
        LHS[i, i] = -(-2 * dy / dx + 3 * dx / dy)
        LHS[i, i - 1] = dy / dx
        LHS[i, i + 1] = dy / dx
        LHS[i, i + Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = -2 * Tbc[i] * dx / dy
                
# Lower part of the domain 
count = 1
for i in range((Ny - 1) * Nx, (Ny * Nx) - 1):

    if i % Nx == 0:
        # Left Hand Side Matrix Coeffcients (BOTTOM LEFT CORNER OF THE DOMAIN)
        LHS[i, i] = -3 * (dy / dx + dx / dy)
        LHS[i, i + 1] = dy / dx
        LHS[i, i - Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = -2 * (Lbc[Ny - 1] * dy / dx + Bbc[0] * dx / dy) 
        
    elif i % Nx == (Nx - 1):
        # Left Hand Side Matrix Coeffcients (BOTTOM RIGHT CORNER OF THE DOMAIN)
        LHS[i, i] = -3 * (dy / dx + dx / dy)
        LHS[i, i - 1] = dy / dx
        LHS[i, i - Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = -2 * (Rbc[Ny - 1] * dy / dx + Bbc[Nx - 1] * dx / dy) 
        
    else:
        # Left Hand Side Matrix Coeffcients
        LHS[i, i] = -(-2 * dy / dx + 3 * dx / dy)
        LHS[i, i - 1] = dy / dx
        LHS[i, i + 1] = dy / dx
        LHS[i, i - Nx] = dx / dy        
        # Right Hand Side Vector values
        RHS[i] = -2 * Bbc[count] * dx / dy
    
    count += 1
    
# Erase counter just for order in variables
del count
        

# Verifying consistency of matrix        
plt.spy(LHS)
plt.savefig('spying.pdf')
plt.show()

# =============================================================================
# 
# =============================================================================        




