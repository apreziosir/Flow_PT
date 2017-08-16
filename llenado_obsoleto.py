#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 08:03:51 2017

@author: apreziosir
"""

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
plt.spy(LHS, size=6)
plt.savefig('spying.pdf')
plt.show()