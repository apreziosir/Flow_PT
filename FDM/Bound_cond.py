#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:28:28 2017
Functions to fill boundary conditions for de FDM for Laplace's equation
@author: apreziosir
"""

import numpy as np

# =============================================================================
# Filling the vector with the boundary conditions that are going to be used for
# the Laplace's problem 
# =============================================================================
    
def fill_tbc(Len, Num, delta, hm, Lambda):
    
    Tbc = np.linspace(0, Len[0], Num[0])    # Defining vector
    k = 2 * np.pi / Lambda                  # Defining k
    
    for i in range(0, len(Tbc)):
        
        # Real function - should work if test function works
        Tbc[i] = hm * np.sin(k * (i * delta[0]))
        
        # Test function - just to test the program under known conditions
#        Tbc[i] = 7 * (i * dx) 
        
        # Test function 2
#        Tbc[i] = i * dx
    
#   Este pedazo imprime la condicion de contorno para comprobar su 
#   funcionamiento 
#    print(Tbc)    
    return Tbc

# =============================================================================
# Fill the Bottom Bbc vector with Neumann Boundary condition
# =============================================================================
    
def fill_bbc_N(Num, delta, q, K):
    
    # Converting q to consistent units (it is given in cm/day and it shall be 
    # converted to m/s)
    q_r = q / (8.64e6)
    
    # It has the sign that has to go in the RHS vector 
    Bbc = -np.ones(int(Num[0])) * (delta[1] * q_r) / K  
    
    return Bbc

# =============================================================================
# Filling left boundary condition - According to top boundary condition
# =============================================================================
    
def fill_lbc(Ly, Ny, Tbc):
    
    dy = np.abs(Ly / (Ny - 1))
    Lbc = np.zeros(Ny - 2)    
    
    for i in range(0, len(Lbc)):
        # Hydrostatic increment of pressure in the left boundary
        Lbc[i] = Tbc + (i + 1) * dy 
#         Constant value of pressure left boundary
        # Lbc[i] = Tbc
        
        # Test function - just to test the program under known conditions
#        Lbc[i] = 8 * dy * (i + 1)
        
        # Test case 2
#        Lbc[i] = (i + 1) * dy
        
    return Lbc

# ==============================================================================
# Filling left boundary condition when there is a V. Neumann imposed with no 
# flow at the boundaries u = 0, v = 0. This routine is used also for the right
# boundary when the gradient of pressure is 0
# ==============================================================================

def fill_lbc_N(Ny):
    
    Lbc = np. zeros(Ny - 2)
    
    return Lbc

# =============================================================================
# Filling right boundary condition - According to top boundary condition
# =============================================================================
    
def fill_rbc(Lx, Ly, Ny, Tbc):
    
    Rbc = np.zeros(Ny - 2)
    dy = np.abs(Ly / Ny)
    
    for i in range(0, len(Rbc)):
        # Hydrostatic increment of pressure in the left boundary
        Rbc[i] = Tbc + (i + 1) * dy 
        # Constant value of pressure left boundary
        # Rbc[i] = Tbc
        
        # Test function - just to test the program under known conditions
#        Rbc[i] = 6 * Lx * (i + 1) * dy + 7 * Lx + 8 * (i + 1) * dy
        
        # Test case 2
#        Rbc[i] = Lx + (i + 1) * dy
        
    return Rbc