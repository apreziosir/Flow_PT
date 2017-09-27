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
    
def fill_tbc(Lx, Nx, hm, Lambda):
    
    Tbc = np.linspace(0, Lx, Nx)            # Defining vector
    k = 2 * np.pi / Lambda                  # Defining k
    dx = np.abs(Lx / (Nx - 1))              # Calculating dx
    
    for i in range(0, len(Tbc)):
        
        # Real function - should work if test function works
        Tbc[i] = hm * np.sin(k * (i * dx))
        
        # Test function - just to test the program under known conditions
#        Tbc[i] = 7 * (i * dx) 
        
        # Test function 2
#        Tbc[i] = i * dx
    
#   Este pedazo imprime la condicion de contorno para comprobar su 
#   funcionamiento 
#    print(Tbc)    
    return Tbc

# =============================================================================
# Filling the bottom boundary condition according to problem
# It must vary (constant, linear)
# =============================================================================

def fill_bbc(Tbc, Nx, Lx, Ly):
    
    Bbc = np.zeros(Nx)
    dx = np.abs(Lx / (Nx - 1))
    Bbc = Tbc - np.abs(Ly)    # Defining vector
    
#   Test case 
#    for i in range(0, len(Bbc)):
        # Test case 1
#        Bbc[i] = 6 * Ly * (i * dx) + 7 * (i * dx) + 8 * Ly
        
        # Test case 2
#        Bbc[i] = Ly + i * dx
        
    return Bbc

# =============================================================================
# Fill the Bottom Bbc vector with Neumann Boundary condition
# =============================================================================
    
def fill_bbc_N(Nx, Dif, q, dy, K):
    
    Bbc = np.ones(Nx) * (2 * Dif * q) / (dy * K)  
    
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