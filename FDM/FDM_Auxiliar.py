#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:24:15 2017

@author: apreziosir
"""

import numpy as np

# =============================================================================
# This function calculates the hm for the Laplace's equation of the particle
# tracking model 
# =============================================================================

def alt_media(U, H, d):
    
    if H/d <= 0.34:        
        hm = 0.28 * (U ** 2 / (2 * 9.81)) * ((H/d) / 0.34) ** (3/8)    
    else:        
        hm = 0.28 * (U ** 2 / (2 * 9.81)) * ((H/d) / 0.34) ** (3/2)   
    return hm 

# =============================================================================
# Filling the vector with the boundary conditions that are going to be used for
# the Laplace's problem 
# =============================================================================
    
def fill_tbc(Lx, Nx, hm, Lambda):
    
    Tbc = np.linspace(0, Lx, Nx + 1)        # Defining vector
    k = 2 * np.pi / Lambda                  # Defining k
    dx = Lx / (Nx - 1)                      # Calculating dx
    
    for i in range(0, len(Tbc)):
        
        # Real function - should work if test function works
#        Tbc[i] = hm * np.sin(k * i * dx)
        
        # Test function - just to test the program under known conditions
        Tbc[i] = 7 * (i * dx) 
        
    return Tbc

# =============================================================================
# Filling the bottom boundary condition according to problem
# It must vary (constant, linear)
# =============================================================================

def fill_bbc(Tbc, Ly):
    
    Bbc = Tbc + Ly    # Defining vector
    
        
    return Bbc

# =============================================================================
# Filling left boundary condition - According to top boundary condition
# =============================================================================
    
def fill_lbc(Ly, Ny, Tbc):
    
    Lbc = np.linspace(0, Ly, Ny)
    dy = Ly / Ny
    
    for i in range(0, len(Lbc)):
        # Hydrostatic increment of pressure in the left boundary
#        Lbc[i] = Tbc + dy / 2 + i * dy 
        # Constant value of pressure left boundary
        # Lbc[i] = Tbc
        
        # Test function - just to test the program under known conditions
        Lbc[i] = -8 * ((dy/2) + i * dy)
        
    return Lbc

# =============================================================================
# Filling right boundary condition - According to top boundary condition
# =============================================================================
    
def fill_rbc(Lx, Ly, Nx, Ny, Tbc):
    
    Rbc = np.linspace(0, Ly, Ny)
    dy = Ly / Ny
    dx = Lx / Nx
    
    for i in range(0, len(Rbc)):
        # Hydrostatic increment of pressure in the left boundary
#        Rbc[i] = Tbc + dy / 2 + i * dy 
        # Constant value of pressure left boundary
        # Rbc[i] = Tbc
        
        # Test function - just to test the program under known conditions
        Rbc[i] = -6 * Lx * ((dy/2) + i * dy) + 7 * ((dx/2) + i * dx)
        
    return Rbc