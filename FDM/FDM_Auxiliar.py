#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:24:15 2017

@author: apreziosir
"""

import numpy as np

# =============================================================================
# Function that calculates nonzero elements using just Nx and Ny
# =============================================================================

def nzero(Nx, Ny):
    
    nel = Nx * Ny
    elext = 2 * Nx + 2 * (Ny - 2)
    elint = nel - elext
    nzeroint = elint * 5 
    nzeroext = elext * 1
    
    nzerotot = nzeroint + nzeroext
    
    print('Number of nodes:')
    print(nel)
    print('Total size of matrix:')
    print(nel ** 2)
    print('Nonzero elements in matrix:')
    print(nzeroint + nzeroext)
    print('Relation nonzero/total:')
    print((nzeroint + nzeroext) / (nel ** 2))
    
    return nzerotot
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
    
    Tbc = np.linspace(0, Lx, Nx)        # Defining vector
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

def fill_bbc(Tbc, Nx, Lx, Ly):
    
    Bbc = np.zeros(Nx)
    dx = Lx / (Nx - 1)
#    Bbc = Tbc + Ly    # Defining vector
    
#   Test case 
    for i in range(0, len(Bbc)):
        
        Bbc[i] = 6 * Ly * (i * dx) + 7 * (i * dx) + 8 * Ly
        
    return Bbc

# =============================================================================
# Filling left boundary condition - According to top boundary condition
# =============================================================================
    
def fill_lbc(Ly, Ny, Tbc):
    
    dy = Ly / (Ny - 1)
    Lbc = np.zeros(Ny - 2)    
    
    for i in range(0, len(Lbc)):
        # Hydrostatic increment of pressure in the left boundary
#        Lbc[i] = Tbc + (i + 1) * dy 
        # Constant value of pressure left boundary
        # Lbc[i] = Tbc
        
        # Test function - just to test the program under known conditions
        Lbc[i] = 8 * dy * (i + 1)
        
    return Lbc

# =============================================================================
# Filling right boundary condition - According to top boundary condition
# =============================================================================
    
def fill_rbc(Lx, Ly, Ny, Tbc):
    
    Rbc = np.zeros(Ny - 2)
    dy = Ly / Ny
    
    for i in range(0, len(Rbc)):
        # Hydrostatic increment of pressure in the left boundary
#        Rbc[i] = Tbc + (i + 1) * dy 
        # Constant value of pressure left boundary
        # Rbc[i] = Tbc
        
        # Test function - just to test the program under known conditions
        Rbc[i] = 6 * Lx * (i + 1) * dy + 7 * Lx + 8 * (i + 1) * dy
        
    return Rbc

# =============================================================================
# Calculate positions of nodes (works for x and y)
# =============================================================================
    
def positions(Lx, Ly, Nx, Ny):
    
    x = np.linspace(0 ,Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    
    xn = np.zeros((Nx * Ny, 3))             # Node positions matrix
    
    xn[:, 0] = np.arange(0, Nx * Ny, 1)
    xn[:, 1] = np.tile(x, Ny)
    xn[:, 2] = np.repeat(y, Nx)
        
        
    return xn

# =============================================================================
# RHS vector construction
# =============================================================================

def RHS_build(Tbc, Bbc, Lbc, Rbc):
    
    
    Nx = len(Tbc)
    Ny = len(Rbc) + 2
    zblock = np.zeros(Nx)
    
    # Esto es lo mismo que decir Nx * Ny
    rhs = np.zeros(Nx * Ny)
    
    # Inicio loop de llenado por bloques
    for i in range(0, Ny):
        
        # bloque superior (Tbc)
        if i == 0:
            rhs[0:Nx] = Tbc
    
        elif i == (Ny - 1):
            rhs[-Nx:] = Bbc
        
        # Bloques inteirores
        else:
            ind = i *Nx
            zblock[0] = Lbc[i - 1]
            zblock[-1] = Rbc[i - 1]
            rhs[ind:ind + Nx] = zblock
            
    return rhs
    
# =============================================================================
# LHS matrix construction
# (Coordinate system storage mode to save space)
# =============================================================================  

def LHS_build(Nx, Ny):
    
    null = nzero(Nx, Ny)
    
    lhs = 0
    return lhs