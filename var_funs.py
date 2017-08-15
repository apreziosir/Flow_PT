#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Functions for script of Laplace's equation for pressure in Particle Tracking 
model
Antonio Preziosi-Ribero / Universidad Nacional de Colombia
August 2017
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
    
    Tbc = np.linspace(0, Lx, Nx)    # Defining vector
    k = 2 * np.pi / Lambda          # Defining k
    dx = Lx / Nx                    # Calculating dx
    
    for i in range(0, len(Tbc)):
        Tbc[i] = hm * np.sin(k * ((dx/2) + i * dx))
        
    return Tbc

# =============================================================================
# Filling the bottom boundary condition according to problem
# It must vary (constant, linear)
# =============================================================================

def fill_bbc(Lx, Nx, hm, Lambda, bbc):
    
    Bbc = np.linspace(0, Lx, Nx)    # Defining vector
    k = 2 * np.pi / Lambda          # Defining k
    dx = Lx / Nx                    # Calculating dx
    
    for i in range(0, len(Bbc)):
        # Sine wave with hydrostatic increment (just water column)
        Bbc[i] = hm * np.sin(k * ((dx/2) + i * dx)) + bbc
        # Constant value
        # Bbc[i] = bbc
        
    return Bbc

# =============================================================================
# Filling left boundary condition - According to top boundary condition
# =============================================================================
    
def fill_lbc(Ly, Ny, Tbc):
    
    Lbc = np.linspace(0, Ly, Ny)
    dy = Ly / Ny
    
    for i in range(0, len(Lbc)):
        # Hydrostatic increment of pressure in the left boundary
        Lbc[i] = Tbc + dy / 2 + i * dy 
        # Constant value of pressure left boundary
        # Lbc[i] = Tbc
        
    return Lbc

# =============================================================================
# Filling right boundary condition - According to top boundary condition
# =============================================================================
    
def fill_rbc(Ly, Ny, Tbc):
    
    Rbc = np.linspace(0, Ly, Ny)
    dy = Ly / Ny
    
    for i in range(0, len(Rbc)):
        # Hydrostatic increment of pressure in the left boundary
        Rbc[i] = Tbc + dy / 2 + i * dy 
        # Constant value of pressure left boundary
        # Rbc[i] = Tbc
        
    return Rbc

# =============================================================================
# Calculate positions of nodes (works for x and y)
# =============================================================================
    
def positions(Lx, Ly, Nx, Ny):
    
    dx = Lx / Nx
    dy = Ly / Ny
    
    xn = np.zeros((Nx * Ny, 3))         # Node positions matrix
    
    for ic in range(0, Nx * Ny):
        xn[ic, 0] = int(ic)                          # Node id
        xn[ic, 1] = dx / 2 + (ic % Nx) * dx     # Node x position
        xn[ic, 2] = dy / 2 + (ic % Ny) * dy     # Node x position
        
        
    return xn
        
