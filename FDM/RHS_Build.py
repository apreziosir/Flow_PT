#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:19:21 2017
Matrix and vector builder scripts for de FDM Laplace's equations
@author: apreziosir
"""

import numpy as np

# =============================================================================
# Internal nodes calculation in a matrix (subroutine for this script)
# =============================================================================

def internos(Nx, Ny):
    
    intern = Nx * Ny - (2 * Nx) - 2 * (Ny - 2)
    
    return intern

# =============================================================================
# RHS vector construction - Works for both Dirichlet and Neumann BC
# It is just a vector!
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
            rhs[0:Nx] = -Tbc
    
        elif i == (Ny - 1):
            rhs[-Nx:] = -Bbc
        
        # Bloques inteirores
        else:
            ind = i *Nx
            zblock[0] = -Lbc[i - 1]
            zblock[-1] = -Rbc[i - 1]
            rhs[ind:ind + Nx] = zblock

#   Este pedazo se descomenta cuando se necesite ver que el vector del lado 
#   derecho está bien construido o no            
#    print(rhs)
    return rhs