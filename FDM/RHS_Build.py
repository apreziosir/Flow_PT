#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:19:21 2017
RHS vector builder script for de FDM Laplace's equation. 
@author: apreziosir
"""

import numpy as np

# ==============================================================================
# RHS vector construction - Works for both Dirichlet and Neumann BC
# It is just a vector! ALL OF THE SIGNS SHALLBE MANAGED IN THE ROUTINE 
# Bound_cond FOR CONSISTENCY IN THE METHOD. dO NOT PLAY AROUND WITH SIGNS
# ==============================================================================

def RHS_build(Tbc, Bbc, Lbc, Rbc):
    
    
    # Defining the Numbers Nx and Ny. They can be imported or estimated from
    # the length of the vectors that sotre the values in the boundary
    Nx = len(Tbc)
    Ny = len(Rbc) + 2
    zblock = np.zeros(Nx)
    
    # Esto es lo mismo que decir Nx * Ny
    rhs = np.zeros(Nx * Ny)
    
    # Inicio loop de llenado por bloques
    for i in range(0, Ny):
        
        # bloque superior (Tbc). Los signos se manejan desde la rutina Bound
        # _cond
        if i == 0:
            rhs[0:Nx] = Tbc
    
        elif i == (Ny - 1):
            rhs[-Nx:] = Bbc
        
        # Bloques inteirores
        else:
            # Los signos (si hay que cambiarlos), se manejan desde la funcion 
            # Bound_Cond. En este caso no hay mucho problema porque son ceros,
            # tener cuidado
            ind = i *Nx
            zblock[0] = Lbc[i - 1]
            zblock[-1] = Rbc[i - 1]
            rhs[ind:ind + Nx] = zblock

#   Este pedazo se descomenta cuando se necesite ver que el vector del lado 
#   derecho está bien construido o no            
    print((rhs, rhs.shape))
            
    return rhs