#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:45:42 2017
Function to build LHS of equation system in a different way. 
Takes into account Dirichlet and Neumann BC for different cases
@author: toni
"""

import numpy as np
import scipy.sparse as scsp

# ==============================================================================
# Building the top boundary of the domain
# ==============================================================================

def build_top(Nx):
    
    # The top row is always a Dirichlet BC, hence it is only built by ones in 
    # the left hand side of the matrix
    Top = np.zeros((3, Nx))
    Top[0,:] = np.ones(Nx)
    Top[1,:] = np.linspace(0, Nx - 1, Nx)
    Top[2,:] = np.linspace(0, Nx - 1, Nx)    
    
    return Top

# ==============================================================================

# ==============================================================================
# General matrix builder - builds all kind of matrices
# ==============================================================================

def gen_build(Nx, Ny, Diff, Neum, N_LR):
    
    # Building top vectors (data, i, j)
    T = build_top(Nx)
    Top_d = T[0,:]
    Top_i = T[1,:]
    Top_j = T[2,:]
    
    # Building bottom nodes vectors (data, i, j)
    Bot_d = 
    Bot_i = 
    Bot_j = 
    
    # Building left nodes vectors (data, i, j)
    Lft_d = 
    Lft_i = 
    Lft_j = 
    
    # Building right nodes vectors (data, i, j)
    Rgh_d = 
    Rgh_i = 
    Rgh_j = 

    # Building internal nodes vectors (data, i, j)    
    Int_d = 
    Int_i =
    Int_j = 
    
    # Concatenating vectors that compose coordinate matrix
    LHS_data = np.concatenate((Top_d, Bot_d, Lft_d, Rgh_d, Int_d), axis=0)
    LHS_i = np.concatenate((Top_i, Bot_i, Lft_i, Rgh_i, Int_i), axis=0)
    LHS_j = np.concatenate((Top_j, Bot_j, Lft_j, Rgh_j, Int_j), axis=0)
    
    # Building coordinate matrix - final step
    lhs = scsp.coo_matrix((LHS_data, (LHS_i, LHS_j)), shape = ((Nx * Ny), 
                       (Nx * Ny)))
    
    return lhs