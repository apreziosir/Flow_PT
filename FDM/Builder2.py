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
# Building the top boundary of the domain - it's always the same snce the BC for
# the top does not change with the type of problem
# ==============================================================================

def build_top(Nx, Ny):
    
    # The top part of the domainis always a Dirichlet BC, hence it is only built 
    # by ones in the left hand side of the matrix. A matrix is defined to take 
    # out all the values needed in the function
    Top = np.zeros((3, Nx))
    Top[0,:] = np.ones(Nx)
    Top[1,:] = np.linspace(Nx * (Ny - 1), (Nx * Ny) - 1, Nx)
    Top[2,:] = np.linspace(Nx * (Ny - 1), (Nx * Ny) - 1, Nx)    
    
    return Top

# ==============================================================================
# Building the bottom BC - It is a Neumann BC always since the flows are imposed
# and are losing, gaining or null fluxes
# ==============================================================================
    
def build_bot(Nx):
    
    # This part builds the bottom of the domain, but it is the top part of the 
    # matrix since the nodes are numbered from bottom to top
    # Defining matrix in order to take out all the values generated in the 
    # function
       
    
    return Bott
# ==============================================================================
# General matrix builder - builds all kind of matrices
# ==============================================================================

def gen_build(Nx, Ny, Lx, Ly, Diff, Neum, N_LR):
    
    # Building top vectors (data, i, j)
    T = build_top(Nx, Ny)
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