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
    
def build_bot(Num, delta, coef):
    
    # This part builds the bottom of the domain, but it is the top part of the 
    # matrix since the nodes are numbered from bottom to top
    # Defining matrix in order to take out all the values generated in the 
    # function
    Bott = np.zeros((3, 2 * Num[0]))
    Bott[0,0:Nx] = np.ones(Num[0]) 
    Bott[0,-Nx:] = -np.ones(Num[0]) 
    Bott[1,:] = np.tile(np.linspace(0, Num[0] - 1, Num[0]), 2)
    Bott[2,:] = np.tile(np.linspace(0 + Num[0], 2 * Num[0] - 1, Num[0]), 2)
    
    
    return Bott

# ==============================================================================
# Building the left BC - Dirichlet case (this is a first approach and has to be 
# modified for modelling purposes)
# ==============================================================================
    
def build_left(Num):
    
    Left = np.zeros((3, Num[1] - 2))    
    Left[0,:] = np.ones(Num[1] - 2)
    Left[1,:] = np.arange(Num[0], Num[0] - 1, Num[0])
    Left[2,:] = np.arange(Num[0], Num[0] - 1, Num[0])    
    
    return Left

# ==============================================================================
# Building the left BC - Neumann case (this is one of the assumptions of our 
# model)
# ==============================================================================
    
def build_left_N(Num):
    
    Left = np.zeros((3, Num[1] - 2))    
    Left[0,:] = np.ones(Num[1] - 2)
    Left[1,:] = np.arange(Num[0], Num[0] - 1, Num[0])
    Left[2,:] = np.arange(Num[0], Num[0] - 1, Num[0])    
    
    return Left

# ==============================================================================
# Building the right BC - Dirichlet case (this is a first approach and has to be 
# modified for modelling purposes)
# ==============================================================================
    
def build_right(Num):
    
    Right = np.zeros((3, Num[1] - 2))    
    Right[0,:] = np.ones(Num[1] - 2)
    Right[1,:] = np.arange(Num[0], Num[0] - 1, Num[0])
    Right[2,:] = np.arange(Num[0], Num[0] - 1, Num[0])    
    
    return Right

# ==============================================================================
# Building the right BC - Neumann case (this is one of the assumptions of our 
# model)
# ==============================================================================
    
def build_right_N(Num):
    
    Right = np.zeros((3, Num[1] - 2))    
    Right[0,:] = np.ones(Num[1] - 2)
    Right[1,:] = np.arange(Num[0], Num[0] - 1, Num[0])
    Right[2,:] = np.arange(Num[0], Num[0] - 1, Num[0])    
    
    return Right

# ==============================================================================
# General matrix builder - builds all kind of matrices
# ==============================================================================

def gen_build(Num, Len, delta, coef, N_LR):
    
    # Building top vectors (data, i, j)
    T = build_top(Num)
    Top_d = T[0,:]
    Top_i = T[1,:]
    Top_j = T[2,:]
    
    # Building bottom nodes vectors (data, i, j)
    B = build_bot(Num, Diff, delta)
    Bot_d = B[0,:]
    Bot_i = B[1,:]
    Bot_j = B[2,:]
    
    # Building left nodes vectors (data, i, j)
    if N_LR == False : L = build_left(Num)
    else : L = build_left_N(Num)
    Lft_d = L[0,:]
    Lft_i = L[1,:]
    Lft_j = L[2,:]
    
    # Building right nodes vectors (data, i, j)
    if N_LR == False : R = build_left(Num)
    else : R = build_left_N(Num)
    Rgh_d = R[0,:]
    Rgh_i = R[1,:]
    Rgh_j = R[2,:]

    # Building internal nodes vectors (data, i, j)    
    I = build_int()
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