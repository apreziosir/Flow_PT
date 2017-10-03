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

def build_top(Num):
    
    # The top part of the domainis always a Dirichlet BC, hence it is only built 
    # by ones in the left hand side of the matrix. A matrix is defined to take 
    # out all the values needed in the function
    Top = np.zeros((3, int(Num[0])))
    Top[0,:] = np.ones(int(Num[0]))
    Top[1,:] = np.linspace(Num[0] * (Num[1] - 1), (Num[0] * Num[1]) - 1, Num[0])
    Top[2,:] = np.linspace(Num[0] * (Num[1] - 1), (Num[0] * Num[1]) - 1, Num[0])    
    
    return Top

# ==============================================================================
# Building the bottom Boundary - It is a Neumann BC always since the flows are 
# imposed and are losing, gaining or null fluxes
# ==============================================================================
    
def build_bot(Num):
    
    # This part builds the bottom of the domain, but it is the top part of the 
    # matrix since the nodes are numbered from bottom to top
    # Defining matrix in order to take out all the values generated in the 
    # function
    Bott = np.zeros((3, 2 * int(Num[0])))
    Bott[0,:] = np.concatenate((np.ones(int(Num[0])), -np.ones(int(Num[0]))), 
                                axis = 0) 
    Bott[1,:] = np.tile(np.linspace(0, Num[0] - 1, Num[0]), 2)
    Bott[2,:] = np.tile(np.linspace(0 + Num[0], 2 * Num[0] - 1, Num[0]), 2)
        
    return Bott

# ==============================================================================
# Building the left Boundary - Dirichlet case (this is a first approach and 
# has to be modified for modelling purposes)
# ==============================================================================
    
def build_left(Num):
    
    Left = np.zeros((3, int(Num[1] - 2)))    
    Left[0,:] = np.ones(int(Num[1] - 2))
    Left[1,:] = np.arange(Num[0], Num[0] - 1, Num[0])
    Left[2,:] = np.arange(Num[0], Num[0] - 1, Num[0])    
    
    return Left

# ==============================================================================
# Building the left Boundary - Neumann case (this is one of the assumptions of 
# our model)
# ==============================================================================
    
def build_left_N(Num):
    
    Left = np.zeros((3, 2 * (int(Num[1] - 2))))    
    Left[0,:] = np.concatenate((np.ones(int(Num[1] - 2)), -np.ones(int(Num[1] - 
                                2))), axis = 0)
    Left[1,:] = np.tile(np.arange(Num[0], Num[0] * (Num[1] - 1), Num[0]), 2)
    Left[2,:] = np.concatenate((np.arange(Num[0], Num[0] * (Num[1] - 1), 
            Num[0]), np.arange(Num[0] + 1, Num[0] * (Num[1] - 1) + 1, Num[0])),
            axis = 0)    
    
    return Left

# ==============================================================================
# Building the right Boundary - Dirichlet case (this is a first approach and 
# has to be modified for modelling purposes)
# ==============================================================================
    
def build_right(Num):
    
    Right = np.zeros((3, int(Num[1] - 2)))    
    Right[0,:] = np.ones(int(Num[1] - 2))
    Right[1,:] = np.arange(Num[0], Num[0] - 1, Num[0])
    Right[2,:] = np.arange(Num[0], Num[0] - 1, Num[0])    
    
    return Right

# ==============================================================================
# Building the right Boundary - Neumann case (this is one of the assumptions 
# of our model)
# ==============================================================================
    
def build_right_N(Num):
    
    Right = np.zeros((3, 2 * (int(Num[1] - 2))))    
    Right[0,:] = np.concatenate((-np.ones(int(Num[1] - 2)), np.ones(int(Num[1] - 
                                2))), axis = 0)
    Right[1,:] = np.tile(np.arange(2 * Num[0] - 1, Num[0] * Num[1] - 1, Num[0]), 
                                   2)
    Right[2,:] = np.concatenate((np.arange(2 * Num[0] - 2, Num[0] * Num[1] - 2, 
            Num[0]), np.arange(2 * Num[0] - 1, Num[0] * Num[1] - 1, Num[0])),
            axis = 0)    
    
    return Right

# ==============================================================================
# Building the internal nodes - Independent of Neumann or Dirichlet BC
# ==============================================================================
def build_int(Num, coef):
    
    # Number of internal nodes in the mesh
    nod_int = int((Num[0] * Num[1]) - (2 * Num[0] + 2 * (Num[1] - 2)))
    
    # Number of internal node coefficients (each node has 5 coefficients)
    Intern = np.zeros((3, 5 * nod_int))
    
    # Arreglo con nodos internos
    a = np.arange(Num[0], Num[0] * (Num[1] - 1), 1)
    b = a[a % Num[0] != 0]
    c = b[(b + 1) % Num[0] != 0]
    print(c)
    del(a, b)
    
    # Array to sum to the i column to obtain the j vector
    d = np.tile((np.array([-Num[0], -1, 0, 1, Num[0]])), nod_int)
    
    Intern[0,:] = np.tile(np.array([coef[2], coef[1], coef[0], coef[1], 
                                   coef[2]]), nod_int)
    Intern[1,:] = np.repeat(c, 5)
    Intern[2,:] = Intern[1,:] + d
    
    return Intern
# ==============================================================================
# General matrix builder - builds all kind of matrices (LHS)
# ==============================================================================

def gen_build(Num, Len, delta, coef, N_LR):
    
    # Building top vectors (data, i, j)
    T = build_top(Num)
    Top_d = T[0,:]
    Top_i = T[1,:]
    Top_j = T[2,:]
    
    # Building bottom nodes vectors (data, i, j)
    B = build_bot(Num)
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
    if N_LR == False : R = build_right(Num)
    else : R = build_right_N(Num)
    Rgh_d = R[0,:]
    Rgh_i = R[1,:]
    Rgh_j = R[2,:]

    # Building internal nodes vectors (data, i, j)    
    I = build_int(Num, coef)
    Int_d = I[0,:]
    Int_i = I[1,:]
    Int_j = I[2,:]
    
    # Concatenating vectors that compose coordinate matrix
    LHS_data = np.concatenate((Top_d, Bot_d, Lft_d, Rgh_d, Int_d), axis=0)
    LHS_i = np.concatenate((Top_i, Bot_i, Lft_i, Rgh_i, Int_i), axis=0)
    LHS_j = np.concatenate((Top_j, Bot_j, Lft_j, Rgh_j, Int_j), axis=0)
    
    # Building coordinate matrix - final step
    lhs = scsp.coo_matrix((LHS_data, (LHS_i, LHS_j)), shape = ((Num[0] * 
                          Num[1]), (Num[0] * Num[1])))
    
    return lhs