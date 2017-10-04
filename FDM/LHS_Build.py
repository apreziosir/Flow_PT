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
    Bott[2,:] = np.tile(np.linspace(Num[0], 2 * Num[0] - 1, Num[0]), 2)
        
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

"""
================================================================================
OLD ROUTINE - IT IS NOT USED CURRENTLY IN THE PROGRAM!!!! DELETE LATER
================================================================================

This part was the previous build of the LSH matrix. It is more expensive, 
inefficient and less modular. However, it was a nice idea
Antonio Preziosi Ribero
20171004

# =============================================================================
# LHS matrix construction - with Dirichlet BC (does not work with Neumann!!!)
# (Coordinate system storage mode to save space)
# =============================================================================  

def LHS_build(Nx, Ny, dx, dy, Dif):
    
    intern = internos(Nx, Ny)
    
    # Calculating matrix coefficients
    cff = coeff(Dif, dx, dy)
    
    # Elementos de la diagonal mayor de la matriz LHS
    Diag_d = np.ones(Nx * Ny) * -1
    Diag_i = np.arange(0, (Nx * Ny), 1)
    Diag_j = np.arange(0, (Nx * Ny), 1)
    
    # Elementos de las diagonales menores de la matriz (cada vector de valores
    # se repite dos veces ya que hay dos de cada digonal)    
    Diag2_d = cff[1] * np.ones(intern * 2)
    Diag3_d = cff[2] * np.ones(intern * 2)
        
    # Coordenadas de elementos, no se llenan sino hasta el loop
    Diag2_i = np.zeros(intern * 2)
    Diag2_j = np.zeros(intern * 2)
    Diag3_i = np.zeros(intern * 2)
    Diag3_j = np.zeros(intern * 2)
    i0 = 0
    i1 = 1
    
    for i in range(Nx, (Nx * (Ny - 1))):
        
        # Probando condiciones para poder hacer el llenado de los vectores de 
        # posiciones de los índices en parte interna. 
#        print('voy en.....')
#        print(i)
        
        if i % Nx != 0 and i % Nx != (Nx - 1):
            
#            print('Entre al IF en i =')
#            print(i)
            # Replace one with coefficient aP in [i, i]
            Diag_d[i] = cff[0]
#            print(Diag_d)
            # Data and coefficients for the off-diagonal elements
            Diag2_i[i0] = i
            Diag2_i[i1] = i
            Diag2_j[i0] = i - 1 
            Diag2_j[i1] = i + 1
#            print(Diag2_j)
            Diag3_i[i0] = i
            Diag3_i[i1] = i
#            print(Diag3_i)
            Diag3_j[i0] = i - Nx
            Diag3_j[i1] = i + Nx
            i0 += 2
            i1 += 2
    
#    print(Diag_d)
    # Ensamblando vectores de datos para la matriz 
    LHS_data = np.concatenate((Diag_d, Diag2_d, Diag3_d), axis=0)
    np.savetxt('dataLHS.csv', LHS_data)
    LHS_i = np.concatenate((Diag_i, Diag2_i, Diag3_i), axis=0)
    np.savetxt('iLHS.csv', LHS_i)
    LHS_j = np.concatenate((Diag_j, Diag2_j, Diag3_j), axis=0)
    np.savetxt('jLHS', LHS_j)
    
    # Ensambalndo matriz en formato coordenado    
    lhs = scsp.coo_matrix((LHS_data, (LHS_i, LHS_j)), shape = ((Nx * Ny), 
                       (Nx * Ny)))
    
    return lhs

# =============================================================================
# LHS matrix construction - with Neumann BC for the pressure in the bottom part
# of the domain. It is a vertical velocity, hence the horizontal velocity 
# equals zero
# (Coordinate system storage mode to save space)
# =============================================================================
    
def LHS_build_N(Nx, Ny, dx, dy, Dif):
    
    # Calculating internal nodes (doesn't include any boundary)
    intern = internos(Nx, Ny)
    
    # Calculating matrix coefficients
    cff = coeff(Dif, dx, dy) 
    
    # Elementos de la diagonal mayor de la matriz LHS
    Diag_d = np.ones(Nx * Ny) * -1
    Diag_i = np.arange(0, (Nx * Ny), 1)
    Diag_j = np.arange(0, (Nx * Ny), 1)
    
    # Elementos de las diagonales menores de la matriz (cada vector de valores
    # se repite dos veces ya que hay dos de cada digonal)    
    # La Diag_3 tiene Nx elementos mas por la BC del fondo del dominio
    Diag2_d = cff[1] * np.ones(intern * 2)
    Diag3_d = cff[2] * np.ones(intern * 2 + Nx)
    
    # Coordenadas de elementos, no se llenan sino hasta el loop. Aplica el 
    # mismo comentario para los arreglos de la Diag3 (i y j)
    Diag2_i = np.zeros(intern * 2)
    Diag2_j = np.zeros(intern * 2)
    Diag3_i = np.zeros(intern * 2 + Nx)
    Diag3_j = np.zeros(intern * 2 + Nx)
    i0 = 0
    i1 = 1 
    
    for i in range(Nx, (Nx * (Ny - 1))):
        
        # Probando condiciones para poder hacer el llenado de los vectores de 
        # posiciones de los índices en parte interna. 
#        print('voy en.....')
#        print(i)
        
        if i % Nx != 0 and i % Nx != (Nx - 1):
            
#            print('Entre al IF en i =')
#            print(i)
            # Replace one with coefficient aP in [i, i]
            Diag_d[i] = cff[0]
#            print(Diag_d)
            # Data and coefficients for the off-diagonal elements
            Diag2_i[i0] = i
            Diag2_i[i1] = i
            Diag2_j[i0] = i - 1 
            Diag2_j[i1] = i + 1
#            print(Diag2_j)
            Diag3_i[i0] = i
            Diag3_i[i1] = i
#            print(Diag3_i)
            Diag3_j[i0] = i - Nx
            Diag3_j[i1] = i + Nx
            i0 += 2
            i1 += 2    
    
    # Esta sección llena el final de los vectores de Diag_3, donde se almacenan
    # los valores y coordenadas de cada una de las filas que hacen la condicion
    # de contorno de Neumann en la matriz
    
    # Impresion antes del llenado de la ultima parte de Diag3
#    print('*--------*')
#    print(Diag3_d.shape)
#    print(Diag3_i.shape)
#    print(Diag3_j.shape)
#    print('*--------*')
#    print(Diag3_d)
#    print(Diag3_i)
#    print(Diag3_j)
#    print('*--------*')
    
    # Corregir este pedazo para poder escribir bien la matriz
    Fn = Nx * (Ny - 1)
    Ln = (Nx * Ny) - 1
    Diag3_d[-Nx:] =  -2 * cff[2]
    Diag3_i[-Nx:] = np.linspace(Fn, Ln, Nx)
    Diag3_j[-Nx:] = np.linspace(Fn - Nx, Ln - Nx, Nx)
    
    
    # Pruebas de impresión en pantalla de la ultima parte de los arreglos que 
    # hacen parte de las diferencias verticales
#    print('*--------*')
#    print(Diag3_d.shape)
#    print(Diag3_i.shape)
#    print(Diag3_j.shape)
#    print('*--------*')
#    print(Diag3_d)
#    print(Diag3_i)
#    print(Diag3_j)
#    print('*--------*')
    
    # Ensamblando vectores de datos para la matriz 
    LHS_data = np.concatenate((Diag_d, Diag2_d, Diag3_d), axis=0)
    np.savetxt('dataLHS.csv', LHS_data)
    LHS_i = np.concatenate((Diag_i, Diag2_i, Diag3_i), axis=0)
    np.savetxt('iLHS.csv', LHS_i)
    LHS_j = np.concatenate((Diag_j, Diag2_j, Diag3_j), axis=0)
    np.savetxt('jLHS', LHS_j)    
    
    # Ensambalndo matriz en formato coordenado    
    lhs = scsp.coo_matrix((LHS_data, (LHS_i, LHS_j)), shape = ((Nx * Ny), 
                       (Nx * Ny)))
    
    return lhs
"""