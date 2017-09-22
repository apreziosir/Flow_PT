#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:19:21 2017
Matrix and vector builder scripts for de FDM Laplace's equations
@author: apreziosir
"""

import numpy as np
import scipy.sparse as scsp

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
    
# =============================================================================
# LHS matrix construction - with Dirichlet BC
# (Coordinate system storage mode to save space)
# =============================================================================  

def LHS_build(Nx, Ny, dx, dy, Dif):
    
    intern = Nx * Ny - (2 * Nx) - 2 * (Ny - 2)
    
    # Calculating matrix coefficients
    aE_aW = (Dif / (dx ** 2)) 
    aN_aS = (Dif / (dy ** 2))
    aP = - 2 * (aE_aW + aN_aS) 
#    print(aP)
    
    # Elementos de la diagonal mayor de la matriz LHS
    Diag_d = np.ones(Nx * Ny)
    Diag_i = np.arange(0, (Nx * Ny), 1)
    Diag_j = np.arange(0, (Nx * Ny), 1)
    
    # Elementos de las diagonales menores de la matriz (cada vector de valores
    # se repite dos veces ya que hay dos de cada digonal)    
    Diag2_d = aE_aW * np.ones(intern * 2)
    Diag3_d = aN_aS * np.ones(intern * 2)
        
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
            Diag_d[i] = aP
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