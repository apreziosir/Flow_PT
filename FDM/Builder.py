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
# Internal nodes calculation in a matrix (subroutine for this script)
# =============================================================================

def internos(Nx, Ny):
    
    intern = Nx * Ny - (2 * Nx) - 2 * (Ny - 2)
    
    return intern

# =============================================================================
# Calculating the matrix coefficients for the LHS (used in both matrix 
# constructions) 
# =============================================================================
    
def coeff(Dif, dx, dy):
    
    # Vector that stores coefficients (0 = ap, 1 = ae_aw, 2 = an_as)
    coeff = np.zeros(3)
    
    # Calculating matrix coefficients
    coeff[1] = (Dif / (dx ** 2)) 
    coeff[2] = (Dif / (dy ** 2))
    coeff[0] = - 2 * (coeff[2] + coeff[1]) 
#    print(coeff)
    
    return coeff

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