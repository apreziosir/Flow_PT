#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 10:24:15 2017
Funciones varias para calcular la ecuación de Laplace en 2 Dimensiones
Se separan funciones de acuerdo a tareas del script principal
@author: apreziosir
"""

import numpy as np
import scipy.sparse as scsp
import os
import matplotlib.pyplot as plt

# =============================================================================
# Function that calculates nonzero elements using just Nx and Ny
# =============================================================================

def nzero(Nx, Ny, Neum):
    
    nel = Nx * Ny
    
    if Neum == False : elext = 2 * Nx + 2 * (Ny - 2)
    else : elext = 5 * Nx + 2 * (Ny - 2)
        
    elint = nel - elext    
    nzeroint = elint * 5 
    nzeroext = elext * 1
    
    nzerotot = nzeroint + nzeroext
    
    print('Number of nodes:')
    print(nel)
    print('Total size of matrix:')
    print(nel ** 2)
    print('Nonzero elements in matrix:')
    print(nzeroint + nzeroext)
    print('Relation nonzero/total:')
    print((nzeroint + nzeroext) / (nel ** 2))
    
    return nzerotot

# =============================================================================
# Calculate positions of nodes (works for x and y)
# =============================================================================
    
def positions(Lx, Ly, Nx, Ny):
    
    x = np.linspace(0 ,Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    
    xn = np.zeros((Nx * Ny, 3))             # Node positions matrix
    
    xn[:, 0] = np.arange(0, Nx * Ny, 1)
    xn[:, 1] = np.tile(x, Ny)
    xn[:, 2] = np.repeat(y, Nx)
        
        
    return xn

# =============================================================================
# Comparison with theoretical functions (just for test cases, not for real 
# case)
# =============================================================================
    
def comp(RTA, Nx, Ny, Lx, Ly):
    
    # Generando espacio para la solucion
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(Ly, 0, Ny)
    X, Y = np.meshgrid(x, y)
    
    # Calculo de matriz respuesta analítica de la primera funcion de prueba
    # 6xy +7x + 8y
    z = 6 * X * Y + 7 * X + 8 * Y
    
    err = z
    
#    # Calculo de la matriz de respuesta analitica de la segunda funcion de 
#    # prueba x + y
#    z = X + Y
#    
#    err = RTA - z
    
    # Imprimiendo en pantalla el valor de la solucion analitica
    print('El valor de la respuesta analitica es... ')
    print(z)
    
    # Plotting the analytical solution
    CS = plt.contourf(X, Y, z)
    cbar = plt.colorbar(CS)
    plt.gca().set_aspect(9, adjustable='box')
    plt.ylim((Ly, 0))
    
    return err