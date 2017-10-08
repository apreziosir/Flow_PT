#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Solver de Laplace con diferencias finitas - Se hace con matrices eficientes
Incluye el calculo de la red de flujo con los valores de presion (gradiente del 
campo escalar de presiones) -ESTO SE DEBE TERMINAR!!!
SUPOSICION FUERTE: El campo de velocidades tiene condiciones no slip en las 
fronteras laterales. Lo mismo sucede en la condicion interior
Created on Tue Sep 12 10:11:10 2017
@author: apreziosir
"""

# Python libraries to make the program work. All of them are standard Python, 
# nothing out of the ordinary
import numpy as np
import scipy.sparse as scsp
import scipy.io 
import matplotlib.pyplot as plt
# Funcions written by author in Python for different purposes
import Exp_cond as EC
import Bound_cond as BC
import FDM_Auxiliar as FDMAux 
import RHS_Build as RHSB
import LHS_Build as LHSB
import Plot_HM as HM

# ==============================================================================
# Input variables - flow conditions
# ==============================================================================

U = 0.15                    # Mean flow velocity in m/s 
H = 0.015                   # Dune height (m)
d = 0.10                    # Mean depth of flow (m)
phi = 0.33                  # Porosity of material (nondimensional)
q = 25                      # Inflow or downflow velocity (+ up / - down)(cm/d)
K = 0.1195                  # Hydraulic conductivity 
N_LR = True                 # Neumann condition in the sides?

# ==============================================================================
# Input variables - domain and bed characteristics
# ==============================================================================

Lx = 6.40           # Length of the flume (m) (considered for numerical model) 
Ly = 0.25           # Depth of bed (m)
Lambda = 0.15       # Wavelength of bedform (m)
Dif = 1.0           # Diffusion coefficient (just for fun - not used in this 
                    # case)

# ==============================================================================
# Numerical model input parameters - modify to refine mesh
# ==============================================================================

Nx = 100                # Nodes in x direction (number)
Ny = 100                 # Nodes in y direction  (number)

# ==============================================================================
# Setting up vectors to carry values into functions with less arguments
# ==============================================================================

# Set up mesh - vectors with lengths, Number of nodes and deltas
Len = np.zeros(2)
Len[0] = Lx
Len[1] = Ly

Num = np.zeros(2)
Num[0] = Nx
Num[1] = Ny

# Calculating delta_x and delta_y for the matrix
delta = FDMAux.dx_dy(Lx, Nx, Ly, Ny)

# Node positions with id for follow up
xn = FDMAux.positions(Lx, Ly, Nx, Ny)

# Calculating coefficients for matrix
coef = FDMAux.coeff(Dif, delta)

# ==============================================================================
# Calculate hm value for the problem assigned, inflow Darcy velocity and 
# estimate the mean free flow velocity value using the Manning formula
# ==============================================================================

# Mean head over bed - formula from Elliott and Brooks (1997)
hm = EC.alt_media(U, H, d)

# Inflow/outflow velocity - from Darcy's law (Bear 1975)
v = EC.inf_vel(phi, q, K)

# Mean free flow velocity from Manning's formula
# ACA SE DEBE PONER LA FORMULA PARA EL CALCULO DE U MEDIA Y NO DEJARLO A LO QUE 
# SE NOS DE LA GANA (VERIFICAR QUE SEA CERCANO A 15 CM/S)

# ==============================================================================
# Creating vectors that store the value of the Boundary condition values for 
# each one of the boundaries. This part takes the BC value as it is, hence the 
# sign should be changed in the construction of the RHS matrix
# ==============================================================================

# Top
Tbc = BC.fill_tbc(Len, Num, delta, hm, Lambda)
# Bottom
Bbc = BC.fill_bbc_N(Num, delta, q, K)
# Left
if N_LR == False : Lbc = BC.fill_lbc(Ly, Ny, Tbc[0])
else : Lbc = BC.fill_lbc_N(Ny)
# Right
if N_LR == False : Rbc = BC.fill_rbc(Lx, Ly, Ny, Tbc[Nx - 1])
else : Rbc = BC.fill_lbc_N(Ny)

# ==============================================================================
# Building RHS of system to solve Laplace Equation 
# (Results vector) - This vector is built the same way if there is a Dirichlet
# or Neumann BC since it is just a vector
# ==============================================================================

RHS = RHSB.RHS_build(Tbc, Bbc, Lbc, Rbc)

# ==============================================================================
# Building LHS matrix to solve the system
# Coordinate system storage - Later transformed to CSR (by Python script)
# ==============================================================================

LHS = LHSB.gen_build(Num, Len, delta, coef, N_LR)
LHS = LHS.tocsr()
scipy.io.mmwrite('matrix_test', LHS)

# Checking matrix construction
plt.spy(LHS, markersize = 2)
plt.show()

# ==============================================================================
# Solving linear system for the pressure field (Laplace's equation)
# ==============================================================================

# Gradiente con jugado - no esta funcionando
#P = scsp.linalg.cg(LHS, RHS)
P = scsp.linalg.spsolve(LHS, RHS)

# ==============================================================================
# Reshaping and plotting solution to check values and different Bc's
# ==============================================================================

RTA = np.reshape(P, (Ny, Nx), order='C')

np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

#print('*---------------------------*')
#print('El valor máximo de RTA es:')
#print(np.amax(RTA))
#
#print('El valor mínimo de RTA es:')
#print(np.amin(RTA))
#
#print('El valor de RTA es... ')
#print(RTA)
#print('*---------------------------*')

# ==============================================================================
# Plotting the pressure field of the problem
# ==============================================================================

HM.Plot_HM(RTA, Len, Num)

# ==============================================================================
# Comparing with test functions values (just for test cases, not real 
# conditions)
# ==============================================================================

#err = comp(RTA, Nx, Ny, Lx, Ly)
##CS5 = plt.contourf(X, Y, err)
##cbar = plt.colorbar(CS5)
##plt.gca().set_aspect(9, adjustable='box')
##plt.ylim((Ly, 0))
##plt.show()
#
## Guardando archivos para comparar en csv
#np.savetxt('sol_numerica.csv', RTA, delimiter=' ', newline='\n')
#np.savetxt('sol_analit.csv', err, delimiter=' ', newline='\n')

# ==============================================================================
# Calculating the velocity field with Darcy's law q = K grad(h). The numpy 
# gradient function is used with the dx and dy taken from the delta array
# ==============================================================================

Vel = np.gradient(RTA, delta[0], delta[1])
u = Vel[0]
v = Vel[1]

# ==============================================================================
# Plotting the velocity fields in each one of the directions that was calculated
# ==============================================================================

HM.Plot_HM(u, Len, Num)
HM.Plot_HM(v, Len, Num)