#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Solver de Laplace con diferencias finitas - Se hace con matrices eficientes
Created on Tue Sep 12 10:11:10 2017
@author: apreziosir
"""

import numpy as np
import scipy.sparse as scsp
import matplotlib.pyplot as plt
from FDM_Auxiliar import alt_media, fill_tbc, fill_bbc, fill_rbc, fill_lbc

# =============================================================================
# Input variables - flow conditions
# =============================================================================

U = 0.35            # Mean flow velocity in m/s
H = 0.05            # Dune height
d = 0.20            # Mean depth of flow

# =============================================================================
# Input variables - domain and bed characteristics
# =============================================================================

Lx = 6.00           # Length of the flume (m) (considered for numerical model) 
Ly = 0.30           # Depth of bed (m)
Lambda = 0.60       # Wavelength of bedform (m)

# =============================================================================
# Numerical model input parameters
# =============================================================================

Nx = 4            # Nodes in x direction (number)
Ny = 4            # Nodes in y direction  (number)

# =============================================================================
# Calculate hm value for the problem assigned and set a vector for the 
# boundary condition values - outside functions
# =============================================================================

hm = alt_media(U, H, d)
Tbc = fill_tbc(Lx, Nx, hm, Lambda)
Bbc = fill_bbc(Tbc, Ly)
#Lbc = fill_lbc(Ly, Ny, Tbc[0])
#Rbc = fill_rbc(Lx, Ly, Nx, Ny, Tbc[Nx - 1])
