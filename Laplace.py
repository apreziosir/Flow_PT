# -*- coding: utf-8 -*-
"""
Script to calculate Laplace equation and then calculate flow for the particle 
tracking model
Antonio Preziosi-Ribero / Universidad Nacional de Colombia
August 2017
"""

import numpy as np
from var_funs import alt_media, fill_tbc, fill_bbc, fill_lbc, fill_rbc

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

Nx = 200            # Elements in x direction (number)
Ny = 10             # Elements in y direction  (number)
bbc = d + Ly

# =============================================================================
# Calculate hm value for the problem assigned and set a vector for the 
# top boundary condition values - outside functions
# =============================================================================

hm = alt_media(U, H, d)
Tbc = fill_tbc(Lx, Nx, hm, Lambda)
Bbc = fill_bbc(Lx, Nx, hm, Lambda, bbc)
Lbc = fill_lbc(Ly, Ny, Tbc[0])
Rbc = fill_rbc(Ly, Ny, Tbc[Nx - 1])




