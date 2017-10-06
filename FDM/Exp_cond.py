#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:23:28 2017
Functions to estimate the values taken out from the lab and from Packman's 
paper. The only condition that is not calculated in this script is the Hm, that 
corresponds to the pressure boundary condition at the top of the domain. 
That value is estimated in the Bound_cond script for making thing easy
@author: apreziosir
"""

# =============================================================================
# This function calculates the hm for the Laplace's equation of the particle
# tracking model 
# =============================================================================

def alt_media(U, H, d):
    
    if H/d <= 0.34:        
        hm = 0.28 * (U ** 2 / (2 * 9.81)) * ((H/d) / 0.34) ** (3/8)    
    else:        
        hm = 0.28 * (U ** 2 / (2 * 9.81)) * ((H/d) / 0.34) ** (3/2)   
    return hm 

# =============================================================================
# Function to calculate the value of the inflow/ outflow velocity
# It uses the Darcy law, porosity and permeability of the material that's 
# been modeled
# =============================================================================
    
def inf_vel(phi, q, K):
    
    # According to inquiries this value can be easily changed
    v = q / K
    
    return v

