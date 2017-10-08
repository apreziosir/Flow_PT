#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 13:33:51 2017
Plotting a 2D field - Function that takes in array and limits to plot in a 
heatmap. 
This function is built since the plotting function is called more than once in 
the LaplaceFDM script. 
@author: toni
"""

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# Plotting heatmap with Matplotlib
# ==============================================================================

def Plot_HM(RTA, Len, Num):
    
    x = np.linspace(0, Len[0], Num[0])
    y = np.linspace(0, Len[1], Num[1])
    X, Y = np.meshgrid(x, y)
    
    CS4 = plt.contourf(X, Y, RTA)
    cbar = plt.colorbar(CS4)
    #CS4.set_clim(vmin=-10000, vmax=1000)
    #plt.clim(-np.amax(RTA),np.amax(RTA))
    plt.gca().set_aspect(20, adjustable='box')
    plt.ylim((Len[1], 0))
    #cbar.Normalize(clip=False)
    #plt.clabel(CS4, fmt='%2.1f', colors='w', fontsize=14)
    plt.show()