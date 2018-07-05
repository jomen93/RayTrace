#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Configuration File

  

@author: ashcat
"""

import numpy as np

''' Screen Definition '''

ScreenType = 1

# 1: Image Plane
# 2: Point camera (not yet!)

############## Image Plane Options ##############

# Screen size (in ligth-years?)
Ssize = 10

# Resolution of the square screen: 
# Number of pixels in each side of the screen (N X N)
N = 64

# Distance to the Black hole 
D = 1000.

# Inclination of the image plane
i = np.pi/3

#################################################


''' Metric Definition '''

Metric = 2

# 1: Minkowski metric 
# 2: Schwarzschild metric 
# 3: Kerr metric (not yet!)

# Parameters in the metric

M = 1.
a = 0.

########### Minkowski metric Options ############


######### Schwarzschild metric Options ##########

# M: Mass of the black hole

############# Kerr metric Options ###############

# M : Mass of the black hole

# a: Angular Momentum Parameter

#################################################


''' Accretion Structue Definition '''

Structure = 2

# Parameters in the Structure

R_in = 3.*M
R_out = 6.

#################################################

# 1: Simple Accretion Disk
# 2: Linear Spectrum Accretion Disk
# 2: Novikov-Thorne Thin Accretion Disk (not yet!)
