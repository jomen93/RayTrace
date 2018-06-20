#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Configuration File

  

@author: ashcat
"""

import numpy as np

############### Screen Definition ###############

ScreenType = 1

# 1: Image Plane
# 2: Point camera (not yet!)


############## Image Plane Options ##############

# Screen size (in ligth-years?)
Ssize = 6

# Resolution of the square screen: 
# Number of pixels in each side of the screen (N X N)
N = 10 

# Distance to the Black hole 
D = 1000.

# Inclination of the image plane
i = np.pi/4

#################################################


############### Metric Definition ###############

Metric = 1

# 1: Minkowski metric 
# 2: Schwarzschild metric 
# 3: Kerr metric (not yet!)

########### Minkowski metric Options ############



######### Schwarzschild metric Options ##########

# Mass of the black hole
M = 1

#################################################



############# Kerr metric Options ###############

# Mass of the black hole
# M = 1

# Angular Momentum Parameter
# a = 0.1 

#################################################