#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Configuration File

  

@author: ashcat
"""

import numpy as np


############## Metric Definition ################

Metric = "Schwarzschild"

# "Minkowski" metric 
# "Schwarzschild" metric 
# "Kerr" metric 

# Parameters in the metric

M = 1. ## Find the units of this mass 
a = 0.7

########### Minkowski metric Options ############

# No options

######### Schwarzschild metric Options ##########

# M: Mass of the black hole

############# Kerr metric Options ###############

# M : Mass of the black hole
# a: Angular Momentum Parameter

#################################################





######### Accretion Structure Definition ########

Structure = 3

# 1: Simple Accretion Disk
# 2: Linear Spectrum Accretion Disk
# 3: Infinite Accretion Disk with a decreasing exponential spectrum
# 4: Novikov-Thorne Thin Accretion Disk (not yet!)

# Parameters in the Structure

R_out = 12.*M

corotation = True


#################################################





############### Screen Definition ###############

ScreenType = 1

# 1: Image Plane
# 2: Point camera (not yet!)

############## Image Plane Options ##############

# Screen size (in units of M)
Ssize = 20.*M


############# Point Camera Options ##############

# Field of Vieew (in radians)
FoV = np.pi/110.


############ General Screen Options #############

# Resolution of the square screen: 
# Number of pixels in each direction of the screen (N X N)
N = 32

# Distance to the Black hole (in units of M)
D = 1000.*M

# Inclination of the image plane
i = np.pi/2.5

#################################################





################## File Output ################## 

# Name of the output file
fileName = Metric + "_" + str(N) + "x" + str(N) 

#################################################


