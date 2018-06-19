#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Main Module

  

@author: ashcat
"""
import numpy as np
from photon import Photon
from initialConditions import initCond 
from metrics import minkowski as m  
from geodesicIntegration import geoInt  
from screen.imagePlane import screen   
 

########## Screen Definition ##########
# Screen size (in ligth-years?)
Ssize = 3

# Resolution of the screen: 
# Number of pixels in each side of the screen (N X N)
N = 6 

# Ranges of the Alpha and Beta coordinates in the image plane given by the
# screen.imagePlane module    
alphaRange, betaRange, numPixels = screen(Ssize, N)

# Distance to the Black hole 
D = 1000.

# Inclination of the image plane
i = np.pi/4

#######################################


rDataArray = np.zeros((numPixels,numPixels))


for j in range(0,numPixels):
    for k in range(0,numPixels):        
        # Define a photon       
        p = Photon(Alpha = alphaRange[k], Beta = betaRange[numPixels-1-j], D = D, i = i)
        
        # Build the initial conditions needed to solve the geodesic equations
        # [t, r, theta, phi, k_t, k_r, k_theta, k_phi] and stores in the variable
        # p.iC of the Photon class
        
        p.initConds(initCond(p.xin, p.kin, m.g(p.xin)))
               
        finalPos, l = geoInt(m.geodesics, p.iC, intStep = 0.1)
        
        p.finalPosition(finalPos)

        print(" Number of Steps to reach the Equatorial Plane:",l)
        print()
        print("Alpha=", k, " Beta=", numPixels-1-j)
        print(p.iC[4], p.fP[4])
        print(p.iC[7], p.fP[7])
        # Store value of the radius of the photon in the equatorial plane
        # in the rDataArray
        
        rDataArray[j, k] = p.fP[1]
        

print(rDataArray)


