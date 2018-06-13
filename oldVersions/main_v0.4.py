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
 

########## Screen Definition ##########

# Number of X and Y pixels
screenX = 100
screenY = 100

if screenX & 1:
    alphaRange = np.arange(-(screenX-1)/2, (screenX-1)/2 +1)
else:
    alphaRange = np.arange(-(screenX)/2, (screenX)/2 +1)

if screenY & 1:
    betaRange = np.arange(-(screenY-1)/2, (screenY-1)/2 +1)
else:
    betaRange = np.arange(-(screenY)/2, (screenY)/2 +1)


print ("Number of Pixels: ", 4*alphaRange[0]*betaRange[0])




# Define a photon       
p = Photon(Alpha = 5., Beta = -5., i=np.pi/4)    

# Initial position and momentum of the photon
# p.xin
# p.kin


# Build the initial conditions needed to solve the geodesic equations
# [t, r, theta, phi, k_t, k_r, k_theta, k_phi] and stores in the variable
# p.iC of the Photon class

p.initConds(initCond(p.xin, p.kin, m.g(p.xin)))


finalPos, j = geoInt(m.geodesics, p.iC)

p.finalPosition(finalPos)


print(" Number of Steps to reach the Equatorial Plane:",j)
print()
print("t  ", p.iC[0], "  ", p.fP[0])
print("r  ", p.iC[1], "  ", p.fP[1])
print("th ", p.iC[2], "  ", p.fP[2])
print("ph ", p.iC[3], "  ", p.fP[3])
print("E ", p.iC[4], "  ", p.fP[4])
print("L ", p.iC[7], "  ", p.fP[7])


