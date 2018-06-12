#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Test Initial Conditions Module

In this module we test the initial conditions obtained form the Photon class
in the photon.py module.

We consider 5 photons in the observer's plane:

      (Alpha,Beta)
1: at (0,0)
2: at (5,5)
3: at (5,-5)
4: at (-5,5)
5: at (-5,-5)

The initial spherical coordinates (r,theta,phi) are obtained using the 
.xin() procedure in the Photon class.

By plotting the location of this photons in the (x,y,z) coordinates, it is
identifed the observer's plane in space

@author: ashcat
"""

import numpy as np
import matplotlib.pyplot as plt
from photon import Photon
        
angle = np.pi/4
# Initial positions of the 5 photons using the Photon class procedures   
xinit = np.array([ Photon(i=angle).xin(),
                  Photon(Alpha=5., Beta=5., i=angle).xin(), 
                  Photon(Alpha=5., Beta=-5., i=angle).xin(), 
                  Photon(Alpha=-5., Beta=5., i=angle).xin(),
                  Photon(Alpha=-5., Beta=-5., i=angle).xin(), 
                  [0,0,0,0] ])
    
print(xinit)    

# Transformation from (r,theta,phi) to cartesian coordinates (x,y,z)     
z = xinit[:,1]*np.cos(xinit[:,2])
x = xinit[:,1]*np.sin(xinit[:,2])*np.cos(xinit[:,3])
y = xinit[:,1]*np.sin(xinit[:,2])*np.sin(xinit[:,3])


# 3D-scatter plot of the points indicating the initial photons 
# in the observer's plane
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z)
ax.set_xlim(0, 90)
ax.set_zlim(0, 90)
ax.set_ylim(0, 90)
plt.show()
