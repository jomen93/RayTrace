#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Test Initial Momentum Module

In this module we test the initial components of the momentum calculated by
the Photon class in the photon.py module.

We consider 5 photons in the observer's plane, located at the points:

      (Alpha, Beta)
1: at (0,0)
2: at (5,5)
3: at (5,-5)
4: at (-5,5)
5: at (-5,-5)

The photons are moving in the direction perpendicular to the plane.
In the observer's plane coordinates, K^\mu = (K0, K0, 0, 0)

We obtain the components of the momentum in spherical coordinates with the
.kin() procedure in the Photon class

@author: ashcat
"""
import numpy as np
from photon import Photon
        
angle = np.pi/4

# Initial positions of the 5 photons using the Photon class procedures   
kinit = np.array([ Photon(i=angle).kin(),
                  Photon(Alpha=5., Beta=5., i=angle).kin(), 
                  Photon(Alpha=5., Beta=-5., i=angle).kin(), 
                  Photon(Alpha=-5., Beta=5., i=angle).kin(),
                  Photon(Alpha=-5., Beta=-5., i=angle).kin()])

# Transformation from (r,theta,phi) to cartesian coordinates (x,y,z)     
print(kinit)
