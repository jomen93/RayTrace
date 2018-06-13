#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Geodesic Integration Module

This module integrtates the geodesic equations using as input:
    
geoEq: Geodesic equations in Hamiltonian form
[dt, dr, dtheta, dphi, dk_t, dk_r, dk_theta, dk_phi]


initConds: Initial conditions for the integration 
[t0, r0, theta0, phi0, k_t0, k_r0, k_theta0, k_phi0]

This module returns

[t, r, theta, phi, k_t, k_r, k_theta, k_phi] : Position of the photon in 
the equatorial plane (theta = pi/2)

j: Number of steps in the integration to arrive to the equatorial plane

@author: ashcat
"""
import numpy as np
from scipy.integrate import odeint


def geoInt(geoEq, initConds):
    geodesics = geoEq
    iC = initConds
    
    # Check if the photon has zero angular momentum. In this case the photon
    # goes directly into the black hole center
    
    if np.abs(iC[7]) < 0.000001:
       return [0,0,0,0,0,0,0,0], 0 
    
    # Tolerance in the location of the equatorial plane (cos †heta < tolerance)
    tolerance = 0.00005
    jmax = 100000
    coords = odeint(geodesics,iC,[0,-0.1])
    
    j = 1
    while np.cos(coords[j,2]) > tolerance :
        coordloop = odeint(geodesics,coords[j],[0,-0.01])
        coords = np.concatenate((coords,[coordloop[1]]))
        j = j + 1

        if j > jmax:
            return [0,0,0,0,0,0,0,0], j
        
    return coords[j-1,:], j