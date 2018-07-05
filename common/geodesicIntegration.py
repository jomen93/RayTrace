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


def geoInt(geoEq, initConds, rEH = 0. ,intStep = 0.1):
    geodesics = geoEq
    iC = initConds

    # Energy and Momentum constants
    E = iC[4]
    L = iC[7]
    # Carter's Constant
    Carter = iC[6]**2 + ((np.cos(iC[2])**2)/(np.sin(iC[2])**2)) * L**2
     
    
    # Tolerance in the location of the equatorial plane (cos †heta < tolerance)
    tolerance = 0.00005
    # Maximum number of steps admitted in the integration
    jmax = 100000
    
    # First iteration (included to define the variable coords, where the solution will be stored)
    coords = odeint(geodesics,iC,[0,-intStep])
    
    j = 1
    while np.cos(coords[j,2]) > tolerance :
        coordloop = odeint(geodesics,coords[j],[0,-intStep])

        # Event Horizon Condition
        if coordloop[1,1] <= rEH + 0.00001 :
            return [0,0,0,0,0,0,0,0], j
        '''
        # Carter's Constant conservation condition
        if coordloop[1,6] > 0 :
            pthNew = np.sqrt(Carter - ((np.cos(coordloop[1,2])**2)/(np.sin(coordloop[1,2])**2)) * L**2 )
        else :
            pthNew = - np.sqrt(Carter - ((np.cos(coordloop[1,2])**2)/(np.sin(coordloop[1,2])**2)) * L**2 )

        percentualDiffth = np.abs(((coordloop[1,6]-pthNew)/coordloop[1,6])*100)

        #print(coordloop[1,6], pthNew, percentualDiffth)
        
        if percentualDiffth > 1 :
            coordloop[1,6] = pthNew

        
        if coordloop[1,4] > 0 :
            prNew = np.sqrt(coordloop[1,5]**2 + coordloop[1,1]**(-2) * coordloop[1,6]**2
                        + coordloop[1,1]**(-2)*(np.sin(coordloop[1,2]))**(-2) *coordloop[1,7]**2)
        else :
            prNew = - np.sqrt(coordloop[1,5]**2 + coordloop[1,1]**(-2) * coordloop[1,6]**2
                        + coordloop[1,1]**(-2)*(np.sin(coordloop[1,2]))**(-2) *coordloop[1,7]**2)

        percentualDiffr = np.abs(((coordloop[1,4]-prNew)/coordloop[1,4])*100)

        print(coordloop[1,4], prNew, percentualDiffr)
        
        
        if np.abs(percentualDiffr) > 1 :
            coordloop[1,4] = prNew
        '''
        # Stores new line of information in coords
        coords = np.concatenate((coords,[coordloop[1]]))
        
        #print(coordloop[1,1])

        '''
        print(((iC[4] + np.sqrt(coordloop[1,5]**2 + coordloop[1,1]**(-2) * coordloop[1,6]**2 
                          + coordloop[1,1]**(-2)*(np.sin(coordloop[1,2]))**(-2) *coordloop[1,7]**2))/iC[4])*100)
        '''

        j = j + 1

        if j > jmax:
            return [0,0,0,0,0,0,0,0], j
        
    return coords[j-1,:], j