#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Image Plane definition Module

This module defines the image plane as screen to receive the photons

Has a procedure called screen to define a square NxN pixels screen
and the Photon class which creates each of the photons that will be traced back
onto the accretion structure

@author: ashcat
"""


import numpy as np

def screen(Ssize, N):
    '''
    Defines a square NxN pixels screen.
    Returns the ranges of the cordinates Alpha and Beta in the image plane and 
    the value of the maximum value of the coordinates
    '''
    
    if N & 1:
        numPixels = N +1
    else:
        numPixels = N 

    alphaRange = np.linspace(-Ssize, Ssize, numPixels)
    betaRange = np.linspace(-Ssize, Ssize, numPixels)
    
    print ("Size of the screen in Pixels:", numPixels, "X", numPixels)
    print ("Number of Pixels: ", (numPixels)*(numPixels))
    
    return alphaRange, betaRange, numPixels

"""
Photon Class

Here we define each of the photons in the observer's plane using the 
Photon class


Input for the created photon:
    
(X,Y): initial coordinates in the image plane
freq: frequency of the photon (needed for the 4-momentum)
D: distance between the image plane and the force center
i: inclination angle


Procedures:

.initConds(initialConds) : Stores the initial values of coordinates and 
                           momentum needed to solve the geodesic equations.  
              
.finalPosition(finalPos): Stores the initial values of coordinates and 
                          momentum needed to solve the geodesic equations.

@author: ashcat
"""

class Photon:
    def __init__(self, Alpha=1., Beta=0., freq=1., D=100., i = np.pi/4):
        '''
        Given the initial coordinates in the image plane (X,Y), the distance D 
        to the force center and inclination angle i, this calculates the 
        initial coordinates in spherical coordinates (r, theta, phi)
        ''' 
        # Initial Cartesian Coordinates in the Image Plane
        self.Alpha = Alpha
        self.Beta = Beta
        self.D = D
        self.i = i
        
        # Transformation from (Alpha, Beta, D) to (r, theta, phi) 
        self.r = np.sqrt(self.Alpha**2 + self.Beta**2 + self.D**2)
        self.theta = np.arccos((self.Beta*np.sin(self.i)
                        +self.D*np.cos(self.i))/self.r)
        self.phi = np.arctan(self.Alpha/(self.D*np.sin(self.i)
                                        - self.Beta*np.cos(self.i)))
        
        '''
        Initial values in spherical coordinates of the photon
        (t=0, r, theta, phi)
        '''   
        self.xin = [0., self.r, self.theta, self.phi]
        
                
        '''
        Given the frequency value k0, this calculates the initial values for 
        the 4-momentum of the photon
        (kt, kr, ktheta, kphi)
        '''  
        self.K0 =  freq
        
        # Initial 4-momentum components
        self.kr =  (self.D/self.r)*self.K0
        
        aux = self.Alpha**2 + (-self.Beta*np.cos(self.i) 
                                + self.D*np.sin(self.i))**2 
        
        self.ktheta = (self.K0/np.sqrt(aux))*(
                - np.cos(self.i) 
                + (self.Beta*np.sin(self.i) + self.D* np.cos(self.i))
                   *(self.D/(self.r**2)))
        
        self.kphi = - self.Alpha*np.sin(self.i)*self.K0/aux
        
        self.kt = np.sqrt(self.kr**2 + self.r**2 * self.ktheta**2 
                          + self.r**2*(np.sin(self.theta))**2 *self.kphi**2)
        
        '''
        Initial values in spherical coordinates for the 4-momentum of 
        the photon(kt, kr, ktheta, kphi)
        '''   
        self.kin = [self.kt, self.kr, self.ktheta, self.kphi]
                
    
    def initConds(self, initialConds):
        '''
        Stores the initial values of coordinates and momentum needed 
        to solve the geodesic equations.
        '''
        self.iC = initialConds
        return 
    
    def finalPosition(self, finalPos):
        '''
        Stores the initial values of coordinates and momentum needed 
        to solve the geodesic equations.
        '''
        self.fP = finalPos
        return 

