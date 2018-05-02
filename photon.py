#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Initial Conditions Module

Input    
(X,Y): initial coordinates in the image plane
freq: frequency of the photon (needed for the 4-momentum)
D: distance between the image plane and the force center
i: inclination angle


Procedures

.initcoords() : returns spherical initial coordinates (t=0, r, theta, phi)
.initk(): returns initial components of the 4-momentum (kt, kr, ktheta, kphi)  

@author: ashcat
"""
import numpy as np


class Photon:
    def __init__(self, X=0, Y=0, freq=1., D=100, i = np.pi/4):
        '''
        Given the initial coordinates in the image plane (X,Y) and the distance 
        to the force center D and inclination angle i, this calculates the 
        initial values for the spherical coordinates (r, theta, phi)
        ''' 
        # Initial Cartesian Coordinates in the Image Plane
        self.X = X
        self.Y = Y
        self.D = D
        self.i = i
        
        # Transformation from (X, Y, D) to (r, theta, phi) 
        self.r = np.sqrt(self.X**2 + self.Y**2 + self.D**2)
        self.theta = np.arccos((self.Y*np.sin(self.i)
                        +self.D*np.cos(self.i))/self.r)
        self.phi = np.arctan(self.X/(self.D*np.sin(self.i)
                                        - self.Y*np.cos(self.i)))
        
        '''
        Given the frequency value k0, this calculates the initial values for 
        the 4-momentum of the photon
        (kt, kr, ktheta, kphi)
        '''  
        self.k0 = - freq
        
        # Initial 4-momentum components
        self.kr = - (self.D/self.r)*self.k0
        aux = np.sqrt(self.X**2 
                      + (self.D*np.sin(self.i) - self.Y*np.cos(self.i))**2)
        self.ktheta = (np.cos(self.i) 
                        - (self.Y*np.sin(self.i) + self.D* np.cos(self.i))
                        *(self.D/(self.r**2)))*self.k0/aux
        self.kphi = self.X*np.sin(self.i)*self.k0/aux
        
        self.kt = np.sqrt(self.kr**2 + self.r**2 * self.ktheta**2 
                          + self.r**2*(np.sin(self.theta))**2*self.kphi**2)
        
    def initcoords(self):
        '''
        Returns the initial values for the spherical coordinates of the photon
        (t=0, r, theta, phi)
        '''   
        return [0., self.r, self.theta, self.phi]


    def initk(self):
        '''
        Returns the initial values for the 4-momentum of the photon
        (kt, kr, ktheta, kphi)
        '''  
        return [self.kt, self.kr, self.ktheta, self.kphi]
        
        
    
print(Photon().initcoords())
print(Photon().initk())