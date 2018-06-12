#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018

Photon Module

Here we define each of the photons in the observer's plane using the 
Photon class


Input for the created photon:
    
(X,Y): initial coordinates in the image plane
freq: frequency of the photon (needed for the 4-momentum)
D: distance between the image plane and the force center
i: inclination angle


Procedures:

.xin() : returns spherical initial coordinates (t=0, r, theta, phi)
.kin(): returns initial components of the 4-momentum (kt, kr, ktheta, kphi)  

@author: ashcat
"""
import numpy as np

class Photon:
    def __init__(self, Alpha=0, Beta=0, freq=1., D=100, i = np.pi/4):
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
        '''
        # Transformation from (Alpha, Beta, D) to (x, y, z) 
        self.x = - self.Alpha*np.cos(self.i) + self.D*np.sin(i)
        self.y = self.Beta
        self.z = self.Alpha*np.sin(self.i) + self.D*np.cos(i)
        '''
                
        # Transformation from (Alpha, Beta, D) to (r, theta, phi) 
        self.r = np.sqrt(self.Alpha**2 + self.Beta**2 + self.D**2)
        self.theta = np.arccos((self.Beta*np.sin(self.i)
                        +self.D*np.cos(self.i))/self.r)
        self.phi = np.arctan(self.Alpha/(self.D*np.sin(self.i)
                                        - self.Beta*np.cos(self.i)))
                
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
        
        self.ktheta = (- np.cos(self.i) 
                        + (self.Beta*np.sin(self.i) + self.D* np.cos(self.i))
                        *(self.D/(self.r**2)))*self.K0/np.sqrt(aux)
        
        self.kphi = - self.Alpha*np.sin(self.i)*self.K0/aux
        
        self.kt = np.sqrt(self.kr**2 + self.r**2 * self.ktheta**2 
                          + self.r**2*(np.sin(self.theta))**2 *self.kphi**2)
                
    def xin(self):
        '''
        Returns the initial values for the spherical coordinates of the photon
        (t=0, r, theta, phi)
        '''   
        return [0., self.r, self.theta, self.phi]


    def kin(self):
        '''
        Returns the initial values for the 4-momentum of the photon
        (kt, kr, ktheta, kphi)
        '''  
        return [self.kt, self.kr, self.ktheta, self.kphi]

