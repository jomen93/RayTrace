#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:51:45 2018

@author: ashcat
"""

# Kerr Metric



import numpy as np

def Sigma(r,theta,a=0.):
    return r**2 + a**2 * (np.cos(theta))**2

def Delta(r, theta, a=0., M=1.):
    return r**2 -2*M*r + a**2

def A(r, theta, a=0., M=1.):
    return (r**2 + a**2)**2 - a**2 * Delta(r, theta, a, M) * (np.sin(theta))**2

def Carter(theta, E, ptheta, Lz, a=0.):
    C = ptheta**2 + Lz**2 * (np.cos(theta)/np.sin(theta))**2 
    - a**2 * E**2 * (np.cos(theta))**2
    return C
    
def R(r, theta, E, ptheta, Lz, a=0., M=1.):
    R = ((r**2 + a**2)*E - a*Lz)**2 
    - Delta(r, theta, a, M)*(Carter(theta, E, ptheta, Lz, a=0.) 
    + (Lz - a*E)**2)
    return R

def Theta(theta, E, ptheta, Lz, a=0.):
    Theta = Carter(theta, E, ptheta, Lz, a=0.) - (-a**2 * E**2 
                      + (Lz/np.sin(theta))**2)*(np.cos(theta))**2
    return Theta