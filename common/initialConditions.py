#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Initial Conditions Module

This module prepares the initial values of coordinates and momentum to
solve the geodesic equations.

Given    
x: initial coordinates calculated for a particular photon
k: initial components of the momentum for a particular photon
g: metric components at the photon location

this module returns

[t, r, theta, phi, k_t, k_r, k_theta, k_phi] :initial conditions needed to solve 
the geodesic equations 

@author: ashcat
"""

def initCond(x, k, g):
    
    # Coordinates and momentum components
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    pt = k[0]
    pr = k[1]
    pth = k[2]
    pphi = k[3]
    
    # Metric components
    gtt = g[0]
    gtphi = g[1]
    grr = g[2]
    gthth = g[3]
    gphph = g[4]
    
    # Lower k
    p_t = gtt*pt + gtphi*pphi
    p_r = grr*pr
    p_th = gthth*pth
    p_phi =gphph*pphi
    
    return [t, r, theta, phi, p_t, p_r, p_th, p_phi]
        