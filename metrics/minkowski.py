#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:51:45 2018

Minkowskian metric

ds^2 = -dt^2 + dr^2 + r^2 d\theta^2 + r^2 \sin^2 (\theta) d\phi^2

@author: ashcat
"""

import numpy as np

def g(x):
    # Coordinates 
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    
    # Metric components
    gtt = -1
    gtphi = 0
    grr = 1
    gthth = r**2
    gphph = r**2 * np.sin(theta)**2
    
    return [gtt, gtphi, grr, gthth, gphph]
    
def geodesics(x, tau):
    '''
    This procedure contains the geodesic equations in Hamiltonian form 
    for the Minkowski metric
    '''
    # Coordinates and momentum components
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    pt = x[4]
    pr = x[5]
    pth = x[6]
    pphi = x[7]
    
    E = - pt
    L = pphi
    
    # Geodesics differential equations 
    dtdtau = E
    drdtau = pr
    dthdtau = pth/r**2
    dphidtau = L/((r**2)*np.sin(theta)**2)
    
    dptdtau = 0
    dprdtau = pth**2/r**3 + L**2/((r**3)*np.sin(theta)**2)
    dpthdtau = (np.cos(theta)/np.sin(theta)**3)*(L**2/r**2)
    dpphidtau = 0
    
    dxdtau = [dtdtau, drdtau, dthdtau, dphidtau, 
              dptdtau, dprdtau, dpthdtau, dpphidtau]
    return dxdtau
