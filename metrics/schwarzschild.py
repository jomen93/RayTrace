#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:51:45 2018

Schwarzschild metric

ds^2 = -(1-2M/r)dt^2 + dr^2 /(1-2M/r) + r^2 dtheta^2 + r^2 sin^2 (theta) dphi^2

- Event horizon at r=2M
- ISCO at r = 6M

@author: ashcat
"""

import numpy as np
import myconfig as cfg


def g(x, M = cfg.M):
    '''
    This procedure contains the Schwarzschild metric components 
    in spherical coordinates
    '''

    # Coordinates 
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    
    # Metric components
    gtt = -(1. - 2.*M/r)
    gtphi = 0
    grr = 1./(1.- 2.*M/r)
    gthth = r**2
    gphph = r**2 * np.sin(theta)**2
    
    return [gtt, gtphi, grr, gthth, gphph]


def rEH(M = cfg.M):
    return 2.*M

def ISCOco(M = cfg.M):
    return 6.*M

def ISCOcounter(M = cfg.M):
    return 6.*M
    
def geodesics(x, tau, M=cfg.M):
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

    # Conserved Quantities
    E = - pt
    L = pphi

    # Geodesics differential equations 
    dtdtau = E*r**2./(r**2 - 2.*M*r)
    drdtau = (1. - 2.*M/r)*pr
    dthdtau = pth/r**2
    dphidtau = L/((r**2)*np.sin(theta)**2)
    
    dptdtau = 0.
    dprdtau = -M*(pr**2/r**2) + pth**2/r**3 + L**2/((r**3)*np.sin(theta)**2) - M*(E**2/(r-2.*M)**2) 
    dpthdtau = (np.cos(theta)/np.sin(theta)**3)*(L**2/r**2)
    dpphidtau = 0.
    

    dxdtau = [dtdtau, drdtau, dthdtau, dphidtau, 
              dptdtau, dprdtau, dpthdtau, dpphidtau]
    return dxdtau
