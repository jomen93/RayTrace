#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:51:45 2018

Kerr metric

ds^2 = -(1-2Mr/Sigma)dt^2 - 4 M a r sin(theta)^ 2 /Sigma dtdphi+ Sigma/Delta dr^2  
       + Sigma dtheta^2 + (r^2 + a^2 + (2M(a^2)r sin(theta)^2)/Sigma) sin(theta)^2 dphi^2

- Event horizon at r=M + sqrt(M^2 - a^2)
- ISCO 

@author: ashcat
"""

import numpy as np
import myconfig as cfg


def g(x, M = cfg.M, a = cfg.a):
    '''
    This procedure contains the Kerr metric components 
    in Boyer-Linquist coordinates
    '''

    # Coordinates 
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    
    # Metric components

    Sigma = r**2 + a**2 * np.cos(theta)**2
    Delta = r**2 - 2.*M*r + a**2

    gtt = -(1. - 2.*M*r/Sigma)
    gtphi = - (2.*M*a*r*np.sin(theta)**2)/Sigma
    grr = Sigma/Delta
    gthth = Sigma
    gphph = (r**2 + a**2 + (2.*M*(a**2)*r*np.sin(theta)**2)/Sigma)*np.sin(theta)**2
    
    return [gtt, gtphi, grr, gthth, gphph]


def rEH(M=cfg.M, a=cfg.a):
    return M + np.sqrt(M**2 - a**2) 

def ISCOco(M=cfg.M, a=cfg.a):
    Z1 = M + np.cbrt(M**2 - a**2)*( np.cbrt(M + a) + np.cbrt(M - a) )
    Z2 = np.sqrt(3*a**2 + Z1**2)
    return 3.*M + Z2 - np.sqrt( (3*M - Z1)*(3*M + Z1 + 2*Z2) )

def ISCOcounter(M=cfg.M, a=cfg.a):
    Z1 = M + np.cbrt(M**2 - a**2)*( np.cbrt(M + a) + np.cbrt(M - a) )
    Z2 = np.sqrt(3*a**2 + Z1**2)
    return 3.*M + Z2 + np.sqrt( (3.*M - Z1)*(3.*M + Z1 + 2.*Z2) )
    
def geodesics(x, tau, M=cfg.M, a=cfg.a):
    '''
    This procedure contains the geodesic equations in Hamiltonian form 
    for the Kerr metric
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

    # Auxiliar Functions
    Sigma = r**2 + a**2 * np.cos(theta)**2
    Delta = r**2 - 2.*M*r + a**2

    W = E*(r**2 + a**2) - a*L 
    partXi = r**2 + (L - a*E)**2 + a**2 *(1 - E**2)*np.cos(theta)**2 + (L**2 *np.cos(theta)**2)/(np.sin(theta)**2)
    Xi = W**2 - Delta*partXi

    dXidE = 2.*W*(r**2 + a**2) + 2.*a*Delta*(L - a*E*np.sin(theta)**2)
    dXidL = -2.*a*W + 2.*a*E*Delta - 2.*L*Delta/(np.sin(theta)**2)

    dXidr = 4.*r*E*W - 2.*(r - M)*partXi - 2.*r*Delta 

    dAdr = (r - M)/Sigma - (r*Delta)/(Sigma**2)
    dBdr = -r/Sigma**2
    dCdr = dXidr/(2.*Delta*Sigma) - (Xi*(r-M))/(Sigma*Delta**2) - r*Xi/(Delta*Sigma**2)

    auxth = (a**2) * np.cos(theta)*np.sin(theta)

    dAdth = Delta*auxth/(Sigma**2)
    dBdth = auxth/(Sigma**2)
    dCdth = ((1-E**2)*auxth + L**2 * np.cos(theta)/(np.sin(theta)**3) )/Sigma + (Xi/(Delta*Sigma**2))*auxth



    # Geodesics differential equations 
    dtdtau = dXidE/(2.*Delta*Sigma)
    drdtau = (Delta/Sigma)*pr
    dthdtau = pth/Sigma
    dphidtau = - dXidL/(2.*Delta*Sigma)
    
    dptdtau = 0.
    dprdtau = -dAdr* pr**2 - dBdr* pth**2 + dCdr 
    dpthdtau = -dAdth* pr**2 - dBdth* pth**2 + dCdth 
    dpphidtau = 0.
    

    dxdtau = [dtdtau, drdtau, dthdtau, dphidtau, 
              dptdtau, dprdtau, dpthdtau, dpphidtau]
    return dxdtau
