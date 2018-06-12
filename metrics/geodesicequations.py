#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:51:45 2018

@author: ashcat
"""

import numpy as np
import metric as g

def geodesics(x):
    # Coordinates and momentum components
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    pt = x[4]
    pr = x[5]
    pth = x[6]
    pphi = x[7]
    
    # Geodesics differential equations 
    dtdtau = g.dAdE/( 2 * g.Delta(r, theta) * g.Sigma(r,theta) )
    drdtau = (g.Delta(r,theta) /g.Sigma(r,theta) ) * pr
    dthdtau = pth/g.Sigma(r, theta) 
    dphidtau = - g.dAdL/( 2 * g.Delta(r, theta) * g.Sigma(r,theta) )
    
    dptdtau = 0
    dprdtau = 0
    dpthdtau = 0
    dpphidtau = 0
    
    dxdtau = [dtdtau, drdtau, dthdtau, dphidtau, 
              dptdtau, dprdtau, dpthdtau, dpphidtau]
    return dxdtau