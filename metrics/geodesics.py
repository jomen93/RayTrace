#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 08:57:20 2018
Geodesic Equations

@author: ashcat
"""
import numpy as np
import metric as g


def geodesics(initcoords, initk):
    r = initcoords[1]
    theta = initcoords[2]
    E = initk[0]
    Lz = initk[3]
    ptheta = initk[2]
    
    Sigma = g.Sigma(r, theta)
    Delta = g.Delta(r,theta)
    A = g.A(r, theta)
    Carter = g.Carter(theta, E, ptheta, Lz)
    R= g.R(r, theta, E, ptheta, Lz, a=0., M=1.)
    Theta = g.Theta(theta, E, ptheta, Lz, a=0.)
    print(R,Theta)
    
    

    
    
geodesics([0., 1., np.pi, 0.], [1., 0., 0., 1.])
