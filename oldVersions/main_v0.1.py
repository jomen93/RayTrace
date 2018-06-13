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
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from photon import Photon
from initialConditions import initCond 
from metrics import minkowski as m       
 
# Define a photon       
p1 = Photon(Alpha = 5.0, Beta = 0.)    
p2 = Photon(Alpha = -5.0, Beta = 5.)    

# Initial position ans momentum of the photon
xin1 = p1.xin() 
kin1 = p1.kin()

xin2 = p2.xin() 
kin2 = p2.kin()


# Build the initial conditions needed to solve the geodesic equations
# [t, r, theta, phi, k_t, k_r, k_theta, k_phi]
coord01 = initCond(xin1, kin1, m.g(xin1))
coord02 = initCond(xin2, kin2, m.g(xin2))
#print(coord01)


# Solve the geodesic Equations
tau = np.linspace(0,-105,105)

coords1 = odeint(m.geodesics,coord01,tau)

coords2 = odeint(m.geodesics,coord02,tau)

'''
print (coords[:,0:4])
print(coords[len(coords)-1,1], "dE =",coords[0,4]-coords[len(coords)-1,4],
             "  dL =", coords[0,7]-coords[len(coords)-1,7])
'''
z1 = coords1[:,1]*np.cos(coords1[:,2])
x1 = coords1[:,1]*np.sin(coords1[:,2])*np.cos(coords1[:,3])
y1 = coords1[:,1]*np.sin(coords1[:,2])*np.sin(coords1[:,3])

z2 = coords2[:,1]*np.cos(coords2[:,2])
x2 = coords2[:,1]*np.sin(coords2[:,2])*np.cos(coords2[:,3])
y2 = coords2[:,1]*np.sin(coords2[:,2])*np.sin(coords2[:,3])


'''
xinit = np.array([ Photon().xin(),Photon(X=5., Y=5.).xin(), 
                  Photon(X=5., Y=-5.).xin(), Photon(X=-5., Y=5.).xin(),
                  Photon(X=-5., Y=-5.).xin(), [0,0,0,0] ])

z = xinit[:,1]*np.cos(xinit[:,2])
x = xinit[:,1]*np.sin(xinit[:,2])*np.cos(xinit[:,3])
y = xinit[:,1]*np.sin(xinit[:,2])*np.sin(xinit[:,3])
'''
print(x1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x1, y1, z1)
ax.plot(x2, y2, z2)
plt.show()
