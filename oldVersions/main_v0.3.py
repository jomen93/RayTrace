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
p = Photon(Alpha = 5., Beta = -5., i=np.pi/4)    

# Initial position and momentum of the photon
xin = p.xin 
kin = p.kin


# Build the initial conditions needed to solve the geodesic equations
# [t, r, theta, phi, k_t, k_r, k_theta, k_phi]

p.initConds(initCond(xin, kin, m.g(xin)))
iC = p.iC

coords = odeint(m.geodesics,iC,[0,-1])

j = 1

while np.cos(coords[j,2]) > 0.00005 :
    coordloop = odeint(m.geodesics,coords[j],[0,-0.01])
    coords = np.concatenate((coords,[coordloop[1]]))
    j = j + 1

print(coords[0,:])
print(coords[j-1,:])
print(j)
print(coords[j,2])
print(np.cos(coords[j,2]))


z = coords[:,1]*np.cos(coords[:,2])
x = coords[:,1]*np.sin(coords[:,2])*np.cos(coords[:,3])
y = coords[:,1]*np.sin(coords[:,2])*np.sin(coords[:,3])




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.set_xlim(0, 90)
ax.set_zlim(0, 90)
ax.set_ylim(0, 90)
plt.show()
