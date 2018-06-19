#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Single Ray Tracer
Input    
(Alpha,Beta): initial coordinates in the image plane
freq: frequency of the photon (needed for the 4-momentum)
D: distance between the image plane and the force center
i: inclination angle of the image plane

 

@author: ashcat
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from photon import Photon
from initialConditions import initCond 
import metrics.schwarzschild as m       


################### INPUT ######################
# Initial coordinates of the photon in the observer's plane
Alpha = 0.
Beta =  0.

#Frequency of the observed photon
freq = 1.

# Location of the image plane
D = 50.
i = np.pi/4

# Final value of the affine parameter in the geodesic
tauf = -500.

# Number of steps in the integration
dtau = 5000


################### NUMERICAL SOLUTION ######################
# Define a photon class     
p = Photon(Alpha, Beta, freq, D, i)    

# Initial position and momentum of the photon
# p.xin
# p.kin

# Build the initial conditions needed to solve the geodesic equations
# [t, r, theta, phi, k_t, k_r, k_theta, k_phi] and stores in the variable
# p.iC of the Photon class
p.initConds(initCond(p.xin, p.kin, m.g(p.xin)))
 
# Solve the geodesic Equations

tau = np.linspace(0, tauf , dtau)

coords = odeint(m.geodesics, p.iC, tau)


################### PLOT OF THE GEODESIC ######################

# Cartesian coordinates for plotting
z = coords[:,1]*np.cos(coords[:,2])
x = coords[:,1]*np.sin(coords[:,2])*np.cos(coords[:,3])
y = coords[:,1]*np.sin(coords[:,2])*np.sin(coords[:,3])

# 3D-Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)

ax.set_title("Photon motion", loc='center')
#ax.set_xlim(-00, 90)
#ax.set_zlim(-00, 90)
#ax.set_ylim(-00, 90)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
