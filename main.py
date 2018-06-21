#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Main Module

  

@author: ashcat
"""
from sys import exit
import numpy as np
import myconfig as cfg
from initialConditions import initCond   
from geodesicIntegration import geoInt  
from writeFits import FITS
from accretionStructures import simpleAccDisk as st

  
# Load the Screen  
if cfg.ScreenType == 1:
    from screen.imagePlane import screen
    from photon import Photon
else:
    print("DEFINE A VALID SCREEN TYPE")
    exit(0)

# Load the Metric    
if cfg.Metric == 1:
    from metrics import minkowski as m
elif cfg.Metric == 2:
    from metrics import schwarzschild as m
else:
    print("DEFINE A VALID METRIC")
    exit(0)


# Ranges of the Alpha and Beta coordinates in the image plane given by the
# screen.imagePlane module    
alphaRange, betaRange, numPixels = screen(cfg.Ssize, cfg.N)

# This array stores the data of the radius at the equatorial plane
# for the received photons
rDataArray = np.zeros((numPixels,numPixels))


    
perc = 0
for j in range(0,numPixels):
    for k in range(0,numPixels): 
        # Define a photon       
        p = Photon(Alpha = alphaRange[k], Beta = betaRange[numPixels-1-j],
                   D = cfg.D, i = cfg.i)
        
        # Build the initial conditions needed to solve the geodesic equations
        # [t, r, theta, phi, k_t, k_r, k_theta, k_phi] and stores in the variable
        # p.iC of the Photon class
        
        p.initConds(initCond(p.xin, p.kin, m.g(p.xin)))
               
        finalPos, l = geoInt(m.geodesics, p.iC, intStep = 0.5)
        
        p.finalPosition(finalPos)

        #print(" Number of Steps to reach the Equatorial Plane:",l)
        #print()
        #print("Alpha=", k, " Beta=", numPixels-1-j)
        
        # Store value of the radius of the photon in the equatorial plane
        # in the rDataArray
        rDataArray[j, k] = p.fP[1]
        
        # Print progress
        perc = perc +1
        print (np.int(perc/(numPixels*numPixels)*100),"% complete ", end="\r")


# Constructs the image using the SimpleDisk structure 
disk = st.SimpleDisk(rDataArray)
diskImage = disk.image

print(diskImage)

# Stores the image in a name.fits file
name = "exampleHD.fits"
imageData = FITS(diskImage, name)
imageData.Write()
imageData.showImage()

