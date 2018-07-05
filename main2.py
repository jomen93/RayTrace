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
from common.initialConditions import initCond   
from common.geodesicIntegration import geoInt  
from common.writeFits import FITS

  
# Load the Screen  
if cfg.ScreenType == 1:
    from screen.imagePlane import screen
    from screen.imagePlane import Photon
else:
    print("DEFINE A VALID SCREEN TYPE")
    exit(0)

from metrics import minkowski as m
from metrics import schwarzschild2 as m2


# Load the Accretion Structure    
if cfg.Structure == 1:
    from accretionStructures import simpleAccDisk as st
elif cfg.Structure == 2:
    from accretionStructures import linearAccDisk as st
else:
    print("DEFINE A VALID ACCRETION STRUCTURE")
    exit(0)


# Ranges of the Alpha and Beta coordinates in the image plane given by the
# screen.imagePlane module    
alphaRange, betaRange, numPixels = screen(cfg.Ssize, cfg.N)

# This array stores the data of the radius at the equatorial plane
# for the received photons
rDataArray = np.zeros((numPixels,numPixels))
rDataArray2 = np.zeros((numPixels,numPixels))


        




 
perc = 0
for j in range(0,numPixels):
    for k in range(0,numPixels): 
        # Define a photon       
        p = Photon(Alpha = alphaRange[k], Beta = betaRange[numPixels-1-j],
                   D = cfg.D, i = cfg.i)
        p2 = Photon(Alpha = alphaRange[k], Beta = betaRange[numPixels-1-j],
                   D = cfg.D, i = cfg.i)
        
        # Build the initial conditions needed to solve the geodesic equations
        # [t, r, theta, phi, k_t, k_r, k_theta, k_phi] and stores in the variable
        # p.iC of the Photon class
        
        p.initConds(initCond(p.xin, p.kin, m.g(p.xin)))
        p2.initConds(initCond(p2.xin, p2.kin, m2.g(p.xin)))


        finalPos, l = geoInt(m.geodesics, p.iC,  m.rEH(), intStep = 0.5)
        finalPos2, l2 = geoInt(m2.geodesics, p2.iC,  m.rEH(), intStep = 0.5)
        


        p.finalPosition(finalPos)
        p2.finalPosition(finalPos2)
        
        #print(m.rEH())
        #print(m2.rEH())
        
        print(p.fP)
        print(p2.fP)
        print()
        '''
        
        # Store value of the radius of the photon in the equatorial plane
        # in the rDataArray
        rDataArray[j, k] = p.fP[1]
        rDataArray2[j, k] = p2.fP[1]
        # Print progress
        perc = perc +1
        '''
        #print (" ", np.int(perc/(numPixels*numPixels)*100),"% complete ", end="\r")



'''
# Constructs the image data using the SimpleDisk structure 
disk = st.Disk(rDataArray, cfg.R_in, cfg.R_out)
diskImage = disk.image

disk2 = st.Disk(rDataArray2, cfg.R_in, cfg.R_out)
diskImage2 = disk2.image
'''


'''
# Stores the image in a name.fits file
name = str(cfg.N) + "x" + str(cfg.N) + ".fits"
imageData = FITS(diskImage, name)
imageData.Write()
imageData.showImage()
'''
'''
print(diskImage)
print()
print(diskImage2)
'''