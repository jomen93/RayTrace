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
from photonSphere import ps 
from common.writeFits import FITS
  
# Load the Screen  
if cfg.ScreenType == 1:
    from screen.imagePlane import screen
    from screen.imagePlane import Photon
    hdrData = {'SCREEN':'Image Plane'}
    hdrData['DISTANCE'] = str(cfg.D)
    hdrData['INCLINAT'] = str(cfg.i)
    screenSize = cfg.Ssize
elif cfg.ScreenType == 2:
    from screen.pointCamera import screen
    from screen.pointCamera import Photon
    hdrData = {'SCREEN':'Point Camera'}
    hdrData['DISTANCE'] = str(cfg.D)
    hdrData['INCLINAT'] = str(cfg.i)
    screenSize = cfg.FoV
else:
    print("DEFINE A VALID SCREEN TYPE")
    exit(0)

# Load the Metric    
if cfg.Metric == "Minkowski":
    from metrics import minkowski as m
    hdrData['METRIC'] = 'Minkowksi'
elif cfg.Metric == "Schwarzschild":
    from metrics import schwarzschild as m
    hdrData['METRIC'] = 'Schwarzschild'
    hdrData['MASS'] = str(cfg.M)
elif cfg.Metric == "Kerr":
    from metrics import kerr as m
    hdrData['METRIC'] = 'Kerr'
    hdrData['MASS'] = str(cfg.M)
    hdrData['ANG_MOM'] = str(cfg.a)
else:
    print("DEFINE A VALID METRIC")
    exit(0)

# Load the Accretion Structure    
if cfg.Structure == 1:
    from accretionStructures import simpleAccDisk as st
    hdrData['ACSTRCTR'] = 'Simplest Accretion Disk'
elif cfg.Structure == 2:
    from accretionStructures import linearAccDisk as st
    hdrData['ACSTRCTR'] = 'Linear Spectrum Accretion Disk'
elif cfg.Structure == 3:
    from accretionStructures import inftyAccDisk as st
    hdrData['ACSTRCTR'] = 'Infinite Accretion Disk with a decreasing exponential spectrum'
else:
    print("DEFINE A VALID ACCRETION STRUCTURE")
    exit(0)


# Ranges of the Alpha and Beta coordinates in the image plane given by the
# screen.imagePlane module    
alphaRange, betaRange, numPixels = screen(screenSize, cfg.N)

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
        #rEH = m.rEH()
               
        finalPos, l = geoInt(m.geodesics, p.iC,  m.rEH(), intStep = 0.5)
        
        p.finalPosition(finalPos)

        # Store value of the radius of the photon in the equatorial plane
        # in the rDataArray
        rDataArray[j, k] = p.fP[1]
        
        # Print progress
        perc = perc +1
        
        print (" ", np.int(perc/(numPixels*numPixels)*100))


# Constructs the image data using the accretion structure 
if cfg.corotation == True:
    R_in = m.ISCOco()
else:
    R_in = m.ISCOcounter()

disk = st.Disk(rDataArray, R_in, cfg.R_out)
diskImage = disk.image

photonsphere= ps.photonSphere()
photonsphere.psPlot()
photonsphere.pixelate()
psData = photonsphere.Fits()

# Composed image with the photon sphere
completeImage = diskImage + psData

# Stores the image in a name.fits file
name = cfg.fileName + ".fits"
imageData = FITS(completeImage, name, hdrData)
imageData.Write()
imageData.showImage()
#####
