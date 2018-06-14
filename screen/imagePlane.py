#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
Image Plane definition Module

This module defines the image plane in terms of pixels


@author: ashcat
"""
import numpy as np

def screen(Ssize, N):
    '''
    Defines a square NxN pixels screen.
    Returns the ranges of the cordinates Alpha and Beta in the image plane and 
    the value of the maximum value of the coordinates
    '''
    
    if N & 1:
        numPixels = N +1
    else:
        numPixels = N 

    alphaRange = np.linspace(-Ssize, Ssize, numPixels)
    betaRange = np.linspace(-Ssize, Ssize, numPixels)
    
    print ("Size of the screen in Pixels:", numPixels, "X", numPixels)
    print ("Number of Pixels: ", (numPixels)*(numPixels))
    
    return alphaRange, betaRange, numPixels

