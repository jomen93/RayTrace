#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 23:58:58 2017
Shadow of a rotating black hole
@author: ashcat
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from PIL import Image
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
import os
import myconfig as cfg


class photonSphere():
    
    def __init__(self):
        # Parameters of the Metric
        self.M = cfg.M          # Mass of the central body in SolarMasses
        self.a = cfg.a          # Spin Parameter (0<=a<=1)
        self.inc = cfg.i      # Inclination angle in radians
        self.N = cfg.N

        # Radius of photon sphere for the corresponding Schwarzschild BH
        self.Rcircle = np.sqrt(27 * self.M**2)

    def psKerr(self):
        # Event Horizon
        self.rH = self.M + np.sqrt(self.M**2 - self.a**2) 
                            
        # Auxiliar Function
        def Delta(r):
                return r**2 - 2*self.M* r + self.a**2
            
        # Angular Momentum
        def xi(r):
            return (self.M*(r**2 - self.a**2) 
            - r*Delta(r))/(self.a*(r - self.M))  
        # Carter Constant
        def eta(r):
            return ((r**3)*(4*self.M*Delta(r) 
            - r*(r-self.M)**2))/(self.a**2 * (r - self.M)**2)   
            
        # Obtaining the Celestial Coordinates
        # Argument in the square root of beta
        def arg(r):
            return eta(r) + self.a**2 * np.cos(self.inc)**2 
            - (xi(r))**2 * (np.cos(self.inc)/np.sin(self.inc))**2
            
        # Root finding of the function arg(r) to define the limits in the plot        
        r0 = newton(arg, self.rH)       
        r1 = newton(arg, 4*self.rH)
                    
        # Range of the parameter r for the plot 
        self.r = np.linspace(r0+0.0000001 ,r1,100000)
            
        # Coordinate alpha
        self.alpha = - xi(self.r) / np.sin(self.inc)
            
        #coordinate beta
        self.beta = np.sqrt(arg(self.r))  


            
    def psPlot(self):
        # Plot of the orbit
        fig, ax = plt.subplots()

        # Background Color
        plt.style.use('dark_background')

        # Axis range and aspect (square)
        axrange = cfg.Ssize
        ax.axis([-axrange,axrange,-axrange,axrange])
        ax.set_aspect('equal')

        # No axis thicks or labels
        ax.set_axis_off()
        
        # Plotting

        if cfg.Metric == "Schwarzschild":
            # Plot of the Schwarzschild's BH shadow for comparison
            t = np.linspace(0, np.pi, 50000)
            ax.plot(self.Rcircle * np.cos(t), self.Rcircle * np.sin(t),'w')
        elif cfg.Metric == "Kerr":
            self.psKerr()
            # Plot the shadow of the rotating BH
            ax.plot(self.alpha,self.beta, 'w')

        #plt.show()
        plt.savefig('tempPlot.jpg', bbox_inches='tight', pad_inches=0)

        ax.clear()





    # Resize the Image (Pixelate it!)
    def pixelate(self, input_file_path='tempPlot.jpg', output_file_path='tempPixel.jpg'):
        image = Image.open(input_file_path)
        image = image.resize(
            (self.N, self.N),
            Image.NEAREST
        )
        image.save(output_file_path)


    def Fits(self):
        '''
        if os.path.isfile('complete.fits'):
            os.remove('complete.fits')
        '''
        # Convert the image into a FITS file
        plt.style.use(astropy_mpl_style)
        image2 = Image.open('tempPixel.jpg')
        xsize, ysize = image2.size

        r, g, b= image2.split()
        r_data = np.array(r.getdata()) # data is now an array of length ysize*xsize
        g_data = np.array(g.getdata())
        b_data = np.array(b.getdata())

        r_data = r_data.reshape(ysize, xsize)
        g_data = g_data.reshape(ysize, xsize)
        b_data = b_data.reshape(ysize, xsize)

        completeData = r_data * g_data * b_data

        if np.amax(completeData)!= 0.:
            completeData = completeData /np.amax(completeData)

        return completeData

        '''
        complete = fits.PrimaryHDU(data=completeData)
        complete.header['LATOBS'] = "32:11:56"
        complete.header['LONGOBS'] = "110:56"
        complete.writeto('complete.fits')
        '''
        #os.remove('tempPlot.jpg')
        #os.remove('tempPixel.jpg')


