#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 21:00:07 2018
FITS writting module

This module writes the .fits file with the created image.

Class FITS has two procedures:
- Write: wirtes the .fits file
- showImage: shows the claculated image of the BH


Input for the class: 
name: Name of the output .fits file
hdrData: matrix data with the calculated image

@author: 
"""

import matplotlib.pyplot as plt
from astropy.io import fits 
import myconfig as cfg
from matplotlib import cm
import os

class FITS:
    def __init__(self,info,name, hdrData):
        self.info = info
        self.name = name
        self.hdrData = hdrData
        self.hdu = fits.PrimaryHDU(data = self.info)
        
    def Write(self):
        if os.path.isfile(self.name):
            os.remove(self.name)
        self.hdu.writeto(self.name)
        self.hdu.header['i'] = '1'

        # Include data of the image in the Header of the .fits file
        for j in self.hdrData.keys():
            fits.setval(self.name, j, value=self.hdrData[j])

    def showImage(self):
        plt.imshow(fits.open(self.name)[0].data,cmap = cm.afmhot)
        plt.xlabel("$x$ [pc]")
        plt.ylabel("$y$ [pc]")
        plt.savefig(str(cfg.N)+"x"+str(cfg.N))
        plt.colorbar()
        plt.show()


        
        
