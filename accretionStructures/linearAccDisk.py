#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 23:16:05 2018

@author: ashcat
"""

import numpy as np 

class Disk:
    def __init__(self, rData, R_in = 3., R_out = 5.):
        self.rData = rData
        self.Shape = np.shape(self.rData)
        self.R_out = R_out*np.ones(self.Shape)
        self.R_in = R_in*np.ones(self.Shape)
        
        self.m = (1.-0.)/(self.R_in - self.R_out)
        self.intens = self.m * (self.rData - self.R_out)
        
        self.image = self.intens*(self.R_in <= self.rData) * (self.rData <= self.R_out)

