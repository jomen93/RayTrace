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
from photonSphere import ps 

gra=ps.photonSphere()
gra.psPlot()
gra.pixelate()
gra.Fits()