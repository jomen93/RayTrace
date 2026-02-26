#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration File - OPTIMIZED EXAMPLE

This is an example configuration file showing all available options
for the optimized parallel version of the ray tracer.

Copy this to myconfig.py and modify as needed.
"""

import numpy as np
from astropy import constants as const
from astropy import units as u


############## Metric Definition ################

Metric = "Kerr"  # Options: "Minkowski", "Schwarzschild", "Kerr"

# Parameters in the metric
Mb = 1.5  # Mass parameter
Gc = const.G*(const.kpc)**-3*const.M_sun*(3.1e-7)**2*u.s**2*(u.kpc**3/(u.M_sun*u.yr**2))
M = float(2*Mb*(Gc)*(const.c.to(u.kpc/u.yr))**-2*(u.M_sun/u.kpc))
a = 0.7   # Spin parameter (0 <= a < M) - only for Kerr

#################################################


######### Accretion Structure Definition ########

Structure = 3

# 1: Simple Accretion Disk
# 2: Linear Spectrum Accretion Disk  
# 3: Infinite Accretion Disk with a decreasing exponential spectrum

# Parameters in the Structure
R_out = 12*M
corotation = True

#################################################


############### Screen Definition ###############

ScreenType = 1  # 1: Image Plane, 2: Point camera

# Image Plane Options
Ssize = 20*M    # Screen size (in units of M)

# Point Camera Options
FoV = np.pi/110.  # Field of View (in radians)

# General Screen Options
N = 64  # Resolution: N x N pixels (use 32 for tests, 256 for high quality)

# Distance to the Black hole (in units of M)
D = 10*M

# Inclination of the image plane (radians)
# np.pi/2 = edge-on, np.pi/4 = 45 degrees
i = np.pi/2.5

#################################################


################## File Output ################## 

fileName = Metric + "_" + str(N) + "x" + str(N) 

#################################################


################## Parallel Options ################## 

# Number of parallel workers for ray tracing
# -1: Use all available CPU cores (recommended for production)
#  1: Serial execution (use for small tests)
#  n: Use n specific cores
n_jobs = -1

# Integration method for solve_ivp
# 'RK45': 4th order Runge-Kutta, good balance (DEFAULT)
# 'RK23': 2nd/3rd order, faster but less accurate
# 'DOP853': 8th order, high precision but slower
# 'Radau': Implicit, good for stiff problems
# 'BDF': Backward differentiation, for stiff problems
# 'LSODA': Automatic stiffness detection
integration_method = 'RK45'

# Integration step size (initial guess)
# Smaller = more accurate but slower
intStep = 0.5

#################################################


################## Advanced Options ################## 

# Use Numba-accelerated geodesics (if available)
# Requires: pip install numba
use_numba = False  # Experimental feature

# Batch size for parallel processing
# Larger = less overhead but more memory
# None = automatic (recommended)
batch_size = None

# Show detailed progress bar
# Requires: pip install tqdm
show_progress = True

#################################################
