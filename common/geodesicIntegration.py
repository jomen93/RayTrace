#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geodesic integration module using scipy.integrate.odeint.

This is the legacy Python implementation kept for reference and fallback.
The production code now uses the C++ implementation for performance.

Author: Johan Mendez
Copyright (c) 2024

Integration Method:
    Uses scipy.integrate.odeint (LSODA) for adaptive integration.
    Integrates backwards in affine parameter until photon either:
    1. Crosses the equatorial plane (theta = pi/2)
    2. Falls into the event horizon
    3. Exceeds maximum step count
    
Note:
    This implementation is kept for reference. For production use,
    the C++ module (raytracer_cpp) provides ~20x better performance.
"""

import numpy as np
from scipy.integrate import odeint


def geoInt(geoEq, initConds, rEH: float = 0.0, intStep: float = 0.05):
    """
    Integrate geodesic equations to equatorial plane.
    
    Legacy Python implementation using scipy.integrate.odeint.
    
    Args:
        geoEq: Geodesic equations function dx/dtau = f(x)
        initConds: Initial state [t, r, theta, phi, p_t, p_r, p_theta, p_phi]
        rEH: Event horizon radius (terminates if crossed)
        intStep: Integration step size
        
    Returns:
        Tuple of (final_state, num_steps) where final_state contains
        the photon state at equatorial crossing, or zeros if failed.
    """
    geodesics = geoEq
    iC = initConds
    
    # Conserved quantities (for reference)
    E = iC[4]   # Energy
    L = iC[7]   # Angular momentum
    
    # Carter's constant (approximately conserved)
    Carter = iC[6]**2 + ((np.cos(iC[2])**2) / (np.sin(iC[2])**2)) * L**2
    
    tolerance = 5.0e-5
    jmax = 100000
    
    # Initial integration step
    coords = odeint(geodesics, iC, [0, -intStep])
    j = 1
    
    # Main integration loop
    while np.cos(coords[j, 2]) > tolerance:
        coordloop = odeint(geodesics, coords[j], [0, -intStep])
        
        # Check horizon crossing
        if coordloop[1, 1] <= rEH + 1e-6:
            return np.zeros(8), j
        
        # Append new point
        coords = np.concatenate((coords, [coordloop[1]]))
        j = j + 1
        
        # Check max steps
        if j > jmax:
            return np.zeros(8), j
    
    return coords[j-1, :], j
