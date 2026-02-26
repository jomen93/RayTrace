#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Infinite accretion disk model with exponential spectrum.

This module implements a simple toy model for an optically thick, 
geometrically thin accretion disk extending from ISCO to infinity.
The emission follows an exponentially decreasing radial profile.

Author: Johan Mendez
Copyright (c) 2024

Physics:
    The model assumes emission intensity I(r) proportional to:
        I(r) ~ exp(-(r - R_in) / scale)
    
    where R_in is the ISCO radius and scale is an arbitrary decay length
    (set to 5M in this implementation).

Note:
    This is a simplified model for visualization purposes. Real accretion
    disks require full radiative transfer calculations.
"""

import numpy as np


class Disk:
    """
    Infinite accretion disk with exponential radial spectrum.
    
    Attributes:
        rData: 2D array of photon radii at equatorial plane
        R_in: Inner disk radius (typically ISCO)
        R_out: Outer disk radius (cutoff for numerical purposes)
        image: Final 2D intensity map
        
    Example:
        >>> radii = raytracer_output  # From C++ module
        >>> disk = Disk(radii, R_in=6.0, R_out=20.0)
        >>> plt.imshow(disk.image, cmap='hot')
    """
    
    def __init__(self, rData: np.ndarray, R_in: float = 6.0, R_out: float = 20.0):
        """
        Initialize accretion disk model.
        
        Args:
            rData: 2D array where each element is the radius r at which
                   the corresponding photon intersected the equatorial plane.
                   Zero values indicate photons that didn't reach the disk.
            R_in: Inner disk edge (ISCO radius in geometric units)
            R_out: Outer disk edge (arbitrary cutoff)
        """
        self.rData = rData
        self.Shape = np.shape(self.rData)
        
        # Broadcast radii to full image size
        self.R_out = R_out * np.ones(self.Shape)
        self.R_in = R_in * np.ones(self.Shape)
        
        # Compute intensity profile: exponential decay from inner edge
        # I(r) ~ exp(-(r - R_in) / 5)
        exponent = (-self.rData + self.R_in) / 5.0
        
        # Clamp to prevent overflow in exp
        exponent = np.clip(exponent, -700, 700)
        self.intens = np.exp(exponent)
        
        # Mask regions outside disk: R_in <= r < infinity
        # Note: we use R_in <= rData (not <= R_out) for "infinite" disk
        self.image = self.intens * (self.R_in <= self.rData)
