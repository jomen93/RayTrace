#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Image plane and photon initialization module.

Defines the observer's image plane and computes initial conditions for
photons in the black hole ray tracing simulation.

Author: Johan Mendez
Copyright (c) 2024

Coordinate Systems:
    Image Plane (Cartesian-like):
        alpha: Horizontal coordinate (x-direction in image plane)
        beta: Vertical coordinate (y-direction in image plane)
        D: Distance from observer to origin (along line of sight)
        
    Spherical (Boyer-Lindquist):
        r: Radial coordinate
        theta: Polar angle from rotation axis
        phi: Azimuthal angle
        
Transformation:
    Photons are traced backwards from the image plane to the black hole,
    hence the momentum vector is oriented towards the black hole.
"""

import numpy as np
from numpy import ndarray
from typing import Tuple


def screen(screen_size: float, N: int) -> Tuple[ndarray, ndarray, int]:
    """
    Define a square N×N pixel screen in the image plane.
    
    Creates coordinate grids for ray tracing. If N is odd, it is rounded
    up to the next even number to ensure symmetric pixelation.
    
    Args:
        screen_size: Half-width of the screen in geometric units (M)
        N: Requested resolution (will be made even if odd)
        
    Returns:
        Tuple of (alpha_range, beta_range, num_pixels) where:
        - alpha_range: Array of horizontal coordinates [-Ssize, Ssize]
        - beta_range: Array of vertical coordinates [-Ssize, Ssize]
        - num_pixels: Actual resolution (N or N+1 if N was odd)
        
    Example:
        >>> alpha, beta, Npix = screen(screen_size=15.0, N=256)
        >>> print(f"Screen: {Npix}x{Npix} pixels")
    """
    # Ensure even number of pixels for symmetry
    if N & 1:
        num_pixels = N + 1
    else:
        num_pixels = N
    
    # Create coordinate grids
    alpha_range = np.linspace(-screen_size, screen_size, num_pixels)
    beta_range = np.linspace(-screen_size, screen_size, num_pixels)
    
    print(f"Screen size: {num_pixels}x{num_pixels} pixels")
    print(f"Total pixels: {num_pixels * num_pixels}")
    
    return alpha_range, beta_range, num_pixels


class Photon:
    """
    Represents a photon in the observer's image plane.
    
    Computes the transformation from image plane coordinates (alpha, beta)
    to spherical coordinates and initializes the 4-momentum for backward
    ray tracing.
    
    Attributes:
        Alpha, Beta: Image plane coordinates
        D: Distance to black hole
        i: Inclination angle (radians)
        r, theta, phi: Spherical position
        kt, kr, ktheta, kphi: Contravariant 4-momentum components
        kin: Initial momentum array [kt, kr, ktheta, kphi]
        
    Example:
        >>> p = Photon(Alpha=1.0, Beta=0.0, D=100.0, i=np.pi/4)
        >>> print(f"Initial position: r={p.r:.2f}, theta={p.theta:.2f}")
    """
    
    def __init__(self, Alpha: float = 1.0, Beta: float = 0.0, 
                 freq: float = 1.0, D: float = 100.0, 
                 i: float = np.pi / 4):
        """
        Initialize photon from image plane coordinates.
        
        Args:
            Alpha: Horizontal image coordinate (M)
            Beta: Vertical image coordinate (M)
            freq: Photon frequency (normalized to 1)
            D: Distance from observer to black hole (M)
            i: Inclination angle (radians, pi/2 = edge-on)
        """
        # Store image plane coordinates
        self.Alpha = Alpha
        self.Beta = Beta
        self.D = D
        self.i = i
        
        # Transform to spherical coordinates
        self.r = np.sqrt(self.Alpha**2 + self.Beta**2 + self.D**2)
        self.theta = np.arccos(
            (self.Beta * np.sin(self.i) + self.D * np.cos(self.i)) / self.r
        )
        self.phi = np.arctan2(
            self.Alpha,
            self.D * np.sin(self.i) - self.Beta * np.cos(self.i)
        )
        
        # Initial position 4-vector [t, r, theta, phi]
        self.xin = [0.0, self.r, self.theta, self.phi]
        
        # Compute initial 4-momentum (contravariant components)
        self.K0 = freq
        self.kr = (self.D / self.r) * self.K0
        
        # Auxiliary quantity for momentum calculation
        aux = (self.Alpha**2 + 
               (-self.Beta * np.cos(self.i) + self.D * np.sin(self.i))**2)
        
        self.ktheta = (self.K0 / np.sqrt(aux)) * (
            -np.cos(self.i) +
            (self.Beta * np.sin(self.i) + self.D * np.cos(self.i)) * 
            (self.D / (self.r**2))
        )
        
        self.kphi = -self.Alpha * np.sin(self.i) * self.K0 / aux
        
        self.kt = np.sqrt(
            self.kr**2 + 
            self.r**2 * self.ktheta**2 +
            self.r**2 * (np.sin(self.theta))**2 * self.kphi**2
        )
        
        # Momentum 4-vector [kt, kr, ktheta, kphi]
        self.kin = [self.kt, self.kr, self.ktheta, self.kphi]
    
    def initConds(self, initialConds: list) -> None:
        """
        Store initial conditions for geodesic integration.
        
        Args:
            initialConds: List [t, r, theta, phi, p_t, p_r, p_theta, p_phi]
                         containing the full phase space initial state
        """
        self.iC = initialConds
    
    def finalPosition(self, finalPos: list) -> None:
        """
        Store final position after geodesic integration.
        
        Args:
            finalPos: List [t, r, theta, phi, ...] at equatorial crossing
        """
        self.fP = finalPos
