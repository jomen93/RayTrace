#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Photon sphere and black hole shadow visualization.

Computes and renders the shadow (silhouette) of a black hole as seen by
a distant observer. For Schwarzschild, this is a perfect circle of radius
sqrt(27)M ≈ 5.196M. For Kerr, the shadow is distorted due to frame dragging.

Author: Johan Mendez
Copyright (c) 2024

References:
    - Bardeen (1973): "Timelike and null geodesics in the Kerr metric"
    - Gralla et al. (2019): "The Kerr null geodesic
"""

import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.optimize import newton

import myconfig as cfg


# Output directory for temporary files
TEMP_DIR = Path('output')
TEMP_DIR.mkdir(parents=True, exist_ok=True)


class photonSphere:
    """
    Computes and renders the black hole photon sphere shadow.
    
    For Schwarzschild metric: Perfectly circular shadow of radius sqrt(27)M.
    For Kerr metric: Asymmetric shadow computed from critical null geodesics.
    
    The shadow represents the boundary between photons that escape to infinity
    and those that are captured by the black hole.
    
    Attributes:
        M: Black hole mass
        a: Kerr spin parameter (0 for Schwarzschild)
        inc: Inclination angle
        N: Image resolution
        Rcircle: Shadow radius for Schwarzschild case
        
    Example:
        >>> ps = photonSphere(N=512)
        >>> ps.psPlot()      # Generate shadow plot
        >>> ps.pixelate()    # Resize to pixel resolution
        >>> shadow_data = ps.Fits()  # Get as numpy array
    """
    
    def __init__(self, N: int = None):
        """
        Initialize photon sphere calculator.
        
        Args:
            N: Image resolution (uses cfg.N if None, rounded to even)
        """
        self.M = cfg.M
        self.a = cfg.a
        self.inc = cfg.i
        
        # Use provided N or config value, ensure even
        if N is None:
            N = cfg.N
        self.N = N if N % 2 == 0 else N + 1
        
        # Schwarzschild shadow radius: sqrt(27) * M
        self.Rcircle = np.sqrt(27 * self.M**2)
    
    def psKerr(self) -> None:
        """
        Compute Kerr shadow boundary using critical null geodesics.
        
        Calculates the celestial coordinates (alpha, beta) of the shadow
        edge by finding the unstable photon orbits at the critical impact
        parameter.
        
        Reference: Bardeen (1973), Gralla et al. (2019)
        """
        # Event horizon radius for Kerr
        self.rH = self.M + np.sqrt(self.M**2 - self.a**2)
        
        # Metric functions
        def Delta(r: float) -> float:
            return r**2 - 2 * self.M * r + self.a**2
        
        # Carter constants for critical orbits
        def xi(r: float) -> float:
            """Angular momentum of unstable photon orbit."""
            return (self.M * (r**2 - self.a**2) - r * Delta(r)) / (self.a * (r - self.M))
        
        def eta(r: float) -> float:
            """Carter constant for unstable photon orbit."""
            return ((r**3) * (4 * self.M * Delta(r) - r * (r - self.M)**2) / 
                   (self.a**2 * (r - self.M)**2))
        
        # Argument of square root for beta coordinate
        def arg(r: float) -> float:
            return (eta(r) + self.a**2 * np.cos(self.inc)**2 - 
                   (xi(r))**2 * (np.cos(self.inc) / np.sin(self.inc))**2)
        
        # Find range of radii that contribute to shadow
        r0 = newton(arg, self.rH)
        r1 = newton(arg, 4 * self.rH)
        
        # Sample critical orbits
        self.r = np.linspace(r0 + 1e-7, r1, 100000)
        
        # Celestial coordinates of shadow edge
        self.alpha = -xi(self.r) / np.sin(self.inc)
        self.beta = np.sqrt(arg(self.r))
    
    def psPlot(self) -> None:
        """
        Render the black hole shadow to an image file.
        
        Creates a plot with black shadow on transparent/white background.
        For Schwarzschild: perfect circle.
        For Kerr: asymmetric shadow from psKerr().
        
        Output saved to output/tempPlot.jpg
        """
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Black background for shadow
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')
        
        # Set axis range
        axrange = cfg.Ssize
        ax.axis([-axrange, axrange, -axrange, axrange])
        ax.set_aspect('equal')
        ax.set_axis_off()
        
        # Render shadow
        if cfg.Metric == "Schwarzschild":
            # Perfect circle for Schwarzschild
            theta = np.linspace(0, 2 * np.pi, 1000)
            ax.fill(self.Rcircle * np.cos(theta), 
                   self.Rcircle * np.sin(theta), 
                   color='black')
            ax.plot(self.Rcircle * np.cos(theta), 
                   self.Rcircle * np.sin(theta), 
                   'k-', linewidth=1)
        elif cfg.Metric == "Kerr":
            # Distorted shadow for Kerr
            self.psKerr()
            ax.fill(self.alpha, self.beta, color='black')
            ax.plot(self.alpha, self.beta, 'k-', linewidth=1)
        
        # Save plot
        output_path = TEMP_DIR / 'tempPlot.jpg'
        plt.savefig(output_path, facecolor='black', edgecolor='none',
                   bbox_inches='tight', pad_inches=0.1)
        plt.close()
    
    def pixelate(self, input_file_path: str = None, 
                 output_file_path: str = None) -> None:
        """
        Resize shadow image to match simulation resolution.
        
        Args:
            input_file_path: Source image (default: tempPlot.jpg)
            output_file_path: Destination (default: tempPixel.jpg)
        """
        if input_file_path is None:
            input_file_path = TEMP_DIR / 'tempPlot.jpg'
        if output_file_path is None:
            output_file_path = TEMP_DIR / 'tempPixel.jpg'
        
        image = Image.open(input_file_path)
        image = image.resize((self.N, self.N), Image.NEAREST)
        image.save(output_file_path)
    
    def Fits(self) -> np.ndarray:
        """
        Convert pixelated shadow to numpy array for composition.
        
        Returns:
            2D array with values in [0, 1] where:
            - 1.0 = Shadow (black hole silhouette)
            - 0.0 = Transparent region (outside shadow)
        """
        image = Image.open(TEMP_DIR / 'tempPixel.jpg')
        
        # Split into RGB channels
        r, g, b = image.split()
        
        # Convert to arrays
        r_data = np.array(r.getdata()).reshape(image.size[1], image.size[0])
        g_data = np.array(g.getdata()).reshape(image.size[1], image.size[0])
        b_data = np.array(b.getdata()).reshape(image.size[1], image.size[0])
        
        # Combine channels (shadow = white in original, becomes 1.0)
        complete_data = r_data * g_data * b_data
        
        # Normalize to [0, 1]
        if np.amax(complete_data) != 0.0:
            complete_data = complete_data / np.amax(complete_data)
        
        return complete_data
