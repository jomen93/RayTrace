#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Configuration file for Black Hole Ray Tracer.

This module defines all simulation parameters for the ray tracing computation.
Values are specified in geometric units where G = c = 1, with length measured
in units of the black hole mass (typically M = 1).

Author: Johan Mendez
Copyright (c) 2024

Example:
    To run a high-resolution Kerr simulation:
    >>> import myconfig as cfg
    >>> print(f"Simulating {cfg.Metric} black hole at {cfg.N}x{cfg.N} resolution")

Note:
    For accurate results, ensure the screen size (Ssize) is large enough to
    capture the photon ring, and the observer distance (D) is sufficiently
    large (D >> M) for the distant observer approximation.
"""

import numpy as np


# =============================================================================
# BLACK HOLE PARAMETERS
# =============================================================================

METRIC: str = "Schwarzschild"
"""Spacetime metric to use. Options: 'Schwarzschild' or 'Kerr'."""

MASS: float = 1.0
"""Black hole mass in geometric units (G = c = 1). Typically set to 1.0."""

SPIN: float = 0.7 * MASS
"""
Kerr spin parameter (dimensionless angular momentum).
Must satisfy 0 <= |a| < M for physical black holes.
For Schwarzschild metric, this value is ignored.
"""


# =============================================================================
# ACCRETION DISK PARAMETERS  
# =============================================================================

STRUCTURE_TYPE: int = 3
"""
Accretion disk emission model.
3 = Infinite disk with exponentially decreasing spectrum
"""

OUTER_RADIUS: float = 20.0 * MASS
"""Outer edge of accretion disk in units of M."""

COROTATING: bool = True
"""
For Kerr metric: True for prograde (corotating) disk, 
False for retrograde (counter-rotating) disk.
"""


# =============================================================================
# OBSERVER/SCREEN PARAMETERS
# =============================================================================

SCREEN_TYPE: int = 1
"""Screen geometry: 1 = Image Plane (distant observer)."""

SCREEN_SIZE: float = 15.0 * MASS
"""
Half-width of the image plane in units of M.
Must be > photon ring radius (~5.2M for Schwarzschild).
"""

RESOLUTION: int = 256
"""
Image resolution in pixels (N x N square image).
Recommended values:
- 256: Fast testing (~45 seconds)
- 512: Standard quality (~3 minutes)  
- 1024: High quality (~12 minutes)
"""

OBSERVER_DISTANCE: float = 100.0 * MASS
"""
Distance from observer to black hole center.
Must satisfy D >> M for the distant observer approximation.
"""

INCLINATION_ANGLE: float = np.pi / 2.1
"""
Viewing angle in radians.
- pi/2 = 90°: Edge-on (equatorial view)
- pi/4 = 45°: Intermediate inclination
- 0: Face-on (pole-on view)
Current value: ~85.7°, nearly edge-on to show disk structure.
"""


# =============================================================================
# OUTPUT PARAMETERS
# =============================================================================

OUTPUT_FILENAME: str = f"{METRIC}_{RESOLUTION}x{RESOLUTION}"
"""Base name for output files (FITS and PNG)."""


# =============================================================================
# BACKWARD COMPATIBILITY ALIASES
# These aliases allow existing code to reference configuration parameters
# using the original shorter names.
# =============================================================================

M = MASS
a = SPIN
N = RESOLUTION
Ssize = SCREEN_SIZE
D = OBSERVER_DISTANCE
i = INCLINATION_ANGLE
Structure = STRUCTURE_TYPE
R_out = OUTER_RADIUS
corotation = COROTATING
ScreenType = SCREEN_TYPE
fileName = OUTPUT_FILENAME
