#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Initial conditions preparation for geodesic integration.

Transforms contravariant momentum components to covariant components
using the metric at the photon's initial position.

Author: Johan Mendez
Copyright (c) 2024

Mathematical Formulation:
    Given contravariant momentum k^μ (components kt, kr, ktheta, kphi)
    and metric g_μν, computes covariant momentum:
    
        p_μ = g_μν k^ν
    
    This lowers the index for Hamiltonian formulation of geodesic equations.
"""

from typing import List


def initCond(x: List[float], k: List[float], g: List[float]) -> List[float]:
    """
    Prepare initial conditions for geodesic integration.
    
    Transforms position and momentum from image plane coordinates to
    the full phase space state vector required by the integrator.
    
    Args:
        x: Position [t, r, theta, phi] in spherical coordinates
        k: Contravariant momentum [kt, kr, ktheta, kphi]
        g: Metric components [gtt, gtphi, grr, gthth, gphph] at position x
        
    Returns:
        State vector [t, r, theta, phi, p_t, p_r, p_theta, p_phi]
        with covariant momentum components
        
    Example:
        >>> x = [0.0, 100.0, np.pi/3, 0.0]  # Position
        >>> k = [1.0, -0.5, 0.0, 0.02]      # Momentum (contravariant)
        >>> g = metric(x)                    # Metric at position
        >>> ic = initCond(x, k, g)
        >>> print(f"Initial state: {ic}")
    """
    # Unpack position
    t, r, theta, phi = x
    
    # Unpack contravariant momentum
    pt, pr, pth, pphi = k
    
    # Unpack metric components
    gtt, gtphi, grr, gthth, gphph = g
    
    # Lower indices: p_μ = g_μν k^ν
    p_t = gtt * pt + gtphi * pphi
    p_r = grr * pr
    p_th = gthth * pth
    p_phi = gphph * pphi
    
    return [t, r, theta, phi, p_t, p_r, p_th, p_phi]
