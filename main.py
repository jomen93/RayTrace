#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Black Hole Ray Tracer - Main Entry Point

A high-performance ray tracing simulation for generating black hole shadow images
using a hybrid C++/Python architecture. The computationally intensive geodesic
integration is performed in optimized C++ code, while Python handles configuration,
visualization, and data I/O.

Author: Johan Mendez
Copyright (c) 2024

Example usage:
    $ ./build_cpp.sh
    $ python main.py

Configuration is read from myconfig.py. Output images are saved to output/.
"""

import os
import sys
from pathlib import Path

import numpy as np

# Import C++ ray tracing module
try:
    import raytracer_cpp
except ImportError as e:
    print(f"Error: C++ module not found: {e}")
    print("Build instructions:")
    print("  1. Ensure pybind11 is installed: pip install pybind11")
    print("  2. Run: ./build_cpp.sh")
    sys.exit(1)

# Project modules
import myconfig as cfg
from accretionStructures.inftyAccDisk import Disk as AccretionDisk
from common.writeFits import FITS
from photonSphere.ps import photonSphere


def ensure_output_directory() -> None:
    """Create output directory if it doesn't exist."""
    Path("output").mkdir(parents=True, exist_ok=True)


def trace_photons_cpp(alpha_range: np.ndarray, beta_range: np.ndarray,
                      distance: float, inclination: float,
                      mass: float, spin: float, num_pixels: int,
                      metric: str) -> np.ndarray:
    """
    Trace photons through the black hole metric using C++ backend.
    
    Args:
        alpha_range: Array of horizontal image plane coordinates
        beta_range: Array of vertical image plane coordinates  
        distance: Distance from observer to black hole (geometric units)
        inclination: Viewing angle in radians (pi/2 = edge-on)
        mass: Black hole mass (geometric units, typically 1.0)
        spin: Kerr spin parameter (0 for Schwarzschild, < 1 for Kerr)
        num_pixels: Image resolution (NxN)
        metric: "Schwarzschild" or "Kerr"
        
    Returns:
        2D array of radii where photons intersect the equatorial plane.
        Zero indicates the photon fell into the black hole or escaped.
    """
    return raytracer_cpp.trace_photons(
        alpha_range.astype(np.float64),
        beta_range.astype(np.float64),
        float(distance),
        float(inclination),
        float(mass),
        float(spin),
        int(num_pixels),
        metric
    )


def compute_isco_radius(mass: float, spin: float, 
                        corotating: bool, metric: str) -> float:
    """
    Compute the Innermost Stable Circular Orbit (ISCO) radius.
    
    For Schwarzschild: ISCO = 6M
    For Kerr: Depends on spin and corotation direction.
    
    Args:
        mass: Black hole mass
        spin: Kerr spin parameter  
        corotating: True for prograde, False for retrograde
        metric: "Schwarzschild" or "Kerr"
        
    Returns:
        ISCO radius in geometric units
    """
    if metric == "Kerr":
        return raytracer_cpp.kerr_r_isco(mass, spin, corotating)
    return raytracer_cpp.schwarzschild_r_isco(mass)


def generate_visualization(radii_data: np.ndarray, isco_radius: float,
                          config, output_name: str) -> None:
    """
    Generate the final black hole image from ray tracing data.
    
    This combines:
    1. Accretion disk emission based on photon radii
    2. Black hole shadow from photon sphere
    3. Proper orientation and intensity mapping
    
    Args:
        radii_data: 2D array of photon radii from C++ module
        isco_radius: Inner edge of accretion disk
        config: Configuration module with simulation parameters
        output_name: Base filename for output images
    """
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    num_pixels = radii_data.shape[0]
    
    # Generate accretion disk image
    disk = AccretionDisk(radii_data, isco_radius, config.R_out)
    
    # Generate photon sphere (black hole shadow)
    sphere = photonSphere(N=num_pixels)
    sphere.psPlot()
    sphere.pixelate()
    
    # Combine disk and shadow
    final_image = disk.image + sphere.Fits()
    
    # Save FITS format for scientific analysis
    fits_path = Path("output") / f"{output_name}.fits"
    header = {
        'AUTHOR': 'Johan Mendez',
        'METRIC': config.Metric,
        'MASS': str(config.M),
        'N_PIXELS': str(num_pixels),
        'INCLINAT': str(config.i),
        'METHOD': 'C++'
    }
    if config.Metric == "Kerr":
        header['ANG_MOM'] = str(config.a)
    
    FITS(final_image, str(fits_path), header).Write()
    
    # Generate PNG visualization
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Load and flip for correct orientation
    data = np.flipud(fits.open(fits_path)[0].data)
    
    im = ax.imshow(data, cmap='afmhot', origin='lower')
    ax.set_title(f"{config.Metric} Black Hole Shadow\n"
                 f"{num_pixels}×{num_pixels} pixels", fontsize=14)
    ax.set_xlabel("x [M]", fontsize=12)
    ax.set_ylabel("y [M]", fontsize=12)
    
    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(im, cax=cax, label='Normalized Intensity')
    
    png_path = Path("output") / f"{output_name}.png"
    plt.savefig(png_path, dpi=200, bbox_inches='tight')
    plt.close()
    
    print(f"  Output: {png_path}")


def main() -> int:
    """
    Main execution entry point.
    
    Performs black hole ray tracing simulation and generates output images.
    
    Returns:
        Exit code (0 for success, 1 for error)
    """
    print("=" * 60)
    print("BLACK HOLE RAY TRACER")
    print("=" * 60)
    
    ensure_output_directory()
    
    # Extract configuration parameters
    mass = cfg.M
    spin = getattr(cfg, 'a', 0.0)
    
    # Ensure even pixel count for numerical stability
    num_pixels = cfg.N if cfg.N % 2 == 0 else cfg.N + 1
    
    # Setup image plane coordinates
    alpha_range = np.linspace(-cfg.Ssize, cfg.Ssize, num_pixels)
    beta_range = np.linspace(-cfg.Ssize, cfg.Ssize, num_pixels)
    total_photons = num_pixels * num_pixels
    
    # Print configuration summary
    print(f"\nConfiguration:")
    print(f"  Metric: {cfg.Metric}")
    print(f"  Resolution: {num_pixels}×{num_pixels}")
    print(f"  Inclination: {np.degrees(cfg.i):.1f}°")
    print(f"  Total photons: {total_photons:,}")
    
    # Perform ray tracing with C++ backend
    print("\nTracing photons...")
    import time
    start_time = time.time()
    
    radii_data = trace_photons_cpp(
        alpha_range, beta_range,
        cfg.D, cfg.i, mass, spin, num_pixels, cfg.Metric
    )
    
    elapsed = time.time() - start_time
    photon_rate = total_photons / elapsed
    
    print(f"✓ Completed in {elapsed:.2f}s ({photon_rate:,.0f} photons/s)")
    
    # Diagnostic statistics
    nonzero_count = np.count_nonzero(radii_data)
    efficiency = 100 * nonzero_count / total_photons
    print(f"\n  Min radius: {np.min(radii_data):.2f} M")
    print(f"  Max radius: {np.max(radii_data):.2f} M")
    print(f"  Non-zero: {nonzero_count}/{total_photons} ({efficiency:.1f}%)")
    
    # Compute ISCO for disk inner edge
    isco_radius = compute_isco_radius(mass, spin, cfg.corotation, cfg.Metric)
    print(f"  ISCO radius: {isco_radius:.2f} M")
    
    # Generate final visualization
    print("\nGenerating image...")
    generate_visualization(radii_data, isco_radius, cfg, cfg.fileName)
    
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
