#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FITS file I/O module for black hole ray tracing.

Provides functionality to write astronomical FITS format images and
render preview visualizations using matplotlib.

Author: Johan Mendez
Copyright (c) 2024

Example:
    >>> from common.writeFits import FITS
    >>> image_data = compute_black_hole_image(...)  # Your computation
    >>> fits_writer = FITS(image_data, "black_hole.fits", header_info)
    >>> fits_writer.Write()
    >>> fits_writer.showImage()
"""

import os
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib import cm

# Default output directory
OUTPUT_DIR = Path('output')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


class FITS:
    """
    Handler for FITS file operations and image visualization.
    
    This class encapsulates writing astronomical FITS format files and
    generating matplotlib visualizations of black hole shadow images.
    
    Attributes:
        info: 2D numpy array containing the image data
        name: Output filename for the FITS file
        hdrData: Dictionary of header metadata
        hdu: Primary HDU object for FITS structure
        
    Example:
        >>> data = np.zeros((512, 512))
        >>> header = {'AUTHOR': 'Johan Mendez', 'METRIC': 'Kerr'}
        >>> fits = FITS(data, "output.fits", header)
        >>> fits.Write()
    """
    
    def __init__(self, info: np.ndarray, name: str, hdrData: Dict[str, str]):
        """
        Initialize FITS writer.
        
        Args:
            info: 2D image data array
            name: Output filename (should include .fits extension)
            hdrData: Dictionary of header key-value pairs
        """
        self.info = info
        self.name = name
        self.hdrData = hdrData
        self.hdu = fits.PrimaryHDU(data=self.info)
    
    def Write(self) -> None:
        """
        Write the FITS file to disk.
        
        Removes existing file if present to avoid conflicts.
        Adds header metadata from hdrData dictionary.
        """
        output_path = Path(self.name)
        
        # Remove existing file to prevent conflicts
        if output_path.exists():
            output_path.unlink()
        
        # Write primary HDU
        self.hdu.writeto(output_path)
        
        # Add custom header entries
        for key, value in self.hdrData.items():
            fits.setval(str(output_path), key, value=value)
    
    def showImage(self) -> None:
        """
        Display the image using matplotlib with proper orientation.
        
        Creates a visualization with:
        - Vertical flip for correct orientation
        - Proper colorbar sizing
        - Astronomical colormap (afmhot)
        - Coordinate labels in geometric units
        
        Saves preview to OUTPUT_DIR/Prueba.png
        """
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Load and orient data
        data = fits.open(self.name)[0].data
        data_flipped = np.flipud(data)
        
        # Render image
        im = ax.imshow(data_flipped, cmap=cm.afmhot, origin='lower')
        ax.set_xlabel("x [M]", fontsize=12)
        ax.set_ylabel("y [M]", fontsize=12)
        
        # Add colorbar with proper sizing
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = plt.colorbar(im, cax=cax, label='Normalized Intensity')
        
        plt.tight_layout()
        
        # Save preview
        preview_path = OUTPUT_DIR / "Prueba.png"
        plt.savefig(preview_path, dpi=150, bbox_inches='tight')
        plt.show()
        plt.close()
