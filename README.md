# Black Hole Ray Tracer

A high-performance ray tracing simulation for generating images of black hole shadows, including the photon ring and accretion disk.

![Example output](output/Schwarzschild_512x512.png)

**Author:** Johan Mendez  
**License:** MIT (2024)

---

## Overview

This project simulates how light bends around black holes using general relativistic ray tracing. It computes null geodesics in Schwarzschild and Kerr spacetimes, producing scientifically accurate visualizations of black hole shadows.

### Key Features

- **Metrics:** Schwarzschild (non-rotating) and Kerr (rotating) black holes
- **Resolution:** Configurable up to 1024×1024 pixels
- **Performance:** C++ backend with ~20× speedup over pure Python
- **Disk Models:** Infinite accretion disk with radial emission profiles
- **Output:** FITS format (astronomical standard) + PNG visualizations

---

## Quick Start

### Prerequisites

```bash
# Python dependencies
pip install numpy scipy matplotlib astropy pybind11 pillow

# macOS: Install OpenMP (optional but recommended)
brew install libomp
```

### Build and Run

```bash
# 1. Compile C++ module
./build_cpp.sh

# 2. Configure simulation
nano myconfig.py

# 3. Run simulation
python main.py
```

Output images are saved to `output/`.

---

## Configuration (myconfig.py)

```python
# Black hole parameters
METRIC = "Schwarzschild"  # or "Kerr"
MASS = 1.0                # Geometric units (G = c = 1)
SPIN = 0.7                # Kerr spin parameter (|a| < M)

# Image parameters
RESOLUTION = 512          # N×N pixels
SCREEN_SIZE = 15.0        # Field of view in units of M
INCLINATION = np.pi / 2.1  # Viewing angle (~85°, nearly edge-on)

# Disk parameters
OUTER_RADIUS = 20.0       # Outer disk edge in M
COROTATING = True         # Prograde disk rotation
```

---

## Architecture

```
┌─────────────┐     ┌──────────────┐     ┌─────────────┐
│  main.py    │────▶│ raytracer_cpp│────▶│  Geodesic   │
│  (Python)   │     │   (C++)      │     │ Integration │
└─────────────┘     └──────────────┘     └─────────────┘
      │                    │                    │
      ▼                    ▼                    ▼
Config & Viz         RK4 Integration      Null Geodesics
FITS/PNG Output      OpenMP Parallel      Kerr/Schwarzschild
```

---

## Performance Benchmarks

Hardware: Apple M1 (8 cores)

| Resolution | Python (legacy) | C++ (optimized) | Speedup |
|------------|----------------|----------------|---------|
| 256×256    | ~15 min        | ~45 s          | 20×     |
| 512×512    | ~60 min        | ~3 min         | 20×     |
| 1024×1024  | ~4 hours       | ~12 min        | 20×     |

---

## File Structure

```
RayTrace/
├── main.py                      # Entry point
├── myconfig.py                  # Configuration
├── build_cpp.sh                 # Build script
├── cpp/
│   ├── raytracer.h             # C++ header
│   ├── raytracer.cpp           # Geodesic integration
│   ├── binding.cpp             # Python bindings
│   └── Makefile                # Build configuration
├── screen/
│   └── imagePlane.py           # Photon initialization
├── common/
│   ├── geodesicIntegration.py  # Legacy Python integrator
│   ├── initialConditions.py    # Metric transformations
│   └── writeFits.py            # FITS I/O
├── photonSphere/
│   └── ps.py                   # Shadow rendering
├── accretionStructures/
│   └── inftyAccDisk.py         # Disk emission models
└── output/                      # Generated images
```

---

## Scientific Background

### Black Hole Shadow

The shadow is the silhouette formed by photons that are captured by the black hole. For Schwarzschild:

$$R_{\text{shadow}} = \sqrt{27}\,M \approx 5.196\,M$$

For Kerr black holes, frame dragging distorts the shadow into an asymmetric shape (see Bardeen 1973).

### Photon Ring

The bright ring around the shadow is the "photon ring" — light that orbited the black hole before escaping. This is the feature observed by the Event Horizon Telescope in M87*.

### Ray Tracing Method

1. **Setup:** Define image plane at distance $D \gg M$
2. **Initial Conditions:** For each pixel, compute photon position and momentum
3. **Integration:** Trace null geodesics backwards in time using RK4
4. **Detection:** Check if photon crosses equatorial plane (hits disk) or horizon
5. **Rendering:** Combine disk emission with black hole shadow

---

## Troubleshooting

### "C++ module not found"

```bash
cd cpp
make clean
make
cp raytracer_cpp.so ..
```

### Compilation errors on macOS

Ensure libomp is installed:
```bash
brew install libomp
export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
./build_cpp.sh
```

---

## References

- Bardeen, J. M. (1973). "Timelike and null geodesics in the Kerr metric." *Black Holes*, 215-239.
- Gralla, S. E., Lupsasca, A., & Marrone, D. P. (2019). "The Kerr null geodesic." *Physical Review D*, 102(12), 124004.
- Event Horizon Telescope Collaboration (2019). "First M87 Event Horizon Telescope Results." *The Astrophysical Journal Letters*, 875(1), L1.

---

## Acknowledgments

Original concept by ashcat (2018).  
C++ optimization and documentation by Johan Mendez (2024).
