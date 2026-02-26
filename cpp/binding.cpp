/**
 * @file binding.cpp
 * @brief Python bindings for C++ ray tracing engine using pybind11.
 * 
 * This module exposes the C++ ray tracer to Python, allowing efficient
 * computation of null geodesics while maintaining a clean Python interface.
 * 
 * @author Johan Mendez
 * @copyright (c) 2024
 * 
 * Example Python usage:
 *     import raytracer_cpp
 *     import numpy as np
 *     
 *     alpha = np.linspace(-15, 15, 256)
 *     beta = np.linspace(-15, 15, 256)
 *     
 *     radii = raytracer_cpp.trace_photons(
 *         alpha, beta, D=100.0, inclination=np.pi/2.1,
 *         M=1.0, a=0.0, num_pixels=256, metric_name="Schwarzschild"
 *     )
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "raytracer.h"

namespace py = pybind11;
using namespace raytracer;

/**
 * @brief Wrapper for trace_photons that accepts NumPy arrays.
 * 
 * Converts NumPy arrays to std::vector, calls the C++ implementation,
 * and returns the result as a NumPy array.
 */
py::array_t<double> trace_photons_py(
    py::array_t<double> alpha_range,
    py::array_t<double> beta_range,
    double D, double inclination, double M, double a,
    int num_pixels, std::string metric_name) {
    
    // Extract data from NumPy arrays
    py::buffer_info alpha_info = alpha_range.request();
    py::buffer_info beta_info = beta_range.request();
    
    double* alpha_ptr = static_cast<double*>(alpha_info.ptr);
    double* beta_ptr = static_cast<double*>(beta_info.ptr);
    
    // Convert to STL containers
    std::vector<double> alpha_vec(alpha_ptr, alpha_ptr + alpha_info.size);
    std::vector<double> beta_vec(beta_ptr, beta_ptr + beta_info.size);
    
    // Call C++ implementation
    std::vector<double> results = trace_photons(
        alpha_vec, beta_vec, D, inclination, M, a, 
        num_pixels, metric_name);
    
    // Return as NumPy array with shape (num_pixels, num_pixels)
    return py::array_t<double>(
        {num_pixels, num_pixels},
        {num_pixels * sizeof(double), sizeof(double)},
        results.data());
}

/**
 * @brief Compute Schwarzschild ISCO radius.
 * @param M Black hole mass
 * @return ISCO radius = 6M
 */
double schwarzschild_r_isco(double M) {
    return 6.0 * M;
}

/**
 * @brief Compute Kerr ISCO radius.
 * 
 * Uses analytical formula from Bardeen et al. (1972).
 * 
 * @param M Black hole mass
 * @param a Kerr spin parameter
 * @param corotating True for prograde, False for retrograde
 * @return ISCO radius
 */
double kerr_r_isco(double M, double a, bool corotating) {
    // Dimensionless spin
    const double a_hat = a / M;
    
    // Auxiliary functions
    const double Z1 = 1.0 + std::cbrt(1.0 - a_hat * a_hat) * 
                      (std::cbrt(1.0 + a_hat) + std::cbrt(1.0 - a_hat));
    const double Z2 = std::sqrt(3.0 * a_hat * a_hat + Z1 * Z1);
    
    if (corotating) {
        return M * (3.0 + Z2 - std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    } else {
        return M * (3.0 + Z2 + std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));
    }
}

PYBIND11_MODULE(raytracer_cpp, m) {
    m.doc() = R"doc(
        C++ Ray Tracer for Black Hole Shadows
        
        High-performance geodesic integration for computing black hole
        shadow images. Uses 4th-order Runge-Kutta integration.
        
        Author: Johan Mendez
        
        Functions:
            trace_photons: Main ray tracing function
            schwarzschild_r_isco: Compute ISCO for Schwarzschild
            kerr_r_isco: Compute ISCO for Kerr
    )doc";
    
    m.def("trace_photons", &trace_photons_py,
          R"doc(
              Trace photons through black hole spacetime.
              
              Args:
                  alpha_range: NumPy array of horizontal image coordinates
                  beta_range: NumPy array of vertical image coordinates
                  D: Observer distance (geometric units)
                  inclination: Viewing angle in radians
                  M: Black hole mass
                  a: Kerr spin parameter (0 for Schwarzschild)
                  num_pixels: Image resolution (NxN)
                  metric_name: "Schwarzschild" or "Kerr"
                  
              Returns:
                  2D NumPy array of radii at equatorial plane
          )doc",
          py::arg("alpha_range"),
          py::arg("beta_range"),
          py::arg("D"),
          py::arg("inclination"),
          py::arg("M"),
          py::arg("a"),
          py::arg("num_pixels"),
          py::arg("metric_name"));
    
    m.def("schwarzschild_r_isco", &schwarzschild_r_isco,
          "Compute ISCO radius for Schwarzschild metric",
          py::arg("M"));
    
    m.def("kerr_r_isco", &kerr_r_isco,
          "Compute ISCO radius for Kerr metric",
          py::arg("M"), py::arg("a"), py::arg("corotating"));
}
