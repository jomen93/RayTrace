/**
 * @file raytracer.h
 * @brief C++ ray tracing engine for black hole shadow simulations.
 * 
 * This header defines the interface for computing null geodesics in Schwarzschild
 * and Kerr spacetimes. The implementation uses 4th-order Runge-Kutta integration
 * with adaptive detection of the equatorial plane crossing.
 * 
 * @author Johan Mendez
 * @copyright (c) 2024
 * 
 * @note Geometric units are used throughout: G = c = 1, lengths in units of M.
 */

#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace raytracer {

// =============================================================================
// CONSTANTS
// =============================================================================

/**
 * @brief Numerical tolerance for equatorial plane detection.
 * 
 * The equatorial plane is defined by theta = pi/2, or equivalently cos(theta) = 0.
 * We consider the plane "crossed" when |cos(theta)| < TOLERANCE.
 */
inline constexpr double TOLERANCE = 5.0e-5;

/**
 * @brief Maximum integration steps before termination.
 * 
 * Prevents infinite loops for photons that neither cross the equator
 * nor fall into the black hole.
 */
inline constexpr int MAX_STEPS = 100000;

// =============================================================================
// TYPE DEFINITIONS
// =============================================================================

/**
 * @brief Photon state vector in 8-dimensional phase space.
 * 
 * Components: [t, r, theta, phi, p_t, p_r, p_theta, p_phi]
 * where the first four are spacetime coordinates and the last four are
 * conjugate momenta (covariant components).
 */
using PhotonState = std::array<double, 8>;

/**
 * @brief Geodesic equation function type.
 * 
 * Computes derivatives dx/dtau given current state x, mass M, and spin a.
 * Output is written to dxdt.
 */
using GeodesicFunc = void(*)(const PhotonState& x, double M, double a, 
                              PhotonState& dxdt);

// =============================================================================
// METRIC FUNCTIONS
// =============================================================================

/**
 * @brief Compute Schwarzschild metric components.
 * 
 * Metric: ds² = -(1-2M/r)dt² + dr²/(1-2M/r) + r²dΩ²
 * 
 * @param[in] x Photon state (position components used)
 * @param[in] M Black hole mass
 * @param[in] a Unused (Schwarzschild is non-rotating)
 * @param[out] gtt,gtr,gtphi,grr,gthth,gphph Metric components
 */
void schwarzschild_metric(const PhotonState& x, double M, double a,
                          double& gtt, double& gtphi, double& grr, 
                          double& gthth, double& gphph);

/**
 * @brief Compute Kerr metric components in Boyer-Lindquist coordinates.
 * 
 * Metric includes cross term g_tphi for frame dragging.
 * 
 * @param[in] x Photon state
 * @param[in] M Black hole mass  
 * @param[in] a Kerr spin parameter (|a| < M)
 * @param[out] gtt,gtr,gtphi,grr,gthth,gphph Metric components
 */
void kerr_metric(const PhotonState& x, double M, double a,
                 double& gtt, double& gtphi, double& grr, 
                 double& gthth, double& gphph);

// =============================================================================
// GEODESIC EQUATIONS
// =============================================================================

/**
 * @brief Schwarzschild geodesic equations in Hamiltonian form.
 * 
 * Computes dx^μ/dτ given current state using conserved quantities
 * E = -p_t and L = p_phi.
 * 
 * @param[in] x Current state
 * @param[in] M Black hole mass
 * @param[in] a Unused (set to 0)
 * @param[out] dxdt State derivatives
 */
void schwarzschild_geodesic(const PhotonState& x, double M, double a, 
                            PhotonState& dxdt);

/**
 * @brief Kerr geodesic equations in Hamiltonian form.
 * 
 * Full Kerr geodesics including frame dragging effects.
 * Uses Carter's constant implicitly through the equations.
 * 
 * @param[in] x Current state
 * @param[in] M Black hole mass
 * @param[in] a Kerr spin parameter
 * @param[out] dxdt State derivatives
 */
void kerr_geodesic(const PhotonState& x, double M, double a, 
                   PhotonState& dxdt);

// =============================================================================
// INTEGRATION
// =============================================================================

/**
 * @brief Perform one RK4 integration step.
 * 
 * Advances state by time step dt using 4th-order Runge-Kutta method.
 * 
 * @param[in] x Current state
 * @param[in] dt Time step (affine parameter increment)
 * @param[in] M Black hole mass
 * @param[in] a Spin parameter
 * @param[in] geodesic Geodesic function to use
 * @param[out] x_next State after time step
 */
void rk4_step(const PhotonState& x, double dt, double M, double a,
              GeodesicFunc geodesic, PhotonState& x_next);

/**
 * @brief Integrate geodesic to equatorial plane.
 * 
 * Main entry point for photon tracing. Integrates backwards in time
 * from the image plane until the photon either:
 * - Crosses the equatorial plane (returns radius r)
 * - Falls into the event horizon (returns 0)
 * - Escapes to infinity (returns 0)
 * 
 * @param[in] init_cond Initial state at image plane
 * @param[in] M Black hole mass
 * @param[in] r_eh Event horizon radius (2M for Schwarzschild)
 * @param[in] dt Integration time step (negative for backwards integration)
 * @param[in] geodesic Geodesic equation function
 * @return Radius at equator, or 0 if photon doesn't reach equator
 */
double integrate_geodesic(const PhotonState& init_cond, double M, 
                          double r_eh, double dt, GeodesicFunc geodesic);

/**
 * @brief Special integration for meridional trajectories (p_phi = 0).
 * 
 * Optimized integration for photons in the meridional plane.
 * These photons have zero angular momentum and move purely in the
 * r-theta plane.
 * 
 * @param[in] init_cond Initial state (must have p_phi ≈ 0)
 * @param[in] M Black hole mass
 * @param[in] r_eh Event horizon radius
 * @param[in] dt Integration time step
 * @param[in] geodesic Geodesic function
 * @return Radius at equator, or 0
 */
double integrate_meridional(const PhotonState& init_cond, double M, 
                            double r_eh, double dt, GeodesicFunc geodesic);

// =============================================================================
// PHOTON INITIALIZATION
// =============================================================================

/**
 * @brief Create photon initial conditions from image plane coordinates.
 * 
 * Transforms (alpha, beta) image plane coordinates to initial state vector
 * using the camera model. Applies the metric at the photon position to
 * convert contravariant momentum to covariant momentum.
 * 
 * @param[in] alpha Horizontal image coordinate
 * @param[in] beta Vertical image coordinate
 * @param[in] D Distance to black hole
 * @param[in] inclination Viewing angle (radians)
 * @param[in] M Black hole mass
 * @return Initial state vector ready for integration
 */
PhotonState create_photon(double alpha, double beta, double D, 
                          double inclination, double M);

// =============================================================================
// MAIN INTERFACE
// =============================================================================

/**
 * @brief Trace all photons for image generation.
 * 
 * Primary interface for Python bindings. Traces N×N photons through the
 * specified metric and returns the radii at the equatorial plane.
 * 
 * @param[in] alpha_range Array of alpha coordinates
 * @param[in] beta_range Array of beta coordinates  
 * @param[in] D Observer distance
 * @param[in] inclination Viewing angle
 * @param[in] M Black hole mass
 * @param[in] a Spin parameter
 * @param[in] num_pixels Image resolution N (N×N grid)
 * @param[in] metric_name "Schwarzschild" or "Kerr"
 * @return Flattened vector of radii (row-major order)
 */
std::vector<double> trace_photons(
    const std::vector<double>& alpha_range,
    const std::vector<double>& beta_range,
    double D, double inclination, double M, double a,
    int num_pixels, const std::string& metric_name);

/**
 * @brief Validate photon state for numerical stability.
 * 
 * Checks for NaN, Inf, and physically impossible values.
 * 
 * @param[in] x State vector to validate
 * @return True if state is valid
 */
bool is_valid_state(const PhotonState& x);

} // namespace raytracer

#endif // RAYTRACER_H
