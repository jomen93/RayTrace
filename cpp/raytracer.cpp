/**
 * @file raytracer.cpp
 * @brief Implementation of black hole ray tracing engine.
 * 
 * This file implements the numerical geodesic integration for computing
 * black hole shadow images. It uses 4th-order Runge-Kutta integration with
 * specialized handling for meridional trajectories.
 * 
 * @author Johan Mendez
 * @copyright (c) 2024
 * 
 * The integration proceeds backwards in affine parameter from the image plane
 * (distant observer) until photons either cross the equatorial plane, fall
 * into the event horizon, or escape to infinity.
 */

#include "raytracer.h"
#include <algorithm>
#include <cstring>
#include <cmath>

// Compatibility macros for non-OpenMP environments
#ifndef _OPENMP
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

namespace raytracer {

// =============================================================================
// METRIC IMPLEMENTATIONS
// =============================================================================

void schwarzschild_metric(const PhotonState& x, double M, double /*a*/,
                          double& gtt, double& gtphi, double& grr, 
                          double& gthth, double& gphph) {
    const double r = x[1];
    const double theta = x[2];
    
    // Metric coefficient: g_tt = -(1 - 2M/r)
    const double one_minus_2M_r = 1.0 - 2.0 * M / r;
    
    gtt = -one_minus_2M_r;
    gtphi = 0.0;  // No frame dragging in Schwarzschild
    grr = 1.0 / one_minus_2M_r;
    gthth = r * r;
    gphph = r * r * std::sin(theta) * std::sin(theta);
}

void kerr_metric(const PhotonState& x, double M, double a,
                 double& gtt, double& gtphi, double& grr, 
                 double& gthth, double& gphph) {
    const double r = x[1];
    const double theta = x[2];
    
    // Boyer-Lindquist auxiliary functions
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);
    const double sigma = r * r + a * a * cos_theta * cos_theta;
    const double delta = r * r - 2.0 * M * r + a * a;
    
    // Kerr metric components with frame dragging
    gtt = -(1.0 - 2.0 * M * r / sigma);
    gtphi = -2.0 * M * a * r * sin_theta * sin_theta / sigma;
    grr = sigma / delta;
    gthth = sigma;
    gphph = (r * r + a * a + 2.0 * M * a * a * r * sin_theta * sin_theta / sigma) 
            * sin_theta * sin_theta;
}

// =============================================================================
// GEODESIC EQUATIONS
// =============================================================================

void schwarzschild_geodesic(const PhotonState& x, double M, double /*a*/, 
                            PhotonState& dxdt) {
    // Unpack state vector
    const double r = x[1];
    const double theta = x[2];
    const double pt = x[4];
    const double pr = x[5];
    const double pth = x[6];
    const double pphi = x[7];
    
    // Conserved quantities
    const double E = -pt;  // Energy (normalized)
    const double L = pphi; // Angular momentum
    
    const double r2 = r * r;
    const double sin_theta = std::sin(theta);
    const double sin_theta2 = sin_theta * sin_theta;
    
    // Geodesic equations from Hamiltonian H = (1/2)g^{μν}p_μp_ν
    const double one_minus_2M_r = 1.0 - 2.0 * M / r;
    
    dxdt[0] = E * r2 / (r2 - 2.0 * M * r);           // dt/dτ
    dxdt[1] = one_minus_2M_r * pr;                    // dr/dτ
    dxdt[2] = pth / r2;                               // dθ/dτ
    dxdt[3] = L / (r2 * sin_theta2);                 // dφ/dτ
    dxdt[4] = 0.0;                                    // dp_t/dτ (conserved)
    dxdt[5] = -M * (pr * pr / r2)                     // dp_r/dτ
              + pth * pth / (r * r2) 
              + L * L / (r * r2 * sin_theta2) 
              - M * (E * E) / ((r - 2.0 * M) * (r - 2.0 * M));
    dxdt[6] = (std::cos(theta) / (sin_theta2 * sin_theta)) * (L * L / r2); // dp_θ/dτ
    dxdt[7] = 0.0;                                    // dp_φ/dτ (conserved)
}

void kerr_geodesic(const PhotonState& x, double M, double a, 
                   PhotonState& dxdt) {
    const double r = x[1];
    const double theta = x[2];
    const double pr = x[5];
    const double pth = x[6];
    const double pphi = x[7];
    
    // Conserved quantities
    const double E = -x[4];
    const double L = pphi;
    
    // Boyer-Lindquist functions
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);
    const double sigma = r * r + a * a * cos_theta * cos_theta;
    const double delta = r * r - 2.0 * M * r + a * a;
    
    // Carter constant components
    const double W = E * (r * r + a * a) - a * L;
    const double sin_theta2 = sin_theta * sin_theta;
    const double cos_theta2 = cos_theta * cos_theta;
    
    const double partXi = r * r + (L - a * E) * (L - a * E) 
                         + a * a * (1.0 - E * E) * cos_theta2
                         + (L * L * cos_theta2) / sin_theta2;
    const double Xi = W * W - delta * partXi;
    
    // Derivatives for Hamilton's equations
    const double dXidE = 2.0 * W * (r * r + a * a) 
                        + 2.0 * a * delta * (L - a * E * sin_theta2);
    const double dXidL = -2.0 * a * W + 2.0 * a * E * delta 
                        - 2.0 * L * delta / sin_theta2;
    const double dXidr = 4.0 * r * E * W - 2.0 * (r - M) * partXi - 2.0 * r * delta;
    
    const double sigma2 = sigma * sigma;
    const double dAdr = (r - M) / sigma - (r * delta) / sigma2;
    const double dBdr = -r / sigma2;
    const double dCdr = dXidr / (2.0 * delta * sigma) 
                       - (Xi * (r - M)) / (sigma * delta * delta)
                       - r * Xi / (delta * sigma2);
    
    const double auxth = a * a * cos_theta * sin_theta;
    const double dAdth = delta * auxth / sigma2;
    const double dBdth = auxth / sigma2;
    const double dCdth = ((1.0 - E * E) * auxth + L * L * cos_theta / (sin_theta2 * sin_theta)) / sigma
                        + (Xi / (delta * sigma2)) * auxth;
    
    // Kerr geodesic equations
    dxdt[0] = dXidE / (2.0 * delta * sigma);
    dxdt[1] = (delta / sigma) * pr;
    dxdt[2] = pth / sigma;
    dxdt[3] = -dXidL / (2.0 * delta * sigma);
    dxdt[4] = 0.0;
    dxdt[5] = -dAdr * pr * pr - dBdr * pth * pth + dCdr;
    dxdt[6] = -dAdth * pr * pr - dBdth * pth * pth + dCdth;
    dxdt[7] = 0.0;
}

// =============================================================================
// NUMERICAL INTEGRATION
// =============================================================================

void rk4_step(const PhotonState& x, double dt, double M, double a,
              GeodesicFunc geodesic, PhotonState& x_next) {
    PhotonState k1, k2, k3, k4, temp;
    
    // Stage 1
    geodesic(x, M, a, k1);
    
    // Stage 2
    for (int i = 0; i < 8; ++i) {
        temp[i] = x[i] + 0.5 * dt * k1[i];
    }
    geodesic(temp, M, a, k2);
    
    // Stage 3
    for (int i = 0; i < 8; ++i) {
        temp[i] = x[i] + 0.5 * dt * k2[i];
    }
    geodesic(temp, M, a, k3);
    
    // Stage 4
    for (int i = 0; i < 8; ++i) {
        temp[i] = x[i] + dt * k3[i];
    }
    geodesic(temp, M, a, k4);
    
    // Combine stages
    for (int i = 0; i < 8; ++i) {
        x_next[i] = x[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

bool is_valid_state(const PhotonState& x) {
    for (int i = 0; i < 8; ++i) {
        if (std::isnan(x[i]) || std::isinf(x[i])) {
            return false;
        }
    }
    // Physical constraints
    if (x[1] < 0) return false;
    if (x[2] < 0 || x[2] > M_PI) return false;
    return true;
}

// =============================================================================
// MERIDIONAL INTEGRATION (p_phi = 0)
// =============================================================================

double integrate_meridional(const PhotonState& init_cond, double M, 
                            double r_eh, double dt, GeodesicFunc geodesic) {
    PhotonState x = init_cond;
    
    // Enforce exact meridional condition
    x[7] = 0.0;
    
    // Clamp initial theta to valid range
    if (x[2] <= 0) x[2] = 1e-10;
    if (x[2] >= M_PI) x[2] = M_PI - 1e-10;
    
    double cos_theta = std::cos(x[2]);
    
    for (int step = 0; step < MAX_STEPS; ++step) {
        PhotonState x_next;
        rk4_step(x, dt, M, 0.0, geodesic, x_next);
        
        // Maintain meridional constraint
        x_next[7] = 0.0;
        
        if (!is_valid_state(x_next)) {
            return 0.0;
        }
        
        // Keep theta in valid range
        if (x_next[2] <= 0) x_next[2] = 1e-10;
        if (x_next[2] >= M_PI) x_next[2] = M_PI - 1e-10;
        
        const double cos_theta_new = std::cos(x_next[2]);
        
        // Check equatorial plane crossing
        if (cos_theta > TOLERANCE && cos_theta_new < TOLERANCE) {
            return x_next[1];
        }
        
        // Check horizon crossing
        if (x_next[1] <= r_eh + 1e-6) {
            return 0.0;
        }
        
        // Check escape without crossing
        if (x_next[1] > 500.0 * M && std::abs(cos_theta_new) > 0.9) {
            return 0.0;
        }
        
        x = x_next;
        cos_theta = cos_theta_new;
    }
    
    return 0.0;
}

// =============================================================================
// GENERAL INTEGRATION
// =============================================================================

double integrate_geodesic(const PhotonState& init_cond, double M, 
                          double r_eh, double dt, GeodesicFunc geodesic) {
    PhotonState x = init_cond;
    
    // Use specialized integration for meridional trajectories
    if (std::abs(x[7]) < 1e-14) {
        return integrate_meridional(init_cond, M, r_eh, dt, geodesic);
    }
    
    // Clamp initial theta
    if (x[2] <= 0) x[2] = 1e-10;
    if (x[2] >= M_PI) x[2] = M_PI - 1e-10;
    
    double cos_theta = std::cos(x[2]);
    
    for (int step = 0; step < MAX_STEPS; ++step) {
        PhotonState x_next;
        rk4_step(x, dt, M, 0.0, geodesic, x_next);
        
        if (!is_valid_state(x_next)) {
            return 0.0;
        }
        
        if (x_next[2] <= 0) x_next[2] = 1e-10;
        if (x_next[2] >= M_PI) x_next[2] = M_PI - 1e-10;
        
        const double cos_theta_new = std::cos(x_next[2]);
        
        // Equatorial plane crossing
        if (cos_theta > TOLERANCE && cos_theta_new < TOLERANCE) {
            return x_next[1];
        }
        
        if (x_next[1] <= r_eh + 1e-6) {
            return 0.0;
        }
        
        if (x_next[1] > 1e6) {
            return 0.0;
        }
        
        x = x_next;
        cos_theta = cos_theta_new;
    }
    
    return 0.0;
}

// =============================================================================
// PHOTON INITIALIZATION
// =============================================================================

PhotonState create_photon(double alpha, double beta, double D, 
                          double inclination, double M) {
    PhotonState x{};
    
    // Transform to spherical coordinates
    const double r = std::sqrt(alpha * alpha + beta * beta + D * D);
    double cos_theta = (beta * std::sin(inclination) + D * std::cos(inclination)) / r;
    
    // Clamp to valid range
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
    
    const double theta = std::acos(cos_theta);
    const double phi = std::atan2(alpha, D * std::sin(inclination) - beta * std::cos(inclination));
    
    // Set position
    x[0] = 0.0;
    x[1] = r;
    x[2] = theta;
    x[3] = phi;
    
    // Compute contravariant momentum components
    const double K0 = 1.0;
    const double kr = (D / r) * K0;
    
    const double aux = alpha * alpha + 
                      (-beta * std::cos(inclination) + D * std::sin(inclination)) * 
                      (-beta * std::cos(inclination) + D * std::sin(inclination));
    
    // Prevent division by zero
    const double safe_aux = std::max(aux, 1e-15);
    const double sqrt_aux = std::sqrt(safe_aux);
    
    const double ktheta = (K0 / sqrt_aux) * 
                          (-std::cos(inclination) + 
                           (beta * std::sin(inclination) + D * std::cos(inclination)) * 
                           (D / (r * r)));
    
    // Set kphi = 0 for meridional photons to avoid numerical issues
    double kphi;
    if (std::abs(alpha) < 1e-12) {
        kphi = 0.0;
    } else {
        kphi = -alpha * std::sin(inclination) * K0 / safe_aux;
    }
    
    const double kt = std::sqrt(kr * kr + r * r * ktheta * ktheta 
                               + r * r * std::sin(theta) * std::sin(theta) * kphi * kphi);
    
    // Convert to covariant components using metric
    const double one_minus_2M_r = 1.0 - 2.0 * M / r;
    const double sin_theta = std::sin(theta);
    
    const double gtt = -one_minus_2M_r;
    const double grr = 1.0 / one_minus_2M_r;
    const double gthth = r * r;
    const double gphph = r * r * sin_theta * sin_theta;
    
    x[4] = gtt * kt;
    x[5] = grr * kr;
    x[6] = gthth * ktheta;
    x[7] = gphph * kphi;
    
    return x;
}

// =============================================================================
// MAIN INTERFACE
// =============================================================================

std::vector<double> trace_photons(
    const std::vector<double>& alpha_range,
    const std::vector<double>& beta_range,
    double D, double inclination, double M, double a,
    int num_pixels, const std::string& metric_name) {
    
    std::vector<double> results(num_pixels * num_pixels, 0.0);
    
    // Select appropriate geodesic function
    GeodesicFunc geodesic;
    double r_eh;
    
    if (metric_name == "Schwarzschild") {
        geodesic = schwarzschild_geodesic;
        r_eh = 2.0 * M;
    } else if (metric_name == "Kerr") {
        geodesic = kerr_geodesic;
        r_eh = M + std::sqrt(M * M - a * a);
    } else {
        geodesic = schwarzschild_geodesic;
        r_eh = 2.0 * M;
    }
    
    // Integration step (negative for backwards integration)
    const double dt = -0.05;
    
    // Trace all photons
    for (int j = 0; j < num_pixels; ++j) {
        for (int k = 0; k < num_pixels; ++k) {
            const double alpha = alpha_range[k];
            const double beta = beta_range[num_pixels - 1 - j];
            
            const PhotonState init_cond = create_photon(alpha, beta, D, inclination, M);
            const double radius = integrate_geodesic(init_cond, M, r_eh, dt, geodesic);
            
            results[j * num_pixels + k] = radius;
        }
    }
    
    return results;
}

} // namespace raytracer
