#ifndef ABLMESOSCALEFORCING_H
#define ABLMESOSCALEFORCING_H

#include <string>
#include "AMReX_REAL.H"
#include "AMReX_Array.H"
#include "AMReX_Vector.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/ABLMesoscaleInput.H"
#include "amr-wind/wind_energy/ABL.H"

extern "C" {
    void dpptrf_(const char *uplo, const int *n, double *AP, int *info);
    void dpptrs_(
        const char *uplo, const int *n, const int *nrhs, double *AP, double *B,
        const int *ldb, int *info);
    void dgemv_(const char* trans, const int* m, const int* n, double* alpha,
        double* A, const int* lda, double* X, const int* incx, double* beta,
        double* Y, const int* incy);
}


namespace amr_wind {

class ABLMesoscaleForcing
{

public:
    ABLMesoscaleForcing(const CFDSim& sim, const std::string identifier);

    void indirectForcingInit();
    void GP_forcingInit();

    static void invertMat(
        const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
        amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im);

protected:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    const ABLMesoscaleInput& m_meso_file;
    const std::string m_identifier;

    int m_idx_time;
    int m_axis{2};

    // inputs
    std::string m_forcing_scheme;  // direct | indirect | Gaussian_process
    amrex::Real m_gain_coeff{0.2}; // relaxation factor for all forcing schemes
    bool m_debug{false};           // print out profiles for debugging
    bool debug() { return m_debug; }

    // planar averaging heights
    amrex::Vector<amrex::Real> m_zht;
    int m_nht;

    //
    // regression vars for indirect forcing
    //
    // int m_ind_polyOrder;
    int m_norm_zmax = 0;
    amrex::Real m_scaleFact =
        1e-3; // to improve conditioning of regression matrix

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> m_im_zTz;

    // - input (or current value of, if forcing_transition is specified)
    // regression weighting profile
    std::string m_forcing_transition;
    amrex::Real m_transition_height; // constant or variable
    bool m_update_transition_height = false;
    amrex::Real m_transition_thickness; // constant

    // - WORKAROUND
    amrex::Vector<amrex::Real> m_transition_height_hist;

    amrex::Vector<amrex::Real> m_weighting_heights; // constant or variable
    amrex::Vector<amrex::Real> m_weighting_values;  // constant or variable

    // - W(z) at planar averaging heights
    amrex::Vector<amrex::Real> m_W;

    //
    // weighting profile manipulation
    //
    void setTransitionWeighting(); // will override weighting profile
    void
    updateWeights(); // interpolate current profile to planar averaging heights

    //
    // forcing profile manipulation
    //
    bool forcingToConstant()
    {
        // - check if m_forcing_transition ends with "Constant"
        if (m_forcing_transition.length() < 8) {
            return false;
        }
        return (
            0 == m_forcing_transition.compare(
                     m_forcing_transition.length() - 8, 8, "Constant"));
    }

    // - slope varies linearly from the actual gradient of the forcing profile
    //   at m_transition_height to 0 at (m_transition_height +
    //   m_transition_thickness)
    void constantForcingTransition(amrex::Vector<amrex::Real>& error);

    // - blend between a lower and upper forcing profile; the blending fraction
    //   is assumed to be equal to W(z), but this is not a requirement; W=1 and
    //   W=0 correspond to the lower and upper profiles, respectively
    void blendForcings(
        const amrex::Vector<amrex::Real> lower,
        const amrex::Vector<amrex::Real> upper,
        amrex::Vector<amrex::Real>& error);

    //
    // Gaussian process (GP) regression
    // 
    // - "observations" are the direct forcing profiles to exactly reproduce
    //   the profiles of the assimilated quantities
    // - "predictions" are the resulting indirect forcing profiles, given
    //   "noisy" obesrvations--i.e., forcing that we don't completely trust
    //   in specified regions
    //
    bool m_update_var_mat = true; // update Sigma_11 at every time step
    bool m_update_covar_mat = false; // update Sigma_12 at every time step
    int m_update_freq = 1; // how often to update (co)variance matrices
    std::string m_covar_func = "sq_exp"; // covariance function (kernel)
    amrex::Real m_length_scale = 1.0; // to scale covariance function
    amrex::Real m_sigma_noise = 1.0; // standard deviation of the base error
    std::string m_spec_err_type = "none"; // or "forcing_variance"

    // lapack input/output data
    std::vector<double> Sigma11; // packed LHS matrix
    std::vector<double> Sigma12; // full RHS matrix

    // - covariance functions
    // hard-coded for now: squared exponential
    inline double GP_covar_func(double z_i, double z_j) {
        double d = (z_i - z_j) / m_length_scale;
        return std::exp(-0.5*d*d);
    }

    // - kernel of observations to predictions (Sigma12), row-major order; 
    //   should only be called once (unless the input data heights change
    //   over time
    void GP_updateSigma12();

    // - kernel of noisy observations (Sigma11) matrix, in upper triangular
    //   form, row-major order; this needs to be updated if the specific error
    //   is time-varying
    void GP_updateSigma11Packed();

    // - calculate the Gaussian process posterior (mean only, we are currently
    //   not using the covariance in any way)
    amrex::Vector<amrex::Real> GP_posteriorMean(
        const amrex::Vector<amrex::Real> y1);

};

} // namespace amr_wind

#endif
