#include "amr-wind/wind_energy/ABLMesoscaleForcing.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_Print.H"
#include "AMReX_ParmParse.H"

// WORKAROUND
#include <fstream>

namespace amr_wind {

ABLMesoscaleForcing::ABLMesoscaleForcing(
    const CFDSim& sim, const std::string identifier)
    : m_time(sim.time()), m_mesh(sim.mesh())
    , m_meso_file(sim.physics_manager().get<amr_wind::ABL>().abl_meso_file())
    , m_identifier(identifier)
{
    amrex::Print() << "Constructing " << identifier << " object" << std::endl;

    amrex::ParmParse pp(identifier);
    pp.query("forcing_scheme", m_forcing_scheme);
    pp.query("control_gain", m_gain_coeff);
    pp.query("debug", m_debug);
    amrex::Print() << "  forcing_scheme : " << m_forcing_scheme << std::endl;
    amrex::Print() << "  control_gain   : " << m_gain_coeff << std::endl;

    if (pp.query("forcing_transition", m_forcing_transition) == 0) {
        amrex::Print() << "  using full profile assimilation by default"
                       << std::endl;
        m_forcing_transition = "none";
    }

    if (amrex::toLower(m_forcing_scheme) == "indirect") {
        if (amrex::toLower(m_forcing_transition) == "none") {
            if (pp.queryarr("weighting_heights", m_weighting_heights) == 1) {
                pp.getarr("weighting_values", m_weighting_values);
                amrex::Print() << "  given weighting profile" << std::endl;
                for (int i = 0; i < m_weighting_heights.size(); ++i) {
                    amrex::Print() << "  " << m_weighting_heights[i] << " "
                                   << m_weighting_values[i] << std::endl;
                }
                AMREX_ALWAYS_ASSERT(
                    m_weighting_heights.size() == m_weighting_values.size());

            } else {
                // default is to have uniform weighting throughout
                amrex::Print() << "  setting default weighting" << std::endl;
                amrex::Real zmin = m_mesh.Geom(0).ProbLo(m_axis);
                amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
                m_weighting_heights = {zmin, zmax};
                m_weighting_values = {1.0, 1.0};
            }
        } else // weightings will be automatically set based on forcing
               // transition
        {
            pp.get(
                "transition_thickness",
                m_transition_thickness); // constant, required
            if (pp.query("constant_transition_height", m_transition_height) ==
                1) {
                // set weighting profile
                setTransitionWeighting();
            } else {
                // expect to read transition_height history in netCDF input file
                // weighting profile will be updated at every step
                m_update_transition_height = true;

                // ***WORKAROUND***
                // after commenting out code to read in transition_heights from
                // the mesoscale input file (see FIXME lines)
                std::string fname;
                pp.get("transition_heights_file", fname);
                std::ifstream datfile(fname);
                AMREX_ALWAYS_ASSERT(datfile.is_open());
                amrex::Real tval, zval;
                int ntimes = 0;
                while (datfile >> tval >> zval) {
                    ntimes++;
                }
                datfile.clear();
                datfile.seekg(0);
                m_transition_height_hist.resize(ntimes);
                amrex::Print() << "WORKAROUND: Reading transition heights from "
                               << fname << std::endl;
                for (int itime = 0; itime < ntimes; itime++) {
                    datfile >> tval >> zval;
                    m_transition_height_hist[itime] = zval;
                    amrex::Print() << tval << " " << zval << std::endl;
                }
                amrex::Print()
                    << "Note: the times in the mesoscale input file "
                    << "must match these " << ntimes << " values" << std::endl;
            }
        }

        if ((pp.query("normalize_by_zmax", m_norm_zmax) == 1) &&
            (m_norm_zmax != 0)) {
            amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
            m_scaleFact = 1.0 / zmax;
            amrex::Print() << "  set scaling factor to " << m_scaleFact
                           << std::endl;
        }
    } // if forcing scheme is "indirect"

    else if (amrex::toLower(m_forcing_scheme) == "gaussian_process") {
#ifdef AMR_WIND_USE_LAPACK
        pp.query("update_var_mat", m_update_var_mat);
        pp.query("update_covar_mat", m_update_var_mat);
        pp.query("update_freq", m_update_freq);
        if (m_update_var_mat || m_update_covar_mat) {
            amrex::Print() << "  update";
            if (m_update_var_mat) amrex::Print() << " Sigma_11";
            if (m_update_covar_mat) amrex::Print() << " Sigma_12";
            amrex::Print() << " every " << m_update_freq << " steps"
                << std::endl;
        }
        pp.query("covariance_function", m_covar_func);
        if (amrex::toLower(m_covar_func) != "sq_exp") {
            amrex::Print() << "  ignoring specified covariance function "
                << m_covar_func
                << ", only squared-exponential kernel is implemented"
                << std::endl;
        }
        pp.query("length_scale", m_length_scale);
        if (m_length_scale <= 0) {
            amrex::Abort("Length scale must be > 0");
        }
        pp.query("sigma_noise", m_sigma_noise);
        pp.query("specified_error", m_spec_err_type);
        if ((amrex::toLower(m_spec_err_type) != "none") &&
            (amrex::toLower(m_spec_err_type) != "forcing_variance")) {
            amrex::Abort("Unrecognized specified_error type");
        }
        amrex::Print() << "  covariance func : " << m_covar_func << std::endl;
        amrex::Print() << "  length scale    : " << m_length_scale << std::endl;
        amrex::Print() << "  base error      : " << m_sigma_noise << std::endl;
        amrex::Print() << "  specified error : " << m_spec_err_type << std::endl;
#else
        amrex::Abort(
            "LAPACK support needed for Gaussian-process profile assimilation.");
#endif
    }

    else if (amrex::toLower(m_forcing_scheme) != "direct") {
        // TOOD: make forcing schemes based on factory approach
        amrex::Print() << "Valid forcing_scheme options: direct indirect Gaussian_process" << std::endl;
        amrex::Abort("Unrecognized forcing_scheme");
    }
}

void ABLMesoscaleForcing::setTransitionWeighting()
{
    amrex::Real zmin = m_mesh.Geom(0).ProbLo(m_axis);
    amrex::Real zmax = m_mesh.Geom(0).ProbHi(m_axis);
    m_weighting_heights = {
        zmin, m_transition_height, m_transition_height + m_transition_thickness,
        zmax};
    m_weighting_values = {1.0, 1.0, 0.0, 0.0};
    amrex::Print() << "setting new weighting profile" << std::endl;
    for (int i = 0; i < m_weighting_heights.size(); ++i) {
        amrex::Print() << "  " << m_weighting_heights[i] << " "
                       << m_weighting_values[i] << std::endl;
    }
}

void ABLMesoscaleForcing::updateWeights()
{
    amrex::Print() << "Updating weights" << std::endl;
    for (int i = 0; i < m_nht; ++i) {
        m_W[i] =
            interp::linear(m_weighting_heights, m_weighting_values, m_zht[i]);
        // amrex::Print() << "  " << m_zht[i] << " " << m_W[i] << std::endl;
    }
}

// TODO: generalize the initialization
void ABLMesoscaleForcing::indirectForcingInit()
{
    if (m_W.empty()) {
        // Will be here for:
        // - full profile assimilation
        // - partial profile assim w/ constant transition height
        // - partial profile assim w/ variable transition height (1st step only)
        amrex::Print() << "Initializing indirect forcing" << std::endl;
        m_W.resize(m_nht);
        updateWeights();
    } else if (amrex::toLower(m_forcing_transition) != "none") {
        // Will be here for:
        // - partial profile assim w/ variable transition height
        amrex::Print() << "Reinitializing indirect forcing" << std::endl;
        updateWeights();
    } else {
        amrex::Print() << "Should not be reinitializing indirect forcing!"
                       << std::endl;
        return;
    }

    amrex::Array2D<amrex::Real, 0, 3, 0, 3> zTz;

    // Generate the matrix Z^T W Z
    for (int irow = 0; irow < 4; irow++) {
        for (int icol = 0; icol < 4; icol++) {

            zTz(irow, icol) = 0.0;

            for (int iht = 0; iht < m_nht; iht++) {
                zTz(irow, icol) =
                    zTz(irow, icol) +
                    m_W[iht] *
                        std::pow(m_zht[iht] * m_scaleFact, (icol + irow));
            }
            // amrex::Print()<< "Z^T W Z ["<<irow<<","<<icol<<"] : " <<
            // zTz(irow,icol) << std::endl;
        }
    }
    // Invert the matrix Z^T W Z
    invertMat(zTz, m_im_zTz);
}

void ABLMesoscaleForcing::GP_forcingInit()
{
    GP_updateSigma11Packed();
    GP_updateSigma12();
}

void ABLMesoscaleForcing::invertMat(
    const amrex::Array2D<amrex::Real, 0, 3, 0, 3>& m,
    amrex::Array2D<amrex::Real, 0, 3, 0, 3>& im)
{
    amrex::Real A2323 = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2);
    amrex::Real A1323 = m(2, 1) * m(3, 3) - m(2, 3) * m(3, 1);
    amrex::Real A1223 = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1);
    amrex::Real A0323 = m(2, 0) * m(3, 3) - m(2, 3) * m(3, 0);
    amrex::Real A0223 = m(2, 0) * m(3, 2) - m(2, 2) * m(3, 0);
    amrex::Real A0123 = m(2, 0) * m(3, 1) - m(2, 1) * m(3, 0);
    amrex::Real A2313 = m(1, 2) * m(3, 3) - m(1, 3) * m(3, 2);
    amrex::Real A1313 = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1);
    amrex::Real A1213 = m(1, 1) * m(3, 2) - m(1, 2) * m(3, 1);
    amrex::Real A2312 = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2);
    amrex::Real A1312 = m(1, 1) * m(2, 3) - m(1, 3) * m(2, 1);
    amrex::Real A1212 = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1);
    amrex::Real A0313 = m(1, 0) * m(3, 3) - m(1, 3) * m(3, 0);
    amrex::Real A0213 = m(1, 0) * m(3, 2) - m(1, 2) * m(3, 0);
    amrex::Real A0312 = m(1, 0) * m(2, 3) - m(1, 3) * m(2, 0);
    amrex::Real A0212 = m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0);
    amrex::Real A0113 = m(1, 0) * m(3, 1) - m(1, 1) * m(3, 0);
    amrex::Real A0112 = m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0);

    amrex::Real det =
        m(0, 0) * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223) -
        m(0, 1) * (m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223) +
        m(0, 2) * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123) -
        m(0, 3) * (m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    det = 1.0 / det;

    im(0, 0) = det * (m(1, 1) * A2323 - m(1, 2) * A1323 + m(1, 3) * A1223);
    im(0, 1) = det * -(m(0, 1) * A2323 - m(0, 2) * A1323 + m(0, 3) * A1223);
    im(0, 2) = det * (m(0, 1) * A2313 - m(0, 2) * A1313 + m(0, 3) * A1213);
    im(0, 3) = det * -(m(0, 1) * A2312 - m(0, 2) * A1312 + m(0, 3) * A1212);
    im(1, 0) = det * -(m(1, 0) * A2323 - m(1, 2) * A0323 + m(1, 3) * A0223);
    im(1, 1) = det * (m(0, 0) * A2323 - m(0, 2) * A0323 + m(0, 3) * A0223);
    im(1, 2) = det * -(m(0, 0) * A2313 - m(0, 2) * A0313 + m(0, 3) * A0213);
    im(1, 3) = det * (m(0, 0) * A2312 - m(0, 2) * A0312 + m(0, 3) * A0212);
    im(2, 0) = det * (m(1, 0) * A1323 - m(1, 1) * A0323 + m(1, 3) * A0123);
    im(2, 1) = det * -(m(0, 0) * A1323 - m(0, 1) * A0323 + m(0, 3) * A0123);
    im(2, 2) = det * (m(0, 0) * A1313 - m(0, 1) * A0313 + m(0, 3) * A0113);
    im(2, 3) = det * -(m(0, 0) * A1312 - m(0, 1) * A0312 + m(0, 3) * A0112);
    im(3, 0) = det * -(m(1, 0) * A1223 - m(1, 1) * A0223 + m(1, 2) * A0123);
    im(3, 1) = det * (m(0, 0) * A1223 - m(0, 1) * A0223 + m(0, 2) * A0123);
    im(3, 2) = det * -(m(0, 0) * A1213 - m(0, 1) * A0213 + m(0, 2) * A0113);
    im(3, 3) = det * (m(0, 0) * A1212 - m(0, 1) * A0212 + m(0, 2) * A0112);
}

void ABLMesoscaleForcing::constantForcingTransition(
    amrex::Vector<amrex::Real>& error)
{
    // based on SOWFA6/src/ABLForcing/drivingForce/drivingForce.C

    // find indices
    int hLevelBlend0 = -1;
    int hLevelBlend1 = -1;
    int hLevelBlendMax = -1;
    for (int iht = 0; iht < m_nht; iht++) {
        if ((hLevelBlend1 < 0) && (m_zht[iht] >= m_transition_height)) {
            hLevelBlend1 = iht;
            hLevelBlend0 = iht - 1;
        } else if (m_zht[iht] >= m_transition_height + m_transition_thickness) {
            hLevelBlendMax = iht;
            break;
        }
    }

    // error checking
    if (hLevelBlend1 < 0) {
        amrex::Print() << "Note: Did not find bottom of transition layer"
                       << std::endl;
        hLevelBlend0 = m_nht - 1;
        hLevelBlend1 = m_nht - 1;
        hLevelBlendMax = m_nht - 1;
    } else if (hLevelBlendMax < 0) {
        amrex::Print() << "Note: Did not find top of transition layer"
                       << std::endl;
        hLevelBlendMax = m_nht - 1;
    }

    amrex::Print() << "Forcing transition to constant"
                   << " from " << m_zht[hLevelBlend1] << " to "
                   << m_zht[hLevelBlendMax] << std::endl;

    // calculate initial slope
    amrex::Real slope0 = (error[hLevelBlend1] - error[hLevelBlend0]) /
                         (m_zht[hLevelBlend1] - m_zht[hLevelBlend0]);
    amrex::Real dslope =
        -slope0 / (m_zht[hLevelBlendMax] - m_zht[hLevelBlend1]);

    // march from hLevelBlend1 (z >= m_transition_height)
    // up to hLevelBlendMax (z >= m_transition_height + m_transition_thickness)
    // as the slope decreases linearly to 0
    for (int iht = hLevelBlend1; iht <= hLevelBlendMax; iht++) {
        amrex::Real slope =
            slope0 + dslope * (m_zht[iht] - m_zht[hLevelBlend1]);
        error[iht] = error[iht - 1] + slope * (m_zht[iht] - m_zht[iht - 1]);
    }

    // set the remaining levels above hLevelBlendMax to the last value
    for (int iht = hLevelBlendMax + 1; iht < m_nht; iht++) {
        error[iht] = error[hLevelBlendMax];
    }
}

void ABLMesoscaleForcing::blendForcings(
    const amrex::Vector<amrex::Real> lower, // W=1
    const amrex::Vector<amrex::Real> upper, // W=0
    amrex::Vector<amrex::Real>& error)
{
    amrex::Print() << "Blending forcings" << std::endl;
    for (int iht = 0; iht < m_nht; iht++) {
        error[iht] = m_W[iht] * lower[iht] + (1.0 - m_W[iht]) * upper[iht];
    }
}

// TODO: covariance function is hard-coded for now; can generalize the
// following updateSigma* functions

void ABLMesoscaleForcing::GP_updateSigma11Packed()
{
    const amrex::Vector<amrex::Real>& x1 = m_meso_file.meso_heights();
    const int n1 = m_meso_file.nheights();
    amrex::Print() << "[GPIPA] Updating Sigma_11 for " << m_identifier
        << std::endl;
    Sigma11.resize(n1*(n1+1)/2);
    int lin_idx = 0;
    double base_err = m_sigma_noise * m_sigma_noise;
    // From https://netlib.org/lapack/double/dpptrf.f
    //   if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j
    for (int j=0; j < n1; j++) {
        for (int i=0; i < j+1; i++) {
            Sigma11[lin_idx] = GP_covar_func(x1[i], x1[j]);
            if (i==j) Sigma11[lin_idx] += base_err;
            lin_idx++;
        }
    }

    // Add specific error

    // Compute Cholesky factorization
    //int info = LAPACKE_dpptrf(LAPACK_COL_MAJOR, 'U', n1, &Sigma11[0]);
    int info;
    const char uplo = 'U';
    dpptrf_(&uplo, &n1, &Sigma11[0], &info);
    if (info != 0) {
        amrex::Print() << "[GPIPA] WARNING: Cholesky factorization returned "
            << info << std::endl;
    }
}

void ABLMesoscaleForcing::GP_updateSigma12()
{
    const amrex::Vector<amrex::Real>& x1 = m_meso_file.meso_heights();
    const int n1 = m_meso_file.nheights();
    amrex::Print() << "[GPIPA] Updating Sigma_12 for " << m_identifier
        << std::endl;
    Sigma12.resize(n1*m_nht);
    int lin_idx = 0;
    // From https://netlib.org/lapack/double/dpptrf.f
    //   if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j
    for (int j=0; j < m_nht; j++) {
        for (int i=0; i < n1; i++) {
            Sigma12[lin_idx++] = GP_covar_func(x1[i], m_zht[j]);
        }
    }
}

amrex::Vector<amrex::Real> ABLMesoscaleForcing::GP_posteriorMean(
    amrex::Vector<amrex::Real> y1)
{
    // TODO: get y1 at meso_heights() if x1 and x2 are different
    const int n1 = m_meso_file.nheights();

    // get mean of input y
    double ymean = 0;
    for (int i=0; i < y1.size(); i++) {
        ymean += y1[i];
    }
    ymean /= y1.size();

    // center input values around 0, get max abs
    double ynorm = 0;
    for (int i=0; i < y1.size(); i++) {
        y1[i] -= ymean;
        if (std::abs(y1[i]) > ynorm) {
            ynorm = std::abs(y1[i]);
        }
    }

    // normalize values
    for (int i=0; i < y1.size(); i++) {
        y1[i] /= ynorm;
    }

    // do regression 
    int info;
    const char uplo = 'U';
    std::vector<double> reg(Sigma12);
    dpptrs_(&uplo, &n1, &m_nht, &Sigma11[0], &reg[0], &n1, &info);
    if (info != 0) {
        amrex::Print() << "[GPIPA] WARNING: SPD solve returned " << info
            << std::endl;
    }

    // evaluate, rescale, and shift values
    // TODO: calculate y2 in place if x1 and x2 are the same
    amrex::Vector<amrex::Real> y2(m_nht, 1.0);
    const char trans = 'T';
    const int incx = 1;
    const int incy = 1;
    // calculate y := alpha*A**T*x + beta*y, where
    //   A == reg
    //   x == y1
    //   alpha == ynorm
    //   beta == ymean
    dgemv_(&trans, &n1, &m_nht, &ynorm, &reg[0], &n1, &y1[0], &incx, &ymean,
        &y2[0], &incy);

    return y2;
}

} // namespace amr_wind

