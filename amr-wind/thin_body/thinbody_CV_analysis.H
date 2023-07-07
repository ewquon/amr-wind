#ifndef THINBODY_CV_ANALYSIS_H
#define THINBODY_CV_ANALYSIS_H

namespace amr_wind::thinbody {

static void sum_momentum_fluxes(
    const FieldRepo& repo,
    int ndim,
    int lev,
    amrex::MFIter& mfi,
    const amrex::Box& bx,
    const amrex::Geometry& geom,
    const amrex::Array4<const amrex::Real>& flux_x_arr,
    const amrex::Array4<const amrex::Real>& flux_y_arr,
    const amrex::Array4<const amrex::Real>& flux_z_arr,
    const amrex::Array4<const amrex::Real>& gp_arr)
{
    auto& flux_sum = repo.get_field("mom_sum");
    auto& xface_mask = repo.get_int_field("mom_xface_mask");
    auto& yface_mask = repo.get_int_field("mom_yface_mask");
    auto& zface_mask = repo.get_int_field("mom_zface_mask");
    const auto& flux_sum_arr = flux_sum(lev).array(mfi);
    const auto& xf_mask_arr = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr = zface_mask(lev).const_array(mfi);

    const auto dx = geom.CellSizeArray();
    amrex::Real cellvol = dx[0];
    for (int i=1; i < ndim; i++)
        cellvol *= dx[i];

    amrex::ParallelFor(bx, ndim,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            flux_sum_arr(i, j, k, n) =
                + xf_mask_arr(i+1,j  ,k  , n) * flux_x_arr(i+1,j  ,k  )
                - xf_mask_arr(i  ,j  ,k  , n) * flux_x_arr(i  ,j  ,k  )
                + yf_mask_arr(i  ,j+1,k  , n) * flux_y_arr(i  ,j+1,k  )
                - yf_mask_arr(i  ,j  ,k  , n) * flux_y_arr(i  ,j  ,k  )
                + zf_mask_arr(i  ,j  ,k+1, n) * flux_z_arr(i  ,j  ,k+1)
                - zf_mask_arr(i  ,j  ,k  , n) * flux_z_arr(i  ,j  ,k  )
                + gp_arr(i,j,k,n)*cellvol;
        });
}
} // namespace amr_wind::thinbody

#endif
