#ifndef THINBODY_CV_ANALYSIS_H
#define THINBODY_CV_ANALYSIS_H

namespace amr_wind::thinbody {

static void sum_momentum_fluxes(
    const FieldRepo& repo,
    int lev,
    amrex::MFIter& mfi,
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& flux_x_arr,
    const amrex::Array4<const amrex::Real>& flux_y_arr,
    const amrex::Array4<const amrex::Real>& flux_z_arr)
{
    auto& flux_sum = repo.get_field("mom_flux_sum");
    auto& xface_mask = repo.get_int_field("mom_xface_mask");
    auto& yface_mask = repo.get_int_field("mom_yface_mask");
    auto& zface_mask = repo.get_int_field("mom_zface_mask");
    const auto& flux_sum_arr = flux_sum(lev).array(mfi);
    const auto& xf_mask_arr = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr = zface_mask(lev).const_array(mfi);

    amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            flux_sum_arr(i, j, k) =
                + xf_mask_arr(i+1,j  ,k  )*flux_x_arr(i+1,j  ,k  )
                - xf_mask_arr(i  ,j  ,k  )*flux_x_arr(i  ,j  ,k  )
                + yf_mask_arr(i  ,j+1,k  )*flux_y_arr(i  ,j+1,k  )
                - yf_mask_arr(i  ,j  ,k  )*flux_y_arr(i  ,j  ,k  )
                + zf_mask_arr(i  ,j  ,k+1)*flux_z_arr(i  ,j  ,k+1)
                - zf_mask_arr(i  ,j  ,k  )*flux_z_arr(i  ,j  ,k  );
        });
}
} // namespace amr_wind::thinbody

#endif
