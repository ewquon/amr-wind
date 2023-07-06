#ifndef THINBODY_CV_ANALYSIS_H
#define THINBODY_CV_ANALYSIS_H

namespace amr_wind::thinbody {

static void sum_momentum_fluxes(
    const FieldRepo& repo,
    int lev,
    const amrex::Box& bx,
    amrex::MFIter& mfi,
    const amrex::Array4<const amrex::Real>& flux_x_arr,
    const amrex::Array4<const amrex::Real>& flux_y_arr,
    const amrex::Array4<const amrex::Real>& flux_z_arr
    //ScratchField& flux_x,
    //ScratchField& flux_y,
    //ScratchField& flux_z
    )
{
    auto& flux_sum = repo.get_field("mom_flux_sum");
    const auto& flux_sum_arr = flux_sum(lev).array(mfi);

    amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            flux_sum_arr(i, j, k) =
                + flux_x_arr(i+1,j  ,k  ) - flux_x_arr(i,j,k)
                + flux_y_arr(i  ,j+1,k  ) - flux_y_arr(i,j,k)
                + flux_z_arr(i  ,j  ,k+1) - flux_z_arr(i,j,k);
        });
}
} // namespace amr_wind::thinbody

#endif
