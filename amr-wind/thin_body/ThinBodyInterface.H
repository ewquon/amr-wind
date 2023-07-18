#ifndef THINBODY_INTERFACE_H
#define THINBODY_INTERFACE_H

namespace amr_wind::thinbody {

/* This linearly interpolates the current face velocity value and subtracts it
 * from adjacent cells such that the face value is zero. Surface faces were
 * previously identified by ActuatorContainer::mark_surface_faces(). This is
 * currently implemented for the ThinBody actuator.
 */
static inline void zero_adjacent_cell_vel(
    const FieldRepo& repo,
    int lev,
    amrex::MFIter& mfi,
    const amrex::Array4<amrex::Real>& vel_arr)
{
    const auto& xface_mask = repo.get_int_field("mom_xface_mask");
    const auto& yface_mask = repo.get_int_field("mom_yface_mask");
    const auto& zface_mask = repo.get_int_field("mom_zface_mask");

    const auto& xf_mask_arr  = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr  = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr  = zface_mask(lev).const_array(mfi);

    //const auto& geom = repo.mesh().Geom();
    //const auto plo = geom[lev].ProbLoArray();
    //const auto dx = geom[lev].CellSizeArray();

    const auto& bx = mfi.tilebox();

    amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (xf_mask_arr(i,j,k)==0) {
                amrex::Real uerr = 0.5*(vel_arr(i,j,k,0)+vel_arr(i-1,j  ,k  ,0));
                vel_arr(i  ,j  ,k  ,0) -= uerr;
                vel_arr(i-1,j  ,k  ,0) -= uerr;
            }
            if (yf_mask_arr(i,j,k)==0) {
                amrex::Real verr = 0.5*(vel_arr(i,j,k,1)+vel_arr(i  ,j-1,k  ,1));
                vel_arr(i  ,j  ,k  ,1) -= verr;
                vel_arr(i  ,j-1,k  ,1) -= verr;
            }
            if (zf_mask_arr(i,j,k)==0) {
                amrex::Real werr = 0.5*(vel_arr(i,j,k,2)+vel_arr(i  ,j  ,k-1,2));
                //amrex::AllPrint() << "Applying interface conditions (lev="<<lev<<")"
                //    << " at ("  << plo[0]+(i+0.5)*dx[0] <<","<< plo[1]+(j+0.5)*dx[1] <<","<< plo[2]+      k*dx[2] <<")" // zface is staggered in k
                //    << " to correct wface=" << werr
                //    << std::endl;
                vel_arr(i  ,j  ,k  ,2) -= werr;
                vel_arr(i  ,j  ,k-1,2) -= werr;
            }
        });
}

static inline void zero_surface_face_vel(
    const FieldRepo& repo,
    int lev,
    amrex::MFIter& mfi,
    const amrex::Array4<amrex::Real>& vel_x_arr,
    const amrex::Array4<amrex::Real>& vel_y_arr,
    const amrex::Array4<amrex::Real>& vel_z_arr)
{
    const auto& xface_mask = repo.get_int_field("mom_xface_mask");
    const auto& yface_mask = repo.get_int_field("mom_yface_mask");
    const auto& zface_mask = repo.get_int_field("mom_zface_mask");

    const auto& xf_mask_arr  = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr  = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr  = zface_mask(lev).const_array(mfi);

    const auto& xbx = mfi.tilebox(amrex::IntVect(1,0,0));
    const auto& ybx = mfi.tilebox(amrex::IntVect(0,1,0));
    const auto& zbx = mfi.tilebox(amrex::IntVect(0,0,1));

    amrex::ParallelFor(xbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            vel_x_arr(i,j,k) *= xf_mask_arr(i,j,k);
        });
    amrex::ParallelFor(ybx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            vel_y_arr(i,j,k) *= yf_mask_arr(i,j,k);
        });
    amrex::ParallelFor(zbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            vel_z_arr(i,j,k) *= zf_mask_arr(i,j,k);
        });
}

/*
static void zero_surface_fluxes(
    const FieldRepo& repo,
    int ndim,
    int lev,
    amrex::MFIter& mfi,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& flux_x_arr,
    const amrex::Array4<amrex::Real>& flux_y_arr,
    const amrex::Array4<amrex::Real>& flux_z_arr)
{
    const auto& xface_mask = repo.get_int_field("mom_xface_mask");
    const auto& yface_mask = repo.get_int_field("mom_yface_mask");
    const auto& zface_mask = repo.get_int_field("mom_zface_mask");
    const auto& grad_p     = repo.get_field("gp"); // cell-centered
    auto& flux_sum         = repo.get_field("mom_sum");

    const auto& xf_mask_arr  = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr  = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr  = zface_mask(lev).const_array(mfi);
    const auto& gp_arr       = grad_p(lev).const_array(mfi);
    const auto& flux_sum_arr = flux_sum(lev).array(mfi);

    const auto&  bx = mfi.tilebox();
    const auto& xbx = mfi.tilebox(amrex::IntVect(1,0,0));
    const auto& ybx = mfi.tilebox(amrex::IntVect(0,1,0));
    const auto& zbx = mfi.tilebox(amrex::IntVect(0,0,1));

    const auto dx = geom.CellSizeArray();
    amrex::Real cellvol = dx[0];
    for (int i=1; i < ndim; i++)
        cellvol *= dx[i];

    amrex::ParallelFor(xbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            flux_x_arr(i,j,k,0) *= xf_mask_arr(i,j,k);
            flux_x_arr(i,j,k,1) *= xf_mask_arr(i,j,k);
            flux_x_arr(i,j,k,2) *= xf_mask_arr(i,j,k);
        });
    amrex::ParallelFor(ybx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            flux_y_arr(i,j,k,0) *= yf_mask_arr(i,j,k);
            flux_y_arr(i,j,k,1) *= yf_mask_arr(i,j,k);
            flux_y_arr(i,j,k,2) *= yf_mask_arr(i,j,k);
        });
    amrex::ParallelFor(zbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            flux_z_arr(i,j,k,0) *= zf_mask_arr(i,j,k);
            flux_z_arr(i,j,k,1) *= zf_mask_arr(i,j,k);
            flux_z_arr(i,j,k,2) *= zf_mask_arr(i,j,k);
        });

    amrex::ParallelFor(bx, ndim,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            flux_sum_arr(i, j, k, n) =
                + flux_x_arr(i+1,j  ,k  ,n) - flux_x_arr(i,j,k,n)
                + flux_y_arr(i  ,j+1,k  ,n) - flux_y_arr(i,j,k,n)
                + flux_z_arr(i  ,j  ,k+1,n) - flux_z_arr(i,j,k,n)
                + gp_arr(i,j,k,n)*cellvol;
        });
}
*/
} // namespace amr_wind::thinbody

#endif
