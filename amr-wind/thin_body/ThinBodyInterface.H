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
    const auto& mask_cell = repo.get_int_field("mask_cell");
    auto& vel_corr = repo.get_field("norm_vel_corr");

    const auto& xf_mask_arr = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr = zface_mask(lev).const_array(mfi);
    const auto& mc_arr = mask_cell(lev).const_array(mfi);

    // We cannot modify vel_arr in place because it needs to be used twice
    // (once per adjacent cell) for each zero-flux face to calculate the
    // velocity error in the normal direction.
    // TODO: Create a temporary MultiFab instead?
    const auto& vel_corr_arr = vel_corr(lev).array(mfi);

    //const auto& geom = repo.mesh().Geom();
    //const auto plo = geom[lev].ProbLoArray();
    //const auto dx = geom[lev].CellSizeArray();

    // Note: We should only set vel_arr(i,j,k,n) and not the ghost cells within
    // this ParallelFor, otherwise FillBoundary will not work as expected.
    const auto& bx = mfi.tilebox();
    amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (xf_mask_arr(i,j,k)==0) {
                // cell is on +x side of thin body, surface is at low face
                vel_corr_arr(i,j,k,0) -= 0.5*(vel_arr(i-1,j,k,0)
                                             +vel_arr(i  ,j,k,0));
            }
            else if (xf_mask_arr(i+1,j,k)==0) {
                // cell is on -x side of thin body, surface is at high face
                vel_corr_arr(i,j,k,0) -= 0.5*(vel_arr(i  ,j,k,0)
                                             +vel_arr(i+1,j,k,0));
            }

            if (yf_mask_arr(i,j,k)==0) {
                // cell is on +y side of thin body, surface is at low face
                vel_corr_arr(i,j,k,1) -= 0.5*(vel_arr(i,j-1,k,1)
                                             +vel_arr(i,j  ,k,1));
            }
            else if (yf_mask_arr(i,j+1,k)==0) {
                // cell is on -y side of thin body, surface is at high face
                vel_corr_arr(i,j,k,1) -= 0.5*(vel_arr(i,j  ,k,1)
                                             +vel_arr(i,j+1,k,1));
            }

            if (zf_mask_arr(i,j,k)==0) {
                // cell is on +z side of thin body, surface is at low face
                vel_corr_arr(i,j,k,2) -= 0.5*(vel_arr(i,j,k-1,2)
                                             +vel_arr(i,j,k  ,2));
                //amrex::AllPrint() << "Applying interface conditions (lev="<<lev<<" at k)"
                //    << " at ("  << plo[0]+(i+0.5)*dx[0] <<","
                //                << plo[1]+(j+0.5)*dx[1] <<","
                //                << plo[2]+      k*dx[2] <<")" // zface is staggered in k
                //    << " to correct wface=" << vel_corr_arr(i,j,k,2)
                //    << std::endl;
            }
            else if (zf_mask_arr(i,j,k+1)==0) {
                // cell is on -z side of thin body, surface is at high face
                vel_corr_arr(i,j,k,2) -= 0.5*(vel_arr(i,j,k  ,2)
                                             +vel_arr(i,j,k+1,2));
                //amrex::AllPrint() << "Applying interface conditions (lev="<<lev<<" at k+1)"
                //    << " at ("  << plo[0]+(i+0.5)*dx[0] <<","<<
                //                   plo[1]+(j+0.5)*dx[1] <<","<<
                //                   plo[2]+(k+1  )*dx[2] <<")" // zface is staggered in k
                //    << " to correct wface=" << vel_corr_arr(i,j,k,2)
                //    << std::endl;
            }
        });
}

/* Sets velocities on masked faces to 0 */
static inline void zero_surface_face_vel(
    const FieldRepo& repo,
    int lev,
    amrex::MFIter& mfi,
    const amrex::Array4<amrex::Real>& umac_arr,
    const amrex::Array4<amrex::Real>& vmac_arr,
    const amrex::Array4<amrex::Real>& wmac_arr)
{
    const auto& xface_mask = repo.get_int_field("mom_xface_mask");
    const auto& yface_mask = repo.get_int_field("mom_yface_mask");
    const auto& zface_mask = repo.get_int_field("mom_zface_mask");

    const auto& xf_mask_arr = xface_mask(lev).const_array(mfi);
    const auto& yf_mask_arr = yface_mask(lev).const_array(mfi);
    const auto& zf_mask_arr = zface_mask(lev).const_array(mfi);

    const auto& xbx = mfi.tilebox(amrex::IntVect(1,0,0));
    const auto& ybx = mfi.tilebox(amrex::IntVect(0,1,0));
    const auto& zbx = mfi.tilebox(amrex::IntVect(0,0,1));

    amrex::ParallelFor(xbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            umac_arr(i,j,k) *= xf_mask_arr(i,j,k);
        });
    amrex::ParallelFor(ybx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            vmac_arr(i,j,k) *= yf_mask_arr(i,j,k);
        });
    amrex::ParallelFor(zbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            //if ((j==0) && (zf_mask_arr(i,j,k) == 0)) {
            //    amrex::AllPrint()<<"zeroed wmac"<<amrex::IntVect(i,j,k)<<"="<<wmac_arr(i,j,k)<<std::endl;
            //}
            wmac_arr(i,j,k) *= zf_mask_arr(i,j,k);
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
