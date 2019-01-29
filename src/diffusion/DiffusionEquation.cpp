#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <DiffusionEquation.H>
#include <diffusion_F.H>
#include <constants.H>

using namespace amrex;

//
// Constructor: 
// We set up everything which doesn't change between timesteps here
//
DiffusionEquation::DiffusionEquation(AmrCore* _amrcore,
                                     Vector<std::unique_ptr<EBFArrayBoxFactory>>* _ebfactory,
                                     Vector<std::unique_ptr<IArrayBox>>& bc_ilo, 
                                     Vector<std::unique_ptr<IArrayBox>>& bc_ihi, 
                                     Vector<std::unique_ptr<IArrayBox>>& bc_jlo, 
                                     Vector<std::unique_ptr<IArrayBox>>& bc_jhi, 
                                     Vector<std::unique_ptr<IArrayBox>>& bc_klo, 
                                     Vector<std::unique_ptr<IArrayBox>>& bc_khi,
                                     int _nghost)
{
    // Get inputs from ParmParse
	readParameters();
    
    if(verbose > 0)
    {
        amrex::Print() << "Constructing DiffusionEquation class" << std::endl;
    }

    // Set AmrCore and ebfactory based on input, fetch some data needed in constructor
    amrcore = _amrcore; 
    ebfactory = _ebfactory;
    nghost = _nghost;
    Vector<Geometry> geom = amrcore->Geom();
    Vector<BoxArray> grids = amrcore->boxArray();
    Vector<DistributionMapping> dmap = amrcore->DistributionMap();
    int max_level = amrcore->maxLevel();
    
    // Whole domain
    Box domain(geom[0].Domain());

    // The boundary conditions need only be set at level 0
    set_diff_bc(bc_lo, bc_hi,
                domain.loVect(), domain.hiVect(),
				&nghost,
                bc_ilo[0]->dataPtr(), bc_ihi[0]->dataPtr(),
                bc_jlo[0]->dataPtr(), bc_jhi[0]->dataPtr(),
                bc_klo[0]->dataPtr(), bc_khi[0]->dataPtr());

    // Resize and reset data
    b.resize(max_level + 1);
    phi.resize(max_level + 1);
    rhs.resize(max_level + 1);
    for(int lev = 0; lev <= max_level; lev++)
    {
        b[lev].resize(3);
        for(int dir = 0; dir < 3; dir++)
        {
            BoxArray edge_ba = grids[lev];
            edge_ba.surroundingNodes(dir);
            b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost));
        }
        phi[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
        rhs[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nghost));
    }

	// Define the matrix.
	LPInfo info;
    matrix.define(geom, grids, dmap, info, GetVecOfConstPtrs(*ebfactory));

    // It is essential that we set MaxOrder to 2 if we want to use the standard
    // phi(i)-phi(i-1) approximation for the gradient at Dirichlet boundaries.
    // The solver's default order is 3 and this uses three points for the gradient.
	matrix.setMaxOrder(2);

	// LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
	matrix.setDomainBC({(LinOpBCType) bc_lo[0], (LinOpBCType) bc_lo[1], (LinOpBCType) bc_lo[2]},
					   {(LinOpBCType) bc_hi[0], (LinOpBCType) bc_hi[1], (LinOpBCType) bc_hi[2]});
}

DiffusionEquation::~DiffusionEquation()
{
}

void DiffusionEquation::readParameters()
{
    ParmParse pp("diffusion");

    pp.query("verbose", verbose);
    pp.query("mg_verbose", mg_verbose);
    pp.query("mg_cg_verbose", mg_cg_verbose);
    pp.query("mg_max_iter", mg_max_iter);
    pp.query("mg_cg_maxiter", mg_cg_maxiter);
    pp.query("mg_max_fmg_iter", mg_max_fmg_iter);
    pp.query("mg_rtol", mg_rtol);
    pp.query("mg_atol", mg_atol);
    pp.query( "bottom_solver_type", bottom_solver_type);
}

void DiffusionEquation::updateInternals(AmrCore* amrcore_in, 
                                        Vector<std::unique_ptr<EBFArrayBoxFactory>>* ebfactory_in)
{
    // This must be implemented when we want dynamic meshing
    //
    amrex::Print() << "ERROR: DiffusionEquation::updateInternals() not yet implemented" << std::endl;
    amrex::Abort(); 
}

//
// Solve the matrix equation
//
void DiffusionEquation::solve(Vector<std::unique_ptr<MultiFab>>& vel, 
                              const Vector<std::unique_ptr<MultiFab>>& ro, 
                              const Vector<std::unique_ptr<MultiFab>>& eta, 
                              Real dt)
{
	BL_PROFILE("DiffusionEquation::solve");

    // Update the coefficients of the matrix going into the solve based on the current state of the
    // simulation. Recall that the relevant matrix is 
    //
    //      alpha a - beta div ( b grad )   <--->   rho - dt div ( eta grad )
    //
    // So the constants and variable coefficients are: 
    //
    //      alpha: 1
    //      beta: dt
    //      a: ro
    //      b: eta
    
    // Set alpha and beta
    matrix.setScalars(1.0, dt);

    // Set the variable coefficients (on faces) to equal the apparent viscosity
    updateCoefficients(eta);

    for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
    {
        // This sets the spatially varying a coefficients
        matrix.setACoeffs(lev, (*ro[lev]));

        // Copy the spatially varying b coefficients into the proper data strutcure
        Vector<const MultiFab*> tmp = GetVecOfConstPtrs(b[lev]);
        std::array<MultiFab const*, AMREX_SPACEDIM> b_tmp;
        for(int dir = 0; dir < 3; dir++)
        {
            b_tmp[dir] = tmp[dir];
        }
        // Set the coefficients
        matrix.setBCoeffs(lev, b_tmp);
    }

    // Loop over the velocity components
    for(int dir = 0; dir < 3; dir++)
    {
        if(verbose > 0)
        {
            amrex::Print() << "Diffusing velocity component " << dir << std::endl;
        }

        for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
        {
            // Set the right hand side to equal ro 
            rhs[lev]->copy(*ro[lev], 0, 0, 1, nghost, nghost);

            // Multiply rhs by vel(dir) to get momentum
            // Note that vel holds the updated velocity:
            //
            //      u_old + dt ( - u grad u + div tau / ro - grad p / rho + gravity )
            //
            MultiFab::Multiply((*rhs[lev]), (*vel[lev]), dir, 0, 1, nghost);

            // By this point we must have filled the Dirichlet values of phi stored in ghost cells
            phi[lev]->copy(*vel[lev], dir, 0, 1, nghost, nghost);
            matrix.setLevelBC(lev, GetVecOfConstPtrs(phi)[lev]);
        }

        MLMG solver(matrix);
        setSolverSettings(solver);
        solver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), mg_rtol, mg_atol);

        for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
        {
            phi[lev]->FillBoundary(amrcore->Geom(lev).periodicity());
            vel[lev]->copy(*phi[lev], 0, dir, 1, nghost, nghost);
        }
    }
}

//
// Computes b = eta at the faces of the scalar cells
//
void DiffusionEquation::updateCoefficients(const Vector<std::unique_ptr<MultiFab>>& eta)
{
	BL_PROFILE("DiffusionEquation::updateCoefficients");

    if(verbose > 0)
    {
        amrex::Print() << "Updating DiffusionEquation coefficients" << std::endl;
    }

	// Directions
	int xdir = 1;
	int ydir = 2;
	int zdir = 3;

    Vector<Geometry> geom = amrcore->Geom();
    for(int lev = 0; lev <= amrcore->finestLevel(); lev++)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(*eta[lev], true); mfi.isValid(); ++mfi)
        {
            // Tileboxes for staggered components
            Box ubx = mfi.tilebox(e_x);
            Box vbx = mfi.tilebox(e_y);
            Box wbx = mfi.tilebox(e_z);

            // X direction
            compute_bcoeff_diff(BL_TO_FORTRAN_BOX(ubx),
                                BL_TO_FORTRAN_ANYD((*(b[lev][0]))[mfi]),
                                BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                &xdir);

            // Y direction
            compute_bcoeff_diff(BL_TO_FORTRAN_BOX(vbx),
                                BL_TO_FORTRAN_ANYD((*(b[lev][1]))[mfi]),
                                BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                &ydir);

            // Z direction
            compute_bcoeff_diff(BL_TO_FORTRAN_BOX(wbx),
                                BL_TO_FORTRAN_ANYD((*(b[lev][2]))[mfi]),
                                BL_TO_FORTRAN_ANYD((*eta[lev])[mfi]),
                                &zdir);
        }
        // TODO: do we need these? 
        b[lev][0]->FillBoundary(geom[lev].periodicity());
        b[lev][1]->FillBoundary(geom[lev].periodicity());
        b[lev][2]->FillBoundary(geom[lev].periodicity());
    }
}

// 
// Set the user-supplied settings for the MLMG solver
// (this must be done every time step, since MLMG is created after updating matrix
//
void DiffusionEquation::setSolverSettings(MLMG& solver)
{
    // The default bottom solver is BiCG
    if(bottom_solver_type == "smoother")
    {
       solver.setBottomSolver(MLMG::BottomSolver::smoother);
    }
    else if(bottom_solver_type == "hypre")
    {
       solver.setBottomSolver(MLMG::BottomSolver::hypre);
    }

    // Maximum iterations for MultiGrid / ConjugateGradients
	solver.setMaxIter(mg_max_iter);
	solver.setMaxFmgIter(mg_max_fmg_iter);
	solver.setCGMaxIter(mg_cg_maxiter);

    // Verbosity for MultiGrid / ConjugateGradients
	solver.setVerbose(mg_verbose);
	solver.setCGVerbose(mg_cg_verbose);

	// This ensures that ghost cells of phi are correctly filled when returned from the solver
	solver.setFinalFillBC(true);
}
