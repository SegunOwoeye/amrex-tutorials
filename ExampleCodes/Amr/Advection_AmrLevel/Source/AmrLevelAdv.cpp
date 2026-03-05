#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuMemory.H>

#include "AmrLevelAdv.H"
#include "Adv_F.H"
#include "Kernels.H"

#include "euler_hllc_muscl_2d.hpp"

using namespace amrex;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9;
int      AmrLevelAdv::do_reflux       = 1;

int      AmrLevelAdv::NUM_STATE       = 4;
int      AmrLevelAdv::NUM_GROW        = 2;

int      AmrLevelAdv::max_phierr_lev  = -1;
int      AmrLevelAdv::max_phigrad_lev = -1;
Real     AmrLevelAdv::gamma           = 1.4;
int      AmrLevelAdv::prob_type       = 0;

Vector<Real> AmrLevelAdv::phierr;
Vector<Real> AmrLevelAdv::phigrad;

int AmrLevelAdv::max_rhograd_lev = -1;
amrex::Vector<amrex::Real> AmrLevelAdv::rho_grad;

#ifdef AMREX_PARTICLES
std::unique_ptr<AmrTracerParticleContainer> AmrLevelAdv::TracerPC =  nullptr;
int AmrLevelAdv::do_tracers                       =  0;
#endif


namespace
{
    // [P0] Enforce rho > 0 and p > 0, then recompute E consistently.
    void enforce_positivity(
        amrex::MultiFab& U,
        const amrex::BoxArray& grids,
        const amrex::DistributionMapping& dmap,
        const int nGrow,
        const amrex::Real gamma)
    {
        constexpr amrex::Real rho_floor = 1e-12_rt;
        constexpr amrex::Real p_floor   = 1e-12_rt;

        for (amrex::MFIter mfi(U); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = amrex::grow(mfi.validbox(), nGrow);
            auto const a = U.array(mfi);

            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                {
                    amrex::Real rho  = a(i,j,0,0);
                    amrex::Real momx = a(i,j,0,1);
                    amrex::Real momy = a(i,j,0,2);
                    amrex::Real E    = a(i,j,0,3);

                    rho = amrex::max(rho, rho_floor);

                    const amrex::Real ux = momx / rho;
                    const amrex::Real uy = momy / rho;

                    const amrex::Real kinetic = 0.5_rt * rho * (ux*ux + uy*uy);
                    amrex::Real p_raw = (gamma - 1.0_rt) * (E - kinetic);

                    if (!(p_raw > p_floor)) {
                        p_raw = p_floor;
                        E = p_raw / (gamma - 1.0_rt) + kinetic;
                    }

                    a(i,j,0,0) = rho;
                    a(i,j,0,1) = rho * ux;
                    a(i,j,0,2) = rho * uy;
                    a(i,j,0,3) = E;
                }
            }
        }
    }
}


/**
 * The basic constructor.
 */
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
                          int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

/**
 * The destructor.
 */
AmrLevelAdv::~AmrLevelAdv ()
{
    delete flux_reg;
}

/**
 * Restart from a checkpoint file.
 */
void
AmrLevelAdv::restart (Amr&          papa,
                      std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

/**
 * Write a checkpoint file.
 */
void
AmrLevelAdv::checkPoint (const std::string& dir,
                         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers && level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

/**
 * Write a plotfile to specified directory.
 */
void
AmrLevelAdv::writePlotFile (const std::string& dir,
                             std::ostream&      os,
                            VisMF::How         how)
{

    AmrLevel::writePlotFile (dir,os,how);

#ifdef AMREX_PARTICLES
    if (do_tracers && level == 0) {
      TracerPC->Checkpoint(dir, "Tracer", true);
    }
#endif
}

/**
 * Define data descriptors.
 */
void
AmrLevelAdv::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);


    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
                           &pc_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];

    if (prob_type == 1) {
        // y-discontinuity: periodic in x, outflow in y
        lo_bc[0] = BCType::int_dir;
        hi_bc[0] = BCType::int_dir;
        lo_bc[1] = BCType::foextrap;
        hi_bc[1] = BCType::foextrap;
    } else if (prob_type == 0) {
        // x-discontinuity: outflow in x, periodic in y
        lo_bc[0] = BCType::foextrap;
        hi_bc[0] = BCType::foextrap;
        lo_bc[1] = BCType::int_dir;
        hi_bc[1] = BCType::int_dir;
    } else {
        // oblique or quadrant: outflow everywhere
        lo_bc[0] = BCType::foextrap;
        hi_bc[0] = BCType::foextrap;
        lo_bc[1] = BCType::foextrap;
        hi_bc[1] = BCType::foextrap;
    }

    BCRec bc(lo_bc, hi_bc);

    StateDescriptor::BndryFunc bndryfunc(nullfill);
    bndryfunc.setRunOnGPU(true);  // I promise the bc function will launch gpu kernels.

    //desc_lst.setComponent(Phi_Type, 0, "phi", bc, bndryfunc);
    desc_lst.setComponent(Phi_Type, 0, "density", bc, bndryfunc);
    desc_lst.setComponent(Phi_Type, 1, "momx", bc, bndryfunc);
    desc_lst.setComponent(Phi_Type, 2, "momy", bc, bndryfunc);
    desc_lst.setComponent(Phi_Type, 3, "energy", bc, bndryfunc);
}

/**
 * Cleanup data descriptors at end of run.
 */
void
AmrLevelAdv::variableCleanUp ()
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif

}

/**
 * Initialize grid data at problem start-up.
 */
void
AmrLevelAdv::initData ()
{
    const Real* dx      = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level "
                       << level << std::endl;
    }

    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx  = mfi.tilebox();      // cells we will fill
        const Box& fbx = S_new[mfi].box();   // full FAB bounds (dlo/dhi must match this)

        const int* lo2  = bx.loVect();
        const int* hi2  = bx.hiVect();

        const int* dlo2 = fbx.loVect();
        const int* dhi2 = fbx.hiVect();

        int lo[3];
        int hi[3];
        int dlo[3];
        int dhi[3];

        lo[0]  = lo2[0];
        lo[1]  = lo2[1];
        lo[2]  = 0;

        hi[0]  = hi2[0];
        hi[1]  = hi2[1];
        hi[2]  = 0;

        dlo[0] = dlo2[0];
        dlo[1] = dlo2[1];
        dlo[2] = 0;

        dhi[0] = dhi2[0];
        dhi[1] = dhi2[1];
        dhi[2] = 0;

        Real* state_ptr = S_new[mfi].dataPtr();

        initdata_(
            &level,
            &cur_time,
            lo,
            hi,
            state_ptr,
            dlo,
            dhi,
            dx,
            prob_lo
        );
    }

#ifdef AMREX_PARTICLES
    init_particles();
#endif

    if (verbose) {
        amrex::Print() << "Done initializing the level "
                       << level << " data " << std::endl;
    }
    enforce_positivity(S_new, grids, dmap, 0, gamma);
}


/**
 * Initialize data on this level from another AmrLevelAdv (during regrid).
 */
void
AmrLevelAdv::init (AmrLevel &old)
{
    AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;

    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);
    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
    enforce_positivity(S_new, grids, dmap, 0, gamma);
}

/**
 * Initialize data on this level after regridding if old level did not previously exist
 */
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
    enforce_positivity(S_new, grids, dmap, 0, gamma);
}


/*
 *Advance grids at this level in time.
 */
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  /*ncycle*/)
{
    // [1] Swap time levels
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);

    // [2] Allocate face-flux MultiFabs for refluxing (x and y, node-centred in each direction)
    const bool need_reflux = (do_reflux && (flux_reg || level < parent->finestLevel()));

    amrex::BoxArray ba_x = amrex::convert(grids, amrex::IntVect::TheDimensionVector(0));
    amrex::BoxArray ba_y = amrex::convert(grids, amrex::IntVect::TheDimensionVector(1));

    MultiFab flux_x0(ba_x, dmap, NUM_STATE, 0);
    MultiFab flux_y0(ba_y, dmap, NUM_STATE, 0);
    MultiFab flux_x1(ba_x, dmap, NUM_STATE, 0);
    MultiFab flux_y1(ba_y, dmap, NUM_STATE, 0);

    MultiFab* pfx0 = need_reflux ? &flux_x0 : nullptr;
    MultiFab* pfy0 = need_reflux ? &flux_y0 : nullptr;
    MultiFab* pfx1 = need_reflux ? &flux_x1 : nullptr;
    MultiFab* pfy1 = need_reflux ? &flux_y1 : nullptr;

    // [3] Fill ghost cells for U^n
    MultiFab U0(grids, dmap, NUM_STATE, NUM_GROW);
    FillPatch(*this, U0, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    U0.FillBoundary(geom.periodicity());

    // [3.0] Clamp ghost cells on fine levels - interpolation across shocks can produce
    //       unphysical values at coarse-fine boundaries
    if (level > 0)
        enforce_positivity(U0, grids, dmap, NUM_GROW, gamma);

    // [3.1] X outflow BCs for U0 - only on patches touching physical boundary
    for (MFIter mfi(U0); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const Ua  = U0.array(mfi);
        const int ilo = vbx.smallEnd(0);
        const int ihi = vbx.bigEnd(0);
        const int jlo = vbx.smallEnd(1);
        const int jhi = vbx.bigEnd(1);

        // [3.1.1] Only fill left ghost if this patch touches the left physical boundary
        const bool at_lo_x = (ilo == geom.Domain().smallEnd(0));
        const bool at_hi_x = (ihi == geom.Domain().bigEnd(0));

        for (int j = jlo; j <= jhi; ++j)
            for (int ig = 1; ig <= NUM_GROW; ++ig)
                for (int n = 0; n < NUM_STATE; ++n) {
                    if (at_lo_x) Ua(ilo - ig, j, 0, n) = Ua(ilo, j, 0, n);
                    if (at_hi_x) Ua(ihi + ig, j, 0, n) = Ua(ihi, j, 0, n);
                }
    }

    // [3.2] Y outflow BCs for U0
    if (!geom.isPeriodic(1))
    {
        for (MFIter mfi(U0); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const Ua  = U0.array(mfi);

            const int ilo = vbx.smallEnd(0);
            const int ihi = vbx.bigEnd(0);
            const int jlo = vbx.smallEnd(1);
            const int jhi = vbx.bigEnd(1);

            const bool at_lo_y = (jlo == geom.Domain().smallEnd(1));
            const bool at_hi_y = (jhi == geom.Domain().bigEnd(1));

            for (int i = ilo; i <= ihi; ++i)
            for (int ig = 1; ig <= NUM_GROW; ++ig)
            for (int n = 0; n < NUM_STATE; ++n)
            {
                if (at_lo_y) Ua(i, jlo-ig, 0, n) = Ua(i, jlo, 0, n);
                if (at_hi_y) Ua(i, jhi+ig, 0, n) = Ua(i, jhi, 0, n);
            }
        }
    }

    enforce_positivity(U0, grids, dmap, NUM_GROW, gamma);

    // [4] Stage 1: U1 = U0 + dt * L(U0)
    MultiFab rhs0(grids, dmap, NUM_STATE, 0);
    euler2d::compute_rhs_2d(U0, rhs0, geom, gamma, pfx0, pfy0);

    MultiFab U1(grids, dmap, NUM_STATE, NUM_GROW);
    U1.setVal(0.0_rt);

    for (MFIter mfi(U1); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto const U0a = U0.const_array(mfi);
        auto const R0a = rhs0.const_array(mfi);
        auto       U1a = U1.array(mfi);
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                for (int n = 0; n < NUM_STATE; ++n)
                    U1a(i,j,0,n) = U0a(i,j,0,n) + dt * R0a(i,j,0,n);
    }

    U1.FillBoundary(geom.periodicity());

    // [4.1] X outflow BCs for U1
    for (MFIter mfi(U1); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const Ua  = U1.array(mfi);
        const int ilo = vbx.smallEnd(0);
        const int ihi = vbx.bigEnd(0);
        const int jlo = vbx.smallEnd(1);
        const int jhi = vbx.bigEnd(1);

        const bool at_lo_x = (ilo == geom.Domain().smallEnd(0));  
        const bool at_hi_x = (ihi == geom.Domain().bigEnd(0));    

        for (int j = jlo; j <= jhi; ++j)
            for (int ig = 1; ig <= NUM_GROW; ++ig)
                for (int n = 0; n < NUM_STATE; ++n) {
                    if (at_lo_x) Ua(ilo - ig, j, 0, n) = Ua(ilo, j, 0, n);
                    if (at_hi_x) Ua(ihi + ig, j, 0, n) = Ua(ihi, j, 0, n); 
                }
    }

    // [4.2] Y outflow BCs for U1
    if (!geom.isPeriodic(1))
    {
        for (MFIter mfi(U1); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            auto const Ua  = U1.array(mfi);

            const int ilo = vbx.smallEnd(0);
            const int ihi = vbx.bigEnd(0);
            const int jlo = vbx.smallEnd(1);
            const int jhi = vbx.bigEnd(1);

            const bool at_lo_y = (jlo == geom.Domain().smallEnd(1));
            const bool at_hi_y = (jhi == geom.Domain().bigEnd(1));

            for (int i = ilo; i <= ihi; ++i)
            for (int ig = 1; ig <= NUM_GROW; ++ig)
            for (int n = 0; n < NUM_STATE; ++n)
            {
                if (at_lo_y) Ua(i, jlo-ig, 0, n) = Ua(i, jlo, 0, n);
                if (at_hi_y) Ua(i, jhi+ig, 0, n) = Ua(i, jhi, 0, n);
            }
        }
    }

    enforce_positivity(U1, grids, dmap, NUM_GROW, gamma);

    // [5] Stage 2: U^{n+1} = 0.5 * (U0 + U1 + dt * L(U1))
    MultiFab rhs1(grids, dmap, NUM_STATE, 0);
    euler2d::compute_rhs_2d(U1, rhs1, geom, gamma, pfx1, pfy1);

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto const U0a = U0.const_array(mfi);
        auto const U1a = U1.const_array(mfi);
        auto const R1a = rhs1.const_array(mfi);
        auto Una = S_new.array(mfi);
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                for (int n = 0; n < NUM_STATE; ++n)
                    Una(i,j,0,n) =
                        0.5_rt * (U0a(i,j,0,n)
                                + U1a(i,j,0,n)
                                + dt * R1a(i,j,0,n));
    }

    enforce_positivity(S_new, grids, dmap, 0, gamma);

    // [6] Accumulate fluxes into FluxRegister
    if (need_reflux)
    {
        const auto dx = geom.CellSizeArray();
        const Real dtw = 0.5_rt * dt;

        const Real area_x = dx[1]; // dy
        const Real area_y = dx[0]; // dx

        const Real scale_x = dtw * area_x;
        const Real scale_y = dtw * area_y;

        if (flux_reg)
        {
            flux_reg->FineAdd(flux_x0, 0, 0, 0, NUM_STATE, scale_x);
            flux_reg->FineAdd(flux_x1, 0, 0, 0, NUM_STATE, scale_x);

            flux_reg->FineAdd(flux_y0, 1, 0, 0, NUM_STATE, scale_y);
            flux_reg->FineAdd(flux_y1, 1, 0, 0, NUM_STATE, scale_y);
        }

        if (level < parent->finestLevel())
        {
            FluxRegister& fr = getFluxReg(level + 1);

            // Always COPY on first call of coarse level
            FluxRegister::FrOp op = FluxRegister::COPY;

            fr.CrseInit(flux_x0, 0, 0, 0, NUM_STATE, -scale_x, op);
            fr.CrseInit(flux_x1, 0, 0, 0, NUM_STATE, -scale_x, FluxRegister::ADD);

            fr.CrseInit(flux_y0, 1, 0, 0, NUM_STATE, -scale_y, op);
            fr.CrseInit(flux_y1, 1, 0, 0, NUM_STATE, -scale_y, FluxRegister::ADD);
        }
    }

    return dt;
}
 

/**
 * Estimate time step.
 */
Real
AmrLevelAdv::estTimeStep (Real)
{
    const Real cur_time = state[Phi_Type].curTime();

    MultiFab U0(grids, dmap, NUM_STATE, NUM_GROW);
    FillPatch(*this, U0, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);

    // [CFL.0] Clamp on fine levels before CFL estimate
    if (level > 0)
        enforce_positivity(U0, grids, dmap, NUM_GROW, gamma);

    Real dt_est = euler2d::compute_dt_cfl_2d(U0, geom, gamma, cfl);

    if (verbose) {
        amrex::Print()
            << "Euler dt level "
            << level
            << " = "
            << dt_est
            << std::endl;
    }

    return dt_est;
}



/**
 * Compute initial time step.
 */
Real
AmrLevelAdv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

/**
 * Compute initial `dt'.
 */
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
                               int                   /*sub_cycle*/,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& /*ref_ratio*/,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

/**
 * Compute new `dt'.
 */
void
AmrLevelAdv::computeNewDt (int                   finest_level,
                           int                   /*sub_cycle*/,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& /*ref_ratio*/,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        AmrLevelAdv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1)
    {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],dt_level[i]);
        }
    }
    else
    {
        //
        // Limit dt's by change_max * old dt
        //
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        }
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

/**
 * Do work after timestep().
 */
void
AmrLevelAdv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

#ifdef AMREX_PARTICLES
    if (TracerPC)
      {
        const int ncycle = parent->nCycle(level);

        if (iteration < ncycle || level == 0)
          {
            int ngrow = (level == 0) ? 0 : iteration;

            TracerPC->Redistribute(level, TracerPC->finestLevel(), ngrow);
          }
      }
#endif
}

/**
 * Do work after regrid().
 */
void
AmrLevelAdv::post_regrid (int lbase, int /*new_finest*/) {
#ifdef AMREX_PARTICLES
  if (TracerPC && level == lbase) {
      TracerPC->Redistribute(lbase);
  }
#else
  amrex::ignore_unused(lbase);
#endif
}

/**
 * Do work after a restart().
 */
void
AmrLevelAdv::post_restart()
{
#ifdef AMREX_PARTICLES
    if (do_tracers && level == 0) {
      BL_ASSERT(TracerPC == 0);
      TracerPC = std::make_unique<AmrTracerParticleContainer>(parent);
      TracerPC->Restart(parent->theRestartFile(), "Tracer");
    }
#endif
}

/**
 * Do work after init().
 */
void
AmrLevelAdv::post_init (Real /*stop_time*/)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

/**
 * Error estimation for regridding.
 */
void
AmrLevelAdv::errorEst (amrex::TagBoxArray& tags,
                       int /*clearval*/,
                       int /*tagval*/,
                       amrex::Real /*time*/,
                       int /*n_error_buf*/,
                       int /*ngrow*/)
{
    using namespace amrex;

    if (level > max_rhograd_lev) return;
    if (rho_grad.empty()) return;

    MultiFab& S_new = get_new_data(Phi_Type);

    MultiFab Stmp;
    const Real cur_time = state[Phi_Type].curTime();
    Stmp.define(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 1);
    FillPatch(*this, Stmp, 1, cur_time, Phi_Type, 0, S_new.nComp());

    MultiFab const& S = Stmp;

    const Real thr = rho_grad[level];
    const int  comp_rho = 0;

    const Real dx = geom.CellSize(0);

    const char setval = TagBox::SET;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();      // valid cells
        auto const s  = S.const_array(mfi);
        auto       t  = tags.array(mfi);

        const int ilo = bx.smallEnd(0);
        const int ihi = bx.bigEnd(0);

        const int ptype = prob_type;

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int im = (i == ilo) ? i : i-1;
        int ip = (i == ihi) ? i : i+1;

        int jlo = bx.smallEnd(1);
        int jhi = bx.bigEnd(1);

        int jm = (j == jlo) ? j : j-1;
        int jp = (j == jhi) ? j : j+1;

        Real drdx = amrex::Math::abs(s(ip,j,k,comp_rho) - s(im,j,k,comp_rho)) / ((ip-im)*dx);
        Real drdy = amrex::Math::abs(s(i,jp,k,comp_rho) - s(i,jm,k,comp_rho)) / ((jp-jm)*geom.CellSize(1));

        Real grad;

        if (ptype == 0) {
            // x-split shock tube
            grad = drdx;
        } else if (ptype == 1) {
            // y-split shock tube
            grad = drdy;
        } else {
            // 2D: use full gradient magnitude
            grad = std::sqrt(drdx*drdx + drdy*drdy);
        }

        if (grad > thr) {
            t(i,j,k) = setval;
        }
    });
    }

    if (verbose && level == 0 && amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::Print() << "[errorEst] level=" << level
                    << " thr=" << thr
                    << " time=" << state[Phi_Type].curTime()
                    << "\n";
    }
}


/**
 * Read parameters from input file.
 */
void
AmrLevelAdv::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");

    pp.query("max_rhograd_lev", max_rhograd_lev);

    if (max_rhograd_lev >= 0) {
        rho_grad.resize(max_rhograd_lev+1, 0.0_rt);
        pp.queryarr("rho_grad", rho_grad, 0, max_rhograd_lev+1);
    }

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);
    pp.query("gamma", gamma);

    pp.query("max_rhograd_lev", max_rhograd_lev);

    ParmParse pprob("prob");
    pprob.query("type", prob_type);

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // This tutorial code only supports Cartesian coordinates.
    if (! gg->IsCartesian()) {
        amrex::Abort("Please set geom.coord_sys = 0");
    }


#ifdef AMREX_PARTICLES
    pp.query("do_tracers", do_tracers);
#endif

    // Read tagging parameters from tagging block in the input file.
    // See Src_nd/Tagging_params.cpp for the function implementation.
    get_tagging_params();
}

void
AmrLevelAdv::reflux ()
{
    BL_ASSERT(level < parent->finestLevel());

    const auto strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type), 1.0, 0, 0, NUM_STATE, geom);

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        auto      end    = amrex::second() - strt;
        ParallelDescriptor::ReduceRealMax(end, IOProc);
        amrex::Print() << "AmrLevelAdv::reflux() at level " << level
                       << " : time = " << end << std::endl;
    }
}


void
AmrLevelAdv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelAdv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);

    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

#ifdef AMREX_PARTICLES
void
AmrLevelAdv::init_particles ()
{
  if (do_tracers && level == 0)
    {
      BL_ASSERT(TracerPC == nullptr);

      TracerPC = std::make_unique<AmrTracerParticleContainer>(parent);

      AmrTracerParticleContainer::ParticleInitData pdata = {{AMREX_D_DECL(0.0, 0.0, 0.0)},{},{},{}};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
    }
}
#endif
