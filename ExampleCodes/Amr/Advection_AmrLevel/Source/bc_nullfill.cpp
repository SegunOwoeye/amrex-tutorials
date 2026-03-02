// bc_nullfill.cpp (or whatever file defines nullfill)

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

// Transmissive in x: copy nearest interior cell into ghost cells.
// Periodic in y: AMReX handles that via FillBoundary(geom.periodicity()).
struct TransmissiveXFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
    {
        const Box dom = geom.Domain();

        const int ilo = dom.smallEnd(0);
        const int ihi = dom.bigEnd(0);

        const int i = iv[0];
        const int j = iv[1];
#if (AMREX_SPACEDIM == 3)
        const int k = iv[2];
#else
        const int k = 0;
#endif

        // Left x boundary ghost
        if (i < ilo) {
            for (int n = 0; n < numcomp; ++n) {
                dest(i,j,k,dcomp+n) = dest(ilo,j,k,dcomp+n);
            }
            return;
        }

        // Right x boundary ghost
        if (i > ihi) {
            for (int n = 0; n < numcomp; ++n) {
                dest(i,j,k,dcomp+n) = dest(ihi,j,k,dcomp+n);
            }
            return;
        }

        // Interior: do nothing
    }
};

void nullfill (Box const& bx, FArrayBox& data,
               const int dcomp, const int numcomp,
               Geometry const& geom, const Real time,
               const Vector<BCRec>& bcr, const int bcomp,
               const int scomp)
{
    GpuBndryFuncFab<TransmissiveXFill> gpu_bndry_func(TransmissiveXFill{});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}