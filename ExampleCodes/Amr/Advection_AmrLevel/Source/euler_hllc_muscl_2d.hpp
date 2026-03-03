#pragma once

// [0] Standard Libraries
#include <algorithm>
#include <cmath>
#include <vector>

// [1] AMReX
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

namespace euler2d
{
    // [0] Component Indices
    enum : int {
        QRHO  = 0,
        QMOMX = 1,
        QMOMY = 2,
        QENER = 3,
        NVAR = 4
    };

    // [1] Small Utilities
    inline amrex::Real clamp_pos(amrex::Real x, amrex::Real eps = 1e-12_rt)
    {
        return std::max(x, eps);
    }

    // [1.1] MC limiter (Monotonized Central)
    inline amrex::Real mc_limiter(amrex::Real a, amrex::Real b)
    {
        if (a * b <= 0.0_rt) return 0.0_rt;

        const amrex::Real s = (a > 0.0_rt) ? 1.0_rt : -1.0_rt;
        const amrex::Real abs_a = std::abs(a);
        const amrex::Real abs_b = std::abs(b);
        const amrex::Real abs_c = 0.5_rt * std::abs(a + b);

        return s * std::min({2.0_rt * abs_a, 2.0_rt * abs_b, abs_c});
    }

    // [2] Primitive Variables
    struct Primitive
    {
        amrex::Real rho;
        amrex::Real ux;
        amrex::Real uy;
        amrex::Real p;
    };

    // [3] Conversions and EoS
    inline Primitive cons_to_prim(amrex::Real rho, amrex::Real momx, amrex::Real momy,
                                 amrex::Real E, amrex::Real gamma)
    {
        Primitive W{};
        W.rho = clamp_pos(rho);
        W.ux  = momx / W.rho;
        W.uy  = momy / W.rho;

        const amrex::Real kinetic = 0.5_rt * (momx*momx + momy*momy) / W.rho;
        const amrex::Real p_raw   = (gamma - 1.0_rt) * (E - kinetic);
        W.p = clamp_pos(p_raw);

        return W;
    }

    inline void prim_to_cons(const Primitive& W, amrex::Real gamma,
                             amrex::Real& rho, amrex::Real& momx, amrex::Real& momy, amrex::Real& E)
    {
        const amrex::Real rr = clamp_pos(W.rho);
        const amrex::Real pp = clamp_pos(W.p);

        rho  = rr;
        momx = rr * W.ux;
        momy = rr * W.uy;

        const amrex::Real kinetic = 0.5_rt * rr * (W.ux*W.ux + W.uy*W.uy);
        E = pp / (gamma - 1.0_rt) + kinetic;
    }

    inline amrex::Real sound_speed(const Primitive& W, amrex::Real gamma)
    {
        return std::sqrt(gamma * clamp_pos(W.p) / clamp_pos(W.rho));
    }

    // [4] Physical Fluxes
    inline void flux_x(const Primitive& W, amrex::Real gamma,
                       amrex::Real& f_rho, amrex::Real& f_momx, amrex::Real& f_momy, amrex::Real& f_E)
    {
        amrex::Real rho, momx, momy, E;
        prim_to_cons(W, gamma, rho, momx, momy, E);

        f_rho = momx;
        f_momx = momx * W.ux + W.p;
        f_momy = momx * W.uy;
        f_E = (E + W.p) * W.ux;
    }

    inline void flux_y(const Primitive& W, amrex::Real gamma,
                       amrex::Real& g_rho, amrex::Real& g_momx, amrex::Real& g_momy, amrex::Real& g_E)
    {
        amrex::Real rho, momx, momy, E;
        prim_to_cons(W, gamma, rho, momx, momy, E);

        g_rho = momy;
        g_momx = momy * W.ux;
        g_momy = momy * W.uy + W.p;
        g_E = (E + W.p) * W.uy;
    }

    // [5] HLLC Riemann Solver (primitive input)
    inline void hllc_flux_x(const Primitive& WL, const Primitive& WR, amrex::Real gamma,
                            amrex::Real& f_rho, amrex::Real& f_momx, amrex::Real& f_momy, amrex::Real& f_E,
                            amrex::Real tol = 1e-12_rt)
    {
        const amrex::Real cL = sound_speed(WL, gamma);
        const amrex::Real cR = sound_speed(WR, gamma);

        const amrex::Real sL = std::min(WL.ux - cL, WR.ux - cR);
        const amrex::Real sR = std::max(WL.ux + cL, WR.ux + cR);

        const amrex::Real num =
            (WR.p - WL.p)
            + WL.rho * WL.ux * (sL - WL.ux)
            - WR.rho * WR.ux * (sR - WR.ux);

        const amrex::Real den_raw =
            WL.rho * (sL - WL.ux)
            - WR.rho * (sR - WR.ux);

        const amrex::Real den =
            (std::abs(den_raw) < tol) ? ((den_raw < 0.0_rt) ? -tol : tol) : den_raw;

        const amrex::Real s_star = num / den;

        amrex::Real rL, mxL, myL, EL;
        amrex::Real rR, mxR, myR, ER;
        prim_to_cons(WL, gamma, rL, mxL, myL, EL);
        prim_to_cons(WR, gamma, rR, mxR, myR, ER);

        amrex::Real FL_r, FL_mx, FL_my, FL_E;
        amrex::Real FR_r, FR_mx, FR_my, FR_E;
        flux_x(WL, gamma, FL_r, FL_mx, FL_my, FL_E);
        flux_x(WR, gamma, FR_r, FR_mx, FR_my, FR_E);

        if (0.0_rt <= sL) {
            f_rho = FL_r; f_momx = FL_mx; f_momy = FL_my; f_E = FL_E;
            return;
        }

        if (sL <= 0.0_rt && 0.0_rt <= s_star)
        {
            const amrex::Real rL_star = WL.rho * (sL - WL.ux) / (sL - s_star);

            const amrex::Real EL_star =
                rL_star * (
                    (EL / WL.rho)
                    + (s_star - WL.ux) * (s_star + WL.p / (WL.rho * (sL - WL.ux)))
                );

            const amrex::Real mxL_star = rL_star * s_star;
            const amrex::Real myL_star = rL_star * WL.uy;

            f_rho = FL_r  + sL * (rL_star  - rL);
            f_momx = FL_mx + sL * (mxL_star - mxL);
            f_momy = FL_my + sL * (myL_star - myL);
            f_E = FL_E  + sL * (EL_star  - EL);
            return;
        }

        if (s_star <= 0.0_rt && 0.0_rt <= sR)
        {
            const amrex::Real rR_star = WR.rho * (sR - WR.ux) / (sR - s_star);

            const amrex::Real ER_star =
                rR_star * (
                    (ER / WR.rho)
                    + (s_star - WR.ux) * (s_star + WR.p / (WR.rho * (sR - WR.ux)))
                );

            const amrex::Real mxR_star = rR_star * s_star;
            const amrex::Real myR_star = rR_star * WR.uy;

            f_rho = FR_r  + sR * (rR_star  - rR);
            f_momx = FR_mx + sR * (mxR_star - mxR);
            f_momy = FR_my + sR * (myR_star - myR);
            f_E = FR_E  + sR * (ER_star  - ER);
            return;
        }


        f_rho = FR_r; f_momx = FR_mx; f_momy = FR_my; f_E = FR_E;
    }


    inline void hllc_flux_y(const Primitive& WL, const Primitive& WR, amrex::Real gamma,
                            amrex::Real& g_rho, amrex::Real& g_momx, amrex::Real& g_momy, amrex::Real& g_E,
                            amrex::Real tol = 1e-12_rt)
    {
        const amrex::Real cL = sound_speed(WL, gamma);
        const amrex::Real cR = sound_speed(WR, gamma);

        const amrex::Real sL = std::min(WL.uy - cL, WR.uy - cR);
        const amrex::Real sR = std::max(WL.uy + cL, WR.uy + cR);

        const amrex::Real num =
            (WR.p - WL.p)
            + WL.rho * WL.uy * (sL - WL.uy)
            - WR.rho * WR.uy * (sR - WR.uy);

        const amrex::Real den_raw =
            WL.rho * (sL - WL.uy)
            - WR.rho * (sR - WR.uy);

        const amrex::Real den =
            (std::abs(den_raw) < tol) ? ((den_raw < 0.0_rt) ? -tol : tol) : den_raw;

        const amrex::Real s_star = num / den;

        amrex::Real rL, mxL, myL, EL;
        amrex::Real rR, mxR, myR, ER;
        prim_to_cons(WL, gamma, rL, mxL, myL, EL);
        prim_to_cons(WR, gamma, rR, mxR, myR, ER);

        amrex::Real GL_r, GL_mx, GL_my, GL_E;
        amrex::Real GR_r, GR_mx, GR_my, GR_E;
        flux_y(WL, gamma, GL_r, GL_mx, GL_my, GL_E);
        flux_y(WR, gamma, GR_r, GR_mx, GR_my, GR_E);

        if (0.0_rt <= sL) {
            g_rho = GL_r; g_momx = GL_mx; g_momy = GL_my; g_E = GL_E;
            return;
        }

        if (sL <= 0.0_rt && 0.0_rt <= s_star)
        {
            const amrex::Real rL_star = WL.rho * (sL - WL.uy) / (sL - s_star);

            const amrex::Real EL_star =
                rL_star * (
                    (EL / WL.rho)
                    + (s_star - WL.uy) * (s_star + WL.p / (WL.rho * (sL - WL.uy)))
                );

            const amrex::Real mxL_star = rL_star * WL.ux;
            const amrex::Real myL_star = rL_star * s_star;

            g_rho = GL_r  + sL * (rL_star  - rL);
            g_momx = GL_mx + sL * (mxL_star - mxL);
            g_momy = GL_my + sL * (myL_star - myL);
            g_E = GL_E  + sL * (EL_star  - EL);
            return;
        }

        
        if (s_star <= 0.0_rt && 0.0_rt <= sR)
        {
            const amrex::Real rR_star = WR.rho * (sR - WR.uy) / (sR - s_star);

            const amrex::Real ER_star =
                rR_star * (
                    (ER / WR.rho)
                    + (s_star - WR.uy) * (s_star + WR.p / (WR.rho * (sR - WR.uy)))
                );

            const amrex::Real mxR_star = rR_star * WR.ux;
            const amrex::Real myR_star = rR_star * s_star;

            g_rho = GR_r  + sR * (rR_star  - rR);
            g_momx = GR_mx + sR * (mxR_star - mxR);
            g_momy = GR_my + sR * (myR_star - myR); 
            g_E = GR_E  + sR * (ER_star  - ER);
            return;
        }


        g_rho = GR_r; g_momx = GR_mx; g_momy = GR_my; g_E = GR_E;
    }

    // [6] Spatial Operator: computes dU/dt = L(U) on valid cells
    inline void compute_rhs_2d(const amrex::MultiFab& U_ng, amrex::MultiFab& rhs,
                               const amrex::Geometry& geom, amrex::Real gamma)
    {
        const auto dx = geom.CellSizeArray();
        const amrex::Real dx0 = dx[0];
        const amrex::Real dx1 = dx[1];

        rhs.setVal(0.0_rt);

        for (amrex::MFIter mfi(rhs); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();

            auto const U = U_ng.const_array(mfi);
            auto R = rhs.array(mfi);

            const int ilo = bx.smallEnd(0);
            const int ihi = bx.bigEnd(0);
            const int jlo = bx.smallEnd(1);
            const int jhi = bx.bigEnd(1);

            const int nxF = (ihi - ilo + 2);
            const int nyF = (jhi - jlo + 1);
            const int nxG = (ihi - ilo + 1);
            const int nyG = (jhi - jlo + 2);

            std::vector<amrex::Real> Fx(nxF * nyF * NVAR, 0.0_rt);
            std::vector<amrex::Real> Gy(nxG * nyG * NVAR, 0.0_rt);

            auto Fx_at = [&](int fi, int fj, int n) -> amrex::Real& {
                return Fx[((fj * nxF) + fi) * NVAR + n];
            };

            auto Gy_at = [&](int fi, int fj, int n) -> amrex::Real& {
                return Gy[((fj * nxG) + fi) * NVAR + n];
            };

            // [6.2] X-direction face fluxes (i+1/2)
            for (int j = jlo; j <= jhi; ++j)
            {
                for (int i = ilo; i <= ihi + 1; ++i)
                {
                    const int iL = i - 1;
                    const int iR = i;

                    auto WLm = cons_to_prim(U(iL-1,j,0,QRHO), U(iL-1,j,0,QMOMX), U(iL-1,j,0,QMOMY), U(iL-1,j,0,QENER), gamma);
                    auto WL = cons_to_prim(U(iL  ,j,0,QRHO), U(iL  ,j,0,QMOMX), U(iL  ,j,0,QMOMY), U(iL  ,j,0,QENER), gamma);
                    auto WLp = cons_to_prim(U(iL+1,j,0,QRHO), U(iL+1,j,0,QMOMX), U(iL+1,j,0,QMOMY), U(iL+1,j,0,QENER), gamma);

                    auto WRm = cons_to_prim(U(iR-1,j,0,QRHO), U(iR-1,j,0,QMOMX), U(iR-1,j,0,QMOMY), U(iR-1,j,0,QENER), gamma);
                    auto WR = cons_to_prim(U(iR  ,j,0,QRHO), U(iR  ,j,0,QMOMX), U(iR  ,j,0,QMOMY), U(iR  ,j,0,QENER), gamma);
                    auto WRp = cons_to_prim(U(iR+1,j,0,QRHO), U(iR+1,j,0,QMOMX), U(iR+1,j,0,QMOMY), U(iR+1,j,0,QENER), gamma);

                    Primitive dWL{};
                    dWL.rho = mc_limiter(WL.rho - WLm.rho, WLp.rho - WL.rho);
                    dWL.ux = mc_limiter(WL.ux  - WLm.ux , WLp.ux  - WL.ux );
                    dWL.uy = mc_limiter(WL.uy  - WLm.uy , WLp.uy  - WL.uy );
                    dWL.p = mc_limiter(WL.p   - WLm.p  , WLp.p   - WL.p  );

                    Primitive dWR{};
                    dWR.rho = mc_limiter(WR.rho - WRm.rho, WRp.rho - WR.rho);
                    dWR.ux = mc_limiter(WR.ux  - WRm.ux , WRp.ux  - WR.ux );
                    dWR.uy = mc_limiter(WR.uy  - WRm.uy , WRp.uy  - WR.uy );
                    dWR.p = mc_limiter(WR.p   - WRm.p  , WRp.p   - WR.p  );

                    Primitive WLf{ WL.rho + 0.5_rt*dWL.rho, WL.ux + 0.5_rt*dWL.ux, WL.uy + 0.5_rt*dWL.uy, WL.p + 0.5_rt*dWL.p };
                    Primitive WRf{ WR.rho - 0.5_rt*dWR.rho, WR.ux - 0.5_rt*dWR.ux, WR.uy - 0.5_rt*dWR.uy, WR.p - 0.5_rt*dWR.p };

                    amrex::Real fr, fmx, fmy, fE;
                    hllc_flux_x(WLf, WRf, gamma, fr, fmx, fmy, fE);

                    const int fi = i - ilo;
                    const int fj = j - jlo;

                    Fx_at(fi, fj, QRHO)  = fr;
                    Fx_at(fi, fj, QMOMX) = fmx;
                    Fx_at(fi, fj, QMOMY) = fmy;
                    Fx_at(fi, fj, QENER) = fE;
                }
            }

            // [6.3] Y-direction face fluxes (j+1/2)
            for (int j = jlo; j <= jhi + 1; ++j)
            {
                for (int i = ilo; i <= ihi; ++i)
                {
                    const int jL = j - 1;
                    const int jR = j;

                    auto WLm = cons_to_prim(U(i,jL-1,0,QRHO), U(i,jL-1,0,QMOMX), U(i,jL-1,0,QMOMY), U(i,jL-1,0,QENER), gamma);
                    auto WL  = cons_to_prim(U(i,jL  ,0,QRHO), U(i,jL  ,0,QMOMX), U(i,jL  ,0,QMOMY), U(i,jL  ,0,QENER), gamma);
                    auto WLp = cons_to_prim(U(i,jL+1,0,QRHO), U(i,jL+1,0,QMOMX), U(i,jL+1,0,QMOMY), U(i,jL+1,0,QENER), gamma);

                    auto WRm = cons_to_prim(U(i,jR-1,0,QRHO), U(i,jR-1,0,QMOMX), U(i,jR-1,0,QMOMY), U(i,jR-1,0,QENER), gamma);
                    auto WR  = cons_to_prim(U(i,jR  ,0,QRHO), U(i,jR  ,0,QMOMX), U(i,jR  ,0,QMOMY), U(i,jR  ,0,QENER), gamma);
                    auto WRp = cons_to_prim(U(i,jR+1,0,QRHO), U(i,jR+1,0,QMOMX), U(i,jR+1,0,QMOMY), U(i,jR+1,0,QENER), gamma);

                    Primitive dWL{};
                    dWL.rho = mc_limiter(WL.rho - WLm.rho, WLp.rho - WL.rho);
                    dWL.ux = mc_limiter(WL.ux - WLm.ux, WLp.ux - WL.ux);
                    dWL.uy = mc_limiter(WL.uy - WLm.uy, WLp.uy - WL.uy);
                    dWL.p = mc_limiter(WL.p - WLm.p, WLp.p - WL.p);

                    Primitive dWR{};
                    dWR.rho = mc_limiter(WR.rho - WRm.rho, WRp.rho - WR.rho);
                    dWR.ux = mc_limiter(WR.ux - WRm.ux, WRp.ux - WR.ux);
                    dWR.uy = mc_limiter(WR.uy - WRm.uy , WRp.uy - WR.uy);
                    dWR.p = mc_limiter(WR.p - WRm.p, WRp.p - WR.p);

                    Primitive WLf{ WL.rho + 0.5_rt*dWL.rho, WL.ux + 0.5_rt*dWL.ux, WL.uy + 0.5_rt*dWL.uy, WL.p + 0.5_rt*dWL.p };
                    Primitive WRf{ WR.rho - 0.5_rt*dWR.rho, WR.ux - 0.5_rt*dWR.ux, WR.uy - 0.5_rt*dWR.uy, WR.p - 0.5_rt*dWR.p };

                    amrex::Real gr, gmx, gmy, gE;
                    hllc_flux_y(WLf, WRf, gamma, gr, gmx, gmy, gE);

                    const int fi = i - ilo;
                    const int fj = j - jlo;

                    Gy_at(fi, fj, QRHO)  = gr;
                    Gy_at(fi, fj, QMOMX) = gmx;
                    Gy_at(fi, fj, QMOMY) = gmy;
                    Gy_at(fi, fj, QENER) = gE;
                }
            }

            // [6.4] Divergence of fluxes
            for (int j = jlo; j <= jhi; ++j)
            {
                for (int i = ilo; i <= ihi; ++i)
                {
                    const int fi = i - ilo;
                    const int fj = j - jlo;

                    for (int n = 0; n < NVAR; ++n)
                    {
                        const amrex::Real dF = (Fx_at(fi+1, fj, n) - Fx_at(fi, fj, n)) / dx0;
                        const amrex::Real dG = (Gy_at(fi, fj+1, n) - Gy_at(fi, fj, n)) / dx1;
                        R(i,j,0,n) = -(dF + dG);
                    }
                }
            }
        }
    }

    // [7] Stable Euler CFL estimate for dt
    inline amrex::Real compute_dt_cfl_2d(const amrex::MultiFab& U_ng,
                                         const amrex::Geometry& geom,
                                         amrex::Real gamma,
                                         amrex::Real cfl)
    {
        const auto dx = geom.CellSizeArray();
        const amrex::Real dx0 = dx[0];
        const amrex::Real dx1 = dx[1];

        amrex::Real dt_min = 1.0e200_rt;
        amrex::Real max_speed = 0.0_rt;

        // Record one bad cell (local) for diagnostics
        int bad_i = 0, bad_j = 0;
        amrex::Real bad_rho = 0.0_rt, bad_p = 0.0_rt, bad_E = 0.0_rt;

        for (amrex::MFIter mfi(U_ng); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            auto const U = U_ng.const_array(mfi);

            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j)
            {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i)
                {
                    const amrex::Real rho  = U(i,j,0,QRHO);
                    const amrex::Real momx = U(i,j,0,QMOMX);
                    const amrex::Real momy = U(i,j,0,QMOMY);
                    const amrex::Real E    = U(i,j,0,QENER);

                    if (rho <= 0.0_rt) {
                        bad_i = i; bad_j = j;
                        bad_rho = rho; bad_E = E; bad_p = 0.0_rt;
                        dt_min = 0.0_rt;
                        continue;
                    }

                    const amrex::Real ux = momx / rho;
                    const amrex::Real uy = momy / rho;

                    const amrex::Real kinetic = 0.5_rt * (momx*momx + momy*momy) / rho;
                    const amrex::Real p_raw= (gamma - 1.0_rt) * (E - kinetic);

                    // [7.1] Positivity floor for CFL 
                    constexpr amrex::Real p_floor = 1e-12_rt;

                    if (p_raw <= 0.0_rt) {
                        bad_i = i; bad_j = j;
                        bad_rho = rho; bad_E = E; bad_p = p_raw;
                    }

                    // Use floored pressure for sound speed
                    const amrex::Real p_use = (p_raw > p_floor) ? p_raw : p_floor;
                    const amrex::Real c= std::sqrt(gamma * p_use / rho);
                    const amrex::Real sx = (std::abs(ux) + c) / dx0;
                    const amrex::Real sy = (std::abs(uy) + c) / dx1;

                    const amrex::Real denom = sx + sy;
                    dt_min = std::min(dt_min, 1.0_rt / denom);

                    const amrex::Real speed = std::max(std::abs(ux) + c, std::abs(uy) + c);
                    max_speed = std::max(max_speed, speed);
                }
            }
        }

        amrex::ParallelDescriptor::ReduceRealMin(dt_min);
        amrex::ParallelDescriptor::ReduceRealMax(max_speed);

        if (dt_min <= 0.0_rt) {
            amrex::Print()
                << "Non-physical state encountered in CFL estimate.\n"
                << "Example bad cell (local indices) i=" << bad_i << " j=" << bad_j << "\n"
                << "rho=" << bad_rho << " p_raw=" << bad_p << " E=" << bad_E << "\n";
            amrex::Abort("Aborting due to non-physical rho or pressure in CFL computation.");
        }

        //amrex::Print() << "Max wavespeed (max(|u|+c, |v|+c)) = " << max_speed << "\n";

        return cfl * dt_min;
    }
}