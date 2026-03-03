#ifndef _Adv_F_H_
#define _Adv_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void initdata_(
        int* level,
        amrex::Real* time,
        const int* lo,
        const int* hi,
        amrex::Real* state,
        const int* dlo,
        const int* dhi,
        const amrex::Real* dx,
        const amrex::Real* prob_lo
    );

#ifdef __cplusplus
}
#endif

#endif