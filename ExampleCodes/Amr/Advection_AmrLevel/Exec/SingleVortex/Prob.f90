! [0] Problem Parameter Module
!     - Stores user-defined primitive states
!     - Handles runtime parameter parsing via ParmParse
!     - Provides primitive to conserved conversion

module prob_params

  ! [0.1] AMReX Dependencies
  use amrex_parmparse_module
  use amrex_fort_module, only: amrex_real

  implicit none

  ! [1] Global Problem Controls
  ! Ensures parameter parsing is performed only once
  logical :: inited = .false.

  ! Problem selector:
  ! 0 -> x-aligned shock tube
  ! 1 -> y-aligned shock tube
  ! 2 -> oblique discontinuity
  ! 3 -> quadrant (Lax–Liu)
  integer :: prob_type = 0

  ! Ratio of specific heats 
  real(amrex_real) :: gamma = 1.4_amrex_real

  ! [2] Discontinuity Geometry Parameters

  ! Reference point of discontinuity
  real(amrex_real) :: x0 = 0.5_amrex_real
  real(amrex_real) :: y0 = 0.5_amrex_real

  ! Normal vector for oblique interface 
  real(amrex_real) :: nx = 1.0_amrex_real
  real(amrex_real) :: ny = 0.0_amrex_real


  ! [3] Two-State Primitive Inputs (Left / Right)
  ! Used for shock tube configurations
  real(amrex_real) :: rhoL=1.0_amrex_real,   uL=0.0_amrex_real, vL=0.0_amrex_real, pL=1.0_amrex_real
  real(amrex_real) :: rhoR=0.125_amrex_real, uR=0.0_amrex_real, vR=0.0_amrex_real, pR=0.1_amrex_real

  
  
  ! [4] Quadrant Primitive Inputs (4-State)
  ! Used for 2D Riemann problems 
  real(amrex_real) :: rhoQ(4), uQ(4), vQ(4), pQ(4)

contains

  ! [5] prob_init()
  !     - Reads problem parameters from inputs file
  !     - Sets defaults if not provided
  !     - Called once at runtime

  subroutine prob_init()
    type(amrex_parmparse) :: pp

    if (inited) return
    inited = .true.

    ! [5.1] Default Lax–Liu Quadrant States
    ! These correspond to the standard 2D test case

    rhoQ = (/ 1.5_amrex_real, 0.5323_amrex_real, 0.138_amrex_real, 0.5323_amrex_real /)
    uQ = (/ 0.0_amrex_real, 1.206_amrex_real, 1.206_amrex_real, 0.0_amrex_real /)
    vQ = (/ 0.0_amrex_real, 0.0_amrex_real, 1.206_amrex_real, 1.206_amrex_real /)
    pQ = (/ 1.5_amrex_real, 0.3_amrex_real, 0.029_amrex_real, 0.3_amrex_real /)

    call amrex_parmparse_build(pp, "prob")

    ! [5.2] Read User Overrides 
    call pp%query("type",  prob_type)
    call pp%query("gamma", gamma)

    call pp%query("x0", x0)
    call pp%query("y0", y0)
    call pp%query("nx", nx)
    call pp%query("ny", ny)

    call pp%query("rhoL", rhoL)
    call pp%query("uL",   uL)
    call pp%query("vL",   vL)
    call pp%query("pL",   pL)

    call pp%query("rhoR", rhoR)
    call pp%query("uR",   uR)
    call pp%query("vR",   vR)
    call pp%query("pR",   pR)

    ! Quadrant overrides
    call pp%query("rhoQ(1)", rhoQ(1))
    call pp%query("rhoQ(2)", rhoQ(2))
    call pp%query("rhoQ(3)", rhoQ(3))
    call pp%query("rhoQ(4)", rhoQ(4))

    call pp%query("uQ(1)", uQ(1))
    call pp%query("uQ(2)", uQ(2))
    call pp%query("uQ(3)", uQ(3))
    call pp%query("uQ(4)", uQ(4))

    call pp%query("vQ(1)", vQ(1))
    call pp%query("vQ(2)", vQ(2))
    call pp%query("vQ(3)", vQ(3))
    call pp%query("vQ(4)", vQ(4))

    call pp%query("pQ(1)", pQ(1))
    call pp%query("pQ(2)", pQ(2))
    call pp%query("pQ(3)", pQ(3))
    call pp%query("pQ(4)", pQ(4))

    call amrex_parmparse_destroy(pp)

  end subroutine prob_init


  ! [6] prim_to_cons()
  !     Converts primitive variables:(rho, u, v, p)
  !     conserved form: (rho, rho u, rho v, E)
  !     Total energy: E = p/(gamma-1) + 1/2 rho (u^2 + v^2)

  subroutine prim_to_cons(rho,u,v,p, Uout)
    real(amrex_real), intent(in)  :: rho,u,v,p
    real(amrex_real), intent(out) :: Uout(4)

    real(amrex_real) :: eint, Etot

    eint = p/(gamma - 1.0_amrex_real)
    Etot = eint + 0.5_amrex_real*rho*(u*u + v*v)

    Uout(1) = rho
    Uout(2) = rho*u
    Uout(3) = rho*v
    Uout(4) = Etot

  end subroutine prim_to_cons

end module prob_params


! [7] initdata()
!     - Sets initial condition on each AMReX tile
!     - Evaluates primitive states based on prob_type
!     - Writes conserved variables into state array

subroutine initdata(level, time, lo, hi, state, dlo, dhi, dx, prob_lo)

  use prob_params
  use amrex_fort_module, only: amrex_real

  implicit none

  ! [7.1] AMReX Loop Metadata

  integer, intent(in) :: level
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: dlo(3), dhi(3)

  real(amrex_real), intent(in) :: time
  real(amrex_real), intent(in) :: dx(3), prob_lo(3)

  ! Conserved variables layout:
  ! state(i,j,k,:) = (rho, rho u, rho v, E)
  real(amrex_real), intent(inout) :: state(dlo(1):dhi(1), &
                                           dlo(2):dhi(2), &
                                           dlo(3):dhi(3), 4)

  integer :: i, j, k, q
  real(amrex_real) :: x, y, sgn
  real(amrex_real) :: Uc(4)

  call prob_init()

  ! 2D solver: k fixed
  k = lo(3)


  ! [7.2] Loop Over Physical Cells
  !       Note: i,j are global indices in AMReX

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        ! Compute cell-centre coordinates
        x = prob_lo(1) + (real(i, amrex_real) + 0.5_amrex_real) * dx(1)
        y = prob_lo(2) + (real(j, amrex_real) + 0.5_amrex_real) * dx(2)

        select case (prob_type)


        ! [7.3] Case 0: x-aligned shock tube
        case (0)
           if (x < x0) then
              call prim_to_cons(rhoL,uL,vL,pL, Uc)
           else
              call prim_to_cons(rhoR,uR,vR,pR, Uc)
           end if

        ! [7.4] Case 1: y-aligned shock tube
        case (1)
           if (y < y0) then
              call prim_to_cons(rhoL,uL,vL,pL, Uc)
           else
              call prim_to_cons(rhoR,uR,vR,pR, Uc)
           end if

        ! [7.5] Case 2: Oblique discontinuity
        ! Interface defined by normal (nx, ny)
        case (2)
           sgn = (x-x0)*nx + (y-y0)*ny
           if (sgn < 0.0_amrex_real) then
              call prim_to_cons(rhoL,uL,vL,pL, Uc)
           else
              call prim_to_cons(rhoR,uR,vR,pR, Uc)
           end if

        ! [7.6] Case 3: 2D Quadrant Riemann 
        case (3)
           if (x > x0 .and. y > y0) then
              q = 1
           else if (x < x0 .and. y > y0) then
              q = 2
           else if (x < x0 .and. y < y0) then
              q = 3
           else
              q = 4
           end if

           call prim_to_cons(rhoQ(q),uQ(q),vQ(q),pQ(q), Uc)

        
        ! [7.7] Default: fallback to left state
        case default
           call prim_to_cons(rhoL,uL,vL,pL, Uc)

        end select

        ! Write conserved state
        state(i,j,k,1) = Uc(1)
        state(i,j,k,2) = Uc(2)
        state(i,j,k,3) = Uc(3)
        state(i,j,k,4) = Uc(4)

     end do
  end do

end subroutine initdata