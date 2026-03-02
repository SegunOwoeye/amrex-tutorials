subroutine initdata(level, time, lo, hi, state, dlo, dhi, dx, prob_lo)

  implicit none

  integer level
  integer lo(3), hi(3)
  integer dlo(3), dhi(3)

  double precision time
  double precision dx(3), prob_lo(3)

  double precision state(dlo(1):dhi(1), &
                         dlo(2):dhi(2), &
                         dlo(3):dhi(3), 4)

  integer i, j, k
  double precision x

  k = lo(3)

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        x = prob_lo(1) + (dble(i - lo(1)) + 0.5d0) * dx(1)

        if (x .lt. 0.5d0) then
           state(i,j,k,1) = 1.0d0
           state(i,j,k,2) = 0.0d0
           state(i,j,k,3) = 0.0d0
           state(i,j,k,4) = 2.5d0
        else
           state(i,j,k,1) = 0.125d0
           state(i,j,k,2) = 0.0d0
           state(i,j,k,3) = 0.0d0
           state(i,j,k,4) = 0.25d0
        endif

     end do
  end do

end subroutine initdata