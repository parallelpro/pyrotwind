!===============================================================================
! subroutine ratext
!
! Rational-function (diagonal) extrapolation to the limit h -> 0, adapted
! from Numerical Recipes Sec. 16.4, routine rzextr ("Richardson
! Extrapolation and the Bulirsch-Stoer Method"). Used by the
! Bulirsch-Stoer step integrators int1zone/int2zone: each call supplies
! one more (step size, estimate) pair from a sequence of increasingly
! fine sub-steps, and ratext extrapolates that sequence to an improved
! estimate at step size zero, along with an error estimate for the
! extrapolation.
!
! Purpose
!
! Extrapolate the sequence of estimates yest built up over successive
! calls (indexed by iest, at squared step size xest) to produce a
! refined estimate yz, and the difference dy between the last two
! extrapolated estimates (used as the error estimate).
!
! Parameters
!
!   Inputs:
!     iest:  index of this call in the extrapolation sequence, starting
!            at 1 and incrementing by 1 each call (one call per
!            successively finer bs sub-step).
!     xest:  the value of the independent variable for this call
!            (squared step size h**2).
!     yest:  array of dependent-variable estimates at xest.
!     nv:    number of elements of yest in use.
!     nuse:  maximum number of previous estimates to use in the
!            extrapolation tableau.
!
!   Outputs:
!     yz:    extrapolated estimates.
!     dy:    difference between the last two extrapolated estimates,
!            used by the caller as an error estimate.
!
! Internal state (the extrapolation tableau x, d) persists across calls
! for a given integration step and is reset by starting a new sequence
! at iest = 1.
!===============================================================================
subroutine ratext(iest, xest, yest, yz, dy, nv, nuse)
    implicit none
    integer, parameter :: imax=11, nmax=15, ncol=7
    ! input
    integer, intent(in) :: iest, nv, nuse
    real(8), intent(in) :: xest
    real(8), intent(in), dimension(nmax) :: yest
    ! output
    real(8), intent(out), dimension(nmax) :: yz, dy
    ! internal (tableau, persists across calls within one extrapolation sequence)
    real(8), dimension(imax), save :: x
    real(8), dimension(nmax,ncol), save :: d
    real(8), dimension(ncol) :: fx
    integer :: j, k, m1
    real(8) :: yy, v, c, b, b1, ddy

    ! same as sr rzextr from numerical recipes, p.566.
    !
    ! save current independent variable.
    x(iest) = xest
    if (iest.eq.1) then
        do j = 1, nv
            yz(j) = yest(j)
            d(j,1) = yest(j)
            dy(j) = yest(j)
        end do
    else
        ! use at most nuse previous members.
        m1 = min(iest,nuse)
        do k = 1, m1-1
            fx(k+1) = x(iest-k)/xest
        end do
        ! evaluate next diagonal in tableau.
        do j = 1, nv
            yy = yest(j)
            v = d(j,1)
            c = yy
            d(j,1) = yy
            do k = 2, m1
                b1 = fx(k)*v
                b = b1 - c
                ! care needed to avoid division by zero.
                if (b.ne.0d0) then
                    b = (c - v)/b
                    ddy = c*b
                    c = b1*b
                else
                    ddy = v
                endif
                v = d(j,k)
                d(j,k) = ddy
                yy = yy + ddy
            end do
            dy(j) = ddy
            yz(j) = yy
        end do
    endif
    return
end
