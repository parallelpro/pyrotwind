!===============================================================================
! subroutine splinc
!
! Computes the second-derivative coefficients of a natural cubic spline
! through n points (x, y), following Numerical Recipes Sec. 3.3
! ("Cubic Spline Interpolation", routine spline) with natural boundary
! conditions (zero second derivative at both endpoints).
!
! The output y2 is not the interpolated function itself: it is the table
! of second derivatives that callers combine with the endpoint values via
! the standard cubic-spline interpolation formula to evaluate y (or its
! derivative) at any point between the tabulated x values. Every splinc
! call in pyrotwind is paired with inline code elsewhere that does this
! combination for a specific quantity (moment of inertia, wind-loss
! structural factor, convective-overturn timescale, etc.) as a function
! of stellar age.
!
! Inputs:
!   x(1:n)  independent variable, strictly increasing (stellar age)
!   y(1:n)  dependent variable tabulated at x
!   n       number of points in use (<= nmod)
! Output:
!   y2(1:n) second-derivative coefficients of the interpolating spline
!===============================================================================
subroutine splinc(x, y, y2, n)
    use params, only: nmod
    implicit none
    !f2py integer, intent(aux) :: nmod
    integer, intent(in) :: n
    real(8), intent(in), dimension(nmod) :: x, y
    real(8), intent(out), dimension(nmod) :: y2
    ! local
    real(8), dimension(nmod) :: u
    real(8) :: sig, p, qn, un
    integer :: i, k

    ! natural spline: zero second derivative at the lower boundary
    y2(1) = 0.0d0
    u(1) = 0.0d0
    do i = 2, n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.0d0
        y2(i) = (sig-1.0d0)/p
        u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    ! natural spline: zero second derivative at the upper boundary
    qn = 0.0d0
    un = 0.0d0
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
    ! back-substitution
    do k = n-1, 1, -1
        y2(k) = y2(k)*y2(k+1)+u(k)
    end do
    return
end
