!===============================================================================
! module bsstep_mod
!
! Generic Bulirsch-Stoer extrapolation driver shared by int1zone/int2zone/
! int3zone, factored out so the retry/extrapolation control flow (nseq,
! imax, nuse, eps, and the ratext-driving loop) exists exactly once instead
! of being duplicated per zone count. The physics -- what happens during one
! trapezoidal sub-step, and how to (re)establish the start-of-step state
! before each retry -- stays with the caller, passed in as two procedures:
!
!   reset(ok)            Called once per retry (ii), before its sequence of
!                        sub-steps. (Re-)establishes any track-quantity
!                        interpolation at t0 that a fresh retry needs. It
!                        is expected to be an internal procedure of the
!                        caller, reaching the caller's own locals (track
!                        arrays, spline coefficients, the physical state
!                        being advanced) via host association -- see
!                        int1zone/int2zone/int3zone. Sets ok = .false. to
!                        abort the whole bs_extrapolate call outright (not
!                        just this retry) -- e.g. int2zone's "envelope i
!                        should never be zero" check, which is a property
!                        of t0 alone and so would fail identically on
!                        every retry. The caller's reset is expected to
!                        have already recorded its own iermsg (via host
!                        association) before returning ok=.false.
!   substep(h1,t1,y,ok)  Called nseq(ii) times per retry, once per
!                        sub-step, advancing the caller's own physical
!                        state (again via host association) from t1-h1 to
!                        t1, and updating y(1:nv) -- the quantity actually
!                        being extrapolated to h->0 (see y0 below). Sets
!                        ok = .false. to abort this call (e.g. a negative
!                        intermediate omega); the caller's substep is
!                        expected to have already recorded its own iermsg
!                        (via host association) before returning ok=.false.
!
! y0(nv) is the value y is reset to at the start of every retry. This
! matters: int1zone extrapolates its running angular momentum itself
! (y0 = [sj0]), while int2zone/int3zone extrapolate zero-based transfer/
! loss sums added back to the start-of-step angular momenta only after
! extrapolation (y0 = 0). These are not interchangeable -- ratext performs
! a *rational* extrapolation, and shifting the extrapolated quantity by a
! constant is not guaranteed to commute with that (unlike for a polynomial
! extrapolation), so which quantity is extrapolated is preserved exactly
! as in the original separate routines, not standardized to one or the
! other.
!
! bs_imax/nseq (2026): extended from the original 11 levels (max sub-step
! count 96) to 15 levels (max 384). The original 11-level sequence still
! comes first unchanged, so any step that already converged within it
! still converges at the exact same retry and gets the exact same
! extrapolated result -- this is a strictly additive change, never a
! behavior change, for every case that isn't hitting the ceiling. It was
! added after finding a real track (an ~11 Gyr subgiant/RGB evolution)
! where int2zone, called from threezoneevol's two-zone fallback, ran out
! of retries and failed to converge (exhausting all 11 levels without the
! extrapolated error dropping below eps) at one specific step -- which,
! since rotwind stops updating everything once a step fails, truncated
! the rest of that track's rotation-rate evolution to zero.
!
! That failure was chased at length before landing here: the same step,
! fed bit-identical time steps and even a bit-identical total moment of
! inertia (temporarily threading the caller's true si through
! threezoneevol instead of its reconstructed sic+sie sum, as an
! experiment) still failed identically, which rules out a reconstruction/
! precision bug in the caller. drevol, solving the same physical step
! through int2zone via a different caller, converges fine. So the two
! callers' trajectories must be diverging from small floating-point-level
! differences building up over ~150 earlier steps (confirmed accumulating
! to a ~0.05% difference in the angular momenta feeding int2zone by the
! failing step) until the failing one lands somewhere int2zone's original
! 11-level ceiling isn't quite enough -- i.e. genuine, if narrow,
! numerical sensitivity in the coupled ODE system at that point, not a
! logic bug. Extending the ceiling gives it the extra refinement it
! needs; the track converges cleanly (iermsg=0) through all ~1500 models
! with this change, with no measurable performance cost (the extra levels
! are only ever reached by steps that would otherwise have failed).
!===============================================================================
module bsstep_mod
    implicit none
    private
    public :: bs_extrapolate, bs_reset_iface, bs_substep_iface
    public :: bs_imax, bs_nuse, bs_nmax, bs_eps

    integer, parameter :: bs_imax = 15, bs_nuse = 7, bs_nmax = 15
    real(8), parameter :: bs_eps = 1.0d-5

    abstract interface
        subroutine bs_reset_iface(ok)
            logical, intent(out) :: ok
        end subroutine bs_reset_iface

        subroutine bs_substep_iface(h1, t1, y, ok)
            real(8), intent(in) :: h1, t1
            real(8), intent(inout) :: y(:)
            logical, intent(out) :: ok
        end subroutine bs_substep_iface
    end interface

contains

    !---------------------------------------------------------------------
    ! Advance y0(nv) across [t0, t0+dtt] via repeated whole-step retries at
    ! increasing sub-step counts (nseq), extrapolating each retry's
    ! end-of-step y to h->0 via ratext, and accepting once the estimated
    ! error (scaled by yscal) drops below eps. See module header for the
    ! reset/substep contract.
    !---------------------------------------------------------------------
    subroutine bs_extrapolate(reset, substep, nv, dtt, t0, y0, yscal, y1, converged)
        procedure(bs_reset_iface) :: reset
        procedure(bs_substep_iface) :: substep
        integer, intent(in) :: nv
        real(8), intent(in) :: dtt, t0
        real(8), intent(in), dimension(nv) :: y0, yscal
        real(8), intent(out), dimension(nv) :: y1
        logical, intent(out) :: converged

        integer, dimension(bs_imax) :: nseq = &
            (/2,4,6,8,12,16,24,32,48,64,96,144,192,256,384/)
        real(8) :: h1, t1, xest, errmax
        real(8), dimension(bs_nmax) :: yest, yout, yerr
        real(8), dimension(nv) :: y
        integer :: ii, jj, kk
        logical :: ok

        converged = .false.
        do ii = 1, bs_imax
            call reset(ok)
            if (.not.ok) return
            h1 = dtt/dble(nseq(ii))
            t1 = t0
            y = y0
            do jj = 1, nseq(ii)
                t1 = t1 + h1
                call substep(h1, t1, y, ok)
                if (.not.ok) return
            end do
            xest = (dtt/dble(nseq(ii)))**2
            yest(1:nv) = y(1:nv)
            call ratext(ii, xest, yest, yout, yerr, nv, bs_nuse)
            errmax = 0.0d0
            do kk = 1, nv
                errmax = dmax1(errmax, dabs(yerr(kk)/yscal(kk)))
            end do
            errmax = errmax/bs_eps
            if (errmax.lt.1.0d0) then
                y1(1:nv) = yout(1:nv)
                converged = .true.
                return
            endif
        end do
    end subroutine bs_extrapolate

end module bsstep_mod
