!===============================================================================
! subroutine int1zone
!
! Bulirsch-Stoer integrator (based on Numerical Recipes bsstep) that
! advances the total angular momentum sj0 -> sj1 of a single rigidly-
! rotating reservoir (solid body, or a track segment with no radiative
! core) across one outer timestep [t0, t0+dtt], under wind-braking torque
! dj/dt = -fc*fstruct*omega*min(omega,omega_crit)**(exw-1).
!
! The retry/extrapolation control flow (sub-step-count sequence, ratext
! extrapolation, eps-based acceptance) is shared with int2zone/int3zone via
! bsstep_mod -- see its header. What's specific to this routine is the
! trapezoidal sub-step physics (reset1/substep1 below), which extrapolates
! the running angular momentum sj1 itself (starting from sj0 each retry),
! obtaining moment of inertia, structural torque factor, centrifugal
! factor, and (if lross) convective-overturn timescale at intermediate
! times by cubic-spline interpolation between the tabulated track points
! j-1 and j (coefficients yi, ystr, ycen, ytau were precomputed by splinc
! in the caller).
!
! If none of the nseq sub-step counts converges to within eps by
! ii = imax, the routine simply falls out of the loop and returns
! whatever sj1 was last computed, without flagging an error via iermsg.
! This matches the routine's original behavior (its Numerical-Recipes-
! style adaptive retry/step-halving path is present in Numerical
! Recipes but is not used here).
!
! Inputs:
!   sage, si, fstruct, staucz, fcen  track arrays (see setup for how the
!                                     latter three are built)
!   yi, ystr, ytau, ycen             splinc second-derivative coefficients
!                                     for si, fstruct, staucz, fcen
!   j           index of the track point ending this step
!   t0, dtt     start time and duration of this step
!   excen, exw  loss-law exponents (see setup)
!   fk2         centrifugal-suppression constant (see setup)
!   sj0         angular momentum at the start of the step
! Inputs/Outputs:
!   lrocrit     the calling track's rossby-cutoff flag (see solidevol):
!               left unchanged unless this step is the one where the
!               critical rossby number is first crossed, in which case
!               it is set .true. here. Local to the caller's track, not
!               module state, so it never leaks between tracks.
! Outputs:
!   sj1         angular momentum at the end of the step
!   iermsg      nonzero on failure (negative intermediate omega)
!===============================================================================
subroutine int1zone(sage, si, fstruct, staucz, yi, ystr, ytau, sj0, sj1, &
                     j, t0, dtt, excen, exw, fcen, fk2, ycen, lrocrit, iermsg)
    use constm
    use params
    use bsstep_mod
    implicit none
    !f2py integer, intent(aux) :: nmod
    ! inputs from the track - cgs or solar units as noted
    ! sage = age (yr), si = cgs moment of inertia total
    ! staucz = convective overturn timescale (sec)
    ! fstruct = rotation-independent terms used in loss law
    ! see sr setup for description of ingredients
    real(8), intent(in), dimension(nmod) :: sage, si, fstruct, staucz, fcen
    real(8), intent(in) :: excen, exw, fk2, t0, dtt, sj0
    integer, intent(in) :: j
    real(8), intent(in), dimension(nmod) :: yi, ystr, ytau, ycen
    ! inout
    logical, intent(inout) :: lrocrit
    ! outputs
    real(8), intent(out) :: sj1
    integer, intent(out) :: iermsg

    ! state shared between reset1/substep1 below via host association:
    ! si0/fstr0/fcen0/stau0/w0/wc0/sw0 are the "start of this sub-step"
    ! track-interpolated/physical values, updated to the "end of this
    ! sub-step" (si1 etc.) values at the close of each substep1 call.
    real(8) :: hh, si0, si1, stau0, stau1, staumin, staumax
    real(8) :: w0, w1, wc0, wc1, ro0, ro1
    real(8) :: fstr0, fstr1, fcen0, fcen1, fcc0, fcc1, fc0, fc1
    real(8) :: sw0, sw1, djdt0, djdt1, sjnew
    real(8), dimension(1) :: y0v, yscalv, y1v
    logical :: converged

    iermsg = 0
    ! tau(cz) can oscillate when close to a fully convective state
    ! to avoid numerical pathologies, restrict the range to within
    ! the endpoint values.
    staumin = min(staucz(j-1), staucz(j))
    staumax = max(staucz(j-1), staucz(j))
    hh = sage(j)-sage(j-1)

    y0v(1) = sj0
    yscalv(1) = sj0
    call bs_extrapolate(reset1, substep1, 1, dtt, t0, y0v, yscalv, y1v, converged)

    if (converged) then
        sj1 = y1v(1)
        ! check whether the critical rossby number was crossed during this step;
        ! if so, interpolate back to the crossing point and flag lrocrit so the
        ! caller stops applying wind loss for the rest of this track. rocrit is
        ! expressed in solar units, so the absolute threshold is rocrit*rosun.
        ! si1/stau1 here are whatever reset1/substep1 last left them at, i.e.
        ! the end of the accepted retry's final sub-step.
        if (2.0d0*cpi*si1/sj1/stau1 .gt. rocrit*rosun) then
            ro1 = 2.0d0*cpi*si1/sj1/stau1
            ro0 = 2.0d0*cpi*si0/sj0/stau0
            sjnew = (rocrit*rosun-ro0)*((sj1-sj0)/(ro1-ro0)) + sj0
            sj1 = sjnew
            lrocrit = .true.
        endif
    endif
    ! if not converged (or aborted early by substep1), sj1 already holds
    ! whatever raw value substep1 last left it at, via host association --
    ! matching the original routine's fallback behavior exactly.
    return

contains

    ! (re-)establish the start-of-substep ("0"-suffixed) state at t0, for
    ! this retry. sj1 also resets to sj0 here, since it doubles as both
    ! the physical running state and the quantity bs_extrapolate is
    ! extrapolating (see bsstep_mod's header on y0). int1zone has no
    ! reset-time failure mode, so ok is always .true. here.
    subroutine reset1(ok)
        logical, intent(out) :: ok
        real(8) :: a, b
        ok = .true.
        sj1 = sj0
        a = (sage(j)-t0)/hh
        b = (t0-sage(j-1))/hh
        si0 = a*si(j-1)+b*si(j)+ &
              ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
        fstr0 = a*fstruct(j-1)+b*fstruct(j)+ &
                ((a**3-a)*ystr(j-1)+(b**3-b)*ystr(j))*(hh**2)/6.0d0
        fcen0 = a*fcen(j-1)+b*fcen(j)+ &
                ((a**3-a)*ycen(j-1)+(b**3-b)*ycen(j))*(hh**2)/6.0d0
        sw0 = sj0/si0
        if (lross) then
            stau0 = a*staucz(j-1)+b*staucz(j)+ &
                    ((a**3-a)*ytau(j-1)+(b**3-b)*ytau(j))*(hh**2)/6.0d0
            stau0 = max(stau0,staumin)
            stau0 = min(stau0,staumax)
            if (iwind.eq.3) then
                w0 = sw0*stau0/soltau
                wc0 = wcrit
            else if (iwind.eq.2) then
                w0 = sw0
                wc0 = wcrit*soltau/stau0
            endif
        else
            w0 = sw0
            wc0 = wcrit
        endif
    end subroutine reset1

    ! advance sj1 (and its mirror y(1)) from t1-h1 to t1 by one trapezoidal
    ! sub-step, exactly as the original inline loop body.
    subroutine substep1(h1, t1, y, ok)
        real(8), intent(in) :: h1, t1
        real(8), intent(inout) :: y(:)
        logical, intent(out) :: ok
        real(8) :: a, b
        ok = .true.
        a = (sage(j)-t1)/hh
        b = (t1-sage(j-1))/hh
        si1 = a*si(j-1)+b*si(j)+ &
              ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
        fstr1 = a*fstruct(j-1)+b*fstruct(j)+ &
                ((a**3-a)*ystr(j-1)+(b**3-b)*ystr(j))*(hh**2)/6.0d0
        fcen1 = a*fcen(j-1)+b*fcen(j)+ &
                ((a**3-a)*ycen(j-1)+(b**3-b)*ycen(j))*(hh**2)/6.0d0
        ! start of step dj/dt
        fcc0 = min(0.5d0, sw0**2*fcen0)
        ! centrifugal-suppression factor currently fixed at 1 (see fc0/fc1 below);
        ! the fk2/fcc0-based form is kept commented as the alternative.
        !fc0 = (fk2/(fk2**2+fcc0)**0.5d0)**excen
        fc0 = 1.0d0
        djdt0 = -fc0*fstr0*sw0*min(w0,wc0)**(exw-1.0d0)
        ! trial corrected estimate for end of step dj/dt
        ! includes both loss and correction for change in i
        sw1 = (sj1 - djdt0*h1)/si1
        if (sw1.lt.0.0d0) then
            iermsg = 3102
            ok = .false.
            return
        endif
        if (lross) then
            stau1 = a*staucz(j-1)+b*staucz(j)+ &
                    ((a**3-a)*ytau(j-1)+(b**3-b)*ytau(j))*(hh**2)/6.0d0
            stau1 = max(stau1,staumin)
            stau1 = min(stau1,staumax)
            if (iwind.eq.3) then
                w1 = sw1*stau1/soltau
                wc1 = wcrit
            else if (iwind.eq.2) then
                w1 = sw1
                wc1 = wcrit*soltau/stau1
            endif
        else
            w1 = sw1
            wc1 = wcrit
        endif
        fcc1 = min(0.5d0, sw1**2*fcen1)
        ! centrifugal-suppression factor currently fixed at 1, see above.
        !fc1 = (fk2/(fk2**2+fcc1)**0.5d0)**excen
        fc1 = 1.0d0
        djdt1 = -fc1*fstr1*sw1*min(w1,wc1)**(exw-1.0d0)
        sj1 = sj1 + 0.5d0*h1*(djdt0+djdt1)
        y(1) = sj1

        ! copy end of step values to start of next step values
        si0 = si1
        fstr0 = fstr1
        fcen0 = fcen1
        stau0 = stau1
        w0 = w1
        sw0 = sj1/si1

        if (sw0.lt.0.0d0) then
            iermsg = 3143
            ok = .false.
            return
        endif
    end subroutine substep1

end subroutine int1zone
