!===============================================================================
! subroutine int2zone
!
! Bulirsch-Stoer integrator (based on Numerical Recipes bsstep) that
! advances the core and envelope angular momenta (sjc0,sje0) ->
! (sjc1,sje1) of a two-zone core/envelope model across one outer
! timestep [t0, t0+dtt]. Three physical processes act during the step:
!
!   1) Wind braking removes angular momentum from the envelope only,
!      at the same rate law as int1zone:
!        dje/dt|wind = -fc*fstruct*omega_e*min(omega_e,omega_crit)**(exw-1)
!   2) Mass transfer between the growing/shrinking convective envelope
!      and the radiative core carries angular momentum with it (the fcz
!      structural factor from setup), at the rotation rate of whichever
!      zone is losing mass to the other:
!        dj_transfer/dt = fcz * omega_source
!   3) If taucouple > 0, core-envelope coupling exchanges angular
!      momentum on that timescale following Denissenkov et al. (2010):
!        dj_coupling/dt = ie*ic/itot * (omega_c - omega_e) / taucouple
!
! Envelope: dje/dt = -dj/dt|wind + dj_transfer/dt + dj_coupling/dt
! Core:     djc/dt = -dj_transfer/dt - dj_coupling/dt
!
! The retry/extrapolation control flow (sub-step-count sequence, ratext
! extrapolation, eps-based acceptance) is shared with int1zone/int3zone via
! bsstep_mod -- see its header. What's specific to this routine is the
! trapezoidal sub-step physics (reset2/substep2 below), which extrapolates
! the zero-based mass-transfer and wind-loss sums (sumtran, sumdj) over the
! step -- not sjc1/sje1 themselves -- and reconstructs sjc1/sje1 from those
! sums only after extrapolation accepts. Track quantities are obtained at
! intermediate times by cubic-spline interpolation between points j-1 and
! j, using coefficients (yi, yie, ycz, ystr, ytau, ycen) precomputed by
! splinc in the caller.
!
! If none of the nseq sub-step counts converges by ii = imax, or an
! intermediate/final angular velocity comes out negative, the routine
! sets lok = .false. and returns without updating sjc1/sje1; the caller
! is responsible for handling non-convergence.
!
! Inputs:
!   sage,si,sie,staucz,fcz,fstruct,fcen   track arrays (setup builds
!                                          fstruct, fcen, fcz)
!   yi,yie,ycz,ystr,ytau,ycen             splinc coefficients for
!                                          si, sie, fcz, fstruct, staucz, fcen
!   j             index of the track point ending this step
!   t0, dtt       start time and duration of this step
!   excen, exw    loss-law exponents (see setup)
!   fcen          centrifugal-term structural factor (see setup)
!   fk2           centrifugal-suppression constant (see setup)
!   taucouple     core-envelope coupling timescale (yr), per-model array;
!                 interpolated across the step exactly like the other
!                 track arrays (fstruct, fcen, ...), so it can vary
!                 along the track
!   ycouple       splinc coefficients for taucouple
! Inputs/Outputs (start-of-step values, restored internally on each
! sub-step-count retry so the caller's copies are unchanged on return):
!   sjc0, sje0    core, envelope angular momentum at the start of the step
! Outputs:
!   sjc1, sje1    core, envelope angular momentum at the end of the step
!   lok           .true. if the step converged
!   iermsg        nonzero on failure: negative core/envelope i or omega
!                 (4118, 4172, 4242, 4267, 4271), or 4299 if every
!                 sub-step count in nseq was exhausted without the
!                 extrapolated error dropping below eps
!
! Bug fix (2026): the original source's sjc1 update read a variable
! spelled djtran11, not djtran1 -- a typo that referenced a never-
! assigned (uninitialized) variable instead of the end-of-substep
! transfer term computed just above it. Confirmed present in the
! original source and traced to real symptoms (surface rotation period
! collapsing toward zero within a few models). Fixed here to read
! djtran1.
!
! Bug fix (2026): the sub-step-count-exhaustion path used to leave
! iermsg at its initial 0 -- indistinguishable from success to a caller
! that only checks iermsg -- silently truncating the rest of a track's
! evolution (rotwind stops updating everything once one step fails).
! int3zone's header already documented and fixed this same gap (its
! 5299 code); fixed here too as 4299, after it was directly observed to
! mask a convergence failure on a real track (threezoneevol's two-zone
! fallback, which reuses int2zone via a reconstructed "sictot" combined
! inertia rather than the track's true total inertia).
!===============================================================================
subroutine int2zone(sage, si, sie, staucz, fcz, fstruct, yi, yie, &
                     ycz, ystr, ytau, sjc0, sje0, sjc1, sje1, j, t0, dtt, &
                     excen, exw, fcen, fk2, taucouple, ycouple, ycen, lok, iermsg)
    use params
    use constm
    use bsstep_mod
    implicit none
    !f2py integer, intent(aux) :: nmod
    ! inputs
    ! inputs from the track - cgs or solar units as noted
    ! sage = age (yr)
    ! si, sie = cgs moment of inertia total, envelope (core = si-sie)
    ! staucz = convective overturn timescale (sec)
    ! fcz, fstruct, fcen: see sr setup for description of ingredients
    real(8), intent(in), dimension(nmod) :: sage, fstruct, si, staucz, fcz, sie, fcen
    real(8), intent(in), dimension(nmod) :: taucouple
    real(8), intent(in) :: fk2, excen, exw, t0, dtt
    real(8), intent(inout) :: sjc0, sje0
    integer, intent(in) :: j
    ! spline interpolation coefficients
    real(8), intent(in), dimension(nmod) :: ycen, yi, ystr, ytau, yie, ycz, ycouple
    ! outputs - cgs units
    real(8), intent(out) :: sjc1, sje1
    logical, intent(out) :: lok
    integer, intent(out) :: iermsg

    ! below this fraction of the total moment of inertia, the core/envelope
    ! split (sic = si-sie) is a small residual of two much larger,
    ! independently spline-interpolated quantities and is not numerically
    ! trustworthy from cubic interpolation alone (see the linear-fallback
    ! block below).
    real(8), parameter :: sicfrac=1.0d-3

    ! state shared between reset2/substep2 below via host association:
    ! the "0"-suffixed quantities are the "start of this sub-step" track-
    ! interpolated/physical values, updated to the "end of this sub-step"
    ! ("1"-suffixed) values at the close of each substep2 call. djcoup0/1
    ! persist across the whole call (not reset per retry), matching the
    ! original routine exactly -- see the note in reset2 below.
    real(8) :: hh
    real(8) :: si0, si1, sie0, sie1, sic0, sic1
    real(8) :: sstr0, sstr1, scen0, scen1, srcz0, srcz1
    real(8) :: stau0, stau1, staumin, staumax
    real(8) :: srczmin, srczmax
    real(8) :: swe0, swc0
    real(8) :: w0, wc0
    real(8) :: djcoup0, djcoup1
    real(8) :: sjcsav, sjesav, tcouple0, tcouple1
    real(8) :: wtest
    real(8), dimension(2) :: y0v, yscalv, y1v
    logical :: converged

    ! initialize iermsg
    iermsg = 0
    ! convergence flag
    lok = .false.
    ! scale error in transport term by envelope j
    yscalv(1) = sje0
    ! tau(cz) can oscillate when close to a fully convective state
    ! to avoid numerical pathologies, restrict the range to within
    ! the endpoint values.
    staumin = min(staucz(j-1), staucz(j))
    staumax = max(staucz(j-1), staucz(j))
    ! fcz is a spline derivative (see setup.f90) and can overshoot sharply
    ! at a non-monotonic feature in sie (e.g. the RGB bump); restrict the
    ! interpolated mass-transfer rate to within the tabulated endpoint
    ! values for the same reason staucz is clamped above.
    srczmin = min(fcz(j-1), fcz(j))
    srczmax = max(fcz(j-1), fcz(j))
    ! scale error in loss rate by initial loss rate
    wtest = sje0/sie(j-1)
    if (lross) then
        w0 = wtest*staumin/soltau
    else
        w0 = wtest
    endif
    yscalv(2) = fstruct(j-1)*wtest*min(w0,wcrit)**(exw-1.0d0)
    ! djcoup0/djcoup1 persist across the whole call (all retries), not
    ! reset per retry -- see reset2 below.
    djcoup0 = 0.0d0
    djcoup1 = 0.0d0
    ! store starting angular momenta
    sjcsav = sjc0
    sjesav = sje0
    hh = sage(j)-sage(j-1)

    y0v = 0.0d0
    call bs_extrapolate(reset2, substep2, 2, dtt, t0, y0v, yscalv, y1v, converged)

    if (converged) then
        sjc0 = sjcsav
        sje0 = sjesav
        sje1 = sje0 + y1v(1) + y1v(2)
        ! compute transfer of angular momentum from core to envelope
        sjc1 = sjc0 - y1v(1)
        lok = .true.
        return
    endif
    ! exhausted every sub-step count in nseq without the extrapolated
    ! error dropping below eps. This used to leave iermsg at its initial
    ! 0 -- indistinguishable from success to a caller that only checks
    ! iermsg -- exactly the silent-failure gap int3zone's own header
    ! already documented and fixed (see its 5299 code). Confirmed to
    ! mask a real convergence failure on a real track (threezoneevol's
    ! two-zone fallback, reusing int2zone via a reconstructed "sictot"),
    ! so fixed here the same way: only set if substep2 didn't already
    ! record a more specific code before returning not-converged.
    if (iermsg.eq.0) iermsg = 4299
    lok = .false.
    return

contains

    ! (re-)establish the start-of-substep ("0"-suffixed) state at t0, for
    ! this retry, restoring sjc0/sje0 from the saved start-of-step values
    ! first (they were mutated by substep2 over the previous retry's
    ! sub-steps). djcoup0/djcoup1 are deliberately NOT reset here: the
    ! original routine only zeroes them once, before its first retry, so
    ! a retry can begin with whatever djcoup0 the previous (non-converged)
    ! retry's final sub-step left behind if this retry's first sub-step
    ! happens to have tcouple0 <= 0 (see substep2) -- preserved exactly
    ! as-is rather than "fixed", since faithfully reproducing the
    ! original's behavior is the point of this refactor.
    subroutine reset2(ok)
        logical, intent(out) :: ok
        real(8) :: a, b
        ok = .true.
        sjc0 = sjcsav
        sje0 = sjesav
        a = (sage(j)-t0)/hh
        b = (t0-sage(j-1))/hh
        ! total moment of inertia i
        si0 = a*si(j-1)+b*si(j)+ &
              ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
        ! envelope moment of inertia ie
        sie0 = a*sie(j-1)+b*sie(j)+ &
               ((a**3-a)*yie(j-1)+(b**3-b)*yie(j))*(hh**2)/6.0d0
        ! core moment of inertia = difference between total and envelope
        sic0 = si0 - sie0
        ! constant and structural terms in the loss prescription here
        sstr0 = a*fstruct(j-1)+b*fstruct(j)+ &
                ((a**3-a)*ystr(j-1)+(b**3-b)*ystr(j))*(hh**2)/6.0d0
        ! centrifugal term
        scen0 = a*fcen(j-1)+b*fcen(j)+ &
                ((a**3-a)*ycen(j-1)+(b**3-b)*ycen(j))*(hh**2)/6.0d0
        ! 2/3 radius of cz base ^2 * dm/dt
        srcz0 = a*fcz(j-1)+b*fcz(j)+ &
                ((a**3-a)*ycz(j-1)+(b**3-b)*ycz(j))*(hh**2)/6.0d0
        srcz0 = max(srcz0, srczmin)
        srcz0 = min(srcz0, srczmax)
        ! core-envelope coupling timescale (yr), interpolated like the other
        ! track arrays, then converted to gyr
        tcouple0 = 1.0d-9*(a*taucouple(j-1)+b*taucouple(j)+ &
                   ((a**3-a)*ycouple(j-1)+(b**3-b)*ycouple(j))*(hh**2)/6.0d0)
        ! surface and core angular velocities
        ! envelope i should never be zero, stop if true
        if (sie0.le.0.0d0) then
            iermsg = 4118
            ok = .false.
            return
        endif
        swe0 = sje0/sie0
        ! convective overturn timescale
        if (lross) then
            stau0 = a*staucz(j-1)+b*staucz(j)+ &
                    ((a**3-a)*ytau(j-1)+(b**3-b)*ytau(j))*(hh**2)/6.0d0
            stau0 = max(stau0,staumin)
            stau0 = min(stau0,staumax)
            if (iwind.eq.3) then
                w0 = swe0*stau0/soltau
                wc0 = wcrit
            else if (iwind.eq.2) then
                w0 = swe0
                wc0 = wcrit*soltau/stau0
            endif
        else
            w0 = swe0
            wc0 = wcrit
        endif
        ! core i can start as zero (or a numerically untrustworthy small
        ! residual, see the sicfrac note above), initialize with envelope
        ! omega if true
        if (sic0.le.0.0d0 .or. sic0.lt.sicfrac*si0) then
            swc0 = swe0
        else
            swc0 = sjc0/sic0
        endif
    end subroutine reset2

    ! advance sje0/sjc0 (and swe0/swc0/...) from t1-h1 to t1 by one
    ! trapezoidal sub-step, accumulating the mass-transfer sum into y(1)
    ! and the wind-loss sum into y(2) -- exactly as the original inline
    ! loop body.
    subroutine substep2(h1, t1, y, ok)
        real(8), intent(in) :: h1, t1
        real(8), intent(inout) :: y(:)
        logical, intent(out) :: ok
        real(8) :: a, b
        real(8) :: swe1, swc1, w1, wc1, fcc0, fcc1, fc0, fc1
        real(8) :: djdt0, djdt1, djtran0, djtran1

        ok = .true.
        a = (sage(j)-t1)/hh
        b = (t1-sage(j-1))/hh
        ! total moment of inertia i
        si1 = a*si(j-1)+b*si(j)+ &
              ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
        ! envelope moment of inertia ie
        sie1 = a*sie(j-1)+b*sie(j)+ &
               ((a**3-a)*yie(j-1)+(b**3-b)*yie(j))*(hh**2)/6.0d0
        ! core moment of inertia = difference between total and envelope
        sic1 = si1 - sie1
        ! although the starting core i can be zero, it should be nonzero
        ! at intermediate times. also fall back to linear interpolation
        ! if sic1 is merely a small (numerically untrustworthy) residual
        ! of si1 and sie1, not just when it is outright non-positive --
        ! see the sicfrac note above.
        if (sic1.le.0.0d0 .or. sic1.lt.sicfrac*si1) then
            ! switch to linear interpolation in moment of inertia
            si1 = si(j-1)+b*(si(j)-si(j-1))
            sie1 = sie(j-1)+b*(sie(j)-sie(j-1))
            sic1 = si1 - sie1
            if (sic1.le.0.0d0) then
                iermsg = 4172
                ok = .false.
                return
            endif
        endif
        ! constant and structural terms in the loss prescription here
        sstr1 = a*fstruct(j-1)+b*fstruct(j)+ &
                ((a**3-a)*ystr(j-1)+(b**3-b)*ystr(j))*(hh**2)/6.0d0
        ! centrifugal term
        scen1 = a*fcen(j-1)+b*fcen(j)+ &
                ((a**3-a)*ycen(j-1)+(b**3-b)*ycen(j))*(hh**2)/6.0d0
        ! radius of cz base, cm
        srcz1 = a*fcz(j-1)+b*fcz(j)+ &
                ((a**3-a)*ycz(j-1)+(b**3-b)*ycz(j))*(hh**2)/6.0d0
        srcz1 = max(srcz1, srczmin)
        srcz1 = min(srcz1, srczmax)
        ! core-envelope coupling timescale (yr), interpolated, then gyr
        tcouple1 = 1.0d-9*(a*taucouple(j-1)+b*taucouple(j)+ &
                   ((a**3-a)*ycouple(j-1)+(b**3-b)*ycouple(j))*(hh**2)/6.0d0)
        ! start of step dj/dt
        fcc0 = min(0.5d0, swe0**2*scen0)
        fc0 = (fk2/(fk2**2+fcc0)**0.5d0)**excen
        djdt0 = -fc0*sstr0*swe0*min(w0,wc0)**(exw-1.0d0)
        ! start of step estimate for moment of inertia of mass
        ! exchanged between envelope and core
        ! two branches: core grows (surface omega for material to core)
        ! or core shrinks (core omega for material to cz)
        if (srcz0.lt.0.0d0) then
            ! shrinking envelope; j+i added to core
            djtran0 = srcz0*h1*swe0
        else
            ! growing envelope - angular momentum lost from core
            ! copy end of step values to start of next step values
            djtran0 = srcz0*h1*swc0
        endif
        if (tcouple0.gt.0.0d0) then
            djcoup0 = sie0*sic0/si0*(swc0-swe0)/tcouple0
        endif
        ! trial corrected estimate for end of step dj/dt
        ! includes both loss and correction for change in i
        swe1 = (sje0 + djtran0 +h1*(djcoup0 - djdt0))/sie1
        if (swe1.lt.0.0d0) then
            print*, 'negative envelope omega, 1'
            ok = .false.
            return
        endif
        ! convective overturn timescale
        if (lross) then
            stau1 = a*staucz(j-1)+b*staucz(j)+ &
                    ((a**3-a)*ytau(j-1)+(b**3-b)*ytau(j))*(hh**2)/6.0d0
            stau1 = max(stau1,staumin)
            stau1 = min(stau1,staumax)
            if (iwind.eq.3) then
                w1 = swe1*stau1/soltau
                wc1 = wcrit
            else if (iwind.eq.2) then
                w1 = swe1
                wc1 = wcrit*soltau/stau1
            endif
        else
            w1 = swe1
            wc1 = wcrit
        endif
        fcc1 = min(0.5d0, swe1**2*scen1)
        fc1 = (fk2/(fk2**2+fcc1)**0.5d0)**excen
        swc1 = (sjc0 - djtran0 - h1*djcoup0)/sic1
        if (swc1.lt.0.0d0) then
            iermsg = 4242
            ok = .false.
            return
        endif
        ! update exchange term (end of timestep value)
        if (srcz0.lt.0.0d0) then
            djtran1 = srcz1*h1*swe1
        else
            djtran1 = srcz1*h1*swc1
        endif
        if (tcouple1.gt.0.0d0) then
            djcoup1 = sie1*sic1/si1*(swc1-swe1)/tcouple1
        endif
        ! update torque estimate (end of timestep value)
        djdt1 = -fc1*sstr1*swe1*min(w1,wc1)**(exw-1.0d0)
        sje1 = sje0+0.5d0*(h1*(djdt0+djdt1 + &
               djcoup0+djcoup1)+djtran0+djtran1)
        swe1 = sje1/sie1
        sjc1 = sjc0 - 0.5d0*(djtran0+djtran1)
        swc1 = sjc1/sic1
        if (swe1.lt.0.0d0) then
            iermsg = 4267
            ok = .false.
            return
        endif
        if (swc1.lt.0.0d0) then
            iermsg = 4271
            ok = .false.
            return
        endif

        ! do not accumulate wind loss if the star is above the critical rossby
        ! number (2*pi/omega_e/tau(cz) > rocrit*rosun), consistent with the
        ! rocrit cutoff applied to the one-zone model in int1zone. rocrit is
        ! expressed in solar units, so the absolute threshold is rocrit*rosun.
        if (2.0d0*cpi/swe1/staucz(j-1) .gt. rocrit*rosun) then
            y(2) = y(2) + 0.0d0
        else
            y(2) = y(2) + 0.5d0*h1*(djdt0+djdt1)
        endif

        y(1) = y(1) + 0.5d0*(djtran0+djtran1)+ &
               0.5d0*h1*(djcoup0+djcoup1)
        si0 = si1
        sie0 = sie1
        sic0 = sic1
        sstr0 = sstr1
        scen0 = scen1
        stau0 = stau1
        w0 = w1
        srcz0 = srcz1
        tcouple0 = tcouple1
        swe0 = swe1
        swc0 = swc1
        sje0 = sje1
        sjc0 = sjc1
    end subroutine substep2

end subroutine int2zone
