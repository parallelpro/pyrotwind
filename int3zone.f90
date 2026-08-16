!===============================================================================
! subroutine int3zone
!
! Bulirsch-Stoer integrator (based on Numerical Recipes bsstep), the
! three-zone generalization of int2zone: advances the core, middle, and
! envelope angular momenta (sjcore0,sjc0,sje0) -> (sjcore1,sjc1,sje1)
! across one outer timestep [t0, t0+dtt]. Zones are chained core-middle-
! envelope; there is no direct core-envelope interaction. Four physical
! processes act during the step:
!
!   1) Wind braking removes angular momentum from the envelope only, at
!      the same rate law as int1zone/int2zone:
!        dje/dt|wind = -fc*fstruct*omega_e*min(omega_e,omega_crit)**(exw-1)
!   2) Mass transfer across the middle/envelope boundary, exactly as in
!      int2zone (fcz, built from -d(i_envelope)/dt in setup):
!        dj_transfer_me/dt = fcz * omega_source
!   3) Mass transfer across the core/middle boundary, built the same way
!      from -d(i_core)/dt (fcm, see setup): dj_transfer_cm/dt = fcm *
!      omega_source. Because the core has only this one boundary, its own
!      moment-of-inertia history is an uncontaminated proxy for it (the
!      middle zone's own di/dt is not used here, since it is affected by
!      both boundaries at once).
!   4) Core-envelope-style torque coupling (Denissenkov et al. 2010)
!      applied independently at each boundary, each on its own
!      timescale:
!        dj_coupling_cm/dt = icore*ic/(icore+ic) * (omega_c - omega_core) / taucouplecm
!        dj_coupling_me/dt = ic*ie/(ic+ie) * (omega_c - omega_e) / taucouple
!
! Core:     djcore/dt = +dj_transfer_cm/dt + dj_coupling_cm/dt
! Middle:   djc/dt    = -dj_transfer_cm/dt - dj_coupling_cm/dt
!                       -dj_transfer_me/dt - dj_coupling_me/dt
! Envelope: dje/dt     = -dj/dt|wind + dj_transfer_me/dt + dj_coupling_me/dt
!
! Note on conservation: unlike int2zone's final update (which folds the
! transfer term into both zones but only folds the coupling term into
! the envelope, not the core -- an asymmetry inherited from the original
! source and not yet resolved there), this routine includes both the
! transfer and coupling terms symmetrically in every zone's final
! update, so total angular momentum core+middle+envelope changes only by
! the wind-loss term over a step, exactly as the equations above imply.
!
! As in int1zone/int2zone, the routine repeats the whole timestep at a
! fixed sequence of sub-step counts nseq(1), nseq(2), ... using the
! trapezoidal rule, and uses ratext to extrapolate the sequence of
! end-of-step sums (core/middle transfer+coupling, middle/envelope
! transfer+coupling, wind loss) to the h -> 0 limit, accepting the
! result once the estimated error drops below eps. Track quantities are
! obtained at intermediate times by cubic-spline interpolation between
! points j-1 and j, using coefficients (yicore, yic, yie, ycm, ycz,
! ystr, ytau, ycen) precomputed by splinc in the caller.
!
! Error scaling: int2zone scales its combined transfer+coupling error
! against the envelope's own angular momentum (much the larger
! reservoir), which was found to under-constrain accuracy right when a
! small reservoir is draining fast (see threezoneevol/drevol notes on
! the core-angular-momentum-drain failure mode). Here each boundary's
! error is instead scaled against the smaller of the two reservoirs it
! moves angular momentum between, so the smaller zone's budget always
! gets the tighter (not the looser) tolerance.
!
! Inputs:
!   sage,sicore,sic,sie,staucz,fcm,fcz,fstruct   track arrays (setup
!                                                 builds fstruct, fcen,
!                                                 fcz, fcm)
!   yicore,yic,yie,ycm,ycz,ystr,ytau,ycen        splinc coefficients for
!                                                 sicore, sic, sie, fcm,
!                                                 fcz, fstruct, staucz,
!                                                 fcen
!   j             index of the track point ending this step
!   t0, dtt       start time and duration of this step
!   excen, exw    loss-law exponents (see setup)
!   fcen          centrifugal-term structural factor (see setup)
!   fk2           centrifugal-suppression constant (see setup)
!   taucouplecm   core-middle coupling timescale (yr), per-model array
!   taucouple     middle-envelope coupling timescale (yr), per-model
!                 array (same quantity/name as int2zone's taucouple)
!   ycouplecm, ycouple   splinc coefficients for the two timescales above
! Inputs/Outputs (start-of-step values, restored internally on each
! sub-step-count retry so the caller's copies are unchanged on return):
!   sjcore0, sjc0, sje0   core, middle, envelope angular momentum at the
!                         start of the step
! Outputs:
!   sjcore1, sjc1, sje1   core, middle, envelope angular momentum at the
!                         end of the step
!   lok           .true. if the step converged
!   iermsg        nonzero on failure (negative zone i or omega; codes in
!                 the 5100-5299 range, distinct from int2zone's 4100s/
!                 4200s so a failing track can be traced to the routine
!                 that raised it)
!===============================================================================
subroutine int3zone(sage, sicore, sic, sie, staucz, fcm, fcz, fstruct, &
                     yicore, yic, yie, ycm, ycz, ystr, ytau, &
                     sjcore0, sjc0, sje0, sjcore1, sjc1, sje1, j, t0, dtt, &
                     excen, exw, fcen, fk2, taucouplecm, taucouple, &
                     ycouplecm, ycouple, ycen, lok, iermsg)
    use params
    use constm
    implicit none
    !f2py integer, intent(aux) :: nmod
    integer, parameter :: nmax=15
    real(8), parameter :: one=1.0d0
    ! inputs from the track - cgs or solar units as noted
    ! sage = age (yr)
    ! sicore, sic, sie = cgs moment of inertia, core/middle/envelope
    ! staucz = convective overturn timescale (sec)
    ! fcm, fcz, fstruct, fcen: see sr setup for description of ingredients
    real(8), intent(in), dimension(nmod) :: sage, fstruct, sicore, sic, &
        sie, staucz, fcm, fcz, fcen
    real(8), intent(in), dimension(nmod) :: taucouplecm, taucouple
    real(8), intent(in) :: fk2, excen, exw, t0, dtt
    real(8), intent(inout) :: sjcore0, sjc0, sje0
    integer, intent(in) :: j
    ! spline interpolation coefficients
    real(8), intent(in), dimension(nmod) :: ycen, yicore, yic, yie, &
        ystr, ytau, ycm, ycz, ycouplecm, ycouple
    ! outputs - cgs units
    real(8), intent(out) :: sjcore1, sjc1, sje1
    logical, intent(out) :: lok
    integer, intent(out) :: iermsg

    ! local
    integer, dimension(11) :: nseq = (/2,4,6,8,12,16,24,32,48,64,96/)
    integer, parameter :: imax=11, nuse=7
    real(8), parameter :: eps=1.0d-5
    ! below this fraction of the core's own inertia at the end of the
    ! step, we are close enough to its formation instant that fcm's
    ! cubic-spline derivative (built globally across all track points,
    ! see setup.f90) is unreliable right at the curvature discontinuity
    ! where sicore turns on -- the same class of overshoot documented for
    ! sic/sie's sicfrac fallback in int2zone, just in the transfer-rate
    ! term instead of the inertia term. below this threshold, srcm falls
    ! back to a simple secant-slope estimate for the whole step instead.
    real(8), parameter :: sicorefrac=1.0d-2

    ! scalars
    real(8) :: t1, hh, h1, a, b
    real(8) :: sicore0v, sicore1v, sic0v, sic1v, sie0v, sie1v
    real(8) :: sstr0, sstr1, scen0, scen1
    real(8) :: srcm0, srcm1, srcmmin, srcmmax
    real(8) :: srcz0, srcz1, srczmin, srczmax
    real(8) :: stau0, stau1, staumin, staumax
    real(8) :: swe0, swe1, swc0, swc1, swcore0, swcore1
    real(8) :: swe1t, swc1t, swcore1t
    real(8) :: w0, w1, wc0, wc1
    real(8) :: fcc0, fcc1, fc0, fc1
    real(8) :: djdt0, djdt1
    real(8) :: djtrancm0, djtrancm1, djtranme0, djtranme1
    real(8) :: djcoupcm0, djcoupcm1, djcoupme0, djcoupme1
    real(8) :: sjcoresav, sjcsav, sjesav
    real(8) :: tcouplecm0, tcouplecm1, tcouple0, tcouple1
    real(8) :: sumtrancm, sumtranme, sumdj, xest, errmax, wtest
    real(8), dimension(3,11) :: ytest
    real(8), dimension(nmax) :: yest, yout, yerr, err
    real(8), dimension(3) :: yscal
    integer :: ii, jj, kk, nv
    logical :: lcmlin
    real(8) :: srcmlin

    iermsg = 0
    lok = .false.
    ! tau(cz) can oscillate when close to a fully convective state; restrict
    ! the range to within the endpoint values, exactly as in int2zone.
    staumin = min(staucz(j-1), staucz(j))
    staumax = max(staucz(j-1), staucz(j))
    ! fcm/fcz are spline derivatives (see setup.f90) and can overshoot
    ! sharply at a non-monotonic feature (e.g. the rgb bump); restrict
    ! each to its own tabulated endpoint values, exactly as int2zone does
    ! for fcz.
    srcmmin = min(fcm(j-1), fcm(j))
    srcmmax = max(fcm(j-1), fcm(j))
    srczmin = min(fcz(j-1), fcz(j))
    srczmax = max(fcz(j-1), fcz(j))
    ! decide once per step (not per sub-step -- the kink is a property of
    ! this step's bracketing track points, not something that should
    ! flip-flop mid-step) whether the core is still close enough to its
    ! formation instant to fall back to a linear estimate of srcm
    lcmlin = (sicore(j-1) .lt. sicorefrac*sicore(j))
    if (lcmlin) then
        srcmlin = -(sicore(j)-sicore(j-1))/(sage(j)-sage(j-1))
    endif
    ! scale error in wind loss rate by initial loss rate, as in int2zone
    wtest = sje0/sie(j-1)
    if (lross) then
        w0 = wtest*staumin/soltau
    else
        w0 = wtest
    endif
    do kk = 1, 11
        ytest(1,kk) = 0.0d0
        ytest(2,kk) = 0.0d0
        ytest(3,kk) = 0.0d0
    end do
    djcoupcm0 = 0.0d0
    djcoupcm1 = 0.0d0
    djcoupme0 = 0.0d0
    djcoupme1 = 0.0d0
    ! store starting angular momenta
    sjcoresav = sjcore0
    sjcsav = sjc0
    sjesav = sje0
    do ii = 1, imax
        sjcore0 = sjcoresav
        sjc0 = sjcsav
        sje0 = sjesav
        h1 = dtt/float(nseq(ii))
        t1 = t0
        hh = sage(j)-sage(j-1)
        ! interpolate between points j and j-1 for structure variables
        a = (sage(j)-t0)/hh
        b = (t0-sage(j-1))/hh
        sicore0v = a*sicore(j-1)+b*sicore(j)+ &
                   ((a**3-a)*yicore(j-1)+(b**3-b)*yicore(j))*(hh**2)/6.0d0
        sic0v = a*sic(j-1)+b*sic(j)+ &
                ((a**3-a)*yic(j-1)+(b**3-b)*yic(j))*(hh**2)/6.0d0
        sie0v = a*sie(j-1)+b*sie(j)+ &
                ((a**3-a)*yie(j-1)+(b**3-b)*yie(j))*(hh**2)/6.0d0
        sstr0 = a*fstruct(j-1)+b*fstruct(j)+ &
                ((a**3-a)*ystr(j-1)+(b**3-b)*ystr(j))*(hh**2)/6.0d0
        scen0 = a*fcen(j-1)+b*fcen(j)+ &
                ((a**3-a)*ycen(j-1)+(b**3-b)*ycen(j))*(hh**2)/6.0d0
        if (lcmlin) then
            srcm0 = srcmlin
        else
            srcm0 = a*fcm(j-1)+b*fcm(j)+ &
                    ((a**3-a)*ycm(j-1)+(b**3-b)*ycm(j))*(hh**2)/6.0d0
            srcm0 = max(srcm0, srcmmin)
            srcm0 = min(srcm0, srcmmax)
        endif
        srcz0 = a*fcz(j-1)+b*fcz(j)+ &
                ((a**3-a)*ycz(j-1)+(b**3-b)*ycz(j))*(hh**2)/6.0d0
        srcz0 = max(srcz0, srczmin)
        srcz0 = min(srcz0, srczmax)
        tcouplecm0 = 1.0d-9*(a*taucouplecm(j-1)+b*taucouplecm(j)+ &
                     ((a**3-a)*ycouplecm(j-1)+(b**3-b)*ycouplecm(j))*(hh**2)/6.0d0)
        tcouple0 = 1.0d-9*(a*taucouple(j-1)+b*taucouple(j)+ &
                   ((a**3-a)*ycouple(j-1)+(b**3-b)*ycouple(j))*(hh**2)/6.0d0)
        ! envelope i should never be zero, stop if true
        if (sie0v.le.0.0d0) then
            iermsg = 5118
            return
        endif
        swe0 = sje0/sie0v
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
        ! middle zone should always be established by the time int3zone is
        ! called (a core zone only forms once a middle zone already
        ! exists), but keep the same defensive fallback pattern as
        ! int2zone's core-i check, just in case.
        if (sic0v.le.0.0d0) then
            swc0 = swe0
        else
            swc0 = sjc0/sic0v
        endif
        ! core i can start as zero right when the core first forms (see
        ! threezoneevol); initialize with the middle zone's omega if so.
        if (sicore0v.le.0.0d0) then
            swcore0 = swc0
        else
            swcore0 = sjcore0/sicore0v
        endif
        ! error scale for the wind-loss sum, as in int2zone
        yscal(3) = fstruct(j-1)*wtest*min(w0,wcrit)**(exw-1.0d0)
        ! error scales for the two transfer+coupling sums: scale against
        ! the smaller of the two reservoirs each boundary moves angular
        ! momentum between (see header note on error scaling). right when
        ! the core first forms, sjcore0 is exactly zero by construction
        ! (see threezoneevol) -- there is no existing core budget to
        ! protect yet, so fall back to scaling against the middle zone's
        ! budget alone rather than dividing the error estimate by zero
        ! (which would read as an infinite error and never converge, no
        ! matter how well-converged the extrapolated value actually is).
        if (sjcore0.gt.0.0d0) then
            yscal(1) = min(sjcore0, sjc0)
        else
            yscal(1) = sjc0
        endif
        yscal(2) = min(sjc0, sje0)
        sumtrancm = 0.0d0
        sumtranme = 0.0d0
        sumdj = 0.0d0
        do jj = 1, nseq(ii)
            t1 = t1 + h1
            a = (sage(j)-t1)/hh
            b = (t1-sage(j-1))/hh
            sicore1v = a*sicore(j-1)+b*sicore(j)+ &
                       ((a**3-a)*yicore(j-1)+(b**3-b)*yicore(j))*(hh**2)/6.0d0
            sic1v = a*sic(j-1)+b*sic(j)+ &
                    ((a**3-a)*yic(j-1)+(b**3-b)*yic(j))*(hh**2)/6.0d0
            sie1v = a*sie(j-1)+b*sie(j)+ &
                    ((a**3-a)*yie(j-1)+(b**3-b)*yie(j))*(hh**2)/6.0d0
            ! core i can legitimately be very small right after formation;
            ! cubic-spline overshoot can carry it slightly negative there.
            ! fall back to linear interpolation in that case (no cancellation
            ! risk here since sicore is a direct input, not a subtraction
            ! residual, so a plain non-positive check is enough).
            if (sicore1v.le.0.0d0) then
                sicore1v = sicore(j-1)+b*(sicore(j)-sicore(j-1))
                if (sicore1v.le.0.0d0) then
                    iermsg = 5171
                    return
                endif
            endif
            if (sic1v.le.0.0d0) then
                sic1v = sic(j-1)+b*(sic(j)-sic(j-1))
                if (sic1v.le.0.0d0) then
                    iermsg = 5172
                    return
                endif
            endif
            sstr1 = a*fstruct(j-1)+b*fstruct(j)+ &
                    ((a**3-a)*ystr(j-1)+(b**3-b)*ystr(j))*(hh**2)/6.0d0
            scen1 = a*fcen(j-1)+b*fcen(j)+ &
                    ((a**3-a)*ycen(j-1)+(b**3-b)*ycen(j))*(hh**2)/6.0d0
            if (lcmlin) then
                srcm1 = srcmlin
            else
                srcm1 = a*fcm(j-1)+b*fcm(j)+ &
                        ((a**3-a)*ycm(j-1)+(b**3-b)*ycm(j))*(hh**2)/6.0d0
                srcm1 = max(srcm1, srcmmin)
                srcm1 = min(srcm1, srcmmax)
            endif
            srcz1 = a*fcz(j-1)+b*fcz(j)+ &
                    ((a**3-a)*ycz(j-1)+(b**3-b)*ycz(j))*(hh**2)/6.0d0
            srcz1 = max(srcz1, srczmin)
            srcz1 = min(srcz1, srczmax)
            tcouplecm1 = 1.0d-9*(a*taucouplecm(j-1)+b*taucouplecm(j)+ &
                         ((a**3-a)*ycouplecm(j-1)+(b**3-b)*ycouplecm(j))*(hh**2)/6.0d0)
            tcouple1 = 1.0d-9*(a*taucouple(j-1)+b*taucouple(j)+ &
                       ((a**3-a)*ycouple(j-1)+(b**3-b)*ycouple(j))*(hh**2)/6.0d0)
            ! start-of-substep wind rate
            fcc0 = min(0.5d0, swe0**2*scen0)
            fc0 = (fk2/(fk2**2+fcc0)**0.5d0)**excen
            djdt0 = -fc0*sstr0*swe0*min(w0,wc0)**(exw-1.0d0)
            ! start-of-substep transfer terms. fcm is built from -d(i_core)/dt
            ! (see setup), so following the same branch convention as fcz:
            ! when negative, use the omega of the zone it was differentiated
            ! from (core); when positive, use the other zone's omega (middle).
            if (srcm0.lt.0.0d0) then
                djtrancm0 = srcm0*h1*swcore0
            else
                djtrancm0 = srcm0*h1*swc0
            endif
            if (srcz0.lt.0.0d0) then
                djtranme0 = srcz0*h1*swe0
            else
                djtranme0 = srcz0*h1*swc0
            endif
            if (tcouplecm0.gt.0.0d0) then
                djcoupcm0 = sicore0v*sic0v/(sicore0v+sic0v)* &
                            (swc0-swcore0)/tcouplecm0
            endif
            if (tcouple0.gt.0.0d0) then
                djcoupme0 = sic0v*sie0v/(sic0v+sie0v)* &
                            (swc0-swe0)/tcouple0
            endif
            ! trial end-of-substep omegas from a forward-euler estimate
            ! using only start-of-substep rates, needed to evaluate
            ! end-of-substep saturation/centrifugal terms and the other
            ! zone's rate for the trapezoidal average below (mirrors
            ! int2zone's trial swe1/swc1 pattern, extended to 3 zones)
            swe1t = (sje0 + djtranme0 + h1*(djcoupme0-djdt0))/sie1v
            if (swe1t.lt.0.0d0) then
                iermsg = 5241
                return
            endif
            swc1t = (sjc0 - djtrancm0 - h1*djcoupcm0 &
                          - djtranme0 - h1*djcoupme0)/sic1v
            if (swc1t.lt.0.0d0) then
                iermsg = 5242
                return
            endif
            swcore1t = (sjcore0 + djtrancm0 + h1*djcoupcm0)/sicore1v
            if (swcore1t.lt.0.0d0) then
                iermsg = 5243
                return
            endif
            if (lross) then
                stau1 = a*staucz(j-1)+b*staucz(j)+ &
                        ((a**3-a)*ytau(j-1)+(b**3-b)*ytau(j))*(hh**2)/6.0d0
                stau1 = max(stau1,staumin)
                stau1 = min(stau1,staumax)
                if (iwind.eq.3) then
                    w1 = swe1t*stau1/soltau
                    wc1 = wcrit
                else if (iwind.eq.2) then
                    w1 = swe1t
                    wc1 = wcrit*soltau/stau1
                endif
            else
                w1 = swe1t
                wc1 = wcrit
            endif
            fcc1 = min(0.5d0, swe1t**2*scen1)
            fc1 = (fk2/(fk2**2+fcc1)**0.5d0)**excen
            djdt1 = -fc1*sstr1*swe1t*min(w1,wc1)**(exw-1.0d0)
            ! end-of-substep transfer terms, using the trial omegas and the
            ! same branch decisions (from srcm0/srcz0) as the start of the
            ! substep, exactly as int2zone does for its single boundary
            if (srcm0.lt.0.0d0) then
                djtrancm1 = srcm1*h1*swcore1t
            else
                djtrancm1 = srcm1*h1*swc1t
            endif
            if (srcz0.lt.0.0d0) then
                djtranme1 = srcz1*h1*swe1t
            else
                djtranme1 = srcz1*h1*swc1t
            endif
            if (tcouplecm1.gt.0.0d0) then
                djcoupcm1 = sicore1v*sic1v/(sicore1v+sic1v)* &
                            (swc1t-swcore1t)/tcouplecm1
            endif
            if (tcouple1.gt.0.0d0) then
                djcoupme1 = sic1v*sie1v/(sic1v+sie1v)* &
                            (swc1t-swe1t)/tcouple1
            endif
            ! final trapezoidal-averaged update: unlike int2zone, the
            ! coupling term is folded symmetrically into every zone it
            ! touches, so core+middle+envelope is conserved except for
            ! the wind-loss term (see header note)
            sje1 = sje0 + 0.5d0*(h1*(djdt0+djdt1+djcoupme0+djcoupme1)+ &
                   djtranme0+djtranme1)
            sjcore1 = sjcore0 + 0.5d0*(djtrancm0+djtrancm1+ &
                      h1*(djcoupcm0+djcoupcm1))
            sjc1 = sjc0 - 0.5d0*(djtrancm0+djtrancm1+djtranme0+djtranme1+ &
                   h1*(djcoupcm0+djcoupcm1+djcoupme0+djcoupme1))
            swe1 = sje1/sie1v
            swc1 = sjc1/sic1v
            swcore1 = sjcore1/sicore1v
            if (swe1.lt.0.0d0) then
                iermsg = 5261
                return
            endif
            if (swc1.lt.0.0d0) then
                iermsg = 5262
                return
            endif
            if (swcore1.lt.0.0d0) then
                iermsg = 5263
                return
            endif

            ! do not accumulate wind loss if the star is above the critical
            ! rossby number, exactly as int2zone
            if (2.0d0*cpi/swe1/staucz(j-1) .gt. rocrit*rosun) then
                sumdj = sumdj + 0.0d0
            else
                sumdj = sumdj + 0.5d0*h1*(djdt0+djdt1)
            endif
            sumtrancm = sumtrancm + 0.5d0*(djtrancm0+djtrancm1)+ &
                        0.5d0*h1*(djcoupcm0+djcoupcm1)
            sumtranme = sumtranme + 0.5d0*(djtranme0+djtranme1)+ &
                        0.5d0*h1*(djcoupme0+djcoupme1)

            ! copy end-of-substep values to start-of-next-substep values
            sicore0v = sicore1v
            sic0v = sic1v
            sie0v = sie1v
            sstr0 = sstr1
            scen0 = scen1
            stau0 = stau1
            w0 = w1
            srcm0 = srcm1
            srcz0 = srcz1
            tcouplecm0 = tcouplecm1
            tcouple0 = tcouple1
            swe0 = swe1
            swc0 = swc1
            swcore0 = swcore1
            sje0 = sje1
            sjc0 = sjc1
            sjcore0 = sjcore1
        end do
        xest = (dtt/float(nseq(ii)))**2
        nv = 3
        yest(1) = sumtrancm
        yest(2) = sumtranme
        yest(3) = sumdj
        call ratext(ii, xest, yest, yout, yerr, nv, nuse)
        errmax = 0.0d0
        do kk = 1, nv
            errmax = dmax1(errmax, dabs(yerr(kk)/yscal(kk)))
            err(kk) = dabs(yerr(kk)/yscal(kk))
        end do
        errmax = errmax/eps
        ytest(1,ii) = yout(1)
        ytest(2,ii) = yout(2)
        ytest(3,ii) = yout(3)
        if (errmax.lt.one) then
            sjcore0 = sjcoresav
            sjc0 = sjcsav
            sje0 = sjesav
            sjcore1 = sjcore0 + yout(1)
            sje1 = sje0 + yout(2) + yout(3)
            sjc1 = sjc0 - yout(1) - yout(2)
            lok = .true.
            return
        endif
    end do
    ! exhausted every sub-step count in nseq without the extrapolated
    ! error dropping below eps. int2zone has this same silent-failure gap
    ! (lok=.false. with iermsg left at its initial 0, indistinguishable
    ! from success to a caller that only checks iermsg); fixed here since
    ! it was directly observed to mask a real convergence failure right at
    ! core formation on a real track.
    iermsg = 5299
    lok = .false.
    return
end
