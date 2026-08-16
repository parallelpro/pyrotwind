!===============================================================================
! subroutine int1zone
!
! Bulirsch-Stoer integrator (based on Numerical Recipes bsstep) that
! advances the total angular momentum sj0 -> sj1 of a single rigidly-
! rotating reservoir (solid body, or a track segment with no radiative
! core) across one outer timestep [t0, t0+dtt], under wind-braking torque
! dj/dt = -fc*fstruct*omega*min(omega,omega_crit)**(exw-1).
!
! Rather than adapting its step size, this routine repeats the whole
! timestep at a fixed sequence of sub-step counts nseq(1), nseq(2), ...
! (a trapezoidal-rule sub-integration each time), and uses ratext to
! extrapolate the sequence of end-of-step estimates to the h -> 0 limit.
! It accepts the extrapolated result once the estimated error (from
! ratext) drops below eps. Moment of inertia, structural torque factor,
! centrifugal factor, and (if lross) convective-overturn timescale are
! all obtained at intermediate times by cubic-spline interpolation
! between the tabulated track points j-1 and j (coefficients yi, ystr,
! ycen, ytau were precomputed by splinc in the caller).
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
    implicit none
    !f2py integer, intent(aux) :: nmod
    integer, parameter :: nmax = 15
    real(8), parameter :: one=1.0d0
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
    ! internal use
    ! angular momentum sj0 = start j, sj1 = end j; t0 = start time, dtt = timestep
    integer, dimension(11) :: nseq = (/2,4,6,8,12,16,24,32,48,64,96/)
    integer, parameter :: imax=11, nuse=7
    real(8), parameter :: eps=1.0d-5
    real(8) :: t1, hh, h1, a, b, si0, si1, stau0, stau1, staumin, staumax
    real(8) :: w0, w1, wc0, wc1, ro0, ro1
    real(8) :: fstr0, fstr1, fcen0, fcen1, fcc0, fcc1, fc0, fc1
    real(8) :: sw0, sw1, djdt0, djdt1, sjnew
    real(8) :: xest, errmax
    real(8), dimension(11) :: ytest
    real(8), dimension(nmax) :: yest, yout, yerr, yscal, err
    integer :: kk, ii, jj, nv

    ! initialize iermsg
    iermsg = 0
    yscal(1) = sj0
    do kk = 1, 11
        ytest(kk) = 0.0d0
    end do
    ! tau(cz) can oscillate when close to a fully convective state
    ! to avoid numerical pathologies, restrict the range to within
    ! the endpoint values.
    staumin = min(staucz(j-1), staucz(j))
    staumax = max(staucz(j-1), staucz(j))
    do ii = 1, imax
        h1 = dtt/float(nseq(ii))
        t1 = t0
        sj1 = sj0
        hh = sage(j)-sage(j-1)
        ! interpolate between points j and j-1 for structure variables
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
        do jj = 1, nseq(ii)
            t1 = t1 + h1
            ! interpolate between points j and j-1 for structure variables
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

            ! copy end of step values to start of next step values
            si0 = si1
            fstr0 = fstr1
            fcen0 = fcen1
            stau0 = stau1
            w0 = w1
            sw0 = sj1/si1

            if (sw0.lt.0.0d0) then
                iermsg = 3143
                return
            endif
        end do
        xest = (dtt/float(nseq(ii)))**2
        nv = 1
        yest(1) = sj1
        call ratext(ii, xest, yest, yout, yerr, nv, nuse)
        errmax = 0.0d0
        errmax = dmax1(errmax, dabs(yerr(1)/yscal(1)))
        err(1) = dabs(yerr(1)/yscal(1))
        errmax = errmax/eps
        ytest(ii) = yout(1)
        if (errmax.lt.one) then
            sj1 = yout(1)
            ! check whether the critical rossby number was crossed during this step;
            ! if so, interpolate back to the crossing point and flag lrocrit so the
            ! caller stops applying wind loss for the rest of this track. rocrit is
            ! expressed in solar units, so the absolute threshold is rocrit*rosun.
            if (2.0d0*cpi*si1/sj1/stau1 .gt. rocrit*rosun) then
                ro1 = 2.0d0*cpi*si1/sj1/stau1
                ro0 = 2.0d0*cpi*si0/sj0/stau0
                sjnew = (rocrit*rosun-ro0)*((sj1-sj0)/(ro1-ro0)) + sj0
                sj1 = sjnew
                lrocrit = .true.
            endif
            return
        endif
    end do
end
