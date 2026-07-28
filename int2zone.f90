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
! As in int1zone, the routine repeats the whole timestep at a fixed
! sequence of sub-step counts nseq(1), nseq(2), ... using the trapezoidal
! rule, and uses ratext to extrapolate the sequence of end-of-step
! (mass-transfer term, wind-loss term) sums to the h -> 0 limit,
! accepting the result once the estimated error drops below eps. Track
! quantities are obtained at intermediate times by cubic-spline
! interpolation between points j-1 and j, using coefficients (yi, yie,
! ycz, ystr, ytau, ycen) precomputed by splinc in the caller.
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
!   iermsg        nonzero on failure (negative core/envelope i or omega)
!
! Bug fix (2026): the original source's sjc1 update read a variable
! spelled djtran11, not djtran1 -- a typo that referenced a never-
! assigned (uninitialized) variable instead of the end-of-substep
! transfer term computed just above it. Confirmed present in the
! original source and traced to real symptoms (surface rotation period
! collapsing toward zero within a few models). Fixed here to read
! djtran1.
!===============================================================================
subroutine int2zone(sage, si, sie, staucz, fcz, fstruct, yi, yie, &
                     ycz, ystr, ytau, sjc0, sje0, sjc1, sje1, j, t0, dtt, &
                     excen, exw, fcen, fk2, taucouple, ycouple, ycen, lok, iermsg)
    use params
    use constm
    implicit none
    !f2py integer, intent(aux) :: nmod
    integer, parameter :: nmax=15
    real(8), parameter :: one=1.0d0
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

    ! local
    integer, dimension(11) :: nseq = (/2,4,6,8,12,16,24,32,48,64,96/)
    integer, parameter :: imax=11, nuse=7
    real(8), parameter :: eps=1.0d-5

    ! scalars
    real(8) :: t1, hh, h1, a, b
    real(8) :: si0, si1, sie0, sie1, sic0, sic1
    real(8) :: sstr0, sstr1, scen0, scen1, srcz0, srcz1
    real(8) :: stau0, stau1, staumin, staumax
    real(8) :: swe0, swe1, swc0, swc1
    real(8) :: w0, w1, wc0, wc1, wtest
    real(8) :: fcc0, fcc1, fc0, fc1
    real(8) :: djdt0, djdt1, djtran0, djtran1, djcoup0, djcoup1
    real(8) :: sjcsav, sjesav, tcouple0, tcouple1
    real(8) :: sumdj, sumtran, xest, errmax
    real(8), dimension(2,11) :: ytest
    real(8), dimension(nmax) :: yest, yout, yerr, err
    real(8), dimension(2) :: yscal
    integer :: ii, jj, kk, nv

    ! initialize iermsg
    iermsg = 0
    ! convergence flag
    lok = .false.
    ! scale error in transport term by envelope j
    yscal(1) = sje0
    ! tau(cz) can oscillate when close to a fully convective state
    ! to avoid numerical pathologies, restrict the range to within
    ! the endpoint values.
    staumin = min(staucz(j-1), staucz(j))
    staumax = max(staucz(j-1), staucz(j))
    ! scale error in loss rate by initial loss rate
    wtest = sje0/sie(j-1)
    if (lross) then
        w0 = wtest*staumin/soltau
    else
        w0 = wtest
    endif
    yscal(2) = fstruct(j-1)*wtest*min(w0,wcrit)**(exw-1.0d0)
    do kk = 1, 11
        ytest(1,kk) = 0.0d0
        ytest(2,kk) = 0.0d0
    end do
    djcoup0 = 0.0d0
    djcoup1 = 0.0d0
    ! store starting angular momenta
    sjcsav = sjc0
    sjesav = sje0
    do ii = 1, imax
        sjc0 = sjcsav
        sje0 = sjesav
        h1 = dtt/float(nseq(ii))
        t1 = t0
        hh = sage(j)-sage(j-1)
        ! interpolate between points j and j-1 for structure variables
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
        ! core-envelope coupling timescale (yr), interpolated like the other
        ! track arrays, then converted to gyr
        tcouple0 = 1.0d-9*(a*taucouple(j-1)+b*taucouple(j)+ &
                   ((a**3-a)*ycouple(j-1)+(b**3-b)*ycouple(j))*(hh**2)/6.0d0)
        ! surface and core angular velocities
        ! envelope i should never be zero, stop if true
        if (sie0.le.0.0d0) then
            iermsg = 4118
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
        ! core i can start as zero,initialize with envelope omega if true
        if (sic0.le.0.0d0) then
            swc0 = swe0
        else
            swc0 = sjc0/sic0
        endif
        ! counter to sum angular momentum loss overall
        sumdj = 0.0d0
        sumtran = 0.0d0
        do jj = 1, nseq(ii)
            t1 = t1 + h1
            ! interpolate between points j and j-1 for structure variables
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
            ! at intermediate times
            if (sic1.le.0.0d0) then
                ! switch to linear interpolation in moment of inertia
                si1 = si(j-1)+b*(si(j)-si(j-1))
                sie1 = sie(j-1)+b*(sie(j)-sie(j-1))
                sic1 = si1 - sie1
                if (sic1.le.0.0d0) then
                    iermsg = 4172
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
                lok = .false.
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
                return
            endif
            if (swc1.lt.0.0d0) then
                iermsg = 4271
                return
            endif

            ! do not accumulate wind loss if the star is above the critical rossby
            ! number (2*pi/omega_e/tau(cz) > rocrit*rosun), consistent with the
            ! rocrit cutoff applied to the one-zone model in int1zone. rocrit is
            ! expressed in solar units, so the absolute threshold is rocrit*rosun.
            if (2.0d0*cpi/swe1/staucz(j-1) .gt. rocrit*rosun) then
                sumdj = sumdj + 0.0d0
            else
                sumdj = sumdj + 0.5d0*h1*(djdt0+djdt1)
            endif

            sumtran = sumtran + 0.5d0*(djtran0+djtran1)+ &
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
        end do
        xest = (dtt/float(nseq(ii)))**2
        nv = 2
        yest(1) = sumtran
        yest(2) = sumdj
        call ratext(ii, xest, yest, yout, yerr, nv, nuse)
        errmax = 0.0d0
        do kk = 1, nv
            errmax = dmax1(errmax, dabs(yerr(kk)/yscal(kk)))
            err(kk) = dabs(yerr(kk)/yscal(kk))
        end do
        errmax = errmax/eps
        ytest(1,ii) = yout(1)
        ytest(2,ii) = yout(2)
        if (errmax.lt.one) then
            sjc0 = sjcsav
            sje0 = sjesav
            sje1 = sje0 + yout(1) + yout(2)
            ! compute transfer of angular momentum from core to envelope
            sjc1 = sjc0 - yout(1)
            lok = .true.
            return
        endif
    end do
    lok = .false.
    return
end
