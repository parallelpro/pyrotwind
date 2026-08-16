!===============================================================================
! subroutine drevol
!
! Angular momentum evolution under the assumption of the Sills et al.
! (2000) j-dot wind-braking law, for the two-zone (core/envelope) model.
! This is the drevol/two-zone counterpart of solidevol: it evolves one
! track through disk-locking and then wind-braking, but tracks the
! envelope and radiative-core angular momenta separately, coupling them
! via mass transfer and (optionally) core-envelope torque coupling. See
! int2zone's header for the governing equations during a single step.
!
! For j = 1..nm:
!   1) j = 1: initialize both zones at the disk-locked rotation rate
!      2*pi/pdisk (rigid rotation while disk-locked).
!   2) While sage < taudisk = tlock/1000 (gyr): stay disk-locked.
!   3) After decoupling, and while the star has no radiative core
!      (sic <= 0): treat it as a single rigid reservoir and integrate j
!      with int1zone, exactly as in solidevol.
!   4) Once a radiative core develops (sic > 0): integrate the envelope
!      and core angular momenta separately with int2zone, subdividing
!      the outer track step into kk sub-steps sized to limit both the
!      wind-loss rate and the core/envelope mass-transfer rate per
!      sub-step (see tmax/dttran/dttran2 below).
!
! If the track has fewer than 2 valid models (nm < 2), nothing can be
! evolved and the routine returns immediately without touching swc/swe/
! sj/sjc/sje/sprot, exactly as in the original multi-track version
! (which simply skipped such a track).
!
! Inputs:
!   nm                     number of valid models on this track
!   sage, si, sie, sic     age (yr), total/envelope/core moment of
!                          inertia (cgs)
!   staucz                 convective overturn timescale (s)
!   excen, exw             loss-law exponents (see setup)
!   srcz, smcz             convection-zone base radius (cgs) and mass
!                          (msun); used only to bound the mass-transfer
!                          sub-stepping limit below
!   fcen, fk2              centrifugal-term factor and suppression
!                          constant (see setup)
!   fstruct                rotation-independent structural factor of
!                          the wind loss law
!   fcz                    core/envelope angular-momentum exchange
!                          structural factor (see setup)
!   taucouple              core-envelope coupling timescale (yr), as a
!                          per-model array so it can vary along the
!                          track; <= 0 at a given model disables
!                          coupling there
! Outputs:
!   swc, swe               core, envelope angular velocity (rad/s)
!   sj, sjc, sje           total, core, envelope angular momentum (cgs)
!   sprot                  surface (envelope) rotation period (days)
!   iermsg                 nonzero on failure (see int1zone/int2zone)
!
! pdisk and tlock (initial rotation period, days, and disk lifetime,
! Myr) are not arguments: they come from the params module (set once
! per run by read_wind_params -- see params.f90).
!===============================================================================
subroutine drevol(nm, sage, si, sie, sic, staucz, excen, exw, &
                   srcz, smcz, fcen, fk2, fstruct, fcz, &
                   taucouple, swc, swe, sj, sjc, sje, sprot, iermsg)
    use constm
    use params
    implicit none
    !f2py integer, intent(aux) :: nmod
    ! inputs
    integer, intent(in) :: nm
    real(8), intent(in), dimension(nmod) :: sage, si, sie, sic, &
        staucz, srcz, smcz, fcen, fstruct, fcz, taucouple
    real(8), intent(in) :: excen, exw, fk2
    ! outputs
    ! output vectors - cgs units
    ! swe = omega(env), swc = omega(core), sj= total j sje = envelope j
    ! sjc = core j sprot = surface period (days)
    real(8), intent(out), dimension(nmod) :: swc, swe, sj, sjc, sje, sprot
    integer, intent(out) :: iermsg
    ! local
    ! spline interpolation vectors
    real(8), dimension(nmod) :: x, y, yi, yie, ystr, ytau, ycz, ycen, ycouple
    logical :: lcore, ldisk, lok, lrocrit
    real(8) :: wc1, wc0, w1, w0
    real(8) :: taudisk, hh, t0, t1, a, b, sii, sje0, sjc0, sj0, sje1, sjc1, sj1
    real(8) :: tmax, dt, dtt, fc0, fc1, djdt0, djdt1, djdt
    real(8) :: dmcz, djtran, dttran, dttran2
    integer :: j, jj, k, kk

    iermsg = 0
    jj = nm
    if (jj.ge.2) then
        ! initialize angular momentum
        swc(1) = 2.0d0*cpi/csecday/pdisk
        swe(1) = swc(1)
        sj(1) = swe(1)*si(1)
        sje(1) = swe(1)*sie(1)
        sjc(1) = swc(1)*sic(1)
        sprot(1) = pdisk
        ! initialize disk and disk lifetime in gyr
        ldisk = .true.
        taudisk = 1.0d-3*tlock
        ! local to this track's int1zone calls below (pre-core phase only);
        ! see solidevol for what this flag means. drevol does not itself
        ! act on it once set, matching prior behavior -- it is only
        ! threaded through here because int1zone's signature requires it.
        lrocrit = .false.
        ! set up spline interpolation in the structure
        ! variables for wind loss, itot and overturn timescale (general case)
        ! also ienv, 2/3 rcz^2 * dmcz/dt (for the decoupled case)
        do j = 1, jj
            x(j) = sage(j)
            y(j) = si(j)
        end do
        ! spline factors for moment of inertia
        call splinc(x, y, yi, jj)
        if (iwind.ne.1) then
            do j = 1, jj
                y(j) = fstruct(j)
            end do
            ! spline factors for wind loss terms
            call splinc(x, y, ystr, jj)
            if (lross) then
                do j = 1, jj
                    y(j) = staucz(j)
                end do
                ! spline factors for overturn timescale
                call splinc(x, y, ytau, jj)
            endif
            ! centrifugal term from pmm wind law
            do j = 1, jj
                y(j) = fcen(j)
            end do
            ! spline factors for wind loss terms
            call splinc(x, y, ycen, jj)
        endif
        ! set up spline interpolation in envelope i, 2/3rcz^2*dmcz/dt
        ! for the decoupled case.
        do j = 1, jj
            y(j) = sie(j)
        end do
        ! spline factors for envelope moment of inertia
        call splinc(x, y, yie, jj)
        do j = 1, jj
            y(j) = fcz(j)
        end do
        ! spline factors for change in cz depth term
        call splinc(x, y, ycz, jj)
        do j = 1, jj
            y(j) = taucouple(j)
        end do
        ! spline factors for core-envelope coupling timescale
        call splinc(x, y, ycouple, jj)
        do j = 2, jj
            if (ldisk) then
                if (sage(j).lt.taudisk) then
                    ! disk locked to initial rotation rate
                    swc(j) = 2.0d0*cpi/csecday/pdisk
                    swe(j) = swc(j)
                    sj(j) = swe(j)*si(j)
                    sje(j) = swe(j)*sie(j)
                    sjc(j) = swc(j)*sic(j)
                    sprot(j) = pdisk
                    cycle
                else
                    ! disk decouples
                    ldisk = .false.
                    ! spline interpolate to total i at decoupling epoch
                    ! interpolate between points j and j-1 for structure variables
                    hh = sage(j)-sage(j-1)
                    t0 = taudisk
                    a = (sage(j)-t0)/hh
                    b = (t0-sage(j-1))/hh
                    sii = a*si(j-1)+b*si(j)+ &
                          ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
                    sj(j) = swe(j-1)*sii
                    swe(j) = sj(j)/si(j)
                    swc(j) = swe(j)
                    sje(j) = swe(j)*sie(j)
                    sjc(j) = swc(j)*sic(j)
                    sprot(j) = 2.0d0*cpi/csecday/swe(j)
                    sje0 = sje(j)
                    sjc0 = sjc(j)
                    sj0 = sje0 + sjc0
                endif
            else
                sje0 = sje(j-1)
                sjc0 = sjc(j-1)
                sj0 = sj(j-1)
                t0 = sage(j-1)
            endif
            t1 = sage(j)
            tmax = t1-t0
            ! evolve from time t0 to time t1 including wind torque if applicable
            if (iwind.ne.1) then
                ! evaluate maximum timestep assuming only the cz spins down
                ! start of timestep rate
                fc0 = (fk2/(fk2**2 + swe(j-1)**2* &
                      fcen(j-1))**0.5d0)**excen
                ! end of timestep rate
                swe(j) = sje(j-1)/sie(j)
                fc1 = (fk2/(fk2**2+swe(j)**2*fcen(j))**0.5d0)**excen
                ! account for saturation
                if (lross) then
                    if (iwind.eq.3) then
                        w0 = swe(j-1)*staucz(j-1)/soltau
                        wc0 = wcrit
                        w1 = swe(j)*staucz(j)/soltau
                        wc1 = wcrit
                    else if (iwind.eq.2) then
                        w0 = swe(j-1)
                        wc0 = wcrit*soltau/staucz(j)
                        w1 = swe(j)
                        wc1 = wcrit*soltau/staucz(j)
                    endif
                else
                    w0 = swe(j-1)
                    wc0 = wcrit
                    w1 = swe(j)
                    wc1 = wcrit
                endif
                djdt0 = fc0*fstruct(j-1)*swe(j-1)* &
                        min(w0,wc0)**(exw-1.0d0)
                djdt1 = fc1*fstruct(j)*swe(j)* &
                        min(w1,wc1)**(exw-1.0d0)
                ! choose the larger of the start or end value to constrain dt
                djdt = max(djdt0,djdt1)
                if (djdt.gt.0.0d0) then
                    tmax = 0.1d0*sje(j-1)/djdt
                else
                    iermsg = 2171
                    return
                endif
            endif
            sj(j) = sj(j-1)
            ! limit timestep based on rate of change of envelope j
            dmcz = smcz(j)-smcz(j-1)
            if (abs(dmcz).gt.1.0d-9) then
                djtran = (srcz(j-1)**2+srcz(j-1)**2)/3.0d0* &
                         dmcz*solm/(sage(j)-sage(j-1))
                ! do not apply condition if rate of change of cz is too small.
                if (djtran.lt.0.0d0) then
                    dttran = abs(0.1d0*sje(j-1)/(djtran*swe(j-1)))
                else if (djtran.gt.0.0d0) then
                    dttran = 0.1d0*sje(j-1)/(djtran*swc(j-1))
                endif
            else
                dttran = t1-t0
            endif
            ! limit timestep based on angular momentum exchange from
            ! envelope to core (or vice versa)
            if (taucouple(j).gt.0.0d0) then
                dttran2 = 0.1d0*taucouple(j)*1.0d-9
            else
                dttran2 = dttran
            endif
            dt = (t1-t0)
            tmax = min(tmax,dttran,dttran2)
            if (dt.gt.tmax) then
                kk = int(dt/tmax)+1
                dtt = dt/float(kk)
            else
                kk = 1
                dtt = dt
            endif
            ! loop for torque calculation
            ! check and see whether there is a radiative core
            if (sic(j).gt.0.0d0) then
                if (sic(j-1).lt.1.0d0) then
                    sje0 = sje(j-1)
                    sjc0 = 0.0d0
                    ! newly developed radiative core.  initialize core angular momentum terms.
                else
                    sje0 = sje(j-1)
                    sjc0 = sjc(j-1)
                endif
                lcore = .true.
            else if (sic(j-1).gt.0.0d0) then
                ! radiative core vanished.  merge angular momentum.
                lcore = .false.
            else
                lcore = .false.
            endif

            if (lcore) then
                ! use a b-s intergrator with spline interpolation to solve for angular momentum loss
                ! across the timestep.  numerical convergence properties are currently hardwired in
                ! bsstep, to be replaced with user-specified parameters.
                ! this sr solves for a 2
                ! zone model instead of a solid body model.
                do k = 1, kk
                    call int2zone(sage, si, sie, staucz, fcz, fstruct, yi, &
                                  yie, ycz, ystr, ytau, sjc0, sje0, sjc1, sje1, j, t0, &
                                  dtt, excen, exw, fcen, fk2, taucouple, ycouple, ycen, lok, iermsg)
                    if (.not.lok) return
                    t0 = t0+dtt
                    sje0 = sje1
                    sjc0 = sjc1
                end do
                sje(j) = sje1
                sjc(j) = sjc1
                sj(j) = sje1+sjc1
                swe(j) = sje(j)/sie(j)
                swc(j) = sjc(j)/sic(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
            else
                do k = 1, kk
                    ! use a b-s intergrator with spline interpolation to solve for angular momentum loss
                    ! across the timestep.  numerical convergence properties are currently hardwired in
                    ! bsstep, to be replaced with user-specified parameters.
                    call int1zone(sage, si, fstruct, staucz, yi, ystr, ytau, sj0, sj1, &
                                  j, t0, dtt, excen, exw, fcen, fk2, ycen, lrocrit, iermsg)
                    t0 = t0+dtt
                    sj0 = sj1
                end do
                sj(j) = sj1
                swe(j) = sj(j)/si(j)
                swc(j) = swe(j)
                sje(j) = swe(j)*sie(j)
                sjc(j) = swc(j)*sic(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
            endif
        end do
    endif
    return
end
