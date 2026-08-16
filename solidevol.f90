!===============================================================================
! subroutine solidevol
!
! Evolves solid-body rotation (single reservoir, no core/envelope split)
! forward in time for one track, applying disk-locking followed by
! wind-braking torque. The loss law has the general form
!   dj/dt = fstruct * omega * min(omega,omega_crit)**(exw-1)
! (the saturated branch caps the omega-dependence at omega_crit once
! the star spins faster than that threshold; for a rossby-number
! scaling the saturation criterion instead compares
! omega*tau(cz)/tau(cz)_sun to omega_crit -- see int1zone). Each step
! is advanced with a Bulirsch-Stoer integrator (int1zone) using cubic-
! spline interpolation of the track structure across the step.
!
! For j = 1..nm:
!   1) j = 1: initialize at the disk-locked rotation rate 2*pi/pdisk.
!   2) While sage < taudisk = tlock/1000 (gyr): stay disk-locked,
!      conserving omega (not j) while i(t) evolves.
!   3) After decoupling: j is conserved as the star contracts (spin-up)
!      until wind loss is turned on, then integrated with int1zone.
!      Each outer track step [t0,t1] may be subdivided into kk
!      sub-steps sized so that the angular momentum lost per sub-step
!      does not exceed ~1% of the current j (see tmax below).
!   4) If iwind = 1 (no loss configured) or, once engaged, the star has
!      fallen below the critical rossby number (lrocrit), j is simply
!      conserved for the rest of the track.
!
! If the track has fewer than 2 valid models (nm < 2), nothing can be
! evolved and the routine returns immediately without touching swe/sj/
! sprot, exactly as in the original multi-track version (which simply
! skipped such a track).
!
! Inputs:
!   nm            number of valid models on this track
!   sage, si      age (yr) and total moment of inertia (cgs) per model
!   staucz        convective overturn timescale (s)
!   excen, exw    loss-law exponents (see setup)
!   fcen, fk2     centrifugal-term factor and suppression constant (see setup)
!   fstruct       rotation-independent structural factor of the loss law
!   taucouple     unused here (core-envelope coupling only applies to
!                 the two-zone model in drevol); per-model array, kept
!                 in the argument list for interface symmetry with drevol.
!
! lrocrit (whether this track has dropped below the critical rossby
! number, after which wind loss is disabled for the rest of the track)
! is local to this routine, initialized fresh at the top of every call
! (i.e. every track) and threaded into int1zone as an inout argument --
! not module state, so it cannot leak from one track's evolution into
! the next track processed in the same run.
!
! pdisk and tlock (initial disk-locked rotation period, days, and disk
! lifetime, Myr) are not arguments: they come from the params module
! (set once per run by read_wind_params -- see params.f90).
! Outputs:
!   swe           envelope==total angular velocity (rad/s)
!   sj            total angular momentum (cgs)
!   sprot         surface rotation period (days)
!   iermsg        0 on success from setup's perspective; this routine
!                 sets it to 1 if no error occurred, or to a nonzero
!                 setup-specific error code on failure (see int1zone)
!===============================================================================
subroutine solidevol(nm, sage, si, staucz, excen, exw, &
                      fcen, fk2, fstruct, taucouple, &
                      swe, sj, sprot, iermsg)
    use constm
    use params
    implicit none
    !f2py integer, intent(aux) :: nmod
    real(8), intent(in), dimension(nmod) :: sage, &
        si, staucz, fstruct, fcen, taucouple
    real(8), intent(in) :: fk2, excen, exw
    ! nm = number of valid models on this track
    integer, intent(in) :: nm
    ! output vectors - cgs units
    ! swe = omega(env), sj = total j, sprot = surface period (days)
    real(8), intent(out), dimension(nmod) :: swe, sj, sprot
    integer, intent(out) :: iermsg
    ! variables used within the routine
    ! spline interpolation vectors
    real(8), dimension(nmod) :: x, y, ycen, yi, ystr, ytau
    ! logicals
    logical :: ldisk, lrocrit
    real(8) :: taudisk, fx, t0, t1, hh, a, b, sii, sj0, sj1
    real(8) :: fc0, fc1, w0, w1, wc0, wc1, djdt0, djdt1, djdt
    real(8) :: tmax, dt, dtt
    integer :: j, jj, k, kk

    iermsg = 0
    jj = nm
    if (jj.ge.2) then
        ! initialize angular momentum
        swe(1) = 2.0d0*cpi/csecday/pdisk
        sj(1) = swe(1)*si(1)
        sprot(1) = pdisk
        ! initialize disk and disk lifetime in gyr
        ldisk = .true.
        taudisk = 1.0d-3*tlock
        ! rossby-number cutoff starts disengaged for this track; the
        ! int1zone call below enables it once (if) the track first drops
        ! below rocrit.
        lrocrit = .false.
        ! initialize torque
        ! the loss law has the general form
        ! dj/dt = fstruct*omega^exw, omega < omega(crit);
        ! dj/dt = fstruct*omega*omega(crit)^(exw-1) otherwise
        ! (n.b. centrifugal term for pmm law, see below)
        ! for a rossby scaling the saturation criterion is
        ! omega * tau(cz)/ tau(cz) sun < omega(crit)(sun)
        ! the code evolves across a single step using a
        ! bulirsch-stoer integrator and spline interpolation
        ! across a step.
        ! begin by setting up spline interpolation in the structure
        ! variables for the wind loss, itot, and overturn timescale.
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
            ! spline factors for loss term independent of omega
            call splinc(x, y, ystr, jj)
            do j = 1, jj
                y(j) = fcen(j)
            end do
            ! spline factor for centrifugal term in pmm wind law
            call splinc(x, y, ycen, jj)
        endif
        if (lross) then
            do j = 1, jj
                y(j) = staucz(j)
            end do
            ! spline factors for overturn timescale
            call splinc(x, y, ytau, jj)
        endif
        ! general loop for angular momentum evolution calculations
        do j = 2, jj
            ! models start with star-disk coupling.  evolve at
            ! fixed omega until disk age reached.  if no disk coupling
            ! is chosen the model will evolve from the starting period.
            if (ldisk) then
                if (sage(j).lt.taudisk) then
                    ! disk locked to initial rotation rate
                    swe(j) = 2.0d0*cpi/csecday/pdisk
                    sj(j) = swe(j)*si(j)
                    sprot(j) = pdisk
                    cycle
                else
                    ! disk decouples
                    ldisk = .false.
                    fx = (taudisk-sage(j-1))/(sage(j)-sage(j-1))
                    t0 = taudisk
                    ! spline interpolate to total i at decoupling epoch
                    ! interpolate between points j and j-1 for structure variables
                    hh = sage(j)-sage(j-1)
                    a = (sage(j)-t0)/hh
                    b = (t0-sage(j-1))/hh
                    sii = a*si(j-1)+b*si(j)+ &
                          ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
                    sj(j) = swe(j-1)*sii
                    swe(j) = sj(j)/si(j)
                    sprot(j) = 2.0d0*cpi/csecday/swe(j)
                    sj0 = sj(j)
                endif
            else
                sj0 = sj(j-1)
                t0 = sage(j-1)
            endif
            t1 = sage(j)

            ! no loss solid body case
            if (iwind.eq.1) then
                sj(j) = sj(j-1)
                swe(j) = sj(j)/si(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
                cycle
            endif

            ! no loss if ro > rocrit*rosun (rocrit is expressed in solar units)
            if (lrocrit) then
                ! check if you drop below rocrit this timestep
                sj(j) = sj(j-1)
                swe(j) = sj(j)/si(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
                if (sprot(j)*csecday/staucz(j) .lt. rocrit*rosun) then
                    lrocrit = .false.
                else
                    cycle
                endif
            endif

            ! evolve from time t0 to time t1 including wind torque.
            ! star conserves angular momentum and spins up
            ! evaluate maximum timestep
            ! limit timestep based on angular momentum loss
            ! start of timestep rate
            ! fc = centrifugal term, pmm
            fc0 = (fk2/(fk2**2 + swe(j-1)**2*fcen(j-1))**0.5d0)**excen
            swe(j) = swe(j-1)*si(j-1)/si(j)
            ! end of timestep rate including structural evolution
            fc1 = 1.0d0
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

            ! saturation accounted for - start of timestep
            djdt0 = fc0*fstruct(j-1)*swe(j-1)* &
                    min(w0,wc0)**(exw-1.0d0)
            ! saturation accounted for - end of timestep
            djdt1 = fc1*fstruct(j)*swe(j)* &
                    min(w1,wc1)**(exw-1.0d0)
            djdt = max(djdt0,djdt1)
            sj(j) = sj(j-1)
            ! limit timestep to a maximum fraction of the angular momentum
            ! removed from the star, assuming the highest rate
            if (djdt.gt.0.0d0) then
                tmax = 0.01d0*sj(j)/djdt
            else
                iermsg = 9160
                return
            endif
            dt = (t1-t0)
            if (dt.gt.tmax) then
                kk = int(dt/tmax)+1
                dtt = dt/float(kk)
            else
                kk = 1
                dtt = dt
            endif
            ! loop for torque calculation
            ! sj0 set above when ldisk flipped, or from the prior model point
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
            sprot(j) = 2.0d0*cpi/csecday/swe(j)
        end do
    endif

    if (iermsg.eq.0) then
        iermsg = 1
    endif
    return
end
