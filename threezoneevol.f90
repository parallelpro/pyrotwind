!===============================================================================
! subroutine threezoneevol
!
! Angular momentum evolution for the three-zone (core/middle/envelope)
! model, the generalization of drevol to a growing helium core. Zones
! are chained core-middle-envelope: the middle zone is the same
! radiative-interior reservoir drevol calls "core", and the (new) core
! zone is a helium core that grows inside it later in the evolution.
! Handles all three zone counts that occur along a track:
!   1) fully convective (sic <= 0): single reservoir, int1zone, exactly
!      as in drevol/solidevol.
!   2) radiative interior present, no helium core yet (sic > 0, sicore
!      <= 0): two-zone model, reusing int2zone unmodified with the
!      middle zone playing its "core" argument. This is done by passing
!      int2zone a combined "total" inertia sic+sie in place of its usual
!      si argument; since int2zone recovers its own core term as
!      si-sie, that difference works out to exactly sic here -- with no
!      cancellation risk despite the subtraction, because sic is a
!      direct track input on this branch, not a residual.
!   3) helium core present (sicore > 0): three-zone model, int3zone (see
!      its header for the governing equations). Per explicit
!      confirmation from the user, a helium core is assumed to never
!      vanish once formed, so unlike the sic-vanishes case below there
!      is no 3-zone -> 2-zone reverse transition to handle.
!
! If the radiative interior itself vanishes after having formed (sic
! drops back to <= 0), the model merges back to a single reservoir via
! int1zone, exactly mirroring drevol's existing "radiative core
! vanished, merge" handling.
!
! Step-size limiting: in addition to the wind-loss-rate limit (as in
! drevol), each mass-transfer boundary gets its own limiter built
! directly from the same fcm/fcz structural factors used inside
! int3zone/int2zone (not from the legacy smcz/dmcz proxy in drevol,
! which is disconnected from fcz and inert whenever smcz is not
! populated with real data). Each boundary's proposed transfer is
! checked against 10% of the smaller of the two reservoirs it moves
! angular momentum between, so a fast-draining core (or middle zone)
! constrains the step just as much as a fast-draining envelope would.
! The one exception is the core/middle limiter during the first step
! after the core forms: sjcore(j-1) is then exactly zero by
! construction, which would make the constraint degenerate, so it is
! skipped for that single step; int3zone's own error codes are the
! safety net for that brief window.
!
! Inputs:
!   nm                          number of valid models on this track
!   sage, sicore, sic, sie      age (yr), core/middle/envelope moment of
!                               inertia (cgs)
!   staucz                      convective overturn timescale (s)
!   excen, exw                  loss-law exponents (see setup)
!   fcen, fk2                   centrifugal-term factor and suppression
!                               constant (see setup)
!   fstruct                     rotation-independent structural factor of
!                               the wind loss law
!   fcm, fcz                    core/middle and middle/envelope angular-
!                               momentum exchange structural factors
!                               (see setup)
!   taucouplecm, taucouple      core-middle and middle-envelope coupling
!                               timescales (yr), per-model arrays; <= 0
!                               at a given model disables that coupling
!                               there
! Outputs:
!   swcore, swc, swe            core, middle, envelope angular velocity
!                               (rad/s)
!   sj, sjcore, sjc, sje        total, core, middle, envelope angular
!                               momentum (cgs)
!   sprot                       surface (envelope) rotation period (days)
!   iermsg                      nonzero on failure (see int1zone/int2zone/
!                               int3zone; 6xxx codes below are raised
!                               directly by this routine)
!
! pdisk and tlock (initial rotation period, days, and disk lifetime,
! Myr) are not arguments: they come from the params module (set once
! per run by read_wind_params -- see params.f90).
!===============================================================================
subroutine threezoneevol(nm, sage, sicore, sic, sie, staucz, excen, exw, &
                          fcen, fk2, fstruct, fcm, fcz, taucouplecm, taucouple, &
                          swcore, swc, swe, sj, sjcore, sjc, sje, sprot, iermsg)
    use constm
    use params
    implicit none
    !f2py integer, intent(aux) :: nmod
    ! inputs
    integer, intent(in) :: nm
    real(8), intent(in), dimension(nmod) :: sage, sicore, sic, sie, &
        staucz, fcen, fstruct, fcm, fcz, taucouplecm, taucouple
    real(8), intent(in) :: excen, exw, fk2
    ! outputs
    real(8), intent(out), dimension(nmod) :: swcore, swc, swe, sj, &
        sjcore, sjc, sje, sprot
    integer, intent(out) :: iermsg
    ! local
    ! spline interpolation vectors
    real(8), dimension(nmod) :: x, y, yi, yicore, yic, yie, ystr, ytau, &
        ycm, ycz, ycen, ycouplecm, ycouple
    ! combined "total" inertia used only to reuse int2zone unmodified on
    ! the two-zone branch (see header), its spline coefficients, and the
    ! full core+middle+envelope total inertia needed by int1zone
    real(8), dimension(nmod) :: sictot, yictot, sitot
    logical :: lthree, ltwo, ldisk, lok, lrocrit
    real(8) :: wc1, wc0, w1, w0
    real(8) :: taudisk, hh, t0, t1, a, b, sii
    real(8) :: sje0, sjc0, sjcore0, sj0
    real(8) :: sje1, sjc1, sjcore1, sj1
    real(8) :: tmax, dt, dtt, fc0, fc1, djdt0, djdt1, djdt
    real(8) :: fcmrate, fczrate, swrel, dttrancm, dttranme, dttran2cm, dttran2me
    real(8) :: wcombined
    integer :: j, jj, k, kk

    iermsg = 0
    jj = nm
    if (jj.ge.2) then
        ! initialize angular momentum
        swcore(1) = 2.0d0*cpi/csecday/pdisk
        swc(1) = swcore(1)
        swe(1) = swcore(1)
        sje(1) = swe(1)*sie(1)
        sjc(1) = swc(1)*sic(1)
        sjcore(1) = swcore(1)*sicore(1)
        sj(1) = sje(1)+sjc(1)+sjcore(1)
        sprot(1) = pdisk
        ldisk = .true.
        taudisk = 1.0d-3*tlock
        ! local to this track's int1zone calls below (fully-convective
        ! phase only); see solidevol for what this flag means.
        ! threezoneevol does not itself act on it once set, matching
        ! drevol's behavior -- it is only threaded through here because
        ! int1zone's signature requires it.
        lrocrit = .false.
        ! set up spline interpolation in the structure variables, exactly
        ! as drevol, plus the extra core-zone quantities
        do j = 1, jj
            x(j) = sage(j)
            sitot(j) = sicore(j)+sic(j)+sie(j)
        end do
        call splinc(x, sitot, yi, jj)
        if (iwind.ne.1) then
            do j = 1, jj
                y(j) = fstruct(j)
            end do
            call splinc(x, y, ystr, jj)
            if (lross) then
                do j = 1, jj
                    y(j) = staucz(j)
                end do
                call splinc(x, y, ytau, jj)
            endif
            do j = 1, jj
                y(j) = fcen(j)
            end do
            call splinc(x, y, ycen, jj)
        endif
        do j = 1, jj
            y(j) = sicore(j)
        end do
        call splinc(x, y, yicore, jj)
        do j = 1, jj
            y(j) = sic(j)
        end do
        call splinc(x, y, yic, jj)
        do j = 1, jj
            y(j) = sie(j)
        end do
        call splinc(x, y, yie, jj)
        do j = 1, jj
            sictot(j) = sic(j)+sie(j)
        end do
        ! sic+sie is a separate quantity from sic on its own (used by
        ! int3zone via yic above), needed only to reuse int2zone
        ! unmodified on the two-zone branch (see header)
        call splinc(x, sictot, yictot, jj)
        do j = 1, jj
            y(j) = fcm(j)
        end do
        call splinc(x, y, ycm, jj)
        do j = 1, jj
            y(j) = fcz(j)
        end do
        call splinc(x, y, ycz, jj)
        do j = 1, jj
            y(j) = taucouplecm(j)
        end do
        call splinc(x, y, ycouplecm, jj)
        do j = 1, jj
            y(j) = taucouple(j)
        end do
        call splinc(x, y, ycouple, jj)

        do j = 2, jj
            if (ldisk) then
                if (sage(j).lt.taudisk) then
                    swcore(j) = 2.0d0*cpi/csecday/pdisk
                    swc(j) = swcore(j)
                    swe(j) = swcore(j)
                    sje(j) = swe(j)*sie(j)
                    sjc(j) = swc(j)*sic(j)
                    sjcore(j) = swcore(j)*sicore(j)
                    sj(j) = sje(j)+sjc(j)+sjcore(j)
                    sprot(j) = pdisk
                    cycle
                else
                    ldisk = .false.
                    hh = sage(j)-sage(j-1)
                    t0 = taudisk
                    a = (sage(j)-t0)/hh
                    b = (t0-sage(j-1))/hh
                    sii = a*(sicore(j-1)+sic(j-1)+sie(j-1))+ &
                          b*(sicore(j)+sic(j)+sie(j))+ &
                          ((a**3-a)*yi(j-1)+(b**3-b)*yi(j))*(hh**2)/6.0d0
                    sj(j) = swe(j-1)*sii
                    swe(j) = sj(j)/(sicore(j)+sic(j)+sie(j))
                    swc(j) = swe(j)
                    swcore(j) = swe(j)
                    sje(j) = swe(j)*sie(j)
                    sjc(j) = swc(j)*sic(j)
                    sjcore(j) = swcore(j)*sicore(j)
                    sprot(j) = 2.0d0*cpi/csecday/swe(j)
                    sje0 = sje(j)
                    sjc0 = sjc(j)
                    sjcore0 = sjcore(j)
                    sj0 = sje0+sjc0+sjcore0
                endif
            else
                sje0 = sje(j-1)
                sjc0 = sjc(j-1)
                sjcore0 = sjcore(j-1)
                sj0 = sj(j-1)
                t0 = sage(j-1)
            endif
            t1 = sage(j)
            tmax = t1-t0
            ! wind-loss-rate step limit, exactly as drevol (envelope-only,
            ! independent of zone count)
            if (iwind.ne.1) then
                fc0 = (fk2/(fk2**2 + swe(j-1)**2* &
                      fcen(j-1))**0.5d0)**excen
                swe(j) = sje(j-1)/sie(j)
                fc1 = (fk2/(fk2**2+swe(j)**2*fcen(j))**0.5d0)**excen
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
                djdt = max(djdt0,djdt1)
                if (djdt.gt.0.0d0) then
                    tmax = 0.1d0*sje(j-1)/djdt
                else
                    iermsg = 6171
                    return
                endif
            endif
            sj(j) = sj(j-1)

            ! middle/envelope transfer-rate step limit, built directly from
            ! fcz (see header); only meaningful once a middle zone exists
            if (sic(j-1).gt.0.0d0) then
                fczrate = max(abs(fcz(j-1)), abs(fcz(j)))
                if (fcz(j-1).lt.0.0d0) then
                    swrel = swe(j-1)
                else
                    swrel = swc(j-1)
                endif
                if (fczrate.gt.0.0d0 .and. swrel.gt.0.0d0) then
                    dttranme = 0.1d0*min(sje(j-1),sjc(j-1))/(fczrate*swrel)
                else
                    dttranme = t1-t0
                endif
            else
                dttranme = t1-t0
            endif
            ! core/middle transfer-rate step limit, only meaningful once a
            ! core zone exists with an established (nonzero) angular
            ! momentum -- see header note on the first-step exception
            if (sicore(j-1).gt.0.0d0 .and. sjcore(j-1).gt.0.0d0) then
                fcmrate = max(abs(fcm(j-1)), abs(fcm(j)))
                if (fcm(j-1).lt.0.0d0) then
                    swrel = swcore(j-1)
                else
                    swrel = swc(j-1)
                endif
                if (fcmrate.gt.0.0d0 .and. swrel.gt.0.0d0) then
                    dttrancm = 0.1d0*min(sjcore(j-1),sjc(j-1))/(fcmrate*swrel)
                else
                    dttrancm = t1-t0
                endif
            else
                dttrancm = t1-t0
            endif
            ! coupling-timescale step limits, one per boundary
            if (taucouple(j).gt.0.0d0) then
                dttran2me = 0.1d0*taucouple(j)*1.0d-9
            else
                dttran2me = t1-t0
            endif
            if (taucouplecm(j).gt.0.0d0) then
                dttran2cm = 0.1d0*taucouplecm(j)*1.0d-9
            else
                dttran2cm = t1-t0
            endif
            tmax = min(tmax,dttranme,dttrancm,dttran2me,dttran2cm)
            dt = (t1-t0)
            if (dt.gt.tmax) then
                kk = int(dt/tmax)+1
                dtt = dt/float(kk)
            else
                kk = 1
                dtt = dt
            endif

            ! decide the zone count in effect at this step
            if (sicore(j).gt.0.0d0) then
                lthree = .true.
                ltwo = .true.
            else if (sic(j).gt.0.0d0) then
                lthree = .false.
                ltwo = .true.
            else
                lthree = .false.
                ltwo = .false.
            endif

            if (lthree) then
                if (sicore(j-1).lt.1.0d0) then
                    ! newly developed helium core: asking int3zone's bs
                    ! integrator to grow a zone from a literal zero
                    ! inertia within a single step is a genuine numerical
                    ! singularity -- its trial/final trapezoidal formula
                    ! divides an angular-momentum estimate by whatever
                    ! (still tiny, rapidly changing) core inertia its own
                    ! sub-stepping has interpolated so far, and any
                    ! attempt to seed a nonzero starting sjcore mismatches
                    ! that tiny interpolated inertia and produces a wildly
                    ! inflated trial omega instead (tried and confirmed).
                    ! avoid the singularity altogether: evolve this one
                    ! transitional step as a combined (middle+core) vs
                    ! envelope two-zone problem with int2zone -- exactly
                    ! the sic(j-1).lt.1.0 branch below, just with the
                    ! about-to-exist core folded into the "core" side via
                    ! sitot (=sicore+sic+sie, already built above for
                    ! int1zone) -- then split the combined result between
                    ! middle and core assuming they end the step
                    ! corotating with each other, since a moment ago they
                    ! were one reservoir.
                    do k = 1, kk
                        call int2zone(sage, sitot, sie, staucz, fcz, fstruct, &
                                      yi, yie, ycz, ystr, ytau, sjc0, sje0, sjc1, &
                                      sje1, j, t0, dtt, excen, exw, fcen, fk2, &
                                      taucouple, ycouple, ycen, lok, iermsg)
                        if (.not.lok) return
                        t0 = t0+dtt
                        sjc0 = sjc1
                        sje0 = sje1
                    end do
                    wcombined = sjc1/(sic(j)+sicore(j))
                    sjc1 = wcombined*sic(j)
                    sjcore1 = wcombined*sicore(j)
                else
                    do k = 1, kk
                        call int3zone(sage, sicore, sic, sie, staucz, fcm, fcz, &
                                      fstruct, yicore, yic, yie, ycm, ycz, ystr, &
                                      ytau, sjcore0, sjc0, sje0, sjcore1, sjc1, &
                                      sje1, j, t0, dtt, excen, exw, fcen, fk2, &
                                      taucouplecm, taucouple, ycouplecm, ycouple, &
                                      ycen, lok, iermsg)
                        if (.not.lok) return
                        t0 = t0+dtt
                        sjcore0 = sjcore1
                        sjc0 = sjc1
                        sje0 = sje1
                    end do
                endif
                sjcore(j) = sjcore1
                sjc(j) = sjc1
                sje(j) = sje1
                sj(j) = sjcore1+sjc1+sje1
                swe(j) = sje(j)/sie(j)
                swc(j) = sjc(j)/sic(j)
                swcore(j) = sjcore(j)/sicore(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
            else if (ltwo) then
                if (sic(j-1).lt.1.0d0) then
                    ! newly developed radiative interior (no helium core
                    ! yet): initialize its angular momentum
                    sjc0 = 0.0d0
                endif
                do k = 1, kk
                    call int2zone(sage, sictot, sie, staucz, fcz, fstruct, &
                                  yictot, yie, ycz, ystr, ytau, sjc0, sje0, sjc1, &
                                  sje1, j, t0, dtt, excen, exw, fcen, fk2, &
                                  taucouple, ycouple, ycen, lok, iermsg)
                    if (.not.lok) return
                    t0 = t0+dtt
                    sjc0 = sjc1
                    sje0 = sje1
                end do
                sjc(j) = sjc1
                sje(j) = sje1
                sjcore(j) = 0.0d0
                sj(j) = sjc1+sje1
                swe(j) = sje(j)/sie(j)
                swc(j) = sjc(j)/sic(j)
                ! no core zone yet; report it corotating with the middle
                ! zone, mirroring how drevol reports swc corotating with
                ! swe during its own single-reservoir phase
                swcore(j) = swc(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
            else
                do k = 1, kk
                    call int1zone(sage, sitot, fstruct, staucz, yi, ystr, ytau, &
                                  sj0, sj1, j, t0, dtt, excen, exw, fcen, &
                                  fk2, ycen, lrocrit, iermsg)
                    t0 = t0+dtt
                    sj0 = sj1
                end do
                sj(j) = sj1
                swe(j) = sj(j)/sitot(j)
                swc(j) = swe(j)
                swcore(j) = swe(j)
                sje(j) = swe(j)*sie(j)
                sjc(j) = swc(j)*sic(j)
                sjcore(j) = swcore(j)*sicore(j)
                sprot(j) = 2.0d0*cpi/csecday/swe(j)
            endif
        end do
    endif
    return
end
