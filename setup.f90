!===============================================================================
! subroutine setup
!
! Angular-momentum evolution under the assumption of the Sills et al.
! (2000) j-dot wind-braking law (as generalized by the pmm-style
! parameterization below). Called once per track (when iwind /= 1) to
! precompute, for every point on the track, the rotation-independent
! factors used by the loss law in solidevol/drevol/int1zone/int2zone:
!
!   fstruct(j) -- structural prefactor of the wind torque law, so that
!                 dj/dt = fstruct * f(omega), with the omega-dependence
!                 (including saturation) applied elsewhere.
!   fcen(j)    -- centrifugal-term factor 0.5*(rsun^3/msun/g)*r^3/m,
!                 used to suppress the torque at fast rotation
!                 (Matt et al. 2012).
!   fcz(j)     -- for two-zone (core/envelope) and three-zone
!                 (core/middle/envelope) models: the rate of
!                 angular-momentum exchange between the envelope and the
!                 zone below it (middle zone, or core if there is no
!                 middle zone) as convection-zone mass is gained or
!                 lost. Formally this should be d(j transferred)/dt =
!                 2/3 r_cz^2 * dm_cz/dt (Denissenkov et al. 2010), but
!                 r_cz and m_cz are not available from the input tracks
!                 here. Instead we use the equivalent full-derivative
!                 form dj/dt = di(source)/dt * omega(source), which for
!                 the envelope moment of inertia i_e works out to
!                 dj/dt = -di_e/dt * omega(source) (the sign flip
!                 keeps the convention consistent with dm_cz/dt: fcz is
!                 evaluated by spline-differentiating -i_e(age)).
!   fcm(j)     -- three-zone models only: the same kind of exchange
!                 term as fcz, but for the boundary between the core
!                 and the middle zone, built the same way from
!                 -i_core(age) (see threezoneevol/int3zone). Zero on
!                 tracks with no core zone (sicore <= 0 throughout).
!
! fk2 is the fixed centrifugal-suppression constant from Matt et al.
! (2012), set here rather than read from the namelist because it is a
! physical constant of that torque-suppression model, not a tunable
! run parameter.
!
! excen, exw (and the local exponents exr, exm, expr, exl used to build
! fstruct) are the power-law exponents of the loss law, selected by
! iwind:
!   iwind = 2 : Kawaler (1988) law, fixed n = 1.5
!   iwind = 3 : generalized pmm-style law in terms of user exponents
!               pmma, pmmb, pmmc, pmmm (see Matt et al. 2012 / pmm12)
!===============================================================================
subroutine setup(nm, sage, sm, smcz, sp, sr, srcz, sie, sicore, &
                  sl, excen, exw, fstruct, fcz, fcm, fcen, fk2)
    use params
    use constm
    implicit none
    !f2py integer, intent(aux) :: nmod
    ! inputs from the track - cgs or solar units as noted
    ! sage = age (yr), sl = l/lsun sr = r/rsun sm = m/msun
    ! srcz = r of cz base (cgs) smcz = mass of surface cz, msun
    ! sie = envelope moment of inertia (cgs)
    ! sicore = core moment of inertia (cgs), three-zone models only;
    !          pass all-zero if the track never has a core zone
    ! sp = atm pressure tau=2/3
    integer, intent(in) :: nm
    real(8), intent(in), dimension(nmod) :: sage, sm, sr, &
        smcz, srcz, sl, sp, sie, sicore
    ! outputs - cgs units. loss law has the form jdot=fstruct*f(omega); for
    ! 2-zone models the material that moves between core and envelope
    ! carries angular momentum accounted for via fcz (see header above).
    ! fcm is the equivalent term for the core/middle boundary in
    ! three-zone models.
    real(8), intent(out), dimension(nmod) :: fstruct, fcz, fcm, fcen
    real(8), intent(out) :: exw, excen, fk2
    ! local: exponents for factors in the loss law
    real(8) :: exr, exm, expr, exl, fwind
    ! local: spline interpolation vectors for d(-i_envelope)/dt and,
    ! for three-zone models, d(-i_core)/dt
    real(8), dimension(nmod) :: ydiedt, xmid, diedt
    real(8), dimension(nmod) :: ydicdt, dicdt
    real(8) :: h, fm
    integer :: j, jjj

    ! constant term for centrifugal torque suppression - see matt et al. 2012
    fk2 = 0.056d0
    do j = 1, nm
        fcen(j) = 0.5d0*(solr**3/solm/cg)* &
                  sr(j)**3/sm(j)
    end do
    if (iwind.eq.2) then
        ! loss law constant for kawaler law (n=1.5)
        fwind = fk*2.036d33*1.452d9**1.5d0*1.0d9*csecyr
        ! set up indices for loss law in kawaler law
        ! exr = exponent for r
        ! exm = exponent for m
        ! exw = exponent for omega
        excen = 0.0d0
        exr = 0.5d0
        exm = -0.5d0
        expr = 0.0d0
        exl = 0.0d0
        exw = 3.0d0
    else if (iwind.eq.3) then
        ! set up indices for loss law in terms of pmm a, b, c,
        ! excen = exponent for centrifugal term, 1/(k2^2+w^2 r^3/gm)^pmmm/2
        excen = pmmm
        exr = 2.0d0+5.0d0*pmmm-4.0d0*pmmm*pmmc
        exm = -pmmm
        ! exl = exponent for l;assume lx~lbol
        ! exp = exponent for p;assume equipartition for spot fields p~b^2
        expr = 2.0d0*pmmm
        exl = 1.0d0 - 2.0d0*pmmm
        exw = 1.0d0+pmma*(1.0d0 - 2.0d0*pmmm)+4.0d0*pmmm*pmmb
        ! loss law constant for pmm law
        fwind = fk*soljdot*1.0d9*csecyr/solw**exw
    endif
    ! product of structural factors used for the angular momentum loss law.
    do j = 1, nm
        fstruct(j) = fwind*(sr(j)**exr)*(sm(j)**exm)* &
                     ((sp(j)/solp)**expr)*(sl(j)**exl)
    end do

    ! core-envelope decoupling is treated in a 2-zone model, where the
    ! transfer of angular momentum between the two is treated as
    ! delta j = omega(source)*2/3*r_cz**2*dm_cz/dt * delta (time).
    ! since r_cz and m_cz are not available here, fcz is instead computed
    ! from the full derivative dj/dt = -di_envelope/dt * omega(source) (see
    ! header above), by spline-differentiating -sie(age).
    if (.not.lsolid) then
        jjj = nm
        do j = 1, jjj
            xmid(j) = sage(j)
            diedt(j) = -sie(j)
        end do
        ! spline factors for -envelope moment of inertia
        call splinc(xmid, diedt, ydiedt, jjj)
        do j = 1, jjj-1
            ! first derivative at a point is
            ! dy/dx(i) = (y(i+1)-y(i))/(x(i+1)-x(i))+
            ! 1/3(x(i+1)-x(i))*y2a(i)-1/6(x(i+1)-x(i))*y2a(i+1)
            ! where y2a is the second derivative in splinc (c.f. numerical recipes 3.3)
            h = xmid(j+1)-xmid(j)
            fm = (diedt(j+1)-diedt(j))/h+ &
                 h*ydiedt(j)/3.0d0-h*ydiedt(j+1)/6.0d0
            fcz(j) = fm
        end do
        ! flip indices for final point
        h = xmid(jjj)-xmid(jjj-1)
        fm = (diedt(jjj)-diedt(jjj-1))/h- &
             h*ydiedt(jjj-1)/6.0d0+h*ydiedt(jjj)/3.0d0
        fcz(jjj) = fm

        ! same construction for the core/middle boundary in three-zone
        ! models, built from -i_core(age) instead of -i_envelope(age).
        ! on a track that never has a core zone (sicore all zero) this
        ! comes out identically zero, so it is always safe to compute.
        do j = 1, jjj
            dicdt(j) = -sicore(j)
        end do
        call splinc(xmid, dicdt, ydicdt, jjj)
        do j = 1, jjj-1
            h = xmid(j+1)-xmid(j)
            fm = (dicdt(j+1)-dicdt(j))/h+ &
                 h*ydicdt(j)/3.0d0-h*ydicdt(j+1)/6.0d0
            fcm(j) = fm
        end do
        h = xmid(jjj)-xmid(jjj-1)
        fm = (dicdt(jjj)-dicdt(jjj-1))/h- &
             h*ydicdt(jjj-1)/6.0d0+h*ydicdt(jjj)/3.0d0
        fcm(jjj) = fm
    endif
    return
end
