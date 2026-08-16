!===============================================================================
! subroutine rotwind
!
! Top-level driver for angular-momentum evolution of one low-mass star
! track (MHP 3/12; namelist-input and single-track conversion 2026).
! This is the entry point called from Python via f2py (module
! "pyrotwind") once per track.
!
! Physical model:
!   1) Rigid rotation within convective regions.
!   2) Star-disk coupling in the protostellar phase: rotation held fixed
!      at period pdisk (days) for a disk lifetime tlock (Myr).
!   3) Solid-body rotation enforced throughout if lsolid = .true.;
!      otherwise a core/envelope model is evolved, with core-envelope
!      coupling as in Denissenkov et al. (2010):
!        dj = (ic-ie)/(ic+ie) * (wc-we);  djc/dt = -dj/taucouple + 2/3 rcz^2 we dmcz
!        dje/dt = -djc/dt - jdot(wind)
!      lthreezone selects which core/envelope model: .false. (default)
!      evolves the ordinary two-zone model (drevol); .true. evolves the
!      three-zone core/middle/envelope model (threezoneevol), for tracks
!      that grow a helium core -- see its header for the governing
!      equations. threezoneevol also handles tracks that never grow a
!      core zone, so lthreezone is safe to leave on across a whole run
!      of mixed tracks.
!   4) Angular-momentum loss law selected by iwind:
!        iwind = 1  no loss
!        iwind = 2  Kawaler (1988) law
!        iwind = 3  generalized pmm-style law (pmm12)
!      (see setup for how the loss-law structural factors are built).
!
! Configuration is a two-step process, split so that reading and
! parsing the namelist file only happens once no matter how many
! tracks are processed:
!   1) Before calling rotwind at all, call constm.init_const(nml_file)
!      and params.read_wind_params(nml_file) once (both are public
!      module subroutines, directly callable from Python via f2py; see
!      rotwind_example.nml for the expected file format). This also
!      sets pdisk and tlock (disk-locked rotation period and disk
!      lifetime), which are applied to every track processed this run.
!   2) Then call rotwind once per track. rotwind itself only records
!      this call's array length (set_nmod) and runs the integration --
!      it no longer touches the namelist file.
!
! Inputs:
!   taucouple  this track's core-envelope (lthreezone=.false.) or
!              middle-envelope (lthreezone=.true.) coupling timescale
!              (yr), as a per-model array so it can vary along the track
!   taucouplecm  this track's core-middle coupling timescale (yr),
!              per-model array; only used when lthreezone=.true.
!   sicore     He core moment of inertia (cgs), per-model array; only
!              used when lthreezone=.true. -- pass all-zero on tracks
!              that never grow one, or if lthreezone=.false.
!   sage,sl,sr,sm,srcz,smcz,si,sie,sic,staucz,sp,sxc,shec,steffl
!              stellar-track arrays, cgs or solar units as noted below
!   nm         number of valid models on this track
!   nnmod      allocated length of the track arrays above
! Outputs:
!   swe,swc,swcore   envelope, middle, core angular velocity (rad/s);
!              swcore is only meaningful when lthreezone=.true.
!   sj,sje,sjc,sjcore   total, envelope, middle, core angular momentum
!              (cgs); sjcore is only meaningful when lthreezone=.true.
!   sprot      surface rotation period (days)
!   faccen     centrifugal-suppression factor actually applied (iwind=3 only)
!   iermsg     nonzero on failure (see solidevol/drevol/threezoneevol)
!===============================================================================
subroutine rotwind(taucouple, taucouplecm, &
                    sage, sl, sr, steffl, smcz, srcz, sxc, si, sie, sic, sicore, &
                    staucz, shec, sp, sm, swe, swc, swcore, sj, sje, sjc, sjcore, &
                    sprot, faccen, nm, nnmod, iermsg)
    use constm
    use params
    implicit none
    integer, intent(in) :: nnmod
    ! inputs from the track - cgs or solar units as noted
    ! sage = age (yr), sl = l/lsun sr = r/rsun sm = m/msun
    ! srcz = r of cz base (cgs) smcz = mass of surface cz, msun
    ! si, sie, sic = cgs moment of inertia total, envelope, middle
    ! sicore = cgs moment of inertia of the he core (three-zone only)
    ! staucz = convective overturn timescale (sec), sp = atm pressure tau=2/3
    ! sxc = central x shec= post-ms he core mass (msun), both used for isochrones
    real(8), intent(in), dimension(nnmod) :: sage, sl, sr, &
        sm, srcz, smcz, si, sie, sic, sicore, staucz, sp, sxc, shec, steffl, &
        taucouple, taucouplecm
    ! nm = number of valid models on this track
    integer, intent(in) :: nm
    ! output vectors - cgs units
    ! swe = omega(env), swc = omega(middle), swcore = omega(he core)
    ! sj = total j, sje = envelope j, sjc = middle j, sjcore = he core j
    ! sprot = surface period (days)
    real(8), intent(out), dimension(nnmod) :: swe, swc, swcore, sj, &
        sje, sjc, sjcore, sprot, faccen
    integer, intent(out) :: iermsg

    ! internal use
    real(8), dimension(nnmod) :: fstruct, fcen, fcz, fcm
    real(8) :: fk2, excen, exw
    real(8) :: bmax, dmdtmax, b_bsol, dmdt
    ! fac_cen is a dead scalar preserved from the original source (see note
    ! below, near its assignment).
    real(8) :: fac_cen
    integer :: k

    ! record this call's array length for the helper routines below (see
    ! params.f90). const/wind_params namelist state is assumed already
    ! initialized by the caller -- see the file header above.
    call set_nmod(nnmod)

    ! set up constants and products of structure variables used in the
    ! angular momentum loss law if angular momentum loss is included
    iermsg = 0
    if (iwind.ne.1) then
        call setup(nm, sage, sm, smcz, sp, sr, srcz, sie, sicore, &
                   sl, excen, exw, fstruct, fcz, fcm, fcen, fk2)
    endif
    if (lsolid) then
        call solidevol(nm, sage, si, staucz, excen, exw, &
                        fcen, fk2, fstruct, taucouple, &
                        swe, sj, sprot, iermsg)
        swc = swe
        swcore = swe
        sjc = 0.0d0
        sjcore = 0.0d0
    else if (lthreezone) then
        call threezoneevol(nm, sage, sicore, sic, sie, staucz, excen, exw, &
                            fcen, fk2, fstruct, fcm, fcz, taucouplecm, taucouple, &
                            swcore, swc, swe, sj, sjcore, sjc, sje, sprot, iermsg)
    else
        call drevol(nm, sage, si, sie, sic, staucz, excen, exw, &
                    srcz, smcz, fcen, fk2, fstruct, fcz, taucouple, &
                    swc, swe, sj, sjc, sje, sprot, iermsg)
        swcore = swc
        sjcore = 0.0d0
    endif

    ! diagnostic output: centrifugal-suppression factor actually applied at
    ! each track point (iwind=3 only), plus the corresponding b-field and
    ! mass-loss-rate scalings relative to solar (b_bsol, dmdt), computed
    ! for reference but not returned (kept from the original source).
    ! note: fac_cen here is a dead scalar, distinct from the faccen(:)
    ! output array below -- this matches the original source, where faccen
    ! is only ever assigned inside the iwind=3 branch of the loop below and
    ! is otherwise left uninitialized on return. preserved as-is.
    fac_cen = 0.0d0
    b_bsol = 0.0d0
    dmdt = 0.0d0
    if (iwind.eq.3) then
        if (wcrit.gt.0.0d0) then
            bmax = (wcrit/solw)**pmmb
            dmdtmax = (wcrit/solw)**pmma
        else
            bmax = 9.99d9
            dmdtmax = 9.99d9
        endif
    endif
    do k = 1, nm
        ! loss law ingredients
        if (iwind.eq.3 .and. swe(k).gt.0.0d0 .and. &
            staucz(k).gt.0.0d0) then
            faccen(k) = (fk2/(fk2**2+swe(k)**2*fcen(k)) &
                        **0.5d0)**excen
            b_bsol = min(bmax, (swe(k)*staucz(k)/solw/soltau)**pmmb)
            b_bsol = b_bsol*(sp(k)/solp)**0.5d0/sr(k)**pmmc
            dmdt = min(dmdtmax, (swe(k)*staucz(k)/ &
                          solw/soltau)**pmma)
            dmdt = dmdt*solmdot*csecyr/solm*sl(k)
        endif
    end do
    return
end
