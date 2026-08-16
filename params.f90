!===============================================================================
! module params
!
! User-tunable angular-momentum-loss-law parameters, plus the track-array
! dimension used to size arrays throughout pyrotwind.
!
! The physics parameters (iwind, lsolid, lross, fk, pmma, pmmb, pmmc,
! pmmm, soljdot, solmdot, wcrit, rocrit) and the disk-locking parameters
! (pdisk, tlock) are set once per run by read_wind_params, which reads
! them from the wind_params namelist group of the run's input file (see
! rotwind_example.nml). Call this once, before any calls to rotwind.
! pdisk and tlock are applied to every track processed this run; the
! core-envelope coupling timescale (taucouple) is still a per-track
! argument to rotwind, since it plausibly does vary by track.
!
! nmod (the number of models allocated for this track's arrays) is NOT
! part of the namelist: it is set once per rotwind call by set_nmod,
! from the actual shape of the stellar-track arrays passed in from
! Python via f2py. pyrotwind processes one track per call, so there is
! no separate "number of tracks" dimension.
!
! lrocrit (whether this track has dropped below the critical rossby
! number, after which wind loss is switched off for the remainder of
! the track) used to live here as run-time state, but that made it a
! single flag shared across every track processed in a run instead of
! being reset per track: read_wind_params (which reset it) only runs
! once per run, not once per rotwind call, so it could leak from one
! track into the next. It is now local to solidevol/drevol, threaded
! explicitly into int1zone as an inout argument instead.
!===============================================================================
module params
    implicit none
    private
    public :: read_wind_params, set_nmod
    public :: fk, pmma, pmmb, pmmc, pmmm, soljdot, solmdot, wcrit, rocrit
    public :: pdisk, tlock
    public :: lsolid, lthreezone, lross
    public :: iwind, nmod

    ! ---- wind law selection and physics parameters (from namelist) ----

    ! iwind chooses the angular-momentum-loss law:
    !   1 = no loss (solid body or 2-zone model conserves j)
    !   2 = kawaler (1988) law, fixed exponent n = 1.5
    !   3 = generalized pmm-style law with user exponents pmma/b/c/m
    integer :: iwind

    ! lsolid = .true.  -> enforce solid-body rotation (see solidevol)
    !        = .false. -> evolve a core/envelope model (see drevol/
    !                     threezoneevol); which of those two is used is
    !                     controlled by lthreezone below
    logical :: lsolid

    ! lthreezone = .true.  -> evolve the three-zone core/middle/envelope
    !                         model (see threezoneevol), which also
    !                         subsumes the ordinary two-zone case for
    !                         tracks that never grow a core zone
    !            = .false. -> evolve the plain two-zone core/envelope
    !                         model (see drevol), as before
    ! only consulted when lsolid = .false.
    logical :: lthreezone

    ! lross = .true. -> scale the saturation threshold by the rossby
    ! number, using the local convective-overturn timescale
    logical :: lross

    ! wind law exponents/normalization, used when iwind = 3.
    real(8) :: fk, pmma, pmmb, pmmc, pmmm

    ! solar calibration values for the wind torque and mass-loss rate.
    real(8) :: soljdot, solmdot

    ! critical angular velocity (or rossby number) for wind saturation.
    real(8) :: wcrit

    ! critical rossby number below which wind loss is disabled, expressed
    ! in solar units: rocrit = (star_prot/star_tau) / (solar_prot/solar_tau).
    ! converted to an absolute threshold via constm's rosun (= solar_prot/
    ! solar_tau) at each comparison site -- see int1zone, int2zone, solidevol.
    real(8) :: rocrit

    ! initial (disk-locked) rotation period (days) and disk lifetime
    ! (myr), applied to every track processed this run.
    real(8) :: pdisk, tlock

    ! ---- track-array dimension (set per call, not from the namelist) ----
    ! nmod = allocated length of this track's arrays (number of models)
    integer :: nmod

contains

    !---------------------------------------------------------------------
    ! Read the wind_params namelist group from nml_file into the
    ! module's wind-law control variables. Call once, before any calls
    ! to rotwind.
    !---------------------------------------------------------------------
    subroutine read_wind_params(nml_file)
        character(*), intent(in) :: nml_file
        integer :: nml_unit

        namelist /wind_params/ lsolid, lthreezone, iwind, lross, fk, pmma, pmmb, &
                                pmmc, pmmm, soljdot, solmdot, wcrit, rocrit, &
                                pdisk, tlock

        open(newunit=nml_unit, file=nml_file, status='old', action='read')
        read(nml_unit, nml=wind_params)
        close(nml_unit)
    end subroutine read_wind_params

    !---------------------------------------------------------------------
    ! Record this track's array length for this call. This comes from
    ! the shape of the arrays passed in from Python/f2py, not from the
    ! namelist, and must be set before any of the other routines run.
    !---------------------------------------------------------------------
    subroutine set_nmod(nmod_in)
        integer, intent(in) :: nmod_in
        nmod = nmod_in
    end subroutine set_nmod

end module params
