!===============================================================================
! module constm
!
! Physical constants (cgs units) and solar calibration values used
! throughout pyrotwind's angular-momentum-evolution routines.
!
! All of them -- including the physical constants that never change
! between runs -- are read once per run by init_const from the const
! namelist group of the run's input file (see rotwind_example.nml).
! Keeping them in the namelist (rather than hardcoding them as compile-
! time parameters) makes it possible to try alternate constants (e.g. a
! more modern G) without recompiling.
!
! cg, cc23, cpi, csecday, csecyr, solm, solr, soll are fixed
! physical/astronomical constants.
!
! solp, soltau are solar calibration values for the reference
! atmospheric pressure and convective-overturn timescale used to
! normalize the angular-momentum loss law (see setup).
!
! solw is the Sun's equatorial angular rotation rate, derived from
! solprot (the Sun's equatorial rotation period, days) once cpi and
! csecday are known -- it is not itself a namelist entry.
!
! rosun is the Sun's Rossby number (rotation period / convective-
! overturn timescale, dimensionless), derived from solprot, csecday,
! and soltau. It is used to convert params.rocrit -- which is
! expressed in solar units, i.e. rocrit = (star_prot/star_tau) /
! (solar_prot/solar_tau) -- to an absolute Rossby-number threshold
! (rocrit*rosun) wherever it is compared against a star's actual
! period/timescale ratio (see int1zone, int2zone, solidevol). Like
! solw, it is derived, not itself a namelist entry.
!
! Callers must invoke init_const before using any of these values.
!===============================================================================
module constm
    implicit none
    private
    public :: init_const
    public :: cg, cc23, cpi, csecday, csecyr, solm, solr, soll, &
              solp, soltau, solprot, solw, rosun

    ! gravitational constant [cm^3 g^-1 s^-2]. note: the value used in
    ! the yrec code (typical default 6.6725d-8) differs slightly from
    ! the modern codata value, but is within 1 sigma.
    real(8) :: cg

    ! 2/3, used in the core/envelope angular-momentum-transfer term
    real(8) :: cc23

    ! pi
    real(8) :: cpi

    ! seconds per day and per (julian) year
    real(8) :: csecday, csecyr

    ! solar mass, radius, luminosity [g, cm, erg/s]
    real(8) :: solm, solr, soll

    ! solar atmospheric pressure and convective-overturn timescale.
    ! these are local values (hp(base)/velocity one pressure scale
    ! height above the base); a global integral of dr/v across the
    ! whole convection zone would be about 3x larger (demarque & kim
    ! 1996). typical values solp = 8.7173d4, soltau = 1.0654d6.
    real(8) :: solp, soltau

    ! solar equatorial rotation period (days). typical value 25.4d0.
    real(8) :: solprot

    ! solar equatorial rotation rate, derived from solprot once cpi,
    ! csecday are known. not a namelist entry.
    real(8) :: solw

    ! solar rossby number (rotation period / convective-overturn
    ! timescale, dimensionless), derived from solprot, csecday, soltau.
    ! not a namelist entry.
    real(8) :: rosun

contains

    !---------------------------------------------------------------------
    ! Read all of the above (except solw and rosun) from the const
    ! namelist group in nml_file, and derive the solar equatorial
    ! rotation rate solw and the solar Rossby number rosun.
    !---------------------------------------------------------------------
    subroutine init_const(nml_file)
        character(*), intent(in) :: nml_file
        integer :: nml_unit

        namelist /const/ cg, cc23, cpi, csecday, csecyr, solm, solr, &
                          soll, solp, soltau, solprot

        open(newunit=nml_unit, file=nml_file, status='old', action='read')
        read(nml_unit, nml=const)
        close(nml_unit)

        ! solar equatorial rotation rate
        solw = 2.0d0*cpi/solprot/csecday

        ! solar rossby number: rotation period (seconds) / convective-
        ! overturn timescale
        rosun = solprot*csecday/soltau
    end subroutine init_const

end module constm
