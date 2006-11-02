!------------------------------------------------------------------------------!
  module atom
!------------------------------------------------------------------------------!
    USE basic_data_types, ONLY: dp
    IMPLICIT NONE
    SAVE
!..General parameters
!..Limits for Nalpha and Lmax
    integer,parameter :: Namax=20
    integer,parameter :: Lamax=3
!..Highest angular momentum quantum number
    integer :: Lmax
!..Number of exponents alpha for each angular momentum
    integer :: Nalpha(0:Lamax)
!..Exponents of the Gaussian basis functions
    REAL(dp) :: alpha(Namax,0:Lamax),alphamin,alphamax
!..Optimization parameters for the exponents
    integer :: alpp(Namax,0:Lamax)
!..Energy gradients of the exponents of the Gaussian basis functions
    REAL(dp) :: galpha(Namax,0:Lamax),g_constr(Namax,0:Lamax)
    REAL(dp) :: g_ekin(Namax,0:Lamax),g_epot(Namax,0:Lamax)
    REAL(dp) :: g_epotadd(Namax,0:Lamax)
    REAL(dp) :: g_exc(Namax,0:Lamax),g_ecoul(Namax,0:Lamax)
!..Number of occupied orbitals for each angular momentum
    integer :: Nocc(0:Lamax)
!..Occupation numbers
    REAL(dp) :: occ(Namax,0:Lamax)
!..Effectiv core charge
    REAL(dp) :: Zval,qgrid,qerror
!..Additional parabolic potential prefactor
    REAL(dp) :: add_pot_factor
!..Parameter for the mixing of old and new density matrix
    REAL(dp) :: dmix,dmix1
!..Required accuracy for Etot in the LDA-SCF calculation
    REAL(dp) :: epsscf
!..Maximum number of SCF steps
    integer :: maxscf
!.. Name of the atom
    character :: atomname*2
!..Exponent optimization (.true.) or total energy calculation (.false.)
    logical :: optexp
!..Required upper limit for gradient convergence
    REAL(dp) :: epsgrad
!..All electron (.true.) or pseudopotential (.false.) calculation
    logical :: allelectron

!------------------------------------------------------------------------------!
  end module atom
!------------------------------------------------------------------------------!
