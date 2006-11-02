!------------------------------------------------------------------------------!
  module rint
!------------------------------------------------------------------------------!
    USE basic_data_types, ONLY: dp
    IMPLICIT NONE
    SAVE
!...Parameters for Gauss-Hermite quadrature
!...Limit for number of points in Gauss-Hermite quadrature
    integer,parameter :: ippmax=2000
!...Number of points in Gauss-Hermite quadrature
    integer           :: ippn
!...Points and weights for Gauss-Hermite quadrature
    REAL(dp)           :: xip(ippmax),ri(ippmax),wi(ippmax)
!------------------------------------------------------------------------------!
  end module rint
!------------------------------------------------------------------------------!


