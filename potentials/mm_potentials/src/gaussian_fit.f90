
!!****h* gfit/gaussian_fit [1.0] * 
!!
!!   NAME
!!     gaussian_fit
!!
!!   FUNCTION
!!     Performs a fit with gaussian of the MM potential energy function
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE Gaussian_fit
  USE kinds, only: dbl
  USE gaussian_fit_types, ONLY: gaussian_fit_p_type
  USE lm_gfit, only: lm_gaussian_fit
  IMPLICIT NONE
  PRIVATE


  PUBLIC :: fit

CONTAINS

  SUBROUTINE fit(pgfs, Fit_Type)
    IMPLICIT NONE
    ! Arguments
    TYPE(gaussian_fit_p_type), DIMENSION(:), POINTER :: pgfs
    INTEGER, INTENT(IN) :: Fit_Type
    ! Local Variables
    INTEGER :: I, stat
    REAL(KIND=dbl) :: ExpFac

    ExpFac = 15.0_dbl
    DO I = 1, SIZE(pgfs)
       stat = 0
       IF (Fit_Type.EQ.1) THEN
          CALL lm_gaussian_fit(pgfs(I)%pgf, stat, ExpFac, Fit_Type)
       ELSEIF (Fit_Type.EQ.3) THEN
          CALL lm_gaussian_fit(pgf=pgfs(I)%pgf, Rstat=stat, Fit_Type=Fit_Type)
       END IF
       If (stat.ne.0) EXIT
    END DO

  END SUBROUTINE FIT

END MODULE gaussian_fit
