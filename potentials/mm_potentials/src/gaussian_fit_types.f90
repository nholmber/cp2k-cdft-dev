MODULE gaussian_fit_types
  USE kinds, ONLY: dbl,&
                   default_string_length
  IMPLICIT NONE
  PRIVATE

!!****s* gfit/Gaussian_fit_type [1.0] *
!!
!!   NAME
!!     Gaussian_fit_types
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Laino Teodoro
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
  TYPE Info_new_fit
     INTEGER         :: Nodf
      REAL(KIND=dbl), DIMENSION(:), POINTER :: ub, lb
  END TYPE Info_new_fit

  TYPE Info_Pot_Type
     REAL(KIND=dbl)  :: Rmax, Rmin
     REAL(KIND=dbl)  :: MaxErr, ChiSq, Eps_fit
     INTEGER         :: Number_of_points
     INTEGER         :: MaxIter
  END TYPE Info_Pot_Type

  TYPE Gaussian_fit_type
     INTEGER :: Number_of_gaussians
     REAL(KIND=dbl)  :: Elp_Radius
     CHARACTER(len=default_string_length) :: IdLabel
     TYPE(Info_Pot_Type)  :: Info
     TYPE(Info_new_fit)   :: Info2
     REAL(KIND=dbl)  :: A0 
     REAL(KIND=dbl), DIMENSION(:), POINTER :: Ak, Gk
  END TYPE Gaussian_fit_type

!!****s* gfit/Gaussian_fit_p_type [1.0] *
!!
!!   NAME
!!     Gaussian_fit_p_type
!!
!!   FUNCTION
!!     represent a pointer to a Gaussian_fit_type, to be able to create arrays
!!     of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - Gaussian_fit_p_type: the pointer to the Gaussian_fit_type
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************

  TYPE Gaussian_fit_p_type
     TYPE(Gaussian_fit_type), POINTER :: pgf
  END TYPE Gaussian_fit_p_type


  PUBLIC :: gaussian_fit_type,&
            gaussian_fit_p_type,&
            Info_Pot_type,&
            Info_new_fit

END MODULE gaussian_fit_types
