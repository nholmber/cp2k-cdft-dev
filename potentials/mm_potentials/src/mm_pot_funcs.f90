!!****h* gfit/mm_pot_funcs [1.0] * 
!!
!!   NAME
!!     mm_pot_funcs
!!
!!   FUNCTION
!!     Expression of potentials to be fitted with gaussians
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE MM_POT_FUNCS
  USE KINDS,      ONLY: dbl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: MMF1, DMMF1, EVALMMF1, EVALDMMF1
  
CONTAINS
!!****h* mm_pot_funcs/mmf1 [1.0] * 
!!
!!   NAME
!!     mmf1
!!
!!   FUNCTION
!!     Evaluates the function (x^4-rc^4)/(x^5-rc^5)
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created
!!
!!   SOURCE
!****************************************************************************
  FUNCTION MMF1(x,rc)   RESULT(MyVal)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), INTENT(IN)  :: x, rc
    ! Local Variables
    REAL(KIND=dbl)  :: MyVal, x2, x4, x5
    REAL(KIND=dbl), SAVE :: rc2, rc4, rc5
    INTEGER         :: Ival
    DATA Ival /0/
    
    IF (Ival.EQ.0) THEN
       rc2  = rc*rc
       rc4  = rc2 * rc2
       rc5  = rc4 * rc
       Ival = 1
    END IF
    x2= x*x
    x4= x2*x2
    x5= x4*x
    MyVal = (x4-rc4)/(x5-rc5)
  END FUNCTION MMF1

  SUBROUTINE EVALMMF1(X,Y,NP,Rc)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), DIMENSION(:), POINTER :: X, Y
    REAL(KIND=dbl), INTENT(IN) :: Rc
    INTEGER, INTENT(IN) :: NP
    ! Local Variables
    INTEGER :: I

    DO I = 1, NP
       Y(I)=MMF1(x(I), Rc)
    END DO
  END SUBROUTINE EVALMMF1

!!****h* mm_pot_funcs/dmmf1 [1.0] * 
!!
!!   NAME
!!     dmmf1
!!
!!   FUNCTION
!!     Evaluates the derivative of the function (x^4-rc^4)/(x^5-rc^5)
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created
!!
!!   SOURCE
!****************************************************************************
  FUNCTION DMMF1(x,rc)   RESULT(MyVal)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), INTENT(IN)  :: x, rc
    ! Local Variables
    REAL(KIND=dbl)  :: MyVal, x2, x4, x5, a, b
    REAL(KIND=dbl), SAVE :: rc2, rc4, rc5
    INTEGER         :: Ival
    DATA Ival /0/
    
    IF (Ival.EQ.0) THEN
       rc2  = rc*rc
       rc4  = rc2 * rc2
       rc5  = rc4 * rc
       Ival = 1
    END IF
    x2= x*x
    x4= x2*x2
    x5= x4*x
    a = (x5-rc5)
    b = (x4-rc4)
    MyVal = (4.0_dbl*x*x2*a-5.0_dbl*x4*b)/(a*a)
  END FUNCTION DMMF1

  SUBROUTINE EVALDMMF1(X,Y,NP,Rc)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), DIMENSION(:), POINTER :: X, Y
    REAL(KIND=dbl), INTENT(IN) :: Rc
    INTEGER, INTENT(IN) :: NP
    ! Local Variables
    INTEGER :: I

    DO I = 1, NP
       Y(I)=DMMF1(x(I), Rc)
    END DO
  END SUBROUTINE EVALDMMF1

END MODULE MM_POT_FUNCS
