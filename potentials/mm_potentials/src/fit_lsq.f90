
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

  FUNCTION ff() RESULT func
    REAL func
    REAL :: 
    CALL norm(pgfs,radius,NG,norm)
  END FUNCTION ff
  SUBROUTINE norm(pgfs,radius,NG,norm)
    IMPLICIT NONE
    INTEGER :: NG
    REAL :: radius(NG),norm
! locals
    REAL :: TT(NG),S_OVERLAP(NG,NG),S_INV(NG,NG),WW(NG)
    INTEGER :: i,j,ig
    REAL r1,r2,sqq
    do i=1,NG
      r1=radius(i)
      do j=1,NG
        r2=radius(j)
        sqq=sqrt(r2 ** 2 + r1 ** 2)
        S_OVERLAP(i,j)= 0.5*sqrt(3.141592654D0) * r1 * r2 erf(sqq*Rc/r1/r2)/sqq
      enddo
    enddo

    do i=1,NG
      TT(i)=func(1)*0.5
      do ig=2,NP-1
        TT(i)=TT(i)+func(ig)*exp(-(x(ig)/radius(i))**2)
      enddo
      TT(i)=TT(i)+0.5*func(NP)*exp(-(x(NP)/radius(i))**2)
    enddo
    TT=TT*dr
    call invert_matrix(S_OVERLAP,S_INV,NG)
    do i=1,NG
      WW(i)=0.d0
      do j=1,NG
        WW(i)=WW(i)+S_INV(i,j)*TT(j)
      enddo
    enddo
    norm=func2_av
    do i=1,NG
      norm=norm-WW(i)*TT(i)
    enddo

  END SUBROUTINE compute_norm

END MODULE gaussian_fit
