!!****h* gfit/output [1.0] * 
!!
!!   NAME
!!     output
!!
!!   FUNCTION
!!     It dumps out the results of the fit.
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
MODULE OUTPUT
  USE gaussian_fit_types,  ONLY : gaussian_fit_p_type,&
                                  gaussian_fit_type
  USE kinds, ONLY: default_string_length,&
                   dbl
CONTAINS
!!****h* gfit/output/punch_results [1.0] * 
!!
!!   NAME
!!     punch_results
!!
!!   FUNCTION
!!     Write the file of the MM Potential Fit
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
  SUBROUTINE PUNCH_RESULTS(pgfs)
    IMPLICIT NONE
    ! Arguments
    TYPE(gaussian_fit_p_type), DIMENSION(:), POINTER :: pgfs
    ! Local Variables
    TYPE(gaussian_fit_type), POINTER :: pgf
    CHARACTER(len=default_string_length) :: OUTFILE
    INTEGER :: J, I

    NULLIFY(pgf)
    OUTFILE="RADIUS"
    DO J = 1, SIZE(pgfs)
       pgf => pgfs(J)%pgf
       WRITE(OUTFILE(7:),'(F4.2)')pgf%Elp_Radius
       WRITE(6,'(2A)')"Punching Information on File: ",OUTFILE(1:11)
       ! Dumping Information on File...
       OPEN(10,file=OUTFILE,status='unknown',form='formatted')
       WRITE(10,'(A)')"&"//TRIM(pgf%IdLabel)
       WRITE(10,'(A,F12.6)')"RADIUS ",pgf%Elp_Radius
       WRITE(10,'(I3,2F7.2,2F15.9)')pgf%Number_of_Gaussians,&
                                    pgf%Info%Rmin,&
                                    pgf%Info%Rmax,&
                                    pgf%Info%ChiSq,&
                                    pgf%Info%MaxErr
       CALL refsor(pgf%Gk, pgf%Ak)
       DO I = pgf%Number_of_Gaussians, 1, -1  
          WRITE(10,'(2X,F15.9,5X,F15.9)')pgf%Ak(I), pgf%Gk(I)
       END DO
       WRITE(10,'(2X,F15.9)')pgf%A0
       WRITE(10,'(A/)')"&END"
       CLOSE(10)
       !
    END DO
  END SUBROUTINE PUNCH_RESULTS
    
!!****h* gfit/output/refsor [1.0] * 
!!
!!   NAME
!!     refsor
!!
!!   FUNCTION
!!     It sorts values in array XDONT and YDONT (comparing XDONT values only)
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
  SUBROUTINE refsor (XDONT, YDONT)
    !  Sorts XDONT into ascending order - Quicksort
    REAL (kind=dbl), DIMENSION (:), POINTER:: XDONT, YDONT

    CALL subsor (XDONT, 1, SIZE (XDONT), YDONT)
    CALL inssor (XDONT, YDONT)
    RETURN
  END SUBROUTINE refsor

!!****h* gfit/output/subsor [1.0] * 
!!
!!   NAME
!!     subsor
!!
!!   FUNCTION
!!     Sorting routines
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
  RECURSIVE SUBROUTINE subsor (XDONT, IDEB1, IFIN1, YDONT)
    !  Sorts XDONT from IDEB1 to IFIN1
    REAL(kind=dbl), DIMENSION (:), POINTER :: XDONT, YDONT
    INTEGER, INTENT (In) :: IDEB1, IFIN1
    INTEGER, PARAMETER :: NINS = 16 ! Max for insertion sort
    INTEGER :: ICRS, IDEB, IDCR, IFIN, IMIL
    REAL(kind=dbl) :: XPIV, XWRK, YPIV, YWRK
    !
    IDEB = IDEB1
    IFIN = IFIN1
    !
    IF ((IFIN - IDEB) > NINS) THEN
       IMIL = (IDEB+IFIN) / 2
       !
       IF (XDONT(IMIL) < XDONT(IDEB)) THEN
          XWRK = XDONT (IDEB)
          YWRK = YDONT (IDEB)
          XDONT (IDEB) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
          YDONT (IDEB) = YDONT (IMIL)
          YDONT (IMIL) = YWRK          
       END IF
       IF (XDONT(IMIL) > XDONT(IFIN)) THEN
          XWRK = XDONT (IFIN)
          YWRK = YDONT (IFIN)
          XDONT (IFIN) = XDONT (IMIL)
          XDONT (IMIL) = XWRK
          YDONT (IFIN) = YDONT (IMIL)
          YDONT (IMIL) = YWRK          
          IF (XDONT(IMIL) < XDONT(IDEB)) THEN
             XWRK = XDONT (IDEB)
             YWRK = YDONT (IDEB)
             XDONT (IDEB) = XDONT (IMIL)
             XDONT (IMIL) = XWRK
             YDONT (IDEB) = YDONT (IMIL)
             YDONT (IMIL) = YWRK             
          END IF
       END IF
       XPIV = XDONT (IMIL)
       YPIV = YDONT (IMIL)
       !
       ICRS = IDEB
       IDCR = IFIN
       ECH2: DO
          DO
             ICRS = ICRS + 1
             IF (ICRS >= IDCR) THEN
                EXIT ECH2
             END IF
             IF (XDONT(ICRS) > XPIV) EXIT
          END DO
          DO
             IF (XDONT(IDCR) <= XPIV) EXIT
             IDCR = IDCR - 1
             IF (ICRS >= IDCR) THEN
                EXIT ECH2
             END IF
          END DO
          !
          XWRK = XDONT (IDCR)
          YWRK = YDONT (IDCR)
          XDONT (IDCR) = XDONT (ICRS)
          XDONT (ICRS) = XWRK
          YDONT (IDCR) = YDONT (ICRS)
          YDONT (ICRS) = YWRK          
       END DO ECH2
       CALL subsor (XDONT, IDEB1, ICRS-1, YDONT)
       CALL subsor (XDONT, IDCR, IFIN1, YDONT)
    END IF
  END SUBROUTINE Subsor

!!****h* gfit/output/Inssor [1.0] * 
!!
!!   NAME
!!     Inssor
!!
!!   FUNCTION
!!     Sorting routine
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
  SUBROUTINE Inssor (XDONT, YDONT)
    !  Sorts XDONT into increasing order (Insertion sort)
    REAL(kind=dbl), DIMENSION (:), INTENT (InOut) :: XDONT, YDONT
    INTEGER :: ICRS, IDCR
    REAL(kind=dbl) :: XWRK, YWRK
    !
    DO ICRS = 2, SIZE (XDONT)
       XWRK = XDONT (ICRS)
       YWRK = YDONT (ICRS)
       IF (XWRK >= XDONT(ICRS-1)) CYCLE
       XDONT (ICRS) = XDONT (ICRS-1)
       YDONT (ICRS) = YDONT (ICRS-1)       
       DO IDCR = ICRS - 2, 1, - 1
          IF (XWRK >= XDONT(IDCR)) EXIT
          XDONT (IDCR+1) = XDONT (IDCR)
          YDONT (IDCR+1) = YDONT (IDCR)
       END DO
       XDONT (IDCR+1) = XWRK
       YDONT (IDCR+1) = YWRK
    END DO
  END SUBROUTINE Inssor

END MODULE OUTPUT
