!!****h* cp2k/input [1.0] * 
!!
!!   NAME
!!     input
!!
!!   FUNCTION
!!     Reading input file input.dat
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
MODULE INPUT
  USE global_variables,    ONLY:  IONAME
  USE kinds,  ONLY:  dbl,&
                     default_string_length
  USE gaussian_fit_types,  ONLY : gaussian_fit_p_type
  IMPLICIT NONE
 
  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN = 'input'
  PUBLIC :: read_input_file

CONTAINS
!!****h* input/read_input_file [1.0] * 
!!
!!   NAME
!!     read_input_file
!!
!!   FUNCTION
!!     This one really reads input file input.dat
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
  SUBROUTINE read_input_file(pgfs, Fit_type)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(OUT) :: Fit_type
    TYPE(gaussian_fit_p_type), DIMENSION(:), POINTER :: pgfs
    ! Local Variables
    INTEGER  :: Number_of_fit, Npoint, stat, Number_of_gaussians
    INTEGER  :: i, MaxIter, nodf
    REAL(KIND=dbl), ALLOCATABLE, DIMENSION(:) :: Radius_to_fit, lb, ub
    REAL(KIND=dbl) :: Rmax, Rmin, Eps_fit
    CHARACTER(len=*), PARAMETER :: routineN = 'read_input_file', &
         routineP = moduleN//':'//routineN
    CHARACTER(len=default_string_length) :: line
    CHARACTER(len=default_string_length), DIMENSION(:), ALLOCATABLE :: IdLabel

    OPEN(10,file=IONAME,status='OLD',form='formatted')
    
    READ(10,*,END=100,ERR=200) Fit_Type
    READ(10,*,END=100,ERR=200) Number_of_fit
    ALLOCATE(Radius_to_fit(Number_of_fit),&
             IdLabel(Number_of_fit),stat=stat)
    IF (stat.NE.0) THEN
       WRITE(6,'(A)')"Error allocating vector in :"//routineN//" ."
       stop 99
    END IF
    READ(10,*,END=100,ERR=200)(Radius_to_fit(i),i=1,Number_of_fit)
    READ(10,'(A)',END=100,ERR=200)line
    DO I = 1, Number_of_fit
       IdLabel(I) = My_TRIM(line)
    END DO
    READ(10,*,END=100,ERR=200) Rmin, Rmax
    READ(10,*,END=100,ERR=200) Npoint, Number_of_gaussians
    READ(10,*,END=100,ERR=200) Eps_fit, MaxIter
    IF (Fit_type.EQ.2.OR.Fit_type.EQ.3) THEN
       READ(10,*,END=100,ERR=200) Nodf
       ALLOCATE(lb(Nodf), ub(Nodf))
       READ(10,*,END=100,ERR=200)(lb(i),i=1,Nodf)
       READ(10,*,END=100,ERR=200)(ub(i),i=1,Nodf)
    END IF

    CLOSE(10)

100 CONTINUE
    CALL build_fit_type(Number_of_fit,&
                        Radius_to_fit,&
                        IdLabel,&
                        Rmax,&
                        Rmin,&
                        Npoint,&
                        Number_of_Gaussians,&
                        Eps_fit,&
                        MaxIter,&
                        Fit_type,&
                        Nodf,&
                        lb,&
                        ub,&
                        pgfs)

    DEALLOCATE(Radius_to_fit,&
               IdLabel,stat=stat)
    IF (Fit_type.EQ.2.OR.Fit_type.EQ.3) THEN
       DEALLOCATE(lb, ub)
    END IF

    IF (stat.NE.0) THEN
       WRITE(6,'(A)')"Error deallocating vector in :"//routineN//" ."
       STOP 99
    END IF
    RETURN

200 WRITE(6,'(A)')"Error reading file: "//IONAME//" in "//routineN//" ."
    STOP 99

  END SUBROUTINE read_input_file

!!****h* input/build_fit_type [1.0] * 
!!
!!   NAME
!!     build_fit_type
!!
!!   FUNCTION
!!     It builds the gaussian_types ones all information have been read
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
  SUBROUTINE build_fit_type( NF, RF, Id, Rmax, Rmin, NP, NG, EMAX, MAXIT,&
       Fit_type, Nodf, lb, ub, pgfs)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NF, NP, NG, MAXIT, Fit_type, Nodf
    CHARACTER(len=default_string_length), DIMENSION(:) :: Id
    REAL(KIND=dbl), INTENT(IN), DIMENSION(:) :: RF, lb, ub
    REAL(KIND=dbl), INTENT(IN) :: Rmax, Rmin, EMAX
    TYPE(gaussian_fit_p_type), DIMENSION(:), POINTER :: pgfs
    ! Local Variables
    INTEGER :: stat, I, J
    CHARACTER(len=*), PARAMETER :: routineN = 'build_fit_type', &
         routineP = moduleN//':'//routineN

    ALLOCATE(pgfs(NF), stat=stat)
    IF (stat.NE.0) THEN
       WRITE(6,'(A)')"Error allocating vector in :"//routineN//" ."
       stop 99
    END IF

    DO I = 1, NF
       ALLOCATE(pgfs(I)%pgf,stat=stat)
       IF (stat.NE.0) THEN
          WRITE(6,'(A)')"Error allocating vector in :"//routineN//" ."
          STOP 99
       END IF
       pgfs(I)%pgf%number_of_gaussians   = NG
       pgfs(I)%pgf%info%number_of_points = NP
       pgfs(I)%pgf%Elp_Radius            = RF(I)
       pgfs(I)%pgf%IdLabel               = Id(I)
       pgfs(I)%pgf%info%Eps_fit          = EMAX
       pgfs(I)%pgf%info%MaxIter          = MAXIT
       pgfs(I)%pgf%info%Rmax             = Rmax
       pgfs(I)%pgf%info%Rmin             = Rmin
       IF (Fit_type.EQ.2.OR.Fit_type.EQ.3) THEN
          pgfs(I)%pgf%info2%nodf            = Nodf
          ALLOCATE( pgfs(I)%pgf%info2%lb(Nodf),&
                    pgfs(I)%pgf%info2%ub(Nodf) )
          DO J = 1, nodf
             pgfs(I)%pgf%info2%lb(j)        = lb(j)
             pgfs(I)%pgf%info2%ub(j)        = ub(j)
          END DO
       END IF

       ALLOCATE(pgfs(I)%pgf%Ak(NG),stat=stat)
       IF (stat.NE.0) THEN
          WRITE(6,'(A)')"Error allocating vector in :"//routineN//" ."
          STOP 99
       END IF       
       ALLOCATE(pgfs(I)%pgf%Gk(NG),stat=stat)
       IF (stat.NE.0) THEN
          WRITE(6,'(A)')"Error allocating vector in :"//routineN//" ."
          STOP 99
       END IF       
       pgfs(I)%pgf%Ak = 0.0_dbl
       pgfs(I)%pgf%Gk = 0.0_dbl
       pgfs(I)%pgf%A0 = 0.0_dbl
    END DO

  END SUBROUTINE build_fit_type
  
!!****h* input/my_trim [1.0] * 
!!
!!   NAME
!!     my_trim
!!
!!   FUNCTION
!!     It trims words from a line
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
  FUNCTION MY_TRIM(LINE) RESULT (ID)
    IMPLICIT NONE
    ! Arguments
    CHARACTER(len=default_string_length), INTENT(INOUT) :: line
    CHARACTER(len=default_string_length) :: ID
    ! Local Variables
    INTEGER :: I
    
    DO WHILE (line(1:1) == " ")
       line(1:) = line(2:)
    END DO
    I = 1
    DO WHILE (line(i:i) .ne. " ")
       i = i+1
    END DO
    ID = line(1:I-1)
    line(1:) = line(i:)
  END FUNCTION MY_TRIM

END MODULE INPUT

