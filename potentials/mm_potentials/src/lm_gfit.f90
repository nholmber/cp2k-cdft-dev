
!!****h* gfit/lm_gfit [1.0] * 
!!
!!   NAME
!!     lm_gfit
!!
!!   FUNCTION
!!     Non-linear Fit with Gaussian Functions. Levenberg-Marquardt Method
!!
!!   NOTES
!!     Based on MINPACK Project. NETLIB.
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE lm_gfit
  USE gaussian_fit_types,       only: gaussian_fit_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: zero,&
                                             one
  USE mm_pot_funcs,                    ONLY: evalmmf1

  IMPLICIT NONE
  PRIVATE

  REAL(KIND=dbl),   PRIVATE, PARAMETER :: rdwarf = 3.834E-20_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: rgiant = 1.304E19_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p1     = 1.000E-1_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p5     = 5.000E-1_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p25    = 2.500E-1_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p75    = 7.500E-1_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p05    = 5.000E-2_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p001   = 1.000E-3_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: p0001  = 1.000E-4_dbl
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: epsmch = EPSILON(0.0_dbl)
  REAL(KIND=dbl),   PRIVATE, PARAMETER :: dwarf  = TINY(0.0_dbl)
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='lm_gfit'
  LOGICAL, PRIVATE, PARAMETER          :: debug_this_module=.true.

  PUBLIC :: lm_gaussian_Fit, CreateData

CONTAINS

!!****f* Lm_gfit/lm_gaussian_fit [1.0] *
!!
!!   NAME
!!    lm_gaussian_fit
!!
!!   FUNCTION
!!     Performs a Gaussian of the QMMM Electrostatic Potential Function
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE Lm_Gaussian_Fit( pgf, Rstat, ExpFac, Fit_type)
    IMPLICIT NONE
    ! Arguments
    TYPE(gaussian_fit_type), POINTER                 :: pgf
    INTEGER, INTENT(INOUT)                           :: Rstat
    INTEGER, INTENT(IN)                              :: Fit_Type
    REAL(KIND=dbl), INTENT(IN), OPTIONAL             :: ExpFac
    ! Local Variables
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'Gaussian_Fit', &
         routineP = moduleN//':'//routineN
    INTEGER :: Ma, Ndata
    INTEGER, DIMENSION(:), POINTER :: Ia, IPvt
    REAL(KIND=dbl) :: Fnorm, Elp_Radius, Rmin, Rmax, LocFac
    REAL(KIND=dbl), DIMENSION(:),   POINTER :: X, Y, A, Work, Fvec
    REAL(KIND=dbl), DIMENSION(:,:), POINTER :: FJac
    INTEGER :: stat, I, II, LdFJac, Lwa, Info
    ! Statements
    LocFac = 1.0_dbl
    IF (PRESENT(ExpFac)) LocFac = ExpFac
    Rstat   = 10
    failure = .false.

    failure=.FALSE.
    NULLIFY(Ia, X, Y, A, Work, FJac, Ipvt)
    ! Allocate Working Spaces
    Ma     = pgf%Number_of_Gaussians*2+1
    Ndata  = pgf%Info%Number_of_Points*NINT(LocFac)
    Lwa    = 5*Ma+Ndata
    LdFJac = Ndata
    ALLOCATE(Ia(Ma), stat=stat)

    ALLOCATE(X(Ndata), stat=stat)

    ALLOCATE(Y(Ndata), stat=stat)

    Elp_Radius = pgf%Elp_Radius*LocFac
    Rmax       = pgf%Info%Rmax*LocFac
    Rmin       = pgf%Info%Rmin*LocFac
    CALL CreateData(X, Y, Ndata, Elp_Radius, Rmin, Rmax)
    ALLOCATE(A(Ma), stat=stat)

    IF (Fit_Type.EQ.1) THEN
       ! Initialize the parameter set
       CALL Initialize_Gauss_Parameters(A, pgf%Number_of_Gaussians,&
                                        X, Y, Ndata, Elp_Radius )
    ELSE
       DO I = 1, pgf%Number_of_Gaussians
          II       = (I-1)*2
          A(II+1)  = pgf%Ak(I)
          A(II+2)  = pgf%Gk(I)
       END DO
    ENDIF
    A(Ma)      = 0.0_dbl
    WRITE (*,'(A)')" STARTING PARAMETERS:"
    WRITE (*,'(2F12.6)')(A(I),I=1,Ma)

    ALLOCATE(Work(Lwa), stat=stat)

    ALLOCATE(FJac(Ndata,Ma), stat=stat)

    ALLOCATE(Fvec(Ndata), stat=stat)

    ALLOCATE(Ipvt(Ma), stat=stat)

    !
    ! This is the Main Fit Loop
    ! 
    CALL LmFit(Ndata, Ma, A, Fvec, X, Y, FJac, LdFJac, pgf%Info%Eps_Fit, Info, IPvt, Work, Lwa,&
               pgf%Info%MaxErr, pgf%Info%ChiSq, pgf%Info%MaxIter )
    FNORM = ENORM(Ndata,FVEC)
    WRITE (*,'(A,F15.9)')"FIT RESULTS ON RADIUS:",pgf%Elp_Radius
    WRITE (*,'(A,F15.9,A,F15.9)')" CHISQ:",FNORM," MAXIMUM ERROR:",pgf%Info%MaxErr
    WRITE (*,'(A)')" OPTIMIZED PARAMETERS:"
    !
    ! Transfer Fit Information into the pgf Type
    ! 
    DO I = 1, pgf%Number_of_Gaussians
       II        = (I-1)*2
       pgf%Ak(I) =     A(II+1)  * LocFac 
       pgf%Gk(I) = ABS(A(II+2)) / LocFac 
       WRITE (*,'(2F12.6)')pgf%Ak(I),  pgf%Gk(I)
    END DO 
    pgf%A0       = A(Ma) * LocFac      
    WRITE (*,'(F12.6)')pgf%A0
    If (Info.ne.0) Rstat = 0
    !
    ! Deallocate all the temporary arrays
    !
    DEALLOCATE(Ia, stat=stat)

    DEALLOCATE(A, stat=stat)

    DEALLOCATE(X, stat=stat)

    DEALLOCATE(Y, stat=stat)

    DEALLOCATE(Work, stat=stat)

    DEALLOCATE(FJac, stat=stat)

    DEALLOCATE(Fvec, stat=stat)

    DEALLOCATE(Ipvt, stat=stat)


  END SUBROUTINE Lm_Gaussian_Fit

!!****f* lm_gfit/evaluate_gaussian [1.0] *
!!
!!   NAME
!!    evaluate_gaussian
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE evaluate_gaussian(xt,yt,m,n,a,fvec,fjac,iflag,MaxErr) 
    ! Arguments
    INTEGER, INTENT(IN) ::  n, m, iflag
    REAL(KIND=dbl), INTENT(INOUT), DIMENSION(m)   :: fvec
    REAL(KIND=dbl), INTENT(INOUT), DIMENSION(m,n) :: fjac
    REAL(KIND=dbl), INTENT(IN) :: a(n), xt(m), yt(m)
    REAL(KIND=dbl), INTENT(INOUT) ::  MaxErr
    ! Local Variables
    INTEGER ::  i, j
    REAL(KIND=dbl) ::  arg,ex,fac,y, Dert

    IF (Iflag.EQ.1) THEN
       MaxErr = -1.0_dbl
       DO j=1,m
          y=0.0_dbl
          ! Sum of Gaussians
          ! Function value 
          !     WRITE(*,*)"POINT : ",j, xt(j)
          DO i=1,n-2,2
             !  WRITE(*,*)"     GAUSSIAN : ",i, a(i), a(i+1)
             arg=xt(j)/a(i+1)
             ex=EXP(-arg**2)
             y=y+a(i)*ex                     ! Function value
          END DO
          ! Constant Term
          y= y + a(n) + Term(n,a,xt,m)
          fvec(j) = yt(j) - y
          MaxErr = MAX(MaxErr, Abs(fvec(j)))
       END DO
    ELSEIF(Iflag.EQ.2) THEN
       ! Sum of Gaussians
       fjac = 0.0_dbl
       LoopOnParameters: DO i=1,n-2,2
          Dert = DTerm(i,n,a,xt,m)
          LoopOnPoints: DO j=1,m
             ! Function derivatives
             arg=xt(j)/a(i+1)
             ex=EXP(-arg**2)
             fac        =   a(i)*ex*2.0_dbl*arg
             fjac(j,i)  = - ex                        ! Derivatives w.r.t. Ak
             fjac(j,i+1)= - fac*arg/a(i+1)  - Dert    ! Derivatives w.r.t. Gk
          END DO LoopOnPoints
       END DO LoopOnParameters
       LoopOnPointsConst: DO j=1,m
          ! Constant Term
          fjac(j,n) = - 1.0_dbl
       END DO LoopOnPointsConst
    END IF 
    
  END SUBROUTINE evaluate_gaussian

  FUNCTION Term(n,a,xt,m) RESULT(MyVal)
    IMPLICIT NONE
    INTEGER :: k, n, m, i
    REAL(KIND=dbl), DIMENSION(:) :: a, xt
    REAL(KIND=dbl) :: g, g1, g2, MyVal
    
    MyVal = 0.0_dbl
    DO i = 3, n-2, 2
       DO k = 1, i-2, 2
          g1 = a(i+1) / a(k+1)
          g2 = 1.0_dbl / g1
          g  = xt(m) * ABS(g1-g2)
          MyVal = MyVal + EXP(-g) 
       END DO
    END DO
  END FUNCTION Term

  FUNCTION DTerm(i,n,a,xt,m) RESULT(MyVal)
    IMPLICIT NONE
    INTEGER :: k, n, m, i
    REAL(KIND=dbl), DIMENSION(:) :: a, xt
    REAL(KIND=dbl) :: g, g1, g2, MyVal
    
    MyVal = 0.0_dbl
    DO k = 1, n-2, 2
       IF (k.NE.i) THEN
          g1 = a(i+1) / a(k+1)
          g2 = 1.0_dbl / g1
          g  = xt(m) * ABS(g1-g2)
          g  = - g * SIGN(g1-g2,.0_dbl) * EXP(-g) * ( 1.0_dbl /a(k+1) + a(k+1)/a(i+1)**2)
          MyVal = MyVal + g
       END IF
    END DO
  END FUNCTION DTerm

!!****f* lm_gfit/CreateData [1.0] *
!!
!!   NAME
!!    CreateData
!!
!!   FUNCTION
!!     Creates a Set of Data Points of the MM Potential In Input
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE CreateData(X, Y, N, Rc, Xmin, Xmax)
    Implicit None
    ! Arguments
    REAL(KIND=dbl), INTENT(IN) :: Rc, Xmin, Xmax
    REAL(KIND=dbl), DIMENSION(:), POINTER :: X, Y
    INTEGER, INTENT(IN) :: N
    ! Local Variables
    INTEGER :: I
    REAL(KIND=dbl) :: Dx

    Dx = (Xmax-Xmin)/REAL(N-1,KIND=dbl)
    X(1) = Xmin
    DO I = 2, N
       X(I) = X(I-1) + Dx
    END DO
    CALL  EVALMMF1(X,Y,N,Rc)
    WRITE(77,'(2F12.6)')(X(I),Y(I),I=1,N)

  END SUBROUTINE CreateData
               

!!****f* lm_gfit/LmFit [1.0] *
!!
!!   NAME
!!     LmFit
!!
!!   FUNCTION
!!     LmderD minimizes the sum of the square of m nonlinear
!!     functions in n variables by a Modification of the 
!!     Levenberger-Marquardt algorithm. 
!!     This version is adapted to GFIT. Rewritten in F90.
!!     The function is supposed to be a sum of gaussian functions.
!!
!!   NOTES
!!     Based on...
!!     MinPack: Jorge More', Burt Garbow, and Ken Hillstrom 
!!              @ Argonne National Laboratory
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE LmFit (M,N,x,fvec,xt,yt,fjac,Ldfjac,tol,Info,Ipvt,wa,Lwa,&
                    MaxErr, ChiSq, MaxIter)
    ! Arguments
    INTEGER, INTENT(IN)  :: M, N, Ldfjac, Lwa, MaxIter
    INTEGER, INTENT(OUT) :: Info
    INTEGER, DIMENSION(:), POINTER :: Ipvt
    REAL(KIND=dbl), INTENT(IN)  ::  Tol
    REAL(KIND=dbl), INTENT(OUT) ::  MaxErr, ChiSq
    REAL(KIND=dbl), POINTER, DIMENSION(:)       :: x, xt, yt, fvec
    REAL(KIND=dbl), DIMENSION(Lwa), INTENT(OUT) :: wa
    REAL(KIND=dbl), POINTER, DIMENSION(:,:)     :: fjac    
    ! Local Variables
    INTEGER :: Maxfev, Mode, Nfev, Njev, Nprint
    REAL(KIND=dbl) :: gtol, xtol, ftol
    REAL(KIND=dbl), PARAMETER :: factor=1.0E2_dbl

    Info = 0
    !
    !     check the input parameters for errors.
    !
    IF (.NOT.((n.LE.0).OR.(m.LT.n).OR.(ldfjac.LT.m).OR.(tol.LT.zero) &
         .OR.(lwa.LT.5*n+m))) THEN
       maxfev = MaxIter
       WRITE(*,*)tol, maxiter
       ftol   = tol
       xtol   = tol
       gtol   = zero
       mode   = 1
       nprint = 1

       CALL LmIterate ( xt=xt,&            ! X data from Rmin to Rmax
                        yt=yt,&            ! Y data of the MM Potential
                        m=m,&              ! Number of Points
                        n=n,&              ! Number of Parameters
                        x=x,&              ! Parameters' value
                        fvec=fvec,&        ! In Output the Error of the Fit 
                        fjac=fjac,&        ! In Output the Jacobian of the Fit Function
                        ldfjac=ldfjac,&    ! Leading dimension of the Jacobian
                        ftol=ftol,&        ! => Tolerance on function fit
                        xtol=xtol,&        ! => Tolerance on parameters value
                        gtol=gtol,&        ! => Tolerance on gradients
                        maxfev=maxfev,&    ! => Maximum number of Iterations
                        diag=wa(1),&       ! Working Array
                        mode=mode,&        ! Working Array
                        factor=factor,&
                        nprint=nprint,&    ! => nprint>0 Print Information on MaxError and ChiSq
                        info=info,&        ! = > Information Output value
                        nfev=nfev,&
                        njev=njev,&
                        ipvt=ipvt,&
                        qtf=wa(n+1),&
                        wa1=wa(2*n+1),&
                        wa2=wa(3*n+1),&
                        wa3=wa(4*n+1),&
                        wa4=wa(5*n+1),&
                        MaxErr=MaxErr,&
                        fnorm=ChiSq )
       IF (info .EQ. 8) info = 4
    END IF
    
    WRITE(*,'(A,I5)')"Fit Ended. Info:",Info
    SELECT CASE(Info)
    CASE(0)
       WRITE(*,'(3X,A)')"Improper input parameters."
    CASE(1)
       WRITE(*,'(3X,A)')"Both actual and predicted relative reductions"//&
                        "in the sum of squares are at most ftol."
    CASE(2)
       WRITE(*,'(3X,A)')"Relative error between two consecutive iterates"//&
                        "is at most xtol."
    CASE(3)
       WRITE(*,'(3X,A)')"Both actual and predicted relative reductions"//&
                        "in the sum of squares are at most ftol.",&
                        "Relative error between two consecutive iterates"//&
                        "is at most xtol."
    CASE(4)
       WRITE(*,'(3X,A)')"The cosine of the angle between fvec and any"//&
                        "column of the jacobian is at most gtol in"//&
                        "absolute value."
    CASE(5)
       WRITE(*,'(3X,A)')"Maximum Number of iteration reached. Please Increase MAXITER."
    CASE(6)
       WRITE(*,'(3X,A)')"Ftol is too small. no further reduction in"//&
                        "the sum of squares is possible."
    CASE(7) 
       WRITE(*,'(3X,A)')"Xtol is too small. no further improvement in"//&
                        "the approximate solution x is possible."
    CASE(8)
       WRITE(*,'(3X,A)')"Gtol is too small. fvec is orthogonal to the"//& 
                        "columns of the jacobian to machine precision." 
    CASE DEFAULT
       WRITE(*,'(3X,A,I5,A)')"Info Value:",Info," No Information Available!"
    END SELECT
  END SUBROUTINE LmFit

!!****f* lm_gfit/Enorm [1.0] *
!!
!!   NAME
!!     Enorm
!!
!!   FUNCTION
!!     given an n-vector x, this function calculates the
!!     euclidean norm of x.
!!
!!     the euclidean norm is computed by accumulating the sum of
!!     squares in three different sums. the sums of squares for the
!!     small and large components are scaled so that no overflows
!!     occur. non-destructive underflows are permitted. underflows
!!     and overflows do not occur in the computation of the unscaled
!!     sum of squares for the intermediate components.
!!     the definitions of small, intermediate and large components
!!     depend on two constants, rdwarf and rgiant. the main
!!     restrictions on these constants are that rdwarf**2 not
!!     underflow and rgiant**2 not overflow. the constants
!!     given here are suitable for every known computer.   
!!
!!   NOTES
!!     Based on...
!!     MinPack: Jorge More', Burt Garbow, and Ken Hillstrom 
!!              @ Argonne National Laboratory
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*************************************************************************

  FUNCTION Enorm(N,x) RESULT(Res)
    ! Arguments
    INTEGER,        INTENT(IN) :: N
    REAL(KIND=dbl), INTENT(IN), DIMENSION(N) :: x
    REAL(KIND=dbl) :: Res
    ! Local Variables
    INTEGER :: i
    REAL(KIND=dbl) :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max

    s1 = Zero
    s2 = Zero
    s3 = Zero
    x1max  = Zero
    x3max  = Zero
    floatn = REAL(N,KIND=dbl)
    agiant = rgiant/floatn
    DO i = 1, n
       xabs = dabs(x(i))
       IF (xabs .LE. rdwarf .OR. xabs .GE. agiant) THEN
          IF (xabs .GT. rdwarf) THEN
             !
             ! sum for large components.
             !
             IF (xabs .GT. x1max) THEN
                s1 = one + s1*(x1max/xabs)**2
                x1max = xabs
             ELSE
                s1 = s1 + (xabs/x1max)**2
             END IF
          ELSE
             !
             ! sum for small components.
             !
             IF (xabs .GT. x3max) THEN
                s3 = one + s3*(x3max/xabs)**2
                x3max = xabs
             ELSE
                IF (xabs .NE. zero) s3 = s3 + (xabs/x3max)**2
             END IF
          ENDIF
       ELSE
          !
          ! sum for intermediate components.
          !
          s2 = s2 + xabs**2
       END IF
    END DO
    !
    !     calculation of norm.
    !
    IF (s1 .NE. zero) THEN 
       Res = x1max*dsqrt(s1+(s2/x1max)/x1max)
    ELSE
       IF (s2 .NE. zero) THEN
          IF (s2 .GE. x3max) Res = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
          IF (s2 .LT. x3max) Res = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
       ELSE
          Res = x3max*dsqrt(s3)
       ENDIF
    ENDIF
  END FUNCTION ENORM

!!****f* lm_gfit/LmIterate [1.0] *
!!
!!   NAME
!!     LmIterate
!!
!!   FUNCTION
!!     Minimize the sum of the squares of m nonlinear functions in n
!!     variables by a modification of the Levenberg-Marquardt algorithm
!!
!!   NOTES
!!     Based on...
!!     MinPack: Jorge More', Burt Garbow, and Ken Hillstrom 
!!              @ Argonne National Laboratory
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE LmIterate(xt, yt, m, n, x, fvec, fjac, ldfjac, ftol, xtol,&
                       gtol, maxfev, diag, mode, factor, nprint, info, &
                       nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4, MaxErr,&
                       fnorm)
    ! Arguments
    INTEGER, INTENT(IN) :: m, n, ldfjac, maxfev, mode, nprint
    INTEGER, INTENT(OUT):: nfev, njev, info
    INTEGER, INTENT(OUT), DIMENSION(N) :: ipvt(n)
    REAL(KIND=dbl), INTENT(IN)  :: ftol, xtol, gtol, factor
    REAL(KIND=dbl), INTENT(OUT) :: MaxErr, Fnorm
    REAL(KIND=dbl), INTENT(OUT),   DIMENSION(M)   :: fvec, wa4
    REAL(KIND=dbl), INTENT(OUT),   DIMENSION(N)   :: wa1, wa2, wa3, &
                                                     diag, qtf
    REAL(KIND=dbl), INTENT(IN ),   DIMENSION(M)   :: xt, yt
    REAL(KIND=dbl), INTENT(INOUT), DIMENSION(N)   :: x
    REAL(KIND=dbl), INTENT(OUT),   DIMENSION(LDFJAC,N) :: fjac
    ! Local Variables
    INTEGER        :: i,iflag,iter,j,l
    REAL(KIND=dbl) :: actred,delta,dirder,fnorm1,gnorm,&
                      par,pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm

    
    info  = 0
    iflag = 0
    nfev  = 0
    njev  = 0
    !
    !     check the input parameters for errors.
    !
    IF (.NOT.((n.LE.0).OR.(m.LT.n).OR.(ldfjac.LT.m).OR.(ftol.LT.zero)&
               .OR.(xtol.LT.zero).OR.(gtol.LT.zero).OR.(maxfev.LE.0)&
               .OR.(factor.LE.zero))) THEN 
       !
       !     evaluate the function at the starting point
       !     and calculate its norm.
       !
       iflag = 1
       CALL evaluate_gaussian(xt,yt,m,n,x,fvec,fjac,iflag,MaxErr)
       nfev  = 1
       fnorm = enorm(m,fvec)
       !
       !     initialize levenberg-marquardt parameter and iteration counter.
       !
       par  = zero
       iter = 1
       !
       !     beginning of the outer loop.
       !
       OuterLoop: DO 
          !
          !        calculate the jacobian matrix.
          !
          iflag = 2
          CALL evaluate_gaussian(xt,yt,m,n,x,fvec,fjac,iflag,MaxErr)
          njev = njev + 1
          !
          !        if requested, call fcn to enable printing of iterates.
          !
          IF (nprint .GT. 0) THEN
             IF (MOD(iter-1,nprint) .EQ. 0) &
                  WRITE(*,*)"ITER: ",iter," MAX ERROR:",MaxErr," CHISQ:",fnorm
          END IF
          !
          !        compute the qr factorization of the jacobian.
          !
          CALL qrfac(m,n,fjac,ldfjac,.TRUE.,ipvt,n,wa1,wa2,wa3)
          !
          !        on the first iteration and if mode is 1, scale according
          !        to the norms of the columns of the initial jacobian.
          !
          IF (iter .EQ. 1) THEN 
             IF (mode .NE. 2) THEN
                DO j = 1, n
                   diag(j) = wa2(j)
                   IF (wa2(j) .EQ. zero) diag(j) = one
                END DO
             END IF
             !
             !        on the first iteration, calculate the norm of the scaled x
             !        and initialize the step bound delta.
             !
             DO j = 1, n
                wa3(j) = diag(j)*x(j)
             END DO
             xnorm = enorm(n,wa3)
             delta = factor*xnorm
             IF (delta .EQ. zero) delta = factor
          END IF
          !
          !        form (q transpose)*fvec and store the first n components in
          !        qtf.
          !
          DO i = 1, m
             wa4(i) = fvec(i)
          END DO

          DO j = 1, n
             IF (fjac(j,j) .NE. zero) THEN
                sum = zero
                DO i = j, m
                   sum = sum + fjac(i,j)*wa4(i)
                END DO
                temp = -sum/fjac(j,j)
                DO i = j, m
                   wa4(i) = wa4(i) + fjac(i,j)*temp
                END DO
             END IF
             fjac(j,j) = wa1(j)
             qtf(j) = wa4(j)
          END DO
          !
          !        compute the norm of the scaled gradient.
          !
          gnorm = zero
          IF (fnorm .NE. zero) THEN
             DO j = 1, n
                l = ipvt(j)
                IF (wa2(l) .NE. zero) THEN
                   sum = zero
                   DO i = 1, j
                      sum = sum + fjac(i,j)*(qtf(i)/fnorm)
                   END DO
                   gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
                END IF
             END DO
          END IF
          !
          !        test for convergence of the gradient norm.
          !
          IF (gnorm .LE. gtol) info = 4
          IF (info .NE. 0) EXIT OuterLoop
          !
          !        rescale if necessary.
          !
          IF (mode .NE. 2) THEN
             DO j = 1, n
                diag(j) = dmax1(diag(j),wa2(j))
             END DO
          END IF
          !
          !        beginning of the inner loop.
          !
          InnerLoop: DO
             !
             !           determine the levenberg-marquardt parameter.
             !
             CALL lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,&
                        wa3,wa4)
             !
             !           store the direction p and x + p. calculate the norm of p.
             !
             DO j = 1, n
                wa1(j) = -wa1(j)
                wa2(j) = x(j) + wa1(j)
                wa3(j) = diag(j)*wa1(j)
             END DO
             pnorm = enorm(n,wa3)
             !
             !           on the first iteration, adjust the initial step bound.
             !
             IF (iter .EQ. 1) delta = dmin1(delta,pnorm)
             !
             !           evaluate the function at x + p and calculate its norm.
             !
             iflag = 1
             CALL evaluate_gaussian(xt,yt,m,n,wa2,wa4,fjac,iflag,MaxErr)
             nfev = nfev + 1
             fnorm1 = enorm(m,wa4)
             !
             !           compute the scaled actual reduction.
             !
             actred = -one
             IF (p1*fnorm1 .LT. fnorm) actred = one - (fnorm1/fnorm)**2
             !
             !           compute the scaled predicted reduction and
             !           the scaled directional derivative.
             !
             DO j = 1, n
                wa3(j) = zero
                l = ipvt(j)
                temp = wa1(l)
                DO i = 1, j
                   wa3(i) = wa3(i) + fjac(i,j)*temp
                END DO
             END DO
             temp1 = enorm(n,wa3)/fnorm
             temp2 = (dsqrt(par)*pnorm)/fnorm
             prered = temp1**2 + temp2**2/p5
             dirder = -(temp1**2 + temp2**2)
             !
             !           compute the ratio of the actual to the predicted
             !           reduction.
             !
             ratio = zero
             IF (prered .NE. zero) ratio = actred/prered
             !
             !           update the step bound.
             !
             IF (ratio .LE. p25) THEN
                IF (actred .GE. zero) temp = p5
                IF (actred .LT. zero) &
                     temp = p5*dirder/(dirder + p5*actred)
                IF (p1*fnorm1 .GE. fnorm .OR. temp .LT. p1) temp = p1
                delta = temp*dmin1(delta,pnorm/p1)
                par = par/temp
             ELSE
                IF (par .EQ. zero .OR. ratio .GE. p75) THEN
                   delta = pnorm/p5
                   par = p5*par
                END IF
             END IF
             !
             !           test for successful iteration.
             !
             IF (ratio .GE. p0001) THEN
                !
                !           successful iteration. update x, fvec, and their norms.
                !
                DO j = 1, n
                   x(j) = wa2(j)
                   wa2(j) = diag(j)*x(j)
                END DO
                DO i = 1, m
                   fvec(i) = wa4(i)
                END DO
                xnorm = enorm(n,wa2)
                fnorm = fnorm1
                iter = iter + 1
             END IF
             !
             !           tests for convergence.
             !
             IF (dabs(actred) .LE. ftol .AND. prered .LE. ftol&
                  .AND. p5*ratio .LE. one) info = 1
             IF (delta .LE. xtol*xnorm) info = 2
             IF (dabs(actred) .LE. ftol .AND. prered .LE. ftol&
                  .AND. p5*ratio .LE. one .AND. info .EQ. 2) info = 3
             IF (info .NE. 0) EXIT OuterLoop
             !
             !           tests for termination and stringent tolerances.
             !
             IF (iter .GE. maxfev) info = 5
             IF (dabs(actred) .LE. epsmch .AND. prered .LE. epsmch&
                  .AND. p5*ratio .LE. one) info = 6
             IF (delta .LE. epsmch*xnorm) info = 7
             IF (gnorm .LE. epsmch) info = 8
             IF (info .NE. 0) EXIT OuterLoop
             !
             !           end of the inner loop. repeat if iteration unsuccessful.
             !
             IF (ratio .GE. p0001) EXIT InnerLoop
          END DO InnerLoop
          !
          !        end of the outer loop.
          !
       END DO OuterLoop
    ENDIF
  END SUBROUTINE LmIterate

!!****f* lm_gfit/LmPar [1.0] *
!!
!!   NAME
!!     LmPar
!!
!!   FUNCTION
!!     Given an m by n matrix a, an n by n nonsingular diagonal
!!     matrix d, an m-vector b, and a positive number delta,
!!     the problem is to determine a value for the parameter
!!     par such that if x solves the system
!!
!!           a*x = b ,     sqrt(par)*d*x = 0 ,
!!
!!     in the least squares sense, and dxnorm is the euclidean
!!     norm of d*x, then either par is zero and
!!
!!           (dxnorm-delta) .le. 0.1*delta ,
!!
!!     or par is positive and
!!
!!           abs(dxnorm-delta) .le. 0.1*delta .
!!
!!     this subroutine completes the solution of the problem
!!     if it is provided with the necessary information from the
!!     qr factorization, with column pivoting, of a. that is, if
!!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!!     columns, and r is an upper triangular matrix with diagonal
!!     elements of nonincreasing magnitude, then lmpar expects
!!     the full upper triangle of r, the permutation matrix p,
!!     and the first n components of (q transpose)*b. on output
!!     lmpar also provides an upper triangular matrix s such that
!!
!!            t   t                   t
!!           p *(a *a + par*d*d)*p = s *s .
!!
!!     s is employed within lmpar and may be of separate interest.
!!
!!   NOTES
!!     Based on...
!!     MinPack: Jorge More', Burt Garbow, and Ken Hillstrom 
!!              @ Argonne National Laboratory
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE LmPar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)
    ! Arguments
    INTEGER, INTENT(IN) :: n, ldr
    INTEGER, INTENT(IN), DIMENSION(N) ::  ipvt
    REAL(KIND=dbl), INTENT(IN)        ::  delta
    REAL(KIND=dbl), INTENT(INOUT)     ::  par
    REAL(KIND=dbl), INTENT(INOUT), DIMENSION(ldr,n) :: r
    REAL(KIND=dbl), INTENT(IN),    DIMENSION(N)  :: diag, qtb
    REAL(KIND=dbl), INTENT(OUT),   DIMENSION(N)  :: x, wa1, wa2, sdiag
    ! Local Variables
    INTEGER                   :: i, iter, j, jm1, jp1, k, l, nsing
    REAL(KIND=dbl)            :: dxnorm, dwarf, fp, gnorm, parc,&
                                 parl, paru, sum, temp

    !
    !     compute and store in x the gauss-newton direction. if the
    !     jacobian is rank-deficient, obtain a least squares solution.
    !
    nsing = n
    DO j = 1, n
       wa1(j) = qtb(j)
       IF (r(j,j) .EQ. zero .AND. nsing .EQ. n) nsing = j - 1
       IF (nsing .LT. n) wa1(j) = zero
    END DO
    IF (nsing .GE. 1) THEN
       DO k = 1, nsing
          j = nsing - k + 1
          wa1(j) = wa1(j)/r(j,j)
          temp = wa1(j)
          jm1 = j - 1
          IF (jm1 .GE. 1) THEN
             DO i = 1, jm1
                wa1(i) = wa1(i) - r(i,j)*temp
             END DO
          END IF
       END DO
    END IF
    DO j = 1, n
       l = ipvt(j)
       x(l) = wa1(j)
    END DO
    !
    !     initialize the iteration counter.
    !     evaluate the function at the origin, and test
    !     for acceptance of the gauss-newton direction.
    !
    iter = 0
    DO j = 1, n
       wa2(j) = diag(j)*x(j)
    END DO
    dxnorm = enorm(n,wa2)
    fp = dxnorm - delta
    IF (fp .GT. p1*delta) THEN
       !
       !     if the jacobian is not rank deficient, the newton
       !     step provides a lower bound, parl, for the zero of
       !     the function. otherwise set this bound to zero.
       !
       parl = zero
       IF (nsing .GE. n) THEN
          DO j = 1, n
             l = ipvt(j)
             wa1(j) = diag(l)*(wa2(l)/dxnorm)
          END DO
          DO j = 1, n
             sum = zero
             jm1 = j - 1
             IF (jm1 .GE. 1) THEN
                DO i = 1, jm1
                   sum = sum + r(i,j)*wa1(i)
                END DO
             ENDIF
             wa1(j) = (wa1(j) - sum)/r(j,j)
          END DO
          temp = enorm(n,wa1)
          parl = ((fp/delta)/temp)/temp
       END IF
       !
       !     calculate an upper bound, paru, for the zero of the function.
       !
       DO j = 1, n
          sum = zero
          DO i = 1, j
             sum = sum + r(i,j)*qtb(i)
          END DO
          l = ipvt(j)
          wa1(j) = sum/diag(l)
       END DO
       gnorm = enorm(n,wa1)
       paru = gnorm/delta
       IF (paru .EQ. zero) paru = dwarf/dmin1(delta,p1)
       !
       !     if the input par lies outside of the interval (parl,paru),
       !     set par to the closer endpoint.
       !
       par = dmax1(par,parl)
       par = dmin1(par,paru)
       IF (par .EQ. zero) par = gnorm/dxnorm
       !
       !     beginning of an iteration.
       !
       LoopOnIter: DO
          iter = iter + 1
          !
          !        evaluate the function at the current value of par.
          !
          IF (par .EQ. zero) par = dmax1(dwarf,p001*paru)
          temp = dsqrt(par)
          DO j = 1, n
             wa1(j) = temp*diag(j)
          END DO
          CALL qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
          DO j = 1, n
             wa2(j) = diag(j)*x(j)
          END DO
          dxnorm = enorm(n,wa2)
          temp = fp
          fp = dxnorm - delta
          !
          !        if the function is small enough, accept the current value
          !        of par. also test for the exceptional cases where parl
          !        is zero or the number of iterations has reached 10.
          !
          IF (dabs(fp) .LE. p1*delta&
               .OR. parl .EQ. zero .AND. fp .LE. temp&
               .AND. temp .LT. zero .OR. iter .EQ. 10) EXIT LoopOnIter
          !
          !        compute the newton correction.
          !
          DO j = 1, n
             l = ipvt(j)
             wa1(j) = diag(l)*(wa2(l)/dxnorm)
          END DO
          DO j = 1, n
             wa1(j) = wa1(j)/sdiag(j)
             temp = wa1(j)
             jp1 = j + 1
             IF (n .GE. jp1) THEN
                DO i = jp1, n
                   wa1(i) = wa1(i) - r(i,j)*temp
                END DO
             END IF
          END DO
          temp = enorm(n,wa1)
          parc = ((fp/delta)/temp)/temp
          !
          !        depending on the sign of the function, update parl or paru.
          !
          IF (fp .GT. zero) parl = dmax1(parl,par)
          IF (fp .LT. zero) paru = dmin1(paru,par)
          !
          !        compute an improved estimate for par.
          !
          par = dmax1(parl,par+parc)
          !
          !        end of an iteration.
          !
       END DO LoopOnIter
    END IF
    !
    !     termination.
    !
    IF (iter .EQ. 0) par = zero
    RETURN
  END SUBROUTINE LmPar

!!****f* lm_gfit/QrSolv [1.0] *
!!
!!   NAME
!!     QrSolv
!!
!!   FUNCTION
!!     Complete Solution of least square problem
!!
!!   NOTES
!!     Based on...
!!     MinPack: Jorge More', Burt Garbow, and Ken Hillstrom 
!!              @ Argonne National Laboratory
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE QrSolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
    ! Arguments
    INTEGER, INTENT(IN) :: n, ldr
    INTEGER, INTENT(IN), DIMENSION(N) :: ipvt
    REAL(KIND=dbl), INTENT(IN),  DIMENSION(N) :: diag, qtb
    REAL(KIND=dbl), INTENT(OUT), DIMENSION(N) :: x, sdiag, wa
    REAL(KIND=dbl), DIMENSION(ldr,n), INTENT(INOUT) :: r
    ! Local Variables
    INTEGER :: i, j, jp1, k, kp1, l, nsing
    REAL(KIND=dbl) :: cos, cotan, qtbpj, sin, sum, tan, temp

    !
    !     copy r and (q transpose)*b to preserve input and initialize s.
    !     in particular, save the diagonal elements of r in x.
    !
    DO j = 1, n
       DO i = j, n
          r(i,j) = r(j,i)
       END DO
       x(j) = r(j,j)
       wa(j) = qtb(j)
    END DO
    !
    !     eliminate the diagonal matrix d using a givens rotation.
    !
    DO j = 1, n
       !
       !        prepare the row of d to be eliminated, locating the
       !        diagonal element using p from the qr factorization.
       !
       l = ipvt(j)
       IF (diag(l) .NE. zero) THEN
          DO k = j, n
             sdiag(k) = zero
          END DO
          sdiag(j) = diag(l)
          !
          !        the transformations to eliminate the row of d
          !        modify only a single element of (q transpose)*b
          !        beyond the first n, which is initially zero.
          !
          qtbpj = zero
          DO k = j, n
             !
             !           determine a givens rotation which eliminates the
             !           appropriate element in the current row of d.
             !
             IF (sdiag(k) .NE. zero) THEN
                IF (dabs(r(k,k)) .LT. dabs(sdiag(k))) THEN
                   cotan = r(k,k)/sdiag(k)
                   sin = p5/dsqrt(p25+p25*cotan**2)
                   cos = sin*cotan
                ELSE
                   tan = sdiag(k)/r(k,k)
                   cos = p5/dsqrt(p25+p25*tan**2)
                   sin = cos*tan
                END IF
                !
                !           compute the modified diagonal element of r and
                !           the modified element of ((q transpose)*b,0).
                !
                r(k,k) = cos*r(k,k) + sin*sdiag(k)
                temp = cos*wa(k) + sin*qtbpj
                qtbpj = -sin*wa(k) + cos*qtbpj
                wa(k) = temp
                !
                !           accumulate the tranformation in the row of s.
                !
                kp1 = k + 1
                IF (n .GE. kp1) THEN
                   DO i = kp1, n
                      temp = cos*r(i,k) + sin*sdiag(i)
                      sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
                      r(i,k) = temp
                   END DO
                END IF
             END IF
          END DO
       END IF
       !
       !        store the diagonal element of s and restore
       !        the corresponding diagonal element of r.
       !
       sdiag(j) = r(j,j)
       r(j,j) = x(j)
    END DO
    !
    !     solve the triangular system for z. if the system is
    !     singular, then obtain a least squares solution.
    !
    nsing = n
    DO j = 1, n
       IF (sdiag(j) .EQ. zero .AND. nsing .EQ. n) nsing = j - 1
       IF (nsing .LT. n) wa(j) = zero
    END DO
    IF (nsing .GE. 1) THEN
       DO k = 1, nsing
          j = nsing - k + 1
          sum = zero
          jp1 = j + 1
          IF (nsing .GE. jp1) THEN
             DO i = jp1, nsing
                sum = sum + r(i,j)*wa(i)
             END DO
          END IF
          wa(j) = (wa(j) - sum)/sdiag(j)
       END DO
    END IF
    !
    !     permute the components of z back to components of x.
    !
    DO j = 1, n
       l = ipvt(j)
       x(l) = wa(j)
    END DO
  END SUBROUTINE QrSolv

!!****f* lm_gfit/QrFac [1.0] *
!!
!!   NAME
!!     QrFac
!!
!!   FUNCTION
!!     Compute qr factorization after row addition
!!
!!   NOTES
!!     Based on...
!!     MinPack: Jorge More', Burt Garbow, and Ken Hillstrom 
!!              @ Argonne National Laboratory
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE QrFac( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)
    ! Arguments
    INTEGER, INTENT(IN)                    :: m, n, lda, lipvt
    INTEGER, INTENT(OUT), DIMENSION(lipvt) :: ipvt
    LOGICAL, INTENT(IN)                    :: pivot
    REAL(KIND=dbl), INTENT(INOUT), DIMENSION(lda,n) :: a
    REAL(KIND=dbl), INTENT(OUT), DIMENSION(n)       :: rdiag, acnorm, wa
    ! Local Variables
    INTEGER :: i, j, jp1, k, kmax, minmn
    REAL(KIND=dbl) :: ajnorm, epsmch, sum, temp


    epsmch = EPSILON(0.0_dbl)
    !
    !     compute the initial column norms and initialize several arrays.
    !
    DO j = 1, n
       acnorm(j) = enorm(m,a(1,j))
       rdiag(j) = acnorm(j)
       wa(j) = rdiag(j)
       IF (pivot) ipvt(j) = j
    END DO
    !
    !     reduce a to r with householder transformations.
    !
    minmn = min0(m,n)
    DO j = 1, minmn
       IF (pivot) THEN
          !
          !        bring the column of largest norm into the pivot position.
          !
          kmax = j
          DO k = j, n
             IF (rdiag(k) .GT. rdiag(kmax)) kmax = k
          END DO
          IF (kmax .NE. j) THEN
             DO i = 1, m
                temp = a(i,j)
                a(i,j) = a(i,kmax)
                a(i,kmax) = temp
             END DO
             rdiag(kmax) = rdiag(j)
             wa(kmax) = wa(j)
             k = ipvt(j)
             ipvt(j) = ipvt(kmax)
             ipvt(kmax) = k
          END IF
       END IF
       !
       !        compute the householder transformation to reduce the
       !        j-th column of a to a multiple of the j-th unit vector.
       !
       ajnorm = enorm(m-j+1,a(j,j))
       IF (ajnorm .NE. zero) THEN
          IF (a(j,j) .LT. zero) ajnorm = -ajnorm
          DO i = j, m
             a(i,j) = a(i,j)/ajnorm
          END DO
          a(j,j) = a(j,j) + one
          !
          !        apply the transformation to the remaining columns
          !        and update the norms.
          !
          jp1 = j + 1
          IF (n .GE. jp1) THEN
             DO k = jp1, n
                sum = zero
                DO i = j, m
                   sum = sum + a(i,j)*a(i,k)
                END DO
                temp = sum/a(j,j)
                DO i = j, m
                   a(i,k) = a(i,k) - temp*a(i,j)
                END DO
                IF (pivot .AND. rdiag(k) .NE. zero) THEN
                   temp = a(j,k)/rdiag(k)
                   rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
                   IF (p05*(rdiag(k)/wa(k))**2 .LE. epsmch) THEN
                      rdiag(k) = enorm(m-j,a(jp1,k))
                      wa(k) = rdiag(k)
                   END IF
                END IF
             END DO
          END IF
       END IF
       rdiag(j) = -ajnorm
    END DO
  END SUBROUTINE QrFac

!!****f* lm_gfit/Initizalize_Gauss_Parameters [1.0] *
!!
!!   NAME
!!     Initialize_Gauss_Parameters
!!
!!   FUNCTION
!!     Initialize the parameters value for the fit
!!
!!   NOTES
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!   SOURCE
!!
!!*** **********************************************************************
  SUBROUTINE Initialize_Gauss_Parameters(A, Nog, X, Y, Ndata, Rc)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), DIMENSION(:),   POINTER :: X, Y, A
    INTEGER, INTENT(IN) :: Ndata, Nog
    REAL(KIND=dbl), INTENT(IN) :: Rc
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'Initialize_Gauss_Parameters', &
         routineP = moduleN//':'//routineN
    REAL(KIND=dbl), DIMENSION(:),   POINTER :: Yloc, Yder
    INTEGER :: stat, Iopt, I, J, IP, IQ, iun
    LOGICAL :: failure
    REAL(KIND=dbl) :: MyVal, Fval

    Failure = .false.
    NULLIFY(Yloc, Yder)
    ALLOCATE(Yloc(Ndata), stat=stat)

    ALLOCATE(Yder(Ndata), stat=stat)

    Yloc = Y
    DO J = 1, Ndata
       Yder(J) = (4.0_dbl*X(J)**3*(X(J)**5-Rc**5)-5.0_dbl*X(J)**4*(X(J)**4-Rc**4))/&
                 (X(J)**5-Rc**5)**2
    END DO

    DO I = 1, Nog
       WRITE(*,'(A,I5)')' Guessing Parameter for Gaussian number=',I
       IP = I  - 1       
       IQ = IP - 1
       IF (IP.GT.0) THEN
          iun = 10 + IP
          open(iun)
          DO J = 1, Ndata           
             Yloc(J) = Yloc(J) - A(IQ*2+1)*EXP(-(X(J)/A(IQ*2+2))**2)
             Yder(J) = Yder(J) + A(IQ*2+1)*EXP(-(X(J)/A(IQ*2+2))**2)&
                                 *2.D0*X(J)/(A(IQ*2+2))**2
             WRITE(iun,*)X(J),Yloc(J)
          END DO
          close(iun)
       END IF

       InnerLoop: DO J = Ndata-1, 2, -1
          MyVal = 0.001_dbl
          IF(ABS(Yloc(J)).GT.Myval) THEN
             Iopt = J
             A(IP*2+2) = ABS(-Yder(Iopt)/(2.d0*X(Iopt)*Yloc(Iopt)))
             A(IP*2+1) = Yloc(Iopt)/Exp(-A(IP*2+2)*X(iopt)**2)
             A(IP*2+2) = SQRT(1.D0/A(IP*2+2))
             
             Fval = A(IP*2+1)* EXP(-(X(1)/A(IP*2+2))**2)
             WRITE(*,*)I,J,A(IP*2+1),A(IP*2+2),Fval, Yloc(1)
             IF (Fval.LT.Yloc(1)) EXIT InnerLoop

          ENDIF
       END DO InnerLoop

       WRITE(*,*)I,Iopt,A(IP*2+1),A(IP*2+2)

    END DO

    DEALLOCATE(Yloc,stat=stat)

    DEALLOCATE(Yder,stat=stat)


  END SUBROUTINE Initialize_Gauss_Parameters

END MODULE lm_gfit
