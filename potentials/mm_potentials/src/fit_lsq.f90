
!!****h* gfit/fit_lsq [1.0] * 
!!
!!   NAME
!!     fit_lsq
!!
!!   FUNCTION
!!     Performs a fit with gaussian of the MM potential energy function
!!     New Method
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Alessandro Laio && Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE Fit_lsq
  USE kinds, ONLY: dbl
  USE mathconstants, ONLY: pi
  USE gaussian_fit_types, ONLY: gaussian_fit_type,&
                                gaussian_fit_p_type
  USE lm_gfit, ONLY: CreateData
  USE erf_fn, ONLY: erf
  USE mm_pot_funcs, ONLY : MMF1
  USE fft, ONLY: rdft

  IMPLICIT NONE
  PRIVATE
  REAL(KIND=dbl), PRIVATE :: func2_av
  INTEGER, DIMENSION(:),          POINTER, PRIVATE  :: IPIV
  REAL(KIND=dbl), DIMENSION(:),   POINTER, PRIVATE  :: TT, WW, func, X, radius
  REAL(KIND=dbl), DIMENSION(:,:), POINTER, PRIVATE  :: S_INV, S_OVERLAP
  REAL(KIND=dbl), PRIVATE :: Rc, Rmax, Rmin, Dr
  INTEGER, PRIVATE :: Ng, Np, ind_int
  LOGICAL, PRIVATE :: tstart
  LOGICAL, PARAMETER, PRIVATE :: debug=.false.

  PUBLIC:: Eval_Opt

CONTAINS

  FUNCTION F2(X) RESULT(MyVal)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), INTENT(IN)  :: x
    REAL(KIND=dbl)  :: MyVal

    MyVal = MMF1(x,Rc) * MMF1(x,Rc)
    
  END FUNCTION F2

  FUNCTION TF(X) RESULT(MyVal)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), INTENT(IN)  :: x
    REAL(KIND=dbl)  :: MyVal

    MyVal = MMF1(x,Rc) * EXP(-(x/radius(ind_int))**2)

  END FUNCTION TF

  SUBROUTINE Fit_lsq_init(pgf)
    IMPLICIT NONE
    TYPE(gaussian_fit_type), POINTER :: pgf

    Rc = pgf%Elp_Radius
    Ng = pgf%Number_of_gaussians
    Np = pgf%Info%Number_of_points
    Rmax = pgf%Info%Rmax
    Rmin = pgf%Info%Rmin
    dr = (Rmax-Rmin)/REAL(Np,kind=dbl)

    NULLIFY(TT, S_OVERLAP, S_INV, IPIV, WW, func, X, RADIUS)

    ALLOCATE( TT(NG),&
              RADIUS(NG),&
              S_OVERLAP(NG,NG),&
              S_INV(NG,NG),&
              IPIV(NG),&
              WW(NG),&
              func(Np),&
              X(Np) )

    CALL CreateData(X, func, Np, Rc, Rmin, Rmax)    
    CALL integrate(F2,Rmin,Rmax,func2_av)

  END SUBROUTINE Fit_lsq_init


  SUBROUTINE Fit_lsq_terminate
    IMPLICIT NONE

    DEALLOCATE( TT,&
                RADIUS,&
                S_OVERLAP,&
                S_INV,&
                IPIV,&
                WW,&
                X,&
                func)
    
  END SUBROUTINE Fit_lsq_terminate


!!****h* fit_lsq/eval_norm [1.0] * 
!!
!!   NAME
!!     eval_norm
!!
!!   FUNCTION
!!     Evaluates the norm of the function to be minimized
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Alessandro Laio && Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE eval_norm(radius,norm)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dbl), DIMENSION(:)     :: radius
    REAL(KIND=dbl), INTENT(OUT)      :: norm
    ! Local Variables
    INTEGER        :: i, j,  info
    REAL(KIND=dbl) :: r1, r2, sqq 


    S_OVERLAP = 0.0_dbl
    S_INV     = 0.0_dbl
    TT        = 0.0_dbl
    WW        = 0.0_dbl
    IPIV      = 0

    DO i=1,NG
       r1=radius(i)
       DO j=1,NG
          r2=radius(j)
          sqq = DSQRT(r2 ** 2 + r1 ** 2)
          S_OVERLAP(i,j) = 0.5_dbl * DSQRT(Pi) * r1 * r2 * erf(sqq*Rmax/(r1*r2)) / sqq
       ENDDO
       S_INV(i,i) = 1.0_dbl
    ENDDO
    
    IF (debug) THEN
       WRITE(*,*)"S_OVERLAP:"
       WRITE(*,'(8f15.9)')S_OVERLAP
    END IF

    CALL integrate_v(TF,Rmin, Rmax, TT)

    IF (debug) THEN
       WRITE(*,*)"TT:"
       WRITE(*,'(8f15.9)')TT
    END IF
    !
    CALL DGESV( Ng, Ng, S_OVERLAP, Ng, IPIV, S_INV, Ng, INFO )
    !
    IF (debug) THEN
       WRITE(*,*)"INFO: ",INFO
    END IF
    IF(INFO.NE.0) THEN
       WRITE(*,*)"ERROR IN DGESV!"
       STOP 9
    END IF
    IF (debug) THEN
       WRITE(*,*)"S_INV:"
       WRITE(*,'(8f30.9)')S_INV
    END IF
    DO i=1,NG
       WW(i)=0.d0
       DO j=1,NG
          WW(i)=WW(i)+S_INV(i,j)*TT(j)
       ENDDO
    ENDDO
    norm=func2_av
    DO i=1,NG
       norm=norm-WW(i)*TT(i)
    ENDDO
    IF (debug) THEN
       WRITE(*,*)"NORM: ", norm
    END IF
    norm = dsqrt(abs(norm))

    IF (debug) THEN
       WRITE(*,*)"WW:"
       WRITE(*,'(8f15.9)')WW
       WRITE(*,*)"norm: ",norm
    END IF

  END SUBROUTINE eval_norm

  SUBROUTINE Eval_Opt(pgfs, Npar)
    IMPLICIT NONE
    ! Arguments
    TYPE(gaussian_fit_p_type), DIMENSION(:), POINTER :: pgfs
    INTEGER, INTENT(IN) :: Npar
    ! Local Variables
    INTEGER :: ItMax, K, I, J
    REAL(KIND=dbl) :: tol, tolbr, fret, optf, chisq
    REAL(KIND=dbl), DIMENSION(:), ALLOCATABLE :: Par, diff, w_fft
    LOGICAL :: ERR
    !
    INTEGER, ALLOCATABLE, dimension(:) :: ip_fft
    
    
    DO K = 1, SIZE(pgfs)
       ! Initialize Optimizer...
       CALL Fit_lsq_init(pgfs(K)%pgf)

       ItMax = 1000
       tol   = 1.d-5
       tolbr = 1.d-5
       ALLOCATE(Par(Npar), diff(0:Np-1), ip_fft(0:Np), w_fft(0:Np*5/4-1) )
       ! Precondition Parameters
       CALL Precond_Opt(Par, pgfs(K)%pgf%info2%lb, pgfs(K)%pgf%info2%ub)
       Err   = .FALSE.
       WRITE(*,*) "INITIAL PARAMETERS: ",Par
       CALL powell(Itmax, Npar, tol, tolbr, par, fret, err)
       WRITE(*,*) "OPTIMIZED PARAMETERS: ",Par
       WRITE(*,*) "OPTIMIZED AMPLITUDES AND RADIUS: "
       !
       CALL build_Radius(Par)
       CALL eval_norm(radius,optf)
       !
       DO I = 1, Ng
          WRITE(*,'(2F12.6)')WW(i), radius(i)
       END DO
       WRITE(*,*)"CHISQ INTEGRATED: ",optf
       chisq = 0.0_dbl
       DO I =1, Np
          optf  = 0.0_dbl
          DO J =1 ,Ng
             optf = optf + WW(J) * EXP(-(X(I)/radius(j))**2)
          END DO
          chisq = chisq + (func(I)-optf)**2
          diff(I-1) = func(I)-optf
          diff(I-1) = diff(I-1) !* (1.-(x(i)/(rmax-1))**6)/(1.-(x(i)/(rmax-1.))**12)
          WRITE(76,'(2F12.6)')X(I),diff(I-1)
       END DO
       WRITE(*,*)"CHISQ SUMMED UP: ",sqrt(chisq)
       ! Terminate optimizer..
       DO I = 1, pgfs(K)%pgf%Number_of_Gaussians
          pgfs(K)%pgf%Ak(I) = WW(i)
          pgfs(K)%pgf%Gk(I) = radius(i)
       END DO
       !
       ! Fourier Transform
       ip_fft(0) = 0
       WRITE(*,*)'fft:',np,ip_fft(0)
       CALL rdft(np, 1, diff, ip_fft, w_fft)
       WRITE(75,'(F15.9)')(diff(I),I=0,np-1)
       !
       CALL  Fit_lsq_terminate
       DEALLOCATE(Par, diff, ip_fft, w_fft)
    END DO

  END SUBROUTINE Eval_Opt

  SUBROUTINE distgen(Npar, par, fret)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Npar
    REAL(KIND=dbl), INTENT(OUT)     :: fret
    REAL(KIND=dbl), DIMENSION(Npar) :: Par
    
    radius = 0.0_dbl
    CALL build_radius(par)
    CALL eval_norm(radius,fret)

  END SUBROUTINE distgen

  SUBROUTINE Precond_Opt(Par,Lb,Ub)
    IMPLICIT NONE
    REAL(KIND=dbl), DIMENSION(:) :: Par, Lb, Ub
    INTEGER :: I, J, NpL
    REAL(KIND=dbl) :: Dp1, Dp2, fret, locfret
    REAL(KIND=dbl), DIMENSION(2) :: LocP

    NpL  = 30
    Dp1  = (Ub(1)-Lb(1))/REAL(NpL,kind=dbl)
    Dp2  = (Ub(2)-Lb(2))/REAL(NpL,kind=dbl)
    fret = HUGE(0.0_dbl) 
    DO I = 1, NpL 
       LocP(1) = lb(1) + REAL((I-1),kind=dbl)*Dp1
       DO J = 1, NpL
          LocP(2) = lb(2) + REAL((J-1),kind=dbl)*Dp2
          CALL distgen(2, LocP, locfret)
          WRITE(88,'(3f15.9)')LocP,locfret
          IF (locfret.LT.fret) THEN
             fret = locfret
             Par  = LocP
          END IF
       END DO
    END DO

  END SUBROUTINE Precond_Opt

  SUBROUTINE build_radius(par)
    IMPLICIT NONE
    REAL(KIND=dbl), DIMENSION(:) :: PAR
    INTEGER :: I

    DO I = 1, Ng
!       Radius(I) = PAR(1) + (REAL(I,KIND=dbl)/Ng)**(1.0_dbl/PAR(2))*Rmax ! Alessandro
       Radius(I) = PAR(1) * PAR(2)**(I-1)  ! Fawzi
    END DO
    IF (debug) THEN
       WRITE(*,*)"Computing gaussian radius: PAR(1):",PAR(1),"  PAR(2):",PAR(2)
       DO i = 1, Ng
          WRITE(*,*)I, radius(i), 1.0_dbl/(radius(i))**2
       END DO
    END IF
    
  END SUBROUTINE build_radius


  SUBROUTINE brent(Npar,ax,bx,cx,tol,fmin,xmin,par,xi)	
    IMPLICIT NONE	
    INTEGER ::  Npar
    REAL(KIND=dbl) :: par(Npar),xi(Npar)
    REAL(KIND=dbl) :: ax,bx,cx,tol,fmin,xmin
    
    INTEGER :: iter
    REAL(KIND=dbl) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,dum
    REAL(KIND=dbl) :: tol1,tol2,u,v,w,x,xm
    
    INTEGER, PARAMETER :: ITMAX=400
    REAL(KIND=dbl), PARAMETER :: CGOLD=.3819660                   
    REAL(KIND=dbl), PARAMETER :: ZEPS=1.0e-8
	
    e=0.0
    a=ax
    b=cx                           
    IF(ax.GT.cx)THEN
       a=cx
       b=ax
    ENDIF
    x=bx
    w=bx
    v=bx
    fx=fmin               !e' fb calcolato da mnbrak
    fv=fmin
    fw=fmin
        
    DO iter=1,ITMAX
       xm=.5*(a+b)                 !valor medio a e b
       tol1=tol*ABS(x)+ZEPS
       tol2=2.*tol1
       IF (ABS(x-xm).LE.(tol2-.5*(b-a))) THEN
          !se |x-xm|<2*tol*|x|-.5*(b-a) uscire       
          xmin=x
          fmin=fx
          WRITE(6,*)'stop brent',iter,'  xmin=',xmin,'  azmin=',fmin
          RETURN
       ENDIF

       IF (ABS(e).GT.tol1) THEN     ! parabolic fit (x,v,w):u
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.0*(q-r)
          IF(q.GT.0.0) p=-p
          q=ABS(q)
          etemp=e
          e=d
          IF(ABS(p).GE.ABS(0.5*q*etemp).OR.p.LE.q*(a-x).OR.p.GE.q*(b-x)) THEN  !parabolic step non acc.
             e=b-x
             IF(x.GE.xm) e=a-x
             d=CGOLD*e                       !step scelto con golden rule
          ELSE        !parabolic step accettato
             d=p/q
             u=x+d
             dum=xm-x
             IF((u-a).LT.tol2.OR.(b-u).LT.tol2) d=SIGN(tol1,dum)
          ENDIF
       ELSE           !niente parabolic fit
          e=b-x
          IF(x.GE.xm) e=a-x
          d=CGOLD*e             !step scelto con golden rule
       ENDIF

       u=x+SIGN(tol1,d)
       IF(ABS(d).GE.tol1) u=x+d
       CALL fline(Npar,par,xi,u,fu)

       IF(fu.LT.fx)THEN
          IF(u.GE.x) a=x 
          IF(u.LT.x) b=x
          CALL shft(v,w,x,u)
          CALL shft(fv,fw,fx,fu)
       ELSE
          IF(u.LT.x) a=u 
          IF(u.GE.x) b=u
          IF(fu.LE.fw.OR.w.EQ.x) THEN
             v=w
             w=u
             fv=fw
             fw=fu
          ELSE IF(fu.LE.fv.OR.v.EQ.x.OR.v.EQ.w) THEN
             v=u
             fv=fu
          ENDIF
       ENDIF
    ENDDO
    WRITE(6,*)'raggiunto MAXITER in brent'
    WRITE(6,*)x,fx
    xmin=x
    fmin=fx
    
    RETURN
  END SUBROUTINE brent

  SUBROUTINE fline(Npar,par,xi,x,fx)    
    IMPLICIT NONE
    
    INTEGER :: Npar
    REAL(KIND=dbl) :: par(Npar),xi(Npar)
    REAL(KIND=dbl) :: x,fx    
    REAL(KIND=dbl) :: parx(Npar)
    INTEGER :: ip
    
    DO ip=1,Npar
       parx(ip)=par(ip)+xi(ip)*x
    ENDDO
    CALL distgen(Npar,parx,fx)
    RETURN
  END SUBROUTINE fline

  SUBROUTINE linmin(Npar,ilin,par,xi,fret,tolbr,err)
    IMPLICIT NONE
    INTEGER :: Npar
    REAL(KIND=dbl) :: par(Npar),xi(Npar),fret,tolbr     
    REAL(KIND=dbl) ::  ax,bx,cx,fa,fb,fc,xmin      
    INTEGER :: i,ilin
    LOGICAL :: err
      
    ilin=ilin+1
    ax=0.d0
    bx=1.d-1
    cx=2.d-1
    CALL  mnbrak(Npar,ax,bx,cx,fa,fb,fc,par,xi,err) 
    IF(err)STOP 'fail mnbrak'
    fret=fb
    CALL  brent(Npar,ax,bx,cx,tolbr,fret,xmin,par,xi)
    DO i=1,Npar
       par(i)=par(i)+xmin*xi(i)
    ENDDO

11  FORMAT(i3,2x,2(f8.4,1x),12(f7.4,1x))

    RETURN
  END SUBROUTINE linmin

  SUBROUTINE mnbrak(Npar,ax,bx,cx,fa,fb,fc,par,xi,err)
    IMPLICIT NONE
    INTEGER  :: Npar
    REAL(KIND=dbl) :: par(Npar),xi(Npar) 
    REAL(KIND=dbl) :: ax,bx,cx,fa,fb,fc
    REAL(KIND=dbl) :: ulim,u,r,q,fu,dum
    INTEGER :: iter
    LOGICAL :: err
    INTEGER       , PARAMETER :: ITMAX=20
    REAL(KIND=dbl), PARAMETER :: GOLD=1.618034d0                   
    REAL(KIND=dbl), PARAMETER :: TINY=1.0d-10
    REAL(KIND=dbl), PARAMETER :: GLIMIT=100.d0
	
    CALL fline(Npar,par,xi,ax,fa)
    CALL fline(Npar,par,xi,bx,fb)
    IF (fb.GT.fa) THEN
       dum=ax
       ax= bx
       bx=dum
       dum=fa
       fa= fb
       fb=dum
    ENDIF
    cx=bx+GOLD*(bx-ax)
    CALL fline(Npar,par,xi,cx,fc)
    DO iter=1,ITMAX
       IF (fc.GT.fb) THEN
          RETURN
       ELSE
          r=(bx-ax)*(fb-fc)
          q=(bx-cx)*(fb-fa)
          dum=q-r
          u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(dum),TINY),dum))
          ulim=bx+GLIMIT*(cx-bx)
          IF((bx-u)*(u-cx).GT.0.d0) THEN       !minimo parab. tra b e c
             CALL fline(Npar,par,xi,u,fu)
             IF (fu.LT.fc) THEN                 !minimo tra b e c
		ax=bx
		bx=u
		fa=fb
		fb=fu
		RETURN
             ELSE IF (fu.GT.fb) THEN              !minimo tra a e u
		cx=u
		fc=fu
		RETURN
             ENDIF
             u=cx+GOLD*(cx-bx)                 !fit parabolico inutile
             CALL fline(Npar,par,xi,u,fu)
          ELSE IF ((cx-u)*(u-ulim).GT.0.d0) THEN 
             CALL fline(Npar,par,xi,u,fu)
             IF (fu.LT.fc) THEN
		bx=cx
                cx=u
                u=cx+GOLD*(cx-bx)
                fb=fc
                fc=fu
	        CALL fline(Npar,par,xi,u,fu)
             ENDIF
          ELSE IF ((u-ulim)*(ulim-cx).GE.0.d0) THEN
	      u=ulim
	      CALL fline(Npar,par,xi,u,fu)
           ELSE
	      u=cx+GOLD*(cx-bx)
	      CALL fline(Npar,par,xi,u,fu)
           ENDIF
           CALL shft(ax,bx,cx,u)
           CALL shft(fa,fb,fc,fu)
        ENDIF
     ENDDO
     WRITE(6,*) 'non riesce a brakettare il minimo'
     WRITE(6,*)fa,fb,fc
     err=.TRUE.     
   END SUBROUTINE mnbrak

   SUBROUTINE shft(a,b,c,d)
     IMPLICIT NONE
     REAL(kind=DBL) :: a,b,c,d
     a=b
     b=c
     c=d
     RETURN
   END SUBROUTINE shft


   SUBROUTINE powell(Itmax,Npar,tol,tolbr,par,fret,err)
     IMPLICIT NONE 
     INTEGER  :: Npar,Itmax
     REAL(KIND=dbl) ::  par(Npar)
     REAL(KIND=dbl) :: fret,tol,tolbr
     
     REAL(KIND=dbl) :: xi(Npar,Npar)
     INTEGER :: i,j,iter,ibig,ilin
     REAL(KIND=dbl) :: pt(Npar),xit(Npar),ptt(Npar)
     REAL(KIND=dbl) :: fp,del,fptt,t
     LOGICAL :: err
     
     tstart=.TRUE.
     ilin=0

     DO i=1,Npar
        pt(i)=par(i)
     ENDDO

     DO i=1,Npar
        DO j=1,Npar
           xi(i,j)=0.d0
           IF(i.EQ.j)xi(i,j)=1.d0
        ENDDO
     ENDDO

     CALL distgen(Npar,par,fret)

11   FORMAT(i3,2x,f8.4,1x,12(f7.4,1x))
     
     iter=0     
1    iter=iter+1 
     fp=fret
     ibig=0
     del=0.d0
     DO i=1,Npar
        DO j=1,Npar
           xit(j)=xi(j,i)
        ENDDO
        CALL linmin(Npar,ilin,par,xit,fret,tolbr,err)
        IF(err)STOP 'fail linmin'
        IF(ABS(fp-fret).GT.del)THEN
           del=ABS(fp-fret)
           ibig=i
        ENDIF
     ENDDO
     WRITE(6,11) ilin,fret,par
     IF(2.d0*ABS(fp-fret).LE.tol*(ABS(fp)+ABS(fret)))THEN

        WRITE(6,*)'powell converged',ABS(fp-fret)/ABS(fret)
        RETURN
     ENDIF
     IF(ilin.GE.itmax)THEN

        WRITE(6,*)'raggiunto maxiter in powell'
        RETURN
     ENDIF
     DO j=1,Npar
        ptt(j)=2.d0*par(j)-pt(j) 
        xit(j)=par(j)-pt(j)
        pt(j)=par(j)
     ENDDO
     CALL distgen(Npar,ptt,fptt)
     IF(fptt.GE.fp) GOTO 1
     t=2.d0*(fp-2.d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
     IF(t.GE.0.d0) GOTO 1
     CALL linmin(Npar,ilin,par,xit,fret,tolbr,err)
     IF(err)THEN
        RETURN
     ENDIF
     DO j=1,Npar
        xi(j,ibig)=xit(j)
     ENDDO
     GOTO 1     
   END SUBROUTINE powell
   
   SUBROUTINE integrate_v(func,a,b,s)
     REAL(KIND=dbl) a,b,func
     REAL(KIND=dbl), DIMENSION(:) :: s
     EXTERNAL :: func     
     
     DO ind_int = 1, Ng
        CALL integrate(func,a,b,s(ind_int))
     END DO

   END SUBROUTINE integrate_v

   SUBROUTINE integrate(func,a,b,s)
     INTEGER, PARAMETER :: JMAX=20
     REAL(KIND=dbl) a,b,func,s
     REAL(KIND=dbl), PARAMETER :: EPS=1.0D-6
     EXTERNAL :: func
     INTEGER :: j
     REAL(KIND=dbl) :: olds
     olds=-1.0E30_dbl
     DO  j=1,JMAX
        CALL trapzd(func,a,b,s,j)
        IF (ABS(s-olds).LT.EPS*ABS(olds)) RETURN
        olds=s
     END DO
     STOP 'too many steps in qtrap'
   END SUBROUTINE integrate

   SUBROUTINE trapzd(func,a,b,s,n)
     INTEGER :: n
     REAL(KIND=dbl) ::  a,b,s,func
     EXTERNAL :: func
     INTEGER :: it,j
     REAL(KIND=dbl) :: del,sum,tnm,x
     IF (n.EQ.1) THEN
        s=0.5*(b-a)*(func(a)+func(b))
     ELSE
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        DO j=1,it
           sum=sum+func(x)
           x=x+del
        END DO
        s=0.5*(s+(b-a)*sum/tnm)
     ENDIF
     RETURN
   END SUBROUTINE trapzd

END MODULE Fit_lsq
