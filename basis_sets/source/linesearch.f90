        SUBROUTINE linesearch(cmat,dalpha,etot)
        USE basic_data_types, ONLY: dp
        USE atom, only: namax,lamax,lmax,alpha,nalpha
        USE upd,  only: nbeta,xi,gbeta
        IMPLICIT NONE
        LOGICAL firstloop,back
        REAL(dp) :: a1,a2,a3,amin
        REAL(dp) :: E1,E2,E3,etot
        REAL(dp) :: ksener(namax,0:lamax)
        REAL(dp) :: cmatsave(namax,namax,0:lamax)
        REAL(dp) :: cmat(namax,namax,0:lamax)
        REAL(dp) :: alphah(namax,0:lamax)
        REAL(dp) :: dalpha(namax,0:lamax)
        REAL(dp) :: d21,d23,f23,f21,s23,s21,amax,amaxmin,lgl,eold,fac,de
        integer i,l,steps
        common /save/ fac

        fac=min(1.d0,fac)
        if (fac.lt.1.d-5) fac=1.d0
        eold=etot
        steps=0
        alphah=alpha
        a1=0.D0
        E1=etot
        back=.false.
        dalpha=dalpha*fac

        amaxmin=1.D10
        do l=0,lmax
          do i=1,nalpha(l)-1
            amax=-(alpha(i,l)-alpha(i+1,l))/(dalpha(i,l)-dalpha(i+1,l))
            if ((amax.gt.0.D0).and.(amax.lt.amaxmin)) amaxmin=amax
          enddo
          amax=-alpha(nalpha(l),l)/dalpha(nalpha(l),l)
          if ((amax.gt.0.D0).and.(amax.lt.amaxmin)) amaxmin=amax
        enddo
        amaxmin=amaxmin
!       write (*,*) 'amax=',amaxmin

        a2=min(1.D0,amaxmin/2.D0)

        firstloop=.false.
!       write (*,*) 'a2=',a2
        alpha=alphah+a2*dalpha
        call wfn_ortho(cmat)
        call lda_scf(cmat,firstloop,ksener,E2)
        cmatsave=cmat

        lgl=0.d0
        xi=a2*xi
        do i=1,nbeta
          lgl=lgl-xi(i)*gbeta(i)
        enddo

        a3=min(lgl/(2*(E2-E1+lgl)),amaxmin/a2)
!       write (*,*) 'a3=',a3*a2
        alpha=alphah+a3*a2*dalpha
        call wfn_ortho(cmat)
        call lda_scf(cmat,firstloop,ksener,E3)

        if (E3.le.E2) then
          etot=E3
        else
          etot=E2
          cmat=cmatsave
          alpha=alphah+a2*dalpha
        endif
        if (etot.gt.eold) then
          write(*,*) 'step not accepted!'
          etot=eold
          alpha=alphah
          call wfn_ortho(cmat)
          fac=fac/10.d0
        else
          de=abs(eold-etot)
!         fac=min(0.2d0,0.5d0/abs(log(de)))
          fac=1.d0
        endif

        call exp_grad (cmat,ksener)

  100   continue

        if (E2.ge.E1) then
!..case1  (E2>=E1)
          steps=steps+1
!..exceptional case E1
          if (steps.gt.5) then
!           write (*,*) 'try backwards !!!'
            steps=0
            a1=0.D0
            E1=etot
            a2=a3
            E2=E3
            alpha=alphah-a3*dalpha
            call wfn_ortho(cmat)
            call lda_scf(cmat,firstloop,ksener,E3)
            if (E3.gt.E1) then
!..caseE1.1  (E3=>E1<=E2)
              a3=-a3
              goto 200
            else
!..caseE1.2  (E3<=E1<=E2)
              back=.true.
              dalpha=-dalpha
              E2=E3
!             write(*,*) 'E1=',E1,' E2=',E2,' a1=',a1,' a2=',a2
              goto 100
            endif
          endif
!..end of exceptional caseE1

          a3=(a1+a2)/2.D0
          alpha=alphah+a3*dalpha
          call wfn_ortho(cmat)
          call lda_scf(cmat,firstloop,ksener,E3)
          if (E3.ge.E1) then
!..case1.1  (E1<=E3<=E2)
            a2=a3
            E2=E3
            goto 100
          else
!..case1.2  (E1=>E3<=E2)
            goto 200
          endif
        else
!..case2  (E1>E2)
          a3=2*a2
          steps=steps+1
!..exceptional case E2
          if (a3.gt.amaxmin) then
            alpha=alphah+a2*dalpha
            dalpha=a2*dalpha
            goto 400
          endif
!..exceptional case E3
          if (steps.gt.3) then
            alpha=alphah+a3*dalpha
            dalpha=a3*dalpha
            goto 400
          endif
!..end of exceptional case E3

          alpha=alphah+a3*dalpha
          call wfn_ortho(cmat)
          call lda_scf(cmat,firstloop,ksener,E3)
          if (E3.lt.E2) then
!..case2.1  (E1=>E2=>E3)
            a1=a2
            E1=E2
            a2=a3
            E2=E3
            goto 100
          else
!..case2.2  (E1=>E3<=E2)
            goto 200
          endif
        endif

!..calculate minimum of parabola through (a1,E1),(a2,E2),(a3,E3)
  200   continue
        d21=a2-a1
        d23=a2-a3
        s21=a2+a1
        s23=a2+a3
        f21=E2-E1
        f23=E2-E3

        amin=(d21*f23*s21-d23*f21*s23)/(2.D0*(d21*f23-f21*d23))
        alpha=alphah+amin*dalpha
        dalpha=amin*dalpha

  400   if (back) amin=-amin
        write(*,*) 'amin=',amin,'    steps=',steps
  401   call wfn_ortho(cmat)
        call lda_scf(cmat,firstloop,ksener,etot)
        call exp_grad (cmat,ksener)

        end
