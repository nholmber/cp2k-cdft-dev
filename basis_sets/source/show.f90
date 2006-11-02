!------------------------------------------------------------------------------!
  subroutine print_quickstep_input(cmat,etot)
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  USE atom
  USE xcfcn
  USE pspot
  IMPLICIT NONE
!in:
  REAL(dp)  :: cmat(namax,namax,0:lamax)
  REAL(dp)  :: etot
!locals:
  REAL(dp)  :: nn(namax,namax,0:lamax)
  integer   :: i,j,k,l
  character :: am(0:4),number(0:9)
!------------------------------------------------------------------------------!

  CALL calcnn(nn)

  CALL print_qs_input(cmat,etot)

!MK  CALL print_q96_input(cmat)

  am(0)='S'
  am(1)='P'
  am(2)='D'
  am(3)='F'
  am(4)='G'
  number(0)='0'
  number(1)='1'
  number(2)='2'
  number(3)='3'
  number(4)='4'
  number(5)='5'
  number(6)='6'
  number(7)='7'
  number(8)='8'
  number(9)='9'

  print*
  print*
  write (*,*) 'QUICKSTEP input: '
  write (*,*) ' '
  write (*,'(a)',advance='NO') ' &'//trim(atomname)//'_'
  if (Zeff<10) then
    write (*,'(i1)',advance='NO') nint(Zeff)
  elseif (Zeff<100) then
    write (*,'(i2)',advance='NO') nint(Zeff)
  else
    write (*,'(i2)',advance='NO') nint(Zeff)
  endif
  write (*,'(a)') 'el_'//trim(xcstring)
  write (*,*) '# Info: '
  if (optexp) then
    write (*,*) '# Optimized basis for ',atomname
  else
    write (*,*) '# Not optimized basis for ',atomname
    write (*,*) '# ^^^'
  endif
  write (*,'(A,F5.1,A)') ' # Nuclear charge is ',Zval,' .'
  if (allelectron) then
    write (*,'(A,F6.2,A)') ' # Allelectron calculation with ',&
      Zval,' electrons'
  else
    write (*,'(A,F6.2,A)') ' # '//trim(ppstring)//&
      ' pseudopotential calculation with ',Zeff,' valence electrons'
  endif
  write (*,*) '# XC functional is ',trim(xcstring)
  write (*,'(a,f11.6,a)') ' # Kohn-Sham energy is ',etot,' Hartree'
  write (*,*) '# '
  write (*,*) '# '
  write (*,*) '# '
  write (*,*) '# '

!.pseudopotential parameters:
  write (*,*) '&POTENTIAL'
  write (*,'(i8,i8,a20)') nint(Zval),nint(Zeff),'Z    Zeff' 
  write (*,*) '   GOEDECKER    '//trim(xcstring)
  write (*,'(i8,a20)') pplmax+1,'LMAX' 
  write (*,'(f14.9,a14)') r_loc,'RC' 
  write (*,'(i8)',advance='NO') EXPnr(0)
  do i=1,EXPnr(0)
    write (*,'(f14.9)',advance='NO') PPC(i)
  enddo
  write (*,'(a)',advance='NO') '     #C' 
  do i=1,EXPnr(0)
    write (*,'(a4)',advance='NO') 'C'//number(i) 
  enddo
  write (*,*) ' '
  do l=0,pplmax
    write (*,'(f14.9)',advance='NO') r_proj(l)
    write (*,'(i4)',advance='NO') KBPROJnr(l)
    do i=1,KBPROJnr(l)
      do j=i,KBPROJnr(l)
        write (*,'(f14.9)',advance='NO') KBV(i,j,l)
      enddo
    enddo
    write (*,'(a)',advance='NO') '     R'//am(l)//'  #PRJ  H('//am(l)//')  '
    do i=1,KBPROJnr(l)
      do j=i,KBPROJnr(l)
        write (*,'(a)',advance='NO') number(i)//number(j)//'  '
      enddo
    enddo
    write (*,*) ' '
  enddo
  write (*,*) '&END'
  write (*,*) ' '

!Basis functions:
  write (*,*) '&BASIS'
  write (*,'(a)',advance='NO') '  '//trim(atomname)//'_'
  if (Zeff<10) then
    write (*,'(i1)',advance='NO') nint(Zeff)
  elseif (Zeff<100) then
    write (*,'(i2)',advance='NO') nint(Zeff)
  else
    write (*,'(i2)',advance='NO') nint(Zeff)
  endif
  write (*,'(a)',advance='NO') 'el_'//trim(xcstring)//'_'
  write (*,'(a)',advance='NO') '('
  if (minval(nocc(0:lmax))+3<10) then
    write (*,'(i1)',advance='NO') minval(nocc(0:lmax))+3
  else
    write (*,'(i2)',advance='NO') minval(nocc(0:lmax))+3
  endif
  do l=0,lmax
    write (*,'(a)',advance='NO') am(l)
  enddo
  write (*,'(a)',advance='NO') '/'
  if (nalpha(0)<10) then
    write (*,'(i1)',advance='NO') nalpha(0)
  else
    write (*,'(i2)',advance='NO') nalpha(0)
  endif
  write (*,'(a)',advance='NO') ')_'
  if (add_pot_factor==0.d0) then
    write (*,'(a)',advance='NO') 'f'
  else
    write (*,'(a)',advance='NO') 'c'
  endif
  write (*,'(a)') ' '
  do l=0,lmax
    write (*,'(a8)',advance='NO') 'EXP '//am(l)
    write (*,'(20F18.12)') (alpha(j,l),j=1,nalpha(l))
    do k=1,nocc(l)+1
      write (*,'(a8)',advance='NO') 'C'//am(l)//number(k)//'  '
      write (*,'(20F18.12)') (cmat(j,k,l),j=1,nalpha(l))
!MK   write (*,'(20F18.12)') (cmat(j,k,l)*sqrt(nn(j,j,l)),j=1,nalpha(l))
    enddo
  enddo

!Polarization functions:
  write (*,'(a)',advance='NO') '  '//trim(atomname)//'_'
  if (Zeff<10) then
    write (*,'(i1)',advance='NO') nint(Zeff)
  elseif (Zeff<100) then
    write (*,'(i2)',advance='NO') nint(Zeff)
  else
    write (*,'(i2)',advance='NO') nint(Zeff)
  endif
  write (*,'(a)',advance='NO') 'el_'//trim(xcstring)//'_'
  write (*,'(a)',advance='NO') '['
  if (nocc(lmax)<10) then
    write (*,'(i1)',advance='NO') nocc(lmax)
  else
    write (*,'(i2)',advance='NO') nocc(lmax)
  endif
  write (*,'(a)',advance='NO') am(lmax+1)
  write (*,'(a)',advance='NO') '/'
  if (nalpha(lmax)<10) then
    write (*,'(i1)',advance='NO') nalpha(lmax)
  else
    write (*,'(i2)',advance='NO') nalpha(lmax)
  endif
  write (*,'(a)',advance='NO') ']_'
  if (add_pot_factor==0.d0) then
    write (*,'(a)',advance='NO') 'f'
  else
    write (*,'(a)',advance='NO') 'c'
  endif
  write (*,'(a)') ' '
  write (*,'(a8)',advance='NO') 'EXP  '
  write (*,'(20F18.12)') (alpha(j,lmax),j=1,nalpha(lmax))
  do k=nocc(lmax)+1,nocc(lmax)+1
    write (*,'(a8)',advance='NO') 'C'//am(lmax+1)//number(1)//'  '
    write (*,'(20F18.12)') (cmat(j,k,lmax),j=1,nalpha(lmax))
!   write (*,'(20F18.12)') (cmat(j,k,lmax)*sqrt(nn(j,j,lmax)&
!                           *alpha(j,lmax)),j=1,nalpha(lmax))
  enddo

  write (*,*) '&END'
  write (*,*) ' '



end subroutine print_quickstep_input




!------------------------------------------------------------------------------!
  subroutine showresults(cmat,ksener,etot)
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  use atom
  use xcfcn
  implicit none
!in:
  REAL(dp)   :: ksener(namax,0:lamax)
  REAL(dp)   :: cmat(namax,namax,0:lamax)
  REAL(dp)   :: etot
!locals:
  integer    :: i,j,k,l
  character  :: am(0:3)
!------------------------------------------------------------------------------!
  am(0)='S'
  am(1)='P'
  am(2)='D'
  am(3)='F'
  write (*,*) ' '
  write (*,*) ' '
  write (*,*) ' '
  write (*,*) 'Final results:'
  write (*,*) ' '
  write (*,'(A,F20.14)') '                Total energy=',etot
  write (*,*) ' '
  write (*,'(A,F20.14)') '    Kohn-Sham EV:'
  write (*,*) ('     l=',i,'     ',i=0,lmax)
  do i=1,nalpha(lmax)
    write (*,'(4F20.14)') (ksener(i,l),l=0,lmax)
  enddo
  write (*,*) ' '

  do l=0,lmax
    write (*,*) 'Exponents for l=',l,':'
    write (*,*) ' '
    write (*,'(20F20.14)') (alpha(j,l),j=1,nalpha(l))
    write (*,*) ' '
    write (*,*) 'Coefficients:'
    write (*,*) ' '
    do k=1,nocc(l)
      write (*,'(20F20.14)') (cmat(j,k,l),j=1,nalpha(l))
    enddo
    write (*,*) ' '
    write (*,*) 'Gradients for l=',l,':'
    write (*,*) ' '
    write (*,'(20F20.14)') (galpha(j,l),j=1,nalpha(l))
    write (*,*) ' '
    write (*,*) ' '
  enddo
  write (*,*) ' '



end subroutine showresults



!------------------------------------------------------------------------------!
  subroutine showorbitals(cmat)
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  use atom
  use rint
  implicit none
!in:
  REAL(dp)   :: amin(0:lamax)
  REAL(dp)   :: cmat(namax,namax,0:lamax)
  REAL(dp)   :: orbit(namax),r,rmax(0:lamax),step(0:lamax)
  integer    :: i,j,l
  REAL(dp), PARAMETER :: PI=3.14159265358979323846264_dp
  REAL(dp):: psi
!------------------------------------------------------------------------------!
!..determine smallest exponent
  do l=0,lmax
    amin(l)=alpha(1,l)
    do i=1,Nalpha(l)
      if (alpha(i,l).lt.amin(l)) then
        amin(l)=alpha(i,l)
      endif
    enddo
    rmax(l)=sqrt(10.D0/amin(l))
    step(l)=rmax(l)/300.D0
  enddo

  OPEN (12,FILE='atom.s1')
  r = 0.0D0
  DO i=1,300
    psi = 0.0D0
    DO j=1,nalpha(0)
      psi = psi + cmat(j,1,0)*exp(-alpha(j,0)*r*r)*&
                  (2**7/PI)**0.25D0*alpha(j,0)**0.75D0
    END DO
    WRITE (12,"(2E20.12)") r,psi
    r = r + step(0)
  END DO
  CLOSE (12)
  WRITE (*,*) 'created file "atom.s1".'

  open(12,FILE='atom.sorbits')
  do r=0,rmax(0),step(0)
    do i=1,nalpha(0)
      orbit(i)=0.D0
      do j=1,nalpha(0)
        orbit(i)=orbit(i)+cmat(j,i,0)*exp(-alpha(j,0)*r*r)*&
          (2**7/PI)**.25D0*alpha(j,0)**.75D0
      enddo
    enddo
    write(12,"(30E20.12)") r,orbit
  enddo
  close(12)
  write (*,*) 'created file "atom.sorbits".'

  open(12,FILE='atom.chi_1s')
  do r=0,rmax(0),step(0)
    do i=1,nalpha(0)
      orbit(i)=0.D0
      do j=i,nalpha(0)
        orbit(i)=orbit(i)+cmat(j,1,0)*exp(-alpha(j,0)*r*r)*&
          (2**7/PI)**.25D0*alpha(j,0)**.75D0
      enddo
    enddo
    write(12,"(30E20.12)") r,orbit
  enddo
  close(12)
  write (*,*) 'created file "atom.chi_s".'

  open(12,FILE='atom.chi_2s')
  do r=0,rmax(0),step(0)
    do i=1,nalpha(0)
      orbit(i)=0.D0
      do j=i,nalpha(0)
        orbit(i)=orbit(i)+cmat(j,2,0)*exp(-alpha(j,0)*r*r)*&
          (2**7/PI)**.25D0*alpha(j,0)**.75D0
      enddo
    enddo
    write(12,"(30E20.12)") r,orbit
  enddo
  close(12)
  write (*,*) 'created file "atom.chi_s".'

  open(12,FILE='atom.chi_3s')
  do r=0,rmax(0),step(0)
    do i=1,nalpha(0)
      orbit(i)=0.D0
      do j=i,nalpha(0)
        orbit(i)=orbit(i)+cmat(j,3,0)*exp(-alpha(j,0)*r*r)*&
        (2**7/PI)**.25D0*alpha(j,0)**.75D0
      enddo
    enddo
    write(12,"(30E20.12)") r,orbit
  enddo
  close(12)
  write (*,*) 'created file "atom.chi_s".'

  if (lmax.ge.1) then
    open(12,FILE='atom.porbits')
    do r=0,rmax(1),step(1)
      do i=1,nalpha(1)
        orbit(i)=0.D0
        do j=1,nalpha(1)
          orbit(i)=orbit(i)+cmat(j,i,1)*r*exp(-alpha(j,1)*r*r)*&
            (2**11/PI/9)**.25D0*alpha(j,0)**1.25D0
        enddo
      enddo
      write(12,"(30E20.12)") r,orbit
    enddo
    close(12)
    write (*,*) 'created file "atom.porbits".'
  endif

  open(12,FILE='atom.chi_p')
  do r=0,rmax(0),step(0)
    do i=1,nalpha(0)
      orbit(i)=0.D0
      do j=i,nalpha(0)
        orbit(i)=orbit(i)+cmat(j,3,1)*r*exp(-alpha(j,1)*r*r)*&
          (2**11/PI/9)**.25D0*alpha(j,0)**1.25D0
      enddo
    enddo
    write(12,"(30E20.12)") r,orbit
  enddo
  close(12)
  write (*,*) 'created file "atom.chi_p".'

  if (lmax.ge.2) then
    open(12,FILE='atom.dorbits')
    do r=0,rmax(2),step(2)
      do i=1,nalpha(2)
        orbit(i)=0.D0
        do j=1,nalpha(2)
          orbit(i)=orbit(i)+cmat(j,i,2)*r*r*exp(-alpha(j,1)*r*r)*&
            (2**15/PI/225)**.25D0*alpha(j,0)**1.75D0
        enddo
      enddo
      write(12,"(30E20.12)") r,orbit
    enddo
    close(12)
    write (*,*) 'created file "atom.dorbits".'
  endif

end subroutine showorbitals



!------------------------------------------------------------------------------!
  subroutine showdensities(pmat,cnn)
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  use atom
  use rint
  implicit none
!in:
  REAL(dp)   :: pmat(*),cnn(*)
!locals:
  integer    :: i
  REAL(dp)   :: rho(ippn),drho(ippn),ddrho(ippn)
  REAL(dp)   :: pi,rhoi,rhotot
!------------------------------------------------------------------------------!
  call charge_density(pmat,cnn,ippn,xip,rho,drho,ddrho)
  pi = 3.14159265358979323846264_dp
  open (12,FILE='atom.densities')
  rhotot= 0.0_dp
  do i=1,ippn
    rhoi = 4.0d0*pi*xip(i)**2*wi(i)*rho(i)
    write (12,"(3ES25.12)") xip(i),rho(i),rhoi
    rhotot= rhotot + rhoi
  enddo
  close(12)
  write (*,*) 'rhotot =',rhotot
  write (*,*) 'created file "atom.densities".'

end subroutine showdensities



!------------------------------------------------------------------------------!
  subroutine showpseudopot
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  use pspot
  implicit none
!locals:
  integer   :: i,l
  REAL(dp)  :: Vps(0:PPlmax),g(0:PPlmax),amin,v,r,rmax,step,erf
!------------------------------------------------------------------------------!
!..Determine smallest exponent
  amin=PPNexp(1,PPlmax)
  do l=0,PPlmax
    do i=1,EXPnr(l)
      if (PPNexp(i,l).lt.amin) then
        amin=PPNexp(i,l)
      endif
    enddo
  enddo

  rmax=sqrt(10.D0/amin)
  step=rmax/100.D0

  open (12,FILE='atom.pspot')
  
  do r=1.D-10,rmax,step
    v=0.D0
    do i=1,ERFnr
      v=v-Zeff/r* PPLerf(i)*ERF(PPNerf(i)*r)
    enddo
    do l=0,PPlmax
      g(l)=0.D0
      do i=1,EXPnr(l)
        g(l)=g(l)+PPLexp(i,l)*r**PPRexp(i,l)*exp(-PPNexp(i,l)*r*r)
      enddo
      Vps(l)=v+g(l)
    enddo
    write (12,*) r,Vps,v,g,0
  enddo

  close(12)
  write (*,*) 'created file "atom.pspot".'

end subroutine showpseudopot



SUBROUTINE print_qs_input(cmat,etot)

  USE basic_data_types, ONLY: dp
  USE atom
  USE xcfcn
  USE pspot

  IMPLICIT NONE

  REAL(dp):: etot
!  REAL(dp):: cmat(namax,namax,0:lamax),nn(namax,namax,0:lamax),zero(namax)
  REAL(dp):: cmat(namax,namax,0:lamax),zero(namax)

  CHARACTER(LEN=240) :: string

  CHARACTER :: am(0:4),number(0:9),fmt*12
  INTEGER   :: i,j,l,llmax

! -----------------------------------------------------------------------------

  zero(:) = 0.0D0

  am(0) = "S"
  am(1) = "P"
  am(2) = "D"
  am(3) = "F"
  am(4) = "G"

  number(0)="0"
  number(1)="1"
  number(2)="2"
  number(3)="3"
  number(4)="4"
  number(5)="5"
  number(6)="6"
  number(7)="7"
  number(8)="8"
  number(9)="9"

  IF (allelectron) THEN
    llmax = 0
  ELSE
    llmax = pplmax + 1
  END IF

  OPEN (12,FILE="qs.input")

!MK  WRITE (12,"(A)",ADVANCE="NO") "&"//TRIM(atomname)//"_"
!MK  IF (Zeff < 10) THEN
!MK    WRITE (12,"(I1)",ADVANCE="NO") NINT(Zeff)
!MK  ELSE
!MK    WRITE (12,"(I2)",ADVANCE="NO") NINT(Zeff)
!MK  END IF
!MK  WRITE (12,"(A)") "el_"//TRIM(xcstring)

  string = ""
  WRITE (UNIT=string,FMT="(A,F6.1)")&
    "# Nuclear charge: ",zval
  CALL compress(string,.FALSE.)
  WRITE (UNIT=12,FMT="(A)") TRIM(string)

  string = ""
  WRITE (UNIT=string,FMT="(A,F6.1)")&
    "# Number of valence electrons: ",zeff
  CALL compress(string,.FALSE.)
  WRITE (UNIT=12,FMT="(A)") TRIM(string)

  WRITE (UNIT=12,FMT="(A)")&
    "# XC functional: "//TRIM(xcstring)

  string = ""
  WRITE (UNIT=string,FMT="(A,I3,A,F15.6,A)")&
    "# Kohn-Sham energy for ",nalpha(0)," Gaussians: ",etot," a.u."
  CALL compress(string,.FALSE.)
  WRITE (UNIT=12,FMT="(A)") TRIM(string)

  IF (ABS(add_pot_factor) > 1.0E-20_dp) THEN
    string = ""
    WRITE (UNIT=string,FMT="(A,F10.3)")&
      "# r_cov: ",SQRT(SQRT(0.125d0/add_pot_factor))
    CALL compress(string,.FALSE.)
    WRITE (UNIT=12,FMT="(A)") TRIM(string)
    WRITE (UNIT=string,FMT="(A,F10.3)")&
      "# r_prb: ",SQRT(SQRT(0.5d0/add_pot_factor))
    CALL compress(string,.FALSE.)
    WRITE (UNIT=12,FMT="(A)") TRIM(string)
  END IF

  IF (.NOT.allelectron) THEN
    WRITE (12,"(/,A)") "&POTENTIAL"
    WRITE (12,"(3I5)") NINT(zval),NINT(zeff),llmax
!   IF (allelectron) THEN
!     WRITE (12,"(T2,A)") "Allelectron"
!   ELSE
!     WRITE (12,"(T2,A)") "Goedecker pseudopotential for "//TRIM(xcstring)
!   END IF
!   WRITE (12,"(I5)") llmax
    WRITE (12,"(F15.8,I5,4F15.8)") r_loc,expnr(0),(ppc(i),i=1,expnr(0))
    DO l=0,llmax-1
      WRITE (12,"(F15.8,I5,4F15.8)")&
        r_proj(l),kbprojnr(l),(kbv(1,j,l),j=1,kbprojnr(l))
      fmt = "(T  ,4F15.8)"
      DO i=2,kbprojnr(l)
        WRITE (fmt(3:4),"(I2)") 15*i + 6
        WRITE (12,fmt) (kbv(i,j,l),j=i,kbprojnr(l))
      END DO
    END DO
    WRITE (12,"(A)") "&END"
  END IF

  CALL print_q96_input(cmat)

!MK  WRITE (12,"(/,A)") "&BASIS"
!MK  WRITE (12,"(A)",ADVANCE="NO") TRIM(atomname)//"_"
!MK  IF (zeff < 10) THEN
!MK    WRITE (12,"(I1)",ADVANCE="NO") NINT(zeff)
!MK  ELSE
!MK    WRITE (12,"(I2)",ADVANCE="NO") NINT(zeff)
!MK  END IF
!MK  WRITE (12,"(A)",ADVANCE="NO") "el_"//TRIM(xcstring)//"_"
!MK  WRITE (12,"(A)",ADVANCE="NO") "("
!MK  IF (MAXVAL(nocc(0:lmax)) + 2 < 10) then
!MK    WRITE (12,"(I1)",ADVANCE="NO") MAXVAL(nocc(0:lmax))
!MK  ELSE
!MK    WRITE (12,"(I2)",ADVANCE="NO") MAXVAL(nocc(0:lmax))
!MK  END IF
!MK  DO l=0,lmax
!MK    WRITE (12,"(A)",ADVANCE="NO") am(l)
!MK  END DO
!MK  WRITE (12,"(A)",ADVANCE="NO") "/"
!MK  IF (nalpha(0) < 10) THEN
!MK    WRITE (12,"(I1)",ADVANCE="NO") nalpha(0)
!MK  ELSE
!MK    WRITE (12,"(I2)",ADVANCE="NO") nalpha(0)
!MK  END IF
!MK  WRITE (12,"(A)",ADVANCE="NO") ")_"
!MK  IF (add_pot_factor == 0.0D0) THEN
!MK    WRITE (12,"(A)") "f"
!MK  ELSE
!MK    WRITE (12,"(A)") "c"
!MK  END IF
!MK  DO l=0,lmax
!MK    WRITE (12,"(T3,A5,(T9,4F18.8))") "EXP "//am(l),(alpha(j,l),j=1,nalpha(l))
!MK    DO k=1,nocc(l)
!MK      WRITE (12,"(T3,A5,(T9,4F18.8))")&
!MK        "C"//am(l)//number(k),(zero(j),j=nalpha(l)+1,nalpha(0)),&
!MK        (cmat(j,k,l),j=1,nalpha(l))
!MK!     WRITE (12,"(T3,A5,(T9,4F18.8))")&
!MK!       "C"//am(l)//number(k),(cmat(j,k,l)*SQRT(nn(j,j,l)),j=1,nalpha(l))
!MK    END DO
!MK    WRITE (12,"(T3,A5,(T9,4F18.8))")&
!MK      "C"//am(l)//number(k),(zero(j),j=1,nalpha(0))
!MK  END DO
!MK
!MK  WRITE (12,"(A)") "&END"

  CLOSE (12)

END SUBROUTINE print_qs_input



SUBROUTINE print_q96_input(cmat)

  USE basic_data_types, ONLY: dp
  USE atom
  USE xcfcn
  USE pspot

  IMPLICIT NONE

  REAL(dp):: cmat(namax,namax,0:lamax),zero(namax),dmat(namax,namax,0:lamax)

  INTEGER   :: i,d,l,n,z

! -----------------------------------------------------------------------------

  zero(:) = 0.0D0
  dmat(:,:,:) = 0.0D0

!MK  OPEN (12,FILE="q96.input")

  WRITE (12,"(A,/,I3)") TRIM(atomname),SUM(nocc(0:lmax))

  z = 0

  DO l=0,lmax
!MK
    WRITE (12,"(5I3)") z,l,l,nalpha(l),nocc(l)
    DO i=1,nalpha(l)
      WRITE (12,"(F20.8,20F15.8)")&
        alpha(i,l),(cmat(i,n,l),n=1,nocc(l)),&
        (0.0D0,n=1,nalpha(l)-3)
    END DO
!MK
    d = nalpha(0) - nalpha(l)
    DO i=1,nalpha(l)
      dmat(i+d,:,l) = cmat(i,:,l)
    END DO
  END DO

  WRITE (12,"(A,/,I3)") TRIM(atomname),SUM(nocc(0:lmax))
  WRITE (12,"(5I3)") z,z,z,z,z

  DO i=1,nalpha(0)
    WRITE (12,"(F20.8,20F15.8)") alpha(i,0),&
                                 ((dmat(i,n,l),n=1,nocc(l)+1),l=0,lmax)
!                                ((dmat(i,n,l),l=0,lmax),n=1,3)
  END DO

!MK  CLOSE (12)

END SUBROUTINE print_q96_input

! *****************************************************************************

  SUBROUTINE compress(string,full)

    ! Purpose: Eliminate multiple space characters in a string.
    !          If full is .TRUE., then all spaces are eliminated.

    ! History: - Creation (23.06.1998,MK)

    ! *************************************************************************

    CHARACTER(LEN=*), INTENT(INOUT)          :: string
    LOGICAL, INTENT(IN)                      :: full

    INTEGER                                  :: i, z
    LOGICAL                                  :: remove_all

    ! -------------------------------------------------------------------------

    !IF (PRESENT(full)) THEN
      remove_all = full
    !ELSE
    !  remove_all = .FALSE.
    !END IF

    z = 1

    DO i=1,LEN_TRIM(string)
      IF ((z == 1).OR.remove_all) THEN
        IF (string(i:i) /= " ") THEN
          string(z:z) = string(i:i)
          z = z + 1
        END IF
      ELSE
        IF ((string(i:i) /= " ").OR.(string(z-1:z-1) /= " ")) THEN
          string(z:z) = string(i:i)
          z = z + 1
        END IF
      END IF
    END DO

    string(z:) = ""

  END SUBROUTINE compress
