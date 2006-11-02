!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  program atomcode
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  use atom
  use rint
  use pspot
  use xcfcn
  implicit none
!.Matrix of coefficients for the wavefunctions
  REAL(dp) :: cmat(namax,namax,0:lamax)
  REAL(dp) :: pmat(namax,namax,0:lamax)
  REAL(dp) :: ksener(namax,0:lamax)
  REAL(dp) :: cnn(namax,namax,0:lamax)
  REAL(dp) :: Ne,etot
  integer :: l,i
  logical :: firstloop
!------------------------------------------------------------------------------!
  write (*,*) ' '
  write (*,*) '       PROGRAM ATOM'
  write (*,*) '       ------------'
  write (*,*) 'reading input file...'

  call input

  write (*,*) 'finished reading input file.'
  write (*,*) ' '
  write (*,*) ' '
  if (optexp) then
    write (*,*) 'Basis Optimization for atom ',atomname
  else
    write (*,*) 'Total energy calculation for atom ',atomname
  endif
  write (*,*) ' '
  Ne=0.D0
  do l=0,Lmax
    do i=1,Nocc(l)
      Ne=Ne+(2*l+1)*occ(i,l)
    enddo
  enddo
  if (allelectron) then
    write (*,'(A,F6.2,A)') ' Allelectron calculation with ',&
      Ne,' electrons.'
    write (*,'(A,F5.1,A)') ' Nuclear charge is ',Zval,' .'
  else
    write (*,'(A,F6.2,A)') ' Pseudopotential calculation with ',Ne,&
      ' valence electrons.'
    call ppinit
    print*,  trim(ppstring)//' pseudopotential' 
    write (*,'(A,F5.1,A)') ' Effective core charge is ',Zeff,' .'
  endif
  write (*,*) ' '

  write (*,*) 'The exchange correlation functional is ',trim(xcstring),' .'
  write (*,*) ' '
  if (add_pot_factor/=0.d0) then
    write (*,*) 'An additinal parabolic potential corresponding'
    write (*,'(A,F10.4)') ' to a covalent radius of ',&
      0.5d0*(2.d0/add_pot_factor)**0.25d0
    write (*,*) 'Bohr is used.'
    write (*,*) ' '
  endif
  write (*,*) ' '



!------------------------------------------------------------------------------!
! EXPONENT OPTIMIZATION OR TOTAL ENERGY CALCULATION

!..calculate the Kohn-Sham energy
  firstloop=.true.
  call lda_scf(cmat,firstloop,ksener,etot)
  firstloop=.false.

!..optimize exponents
  if (optexp) call exp_opt(cmat,ksener,etot)

!..final total energy
  write (*,*) 'Final total energy is:  ',etot

!..gradient test
!  call test_grad(cmat)
   call test_grad_direct(cmat)

!..derivative of the wavefunctions w.r.t. occupation number of HOAO
!  write (*,*) ' '
!  write (*,*) 'LUAO calculation:'
!  call numd_wfn(cmat)
!  write (*,*) ' '

!..transfer pseudopotential to Kleinman-Bylander form
  if ((.NOT.allelectron).AND.(pptype.ne.5)) call kbppfit(cmat)

!------------------------------------------------------------------------------!
! OUTPUT SECTION

  call showresults(cmat,ksener,etot)
  write (*,*) ' '
  write (*,*) 'created file "atom.output".'

  call showorbitals(cmat)
  if (.not.allelectron) call showpseudopot
  call print_quickstep_input(cmat,etot)

  call calcnn(cnn)
  call denmat(pmat,cmat)
  call showdensities(pmat,cnn)

  write (*,*) ' '
  write (*,*) 'ATOM finished.'

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  end program atomcode
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
