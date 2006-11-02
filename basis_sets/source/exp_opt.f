C     ==================================================================
        subroutine exp_opt(cmat,ksener,etot)
C     ==================================================================
        use atom
        use upd
        implicit none
        real*8 cmat(namax,namax,0:lamax)
        real*8 ksener(namax,0:lamax)
        real*8 dalpha(namax,0:lamax)
        real*8 etoth,etot,sumg,epsscforig,epslimit
        logical firstloop,check,converged,normal
        integer i,j,k,l,iter,iiter,idis,step,who
C     ------------------------------------------------------------------
c
c       call test_grad(cmat)
c       call test_grad_direct(cmat)
        firstloop=.false.
	epsscforig=epsscf
        epslimit=epsgrad
	call exp_grad (cmat,ksener)
	call alpha2beta
	gbetaold=gbeta
	betaold=beta
        normal=.true.
c
C     ------------------------------------------------------------------
c	Optimization loop
c
        step=8
	do iiter=0,10000
	  call hessunit
	  do iter=1,step
c
c..check gradient convergence
	    sumg=0.D0
	    do k=1,nbeta
	      sumg=sumg+gbeta(k)**2
	    enddo
            sumg=sqrt(sumg)
	    if (sumg.lt.epsgrad) then
	      write  (*,*) 'gradient converged !'
	      write  (*,*) 'checking Hessian:'
	      call matcheck(hessinv,namax*(lamax+1),nbeta,check)
	      goto 101
	      write  (*,*) 'Hessian not positive definite.'
	    endif
c
	    write(*,*) step*iiter+iter,'New direction:'
            if (normal) then
              call updalpha_eps(epslimit,converged,who,dalpha)
              print*,' ',who
              if (converged) then
                print*,'converged  epslimit=',epslimit
                epslimit=epslimit*1.d-1
              else
c
c...........minimization step
                etoth=etot
	        call linesearch(cmat,dalpha,etot)
	        epsscf=max(epsscforig,abs(etoth-etot)/1.D2)
c...........try diis step
c               idis=20*iiter+iter
c               call expdiis(cmat,etot,idis)
              endif
            else
              call updalpha_allin1(normal)
              call wfn_ortho(cmat)
              call lda_scf(cmat,firstloop,ksener,etot)
              call exp_grad (cmat,ksener)
            endif
c
c           call test_grad(cmat)
c           call test_grad_direct(cmat)
c
c
c..update and check Hessian
	    call updhess
	    call matcheck(hessinv,namax*(lamax+1),nbeta,check)
	    if (.not.check) then
	      write(*,*) 'Hessian not positive definite !! ...',
     &			'corrected Hessian.'
	    endif
c
C     ------------------------------------------------------------------
c..output
c	    write (*,*) '                              ',
c     &			'Final total energy:  E=',etot
	    epsscf=max(epsscforig,abs(etoth-etot)/1.D2)
c	    write (*,*) ' '
c
	    open (33,FILE='atom.exp',position="REWIND")
	    write  (33,*) 'alpha:  '
	    do k=0,lmax
	      write (33,'(4F20.10)') (alpha(i,k),i=1,nalpha(k))
	    enddo
	    write  (33,*) 'galpha:  '
	    do k=0,lmax
	      write (33,'(4F20.10)') (galpha(i,k),i=1,nalpha(k))
	    enddo
	    write (33,*) 'Maximum gradient component:',sumg
	    write (33,*) ' ' 
	    write (33,*) step*iiter+iter,
     &                   '                        Etot=',etot
	    write (33,*) '  '
	    close(33)
c
c
C     ------------------------------------------------------------------
c..output
	  write  (*,*) '  alpha:  '
	  do k=0,lmax
	    write (*,'(4F20.10)') (alpha(i,k),i=1,nalpha(k))
	  enddo
c	  write  (*,*) 'galpha:  '
c	  do k=0,lmax
c	    write (*,'(4F20.10)') (galpha(i,k),i=1,nalpha(k))
c	  enddo
          write(*,*) '  gradient modulus=',sumg
	  write (*,*) ' ' 
	  write (*,*) step*iiter+iter,
     &                '                        Etot=',etot
	  write  (*,*) '  '
c
	  enddo
	enddo
C     ------------------------------------------------------------------
c
101	continue
c..output
	    write (*,*) '                              ',
     &			'Final total energy:  E=',etot
	    epsscf=max(epsscforig,abs(etoth-etot)/1.0D2)
	    write (*,*) ' '
c
	    open (33,FILE='atom.exp',position="REWIND")
	    write  (33,*) 'alpha:  '
	    do k=0,lmax
	      write (33,'(4F20.10)') (alpha(i,k),i=1,nalpha(k))
	    enddo
	    write  (33,*) 'galpha:  '
	    do k=0,lmax
	      write (33,'(4F20.10)') (galpha(i,k),i=1,nalpha(k))
	    enddo
	    write (33,*) 'Maximum gradient component:',sumg
	    write (33,*) ' ' 
	    write (33,*) step*iiter+iter,'Etot=',etot
	    write  (33,*) '  '
	    close(33)
c
	epsscf=epsscforig
c
	return
	end
