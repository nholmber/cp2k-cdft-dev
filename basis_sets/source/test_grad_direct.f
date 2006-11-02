C     ==================================================================
      subroutine test_grad_direct(cmat)
C     ==================================================================
      use atom
      use upd
      use energies
      implicit none
      real*8 cmat(namax,namax,0:lamax)
      real*8 ksener(namax,0:lamax)
      real*8 etoth,etot,sumg,epsscforig,delta
      logical firstloop,check
      integer i,j,k,l,iter,iiter,idis,dj,dl
      real*8 dalpha(namax,0:lamax),ldl
      REAL*8 alphasave(namax,0:lamax)
      REAL*8 ng_ekin,ng_epot,ng_epotadd
      REAL*8 ng_exc,ng_ecoul
      REAL*8 ngalpha,ng_constr
      REAL*8 ag_ekin,ag_epot
      REAL*8 ag_exc,ag_ecoul
      REAL*8 agalpha,ag_constr
C     ------------------------------------------------------------------
c
      alphasave=alpha
c
      call lda_scf(cmat,.false.,ksener,etot)
c
      write  (*,*) 'gradient test:'
      write  (*,*) ' '
      write  (*,*) 'energies:    etot =',etot
      write  (*,*) '             esum =',ekin+epot+epotadd+exc+ecb
      write  (*,*) '             ekin =',ekin
      write  (*,*) '             epot =',epot
      write  (*,*) '             epadd=',epotadd
      write  (*,*) '             exc  =',exc
      write  (*,*) '             ecoul=',ecb
      write  (*,*) ' '
c
      call exp_grad (cmat,ksener)
      delta=1.d-5
c
      do l=0,lmax
        do i=1,nalpha(l)

          alpha(i,l)=alpha(i,l)-0.5*delta
          call calc_energy(cmat,etot)
          ngalpha=-etot/delta
          ng_ekin=-ekin/delta
          ng_epot=-epot/delta
          ng_epotadd=-epotadd/delta
          ng_exc=-exc/delta
          ng_ecoul=-ecb/delta
c
          alpha(i,l)=alpha(i,l)+delta
          call calc_energy(cmat,etot)
          ngalpha=ngalpha+etot/delta
          ng_ekin=ng_ekin+ekin/delta
          ng_epotadd=ng_epotadd+epotadd/delta
          ng_epot=ng_epot+epot/delta
          ng_exc=ng_exc+exc/delta
          ng_ecoul=ng_ecoul+ecb/delta
          alpha(i,l)=alpha(i,l)-0.5d0*delta
c
          write  (*,*) 'gradients:'
          write  (*,*) 'l=',l,'  i=',i,' analytical/numerical'
          write  (*,*) ' '
          write  (*,'(a,f12.5,a,f12.5)') '  g_ekin =',
     &                                     g_ekin(i,l),' ',ng_ekin
          write  (*,'(a,f12.5,a,f12.5)') '  g_epot =',
     &                                     g_epot(i,l),' ',ng_epot
          write  (*,'(a,f12.5,a,f12.5)') '  g_epadd=',
     &                                     g_epotadd(i,l),' ',ng_epotadd
          write  (*,'(a,f12.5,a,f12.5)') '  g_exc  =',
     &                                     g_exc(i,l),' ',ng_exc
          write  (*,'(a,f12.5,a,f12.5)') '  g_ecoul=',
     &                                     g_ecoul(i,l),' ',ng_ecoul
          write  (*,'(a,f12.5,a,f12.5)') '  g_const=',
     &                                     g_constr(i,l),' ',0.0d0
          write  (*,'(a,f12.5,a,f12.5)') '  g_etot =',
     &                       galpha(i,l)-g_constr(i,l),' ',ngalpha
          write  (*,*) ' '

        enddo
      enddo
c
      alpha=alphasave
      call calc_energy(cmat,etot)
c	    
      return
      end
