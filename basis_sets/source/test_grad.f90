      SUBROUTINE test_grad(cmat)

      USE basic_data_types, ONLY: dp
      USE atom
      USE upd
      USE energies
      IMPLICIT NONE
      REAL(dp) :: cmat(namax,namax,0:lamax)
      REAL(dp) :: ksener(namax,0:lamax)
      REAL(dp) :: etot,delta
      INTEGER i,l
      REAL(dp) :: dalpha(namax,0:lamax),ldl
      REAL(dp) :: alphasave(namax,0:lamax)
      REAL(dp) :: ngalpha
      REAL(dp) :: agalpha

!     ------------------------------------------------------------------

      call lda_scf(cmat,.false.,ksener,etot)

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

      call exp_grad (cmat,ksener)
      call alpha2beta
      call beta2dalpha(dalpha)
      ldl=sqrt(sum(dalpha**2))

      agalpha=0.d0
      do l=0,lmax
        do i=1,nalpha(l)
            agalpha=agalpha-galpha(i,l)*dalpha(i,l)/ldl
        enddo
      enddo

      delta=1.d-5

      alphasave=alpha
      alpha=alpha+0.5*delta*dalpha/ldl
      call lda_scf(cmat,.false.,ksener,etot)
      ngalpha=-etot/delta

      alpha=alpha-delta*dalpha/ldl
      call lda_scf(cmat,.false.,ksener,etot)
      ngalpha=ngalpha+etot/delta

      write  (*,*) 'analytical gradient:    g_etot =',agalpha
      write  (*,*) 'numerical  gradient:    g_etot =',ngalpha
      write  (*,*) ' '

      alpha=alphasave

      END
