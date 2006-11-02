      subroutine calc_energy(cmat,etot)

      USE basic_data_types, ONLY: dp
      USE atom
      USE rint
      USE pspot
      USE energies
      IMPLICIT NONE
      REAL(dp) cmat(namax,namax,0:lamax)
!..locals
      REAL(dp) smat(namax,namax,0:lamax)
      REAL(dp) pmat(namax,namax,0:lamax)
      REAL(dp) tmat(namax,namax,0:lamax)
      REAL(dp) umat(namax,namax,0:lamax)
      REAL(dp) uaddmat(namax,namax,0:lamax)
      REAL(dp) h1mat(namax,namax,0:lamax)
      REAL(dp) h2mat(namax,namax,0:lamax)
      REAL(dp) vxcmat(namax,namax,0:lamax)
      REAL(dp) excmat(namax,namax,0:lamax)
      REAL(dp) ksmat(namax,namax,0:lamax)
      REAL(dp) cnn(namax,namax,0:lamax)
      REAL(dp) eone,etot,trace
      REAL(dp) f,cost,sint,sint2,t,w,x,pi
      INTEGER k,l,i,j,n

!     ------------------------------------------------------------------

!..determine smallest and biggest exponent
        alphamin=alpha(1,lmax)
        alphamax=alpha(1,lmax)
        do l=0,lmax
          do i=1,Nalpha(l)
            if (alpha(i,l).lt.alphamin) alphamin=alpha(i,l)
            if (alpha(i,l).gt.alphamax) alphamax=alpha(i,l)
          enddo
        enddo
!..calculate ippn basis points for interpolation
!      c=1.015D0
!      al=dlog(c)
!      ippn1=dlog(60800.0d0*zeff)/al
!      ippn=max(ippn,ippn1)
!      xip(1)=.000625d0/max(1.d0,zeff)
!      do i=2,ippn
!       xip(i)=xip(i-1)*c
!      enddo
!MK
!     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
!     *** u [-1,+1] -> r [0,infinity] => r = ln(2/(1 - u))/ln(2)            ***

      n = ippn
      pi = 3.14159265358979323846264D0
      f = pi/REAL(n + 1,dp)

      DO i=1,n
        t = REAL(i,dp)*f
        cost = DCOS(t)
        sint = DSIN(t)
        sint2 = sint**2
        x = REAL(2*i - n - 1,dp)/REAL(n + 1,dp) -&
            2.0D0*(1.0D0 + 2.0D0*sint2/3.0D0)*cost*sint/pi
        w = 16.0D0*sint2**2/REAL(3*(n + 1),dp)
        ri(i) = DLOG(2.0D0/(1.0D0 - x))/DLOG(2.0D0)
        wi(i) = w/(DLOG(2.0d0)*(1.0D0 - x))
      END DO

      DO i=1,n
!MK        write(*,"(I6,3F20.12)") i,xip(i),ri(i),wi(i)
        xip(i) = ri(i)
      END DO
!MK
!
!..Normalization factors
      call calcnn(cnn)
!..Overlap matrix
      call overlap(smat)
!..Kinetic energy matrix
      call kinetic_int(tmat,smat)
!..Potential energy matrix
      if (allelectron) then
        call nuc_att(umat,smat)
        do k=0,lmax
          do i=1,nalpha(k)
            do j=1,nalpha(k)
                umat(i,j,k)=-zval*umat(i,j,k)
             enddo
          enddo
        enddo
      else
        call pseudopot(umat,cnn)
      endif
!..additional parabolic potential
      call v_add(smat,uaddmat)
!..one electron part
        do k=0,lmax
          do i=1,nalpha(k)
            do j=1,nalpha(k)
              h1mat(i,j,k)=tmat(i,j,k)+umat(i,j,k) +&
                           add_pot_factor*uaddmat(i,j,k)
            enddo
          enddo
        enddo
!..Density matrix
      call denmat(pmat,cmat)
!
!..calculate Kohn-Sham matrix
        call ks_matrix(ksmat,smat,pmat,h1mat,h2mat,vxcmat,excmat,cnn)
!
!
!..energies
        ekin=trace(pmat,tmat)
        epot=trace(pmat,umat)
        epotadd=add_pot_factor*trace(pmat,uaddmat)
        eone=trace(pmat,h1mat)
        ecb=0.5D0*trace(pmat,h2mat)
        exc=trace(pmat,excmat)
        etot=eone+ecb+exc
!
      END
