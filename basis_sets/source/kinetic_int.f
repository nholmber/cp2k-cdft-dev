      subroutine kinetic_int(t,s)
      use atom
      implicit none
      real*8 t(namax,namax,0:lamax),s(namax,namax,0:lamax)
      real*8 p,q
      integer l,i,j
c
      do l=0,lmax
        do i=1,nalpha(l)
          p=alpha(i,l)
          do j=i,nalpha(l)
            q=alpha(j,l)
            t(i,j,l)=(2.D0*l+3.D0)*(p*q)/(p+q)*s(i,j,l)
            t(j,i,l)=t(i,j,l)
          enddo
        enddo
      enddo
      return
      end
