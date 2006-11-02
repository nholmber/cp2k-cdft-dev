      subroutine nuc_attd(s,sd,ud)
      use atom
      implicit none
      real*8 ud(namax,namax,0:lamax)
      real*8 s(namax,namax,0:lamax)
      real*8 sd(namax,namax,0:lamax)
      real*8 tpi,facul,dblfac,p,q,pi
      integer i,j,l,k
      pi=3.14159265358979323846264D0
      tpi=2.0D0*pi
      do l=0,lmax
        do k=1,nalpha(l)
          p=alpha(k,l)
          do j=1,nalpha(l)
            q=alpha(j,l)
            ud(k,j,l)=facul(l)/dblfac(2*l+1)
     &			*2**(l+1.5D0)
     &			*sqrt((p+q)/tpi)
     &        		* (0.5D0/(p+q)*s(k,j,l)+sd(k,j,l))
          enddo
        enddo
      enddo
c
      return
      end
c
