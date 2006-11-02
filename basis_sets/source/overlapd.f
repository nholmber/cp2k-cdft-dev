      subroutine overlapd(s,sd)
      use atom
      implicit none
      real*8 s(namax,namax,0:lamax)
      real*8 sd(namax,namax,0:lamax)
      real*8 sint,p,q
      integer l,k,j
c
      do l=0,lmax
        do k=1,nalpha(l)
          p=alpha(k,l)
          do j=1,nalpha(l)
            q=alpha(j,l)
            sd(k,j,l)=(2*l+3)/4.D0*(p-q)/q/(p+q)*s(k,j,l)
          enddo
        enddo
      enddo
c
      return
      end
