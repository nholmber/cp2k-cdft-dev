      subroutine kinetic_intd(s,sd,td)
      use atom
      implicit none
      real*8 td(namax,namax,0:lamax)
      real*8 s(namax,namax,0:lamax)
      real*8 sd(namax,namax,0:lamax)
      real*8 p,q
      integer l,k,j
c
      do l=0,lmax
        do k=1,nalpha(l)
          q=alpha(k,l)
          do j=1,nalpha(l)
            p=alpha(j,l)
            td(k,j,l)=(2*l+3)/(p+q)*(q*q/(p+q)*s(k,j,l)+p*q*sd(k,j,l))
          enddo
        enddo
      enddo
c
      return
      end
