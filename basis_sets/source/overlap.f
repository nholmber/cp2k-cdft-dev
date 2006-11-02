      subroutine overlap(s)
      use atom
      implicit none
      real*8 s(namax,namax,0:lamax)
      real*8 sint,p,q
      integer l,i,j
c
      do l=0,lmax
        do i=1,nalpha(l)
          p=alpha(i,l)
          do j=i,nalpha(l)
            q=alpha(j,l)
            s(i,j,l)=sint(l,p,q)
            s(j,i,l)=s(i,j,l)
          enddo
        enddo
      enddo
c
      return
      end
