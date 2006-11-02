      function trace(a,b)
      use atom
      implicit none
      real*8 a(namax,namax,0:lamax),b(namax,namax,0:lamax)
      real*8 trace,t1
      integer l,i,j
c
      trace=0.d0
      do l=0,lmax
        t1=0.d0
        do i=1,nalpha(l)
          do j=1,nalpha(l)
            t1=t1+a(j,i,l)*b(j,i,l)
          enddo
        enddo
        trace=trace+(2.d0*l+1.d0)*t1
      enddo
c
      return
      end
