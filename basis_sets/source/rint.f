      real*8 function frint(nr,r,w,f,l,factor)
      implicit none
      integer nr,l,i
      real*8 r(nr),w(nr),f(nr),factor
      frint=0.0d0
      do i=1,nr
        frint=frint+w(i)*f(i)*r(i)**(2*l+2)
      enddo
	frint=frint*factor
      return
      end
