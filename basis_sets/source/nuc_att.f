      subroutine nuc_att(u,s)
      use atom
      implicit none
      real*8 u(namax,namax,0:lamax)
      real*8 s(namax,namax,0:lamax)
      real*8 tpi,facul,dblfac,p,q,pi
      integer i,j,l
      pi=3.14159265358979323846264D0
      tpi=2.0D0*pi
      do l=0,lmax
        do i=1,nalpha(l)
          p=alpha(i,l)
          do j=i,nalpha(l)
            q=alpha(j,l)
            u(i,j,l)=facul(l)/dblfac(2*l+1)*sqrt((p+q)/tpi)*
     *        2.d0**(l+1.5D0)*s(i,j,l)
            u(j,i,l)=u(i,j,l)
          enddo
        enddo
      enddo
c
      return
      end
c
