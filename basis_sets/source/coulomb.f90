      SUBROUTINE coulomb(nr,r,cpot,smat,pmat)
      USE basic_data_types, ONLY: dp
      USE atom
      IMPLICIT NONE
      INTEGER nr
      REAL(dp) r(nr),cpot(nr)
      REAL(dp) smat(namax,namax,0:lamax)
      REAL(dp) pmat(namax,namax,0:lamax)
      REAL(dp) gam1,gam2,dblfac
      REAL(dp) p,q,spq,fs,f1,f2,z,pi,control(3)
      INTEGER i,j,l,k
      parameter(pi=3.14159265358979323846264_dp)

      control(1)=0
      control(2)=0
      control(3)=0
      do i=1,nr
        cpot(i)=0.0_dp
      enddo
      do l=0,lmax
        do i=1,nalpha(l)
          p=alpha(i,l)
          do j=1,nalpha(l)
            q=alpha(j,l)
            spq=sqrt(p+q)
            fs=(2*l+1)*smat(i,j,l)*pmat(i,j,l)
            f1=2.**(l+1)/dsqrt(pi)/dblfac(2*l+1)
            f2=f1*spq
            do k=1,nr
              z=spq*r(k)
              cpot(k)=cpot(k)+fs*(1.0_dp/r(k)-f1*gam1(l,z)/r(k)+f2*gam2(l,z))
!             control(1)=fs/r(k)
!             control(2)=-fs*f1*gam1(l,z)/r(k)
!             control(3)=fs*f2*gam2(l,z)
!             write (*,*) fs,f1,f2,r(k)
!             write (*,*) control
!             write (*,*) l,z,gam1(l,z),gam2(l,z)
!             write (*,"(3F20.10)") r(k),gam1(l,z),gam2(l,z)
            enddo
          enddo
        enddo
      enddo
      end
