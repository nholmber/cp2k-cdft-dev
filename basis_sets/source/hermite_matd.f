      subroutine hermite_matd(smat,pmat,cnn,cpotdmat,vxcdmat)
      use atom
      use rint
      use pspot
      implicit none
      real*8 smat(namax,namax,0:lamax)
      real*8 pmat(namax,namax,0:lamax)
      real*8 cnn(namax,namax,0:lamax)
      real*8 cpotdmat(namax,namax,0:lamax)
      real*8 vxcdmat(namax,namax,0:lamax)
      real*8 rho(ippn),drho(ippn),ddrho(ippn),scrat(4*ippn)
      real*8 cpot(ippn),vxc(ippn),exc(ippn)
      real*8 ddcpot(ippn),ddvxc(ippn)
      real*8 cpoti(ippn),vxci(ippn)
      real*8 p,q,pq,v,r,curv2,ee
      integer l,i,j,num,k,ierr
c
c..calculate potentials and provide spline representation     	   
      call coulomb(ippn,xip,cpot,smat,pmat)
      call charge_density(pmat,cnn,ippn,xip,rho,drho,ddrho)
      call evxc(ippn,xip,rho,drho,ddrho,vxc,exc,scrat)
      do i=1,ippn
	cpot(i)=cpot(i)*xip(i)**2
	vxc(i)=vxc(i)*xip(i)**2
      enddo
c
      do l=0,lmax
        do i=1,nalpha(l)
          p=alpha(i,l)
          do j=1,nalpha(l)
            q=alpha(j,l)
	    do num=1,ippn
             ee=exp(-(p+q)*xip(num)**2)*cnn(i,j,l)*xip(num)**(2*l+2)
              vxci(num)=vxc(num)*ee
              cpoti(num)=cpot(num)*ee
	    enddo
c
c
c            call dptnq(xip,cpoti,ippn,v,ee)
            v = 0.0d0
            do num=1,ippn
              v = v + wi(num)*cpoti(num)
            end do
            cpotdmat(i,j,l)=v
            cpotdmat(j,i,l)=cpotdmat(i,j,l)
c            call dptnq(xip,vxci,ippn,v,ee)
            v = 0.0d0
            do num=1,ippn
              v = v + wi(num)*vxci(num)
            end do
            vxcdmat(i,j,l)=v
            vxcdmat(j,i,l)=vxcdmat(i,j,l)
          enddo
        enddo
      enddo
c
      return
      end
