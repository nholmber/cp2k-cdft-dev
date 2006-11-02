      subroutine hermite_mat(smat,pmat,cnn,h2mat,excmat,vxcmat)
      use atom
      use rint
      use pspot
      implicit none
      real*8 smat(namax,namax,0:lamax)
      real*8 pmat(namax,namax,0:lamax)
      real*8 cnn(namax,namax,0:lamax)
      real*8 h2mat(namax,namax,0:lamax)
      real*8 excmat(namax,namax,0:lamax)
      real*8 vxcmat(namax,namax,0:lamax)
      real*8 rho(ippn),drho(ippn),ddrho(ippn),scrat(4*ippn)
      real*8 cpot(ippn),vxc(ippn),exc(ippn)
      real*8 ddcpot(ippn),ddvxc(ippn),ddexc(ippn)
      real*8 cpoti(ippn),vxci(ippn),exci(ippn)
      real*8 p,q,pq,v,r,curv2,ee,pi
      parameter (PI=3.14159265358979323846264D0)
      integer l,i,j,num,k,ierr
c
c..calculate potentials and provide spline representation     	   
      call coulomb(ippn,xip,cpot,smat,pmat)
      call charge_density(pmat,cnn,ippn,xip,rho,drho,ddrho)
      qgrid = 0.0d0
      do i=1,ippn
        qgrid = qgrid + 4.0d0*pi*xip(i)**2*wi(i)*rho(i)
      end do
      call evxc(ippn,xip,rho,drho,ddrho,vxc,exc,scrat)
c
      do l=0,lmax
        do i=1,nalpha(l)
          p=alpha(i,l)
          do j=i,nalpha(l)
            q=alpha(j,l)
	    do num=1,ippn
             ee=exp(-(p+q)*xip(num)**2)*cnn(i,j,l)*xip(num)**(2*l+2)
              vxci(num)=vxc(num)*ee
              exci(num)=exc(num)*ee
              cpoti(num)=cpot(num)*ee
	    enddo
c
c            call dptnq(xip,cpoti,ippn,v,ee)
            v = 0.0d0
            do num=1,ippn
              v = v + wi(num)*cpoti(num)
            end do
            h2mat(i,j,l)=v
            h2mat(j,i,l)=h2mat(i,j,l)
c            call dptnq(xip,exci,ippn,v,ee)
            v = 0.0d0
            do num=1,ippn
              v = v + wi(num)*exci(num)
            end do
            excmat(i,j,l)=v
            excmat(j,i,l)=excmat(i,j,l)
c            call dptnq(xip,vxci,ippn,v,ee)
            v = 0.0d0
            do num=1,ippn
              v = v + wi(num)*vxci(num)
            end do
            vxcmat(i,j,l)=v
            vxcmat(j,i,l)=vxcmat(i,j,l)
          enddo
        enddo
      enddo
c
      return
      end
