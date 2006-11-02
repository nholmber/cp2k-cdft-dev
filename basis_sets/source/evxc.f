      subroutine evxc(nr,r,rho,drho,ddrho,vxc,exc,v2)
      use atom
      implicit real*8(a-h,o-z)
      integer nr
      dimension r(nr),rho(nr),drho(nr),ddrho(nr),exc(nr),vxc(nr),
     *          v2(nr,*)
      do k=1,nr
	if (rho(k).gt.1D-30) then
          call xc(rho(k),ex,ec,vx,vc)
          dr2=drho(k)*drho(k)
          call gcxc(rho(k),dr2,sx,sc,v1x,v2x,v1c,v2c)
          exc(k)=ex+ec+(sx+sc)/rho(k)
          vxc(k)=vx+vc+v1x+v1c
          v2(k,1)=r(k)
          v2(k,2)=(v2x+v2c)*drho(k)/r(k)
	else
	  exc(k)=0.D0
          vxc(k)=0.D0
          v2(k,1)=r(k)
          v2(k,2)=0.D0
	endif
      enddo
      vv2=dasum(nr,v2(1,2),1)
      if(vv2.gt.1.d-14) then
        call curv1(nr,v2(1,1),v2(1,2),0.D0,0.D0,3,v2(1,3),v2(1,4),
     *             0.D0,ierr)
        do k=1,nr
          v2(k,4)=curvd(v2(k,1),nr,v2(1,1),v2(1,2),v2(1,3),0.d0)
        enddo
        do k=1,nr
          vxc(k)=vxc(k)-3.d0*v2(k,2)-r(k)*v2(k,4)
        enddo
      endif
      return
      end
