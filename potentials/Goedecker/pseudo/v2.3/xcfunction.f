      SUBROUTINE XCFUNCTION(EXC,RHOE,V,VTMP,GRAD)
c     input:  rho=RHO, GRAD=grad(rho)
c     output: exc=energy density: energy = \int EXC dr
c             v=dEXC/dRHO, vtmp=dEXC/dGRAD
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
C     ==--------------------------------------------------------------==

c     using hutter's routine
      exc=0.0d0
      v=0.0d0
      rhoe=rhoe
      call xc(rhoe,ex,ec,vx,vc)
      exc=ex+ec
      v=vx+vc
      Vtmp=0.d0
      if (rhoe.gt.1.0d-18) then
         dr2=grad*grad
         call gcxc(rhoe,dr2,sx,sc,v1x,v2x,v1c,v2c)
         exc=exc+(sx+sc)/rhoe
         v=v+v1x+v1c
         VTMP= (V2X+V2C)*abs(grad)
      endif
      return
      end
