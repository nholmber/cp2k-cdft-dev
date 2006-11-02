      FUNCTION gam1(l,z)
      USE basic_data_types, ONLY: dp
      IMPLICIT NONE
      REAL(dp) :: z,zz,ez,ef,pz,gam1,dblfac
      integer l
      REAL(dp), PARAMETER :: pi=3.14159265358979323846264_dp
      REAL(dp) :: derfc
      zz=z*z
      ez=z*exp(-zz)
      ef=sqrt(pi)*derfc(z)*dblfac(2*l+1)/2.d0**(l+1)
      if(l.eq.0) then
        pz=1.d0
      elseif(l.eq.1) then
        pz=1.5d0+zz
      elseif(l.eq.2) then
        pz=3.75d0+2.5d0*zz+zz*zz
      elseif(l.eq.3) then
        pz=13.125d0+8.75d0*zz+3.5d0*zz*zz+zz*zz*zz
      else
        write(*,*) ' function gam1 not implemented for l=',l
        stop
      endif
      gam1=pz*ez+ef
      END

      FUNCTION gam2(l,z)
      USE basic_data_types, ONLY: dp
      IMPLICIT NONE
      real(dp) :: z,zz,ez,pz,gam2
      integer l
      zz=z*z
      ez=exp(-zz)
      if(l.eq.0) then
        pz=1.d0
      elseif(l.eq.1) then
        pz=1.d0+zz
      elseif(l.eq.2) then
        pz=2.d0+2.d0*zz+zz*zz
      elseif(l.eq.3) then
        pz=6.d0+6.d0*zz+3.d0*zz*zz+zz*zz*zz
      else
        write(*,*) ' function gam2 not implemented for l=',l
        stop
      endif
      gam2=pz*ez
      END
