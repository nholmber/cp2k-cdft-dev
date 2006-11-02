      FUNCTION facul(n)
      USE basic_data_types, ONLY: dp
      IMPLICIT NONE
      real(dp) facul,facl(0:7)
      INTEGER n,i
      DATA facl / 1.0_dp,1.0_dp,2.0_dp,6.0_dp,24.0_dp,120.0_dp,720.0_dp,&
                  5040.0_dp/
      if(n.lt.0) then
        write(*,*) ' factorial of ',n,' not defined'
        stop 'facul'
      elseif(n.lt.8) then
        facul=facl(n)
      else
        facul=1.0_dp
        do i=2,n
          facul=facul*REAL(i,KIND=dp)
        enddo
      endif
      end

      FUNCTION dblfac(n)
      USE basic_data_types, ONLY: dp
      IMPLICIT NONE
      real(dp) dblfac,dfacl(-1:10)
      integer n,i
      data dfacl / 1.0_dp,1.0_dp,1.0_dp,2.0_dp,3.0_dp,8.0_dp,15.0_dp,48.0_dp,&
                   105.0_dp,384.0_dp,945.0_dp,3840.0_dp/
      if(n.lt.-1) then
        write(*,*) ' double factorial of ',n,' not defined'
        stop 'dblfac'
      elseif(n.lt.11) then
        dblfac=dfacl(n)
      elseif(mod(n,2).eq.0) then
        dblfac=1.0_dp
        do i=2,n,2
          dblfac=dblfac*REAL(i,KIND=dp)
        enddo
      else
        dblfac=1.0_dp
        do i=1,n,2
          dblfac=dblfac*REAL(i,KIND=dp)
        enddo
      endif
      END
