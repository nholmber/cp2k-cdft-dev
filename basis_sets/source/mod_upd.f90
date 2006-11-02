!------------------------------------------------------------------------------!
  module upd
!------------------------------------------------------------------------------!
    USE basic_data_types, ONLY: dp
    USE atom, ONLY:namax,lamax
    IMPLICIT NONE
    PRIVATE
    SAVE
    REAL(dp),public   :: beta(namax*(lamax+1)),&
                        betaold(namax*(lamax+1))
    REAL(dp),public   :: gbeta(namax*(lamax+1)),&
                        gbetaold(namax*(lamax+1))
    REAL(dp),public   :: hessinv(namax*(lamax+1),namax*(lamax+1))
    REAL(dp),public   :: xi(namax*(lamax+1))
    integer,public   :: nbeta
!------------------------------------------------------------------------------!
  end module upd
!------------------------------------------------------------------------------!

 
