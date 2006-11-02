      FUNCTION sint(l,p,q)
      USE basic_data_types, ONLY: dp
      IMPLICIT none
      REAL(dp) :: sint,p,q
      INTEGER l
      sint=(2.d0*sqrt(p*q)/(p+q))**(l+1.5d0)
      END
