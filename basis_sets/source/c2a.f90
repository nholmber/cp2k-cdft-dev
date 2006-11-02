      SUBROUTINE C2A
      USE basic_data_types,  ONLY: dp
      USE pspot
      IMPLICIT NONE
      REAL(dp) C(6),Q(6,6),S(6,6),alpha(3)
      INTEGER i,l,ipvt(6),info

      DO l=0,PPlmax
         DO i=1,3
            alpha(i)=PPNexp(i,l)
         ENDDO
         DO i=1,6
            C(i)=PPLexp(i,l)
         ENDDO
!..Calculate S,Q
         CALL CALCS(alpha,S)
         CALL CALCQ(S,Q)
!..Calculate A's from C's
         CALL dgetrf(6,6,q,6,ipvt,info)
         CALL dgetrs("N",6,1,q,6,ipvt,c,6,info)

         DO i=1,6
            PPLexp(i,l)=-C(i)
         ENDDO
      ENDDO
      END

      SUBROUTINE CALCQ(S,Q)
      USE basic_data_types,  ONLY: dp
      IMPLICIT NONE
      REAL(dp) S(6,6),Q(6,6)
      INTEGER i,k,l

!     Calculate Q from S
      DO i=1,6
         DO l=1,6
            IF (i.GT.l) THEN
               Q(i,l)=0.0_dp
            ELSE IF (i.EQ.l) THEN
               Q(i,l)=S(i,l)
               DO k=1,i-1
                  Q(i,l)=Q(i,l)-Q(k,i)**2.0_dp
               ENDDO
               Q(i,l)=Q(i,l)**0.5_dp
            ELSE
               Q(i,l)=S(i,l)
               DO k=1,i-1
                  Q(i,l)=Q(i,l)-Q(k,i)*Q(k,l)
               ENDDO
               Q(i,l)=Q(i,l)/Q(i,i)
            ENDIF
         ENDDO
      ENDDO
      END

      SUBROUTINE CALCS(alpha,S)
      USE basic_data_types,  ONLY: dp
      IMPLICIT NONE
      REAL(dp) alpha(3),alfa(6),S(6,6),ROOTPI,PI
      INTEGER m,n
      PARAMETER (PI=3.14159265358979323846264_dp)

      ROOTPI=DSQRT(PI)
!     Calculate S

      DO m=1,3
         alfa(m)=alpha(m)
         alfa(m+3)=alpha(m)
      ENDDO

      DO m=1,6
         DO n=1,6
            IF (m.le.3) THEN
               IF (n.le.3) THEN
                  S(m,n)=ROOTPI/4.0_dp*(alfa(m)+alfa(n))**(-3.0_dp/2.0_dp)
               ELSE
                 S(m,n)=3.0_dp/8.0_dp*ROOTPI*(alfa(m)+alfa(n))**(-5.0_dp/2.0_dp)
               END IF
            ELSE
               IF (n.le.3) THEN
                 S(m,n)=3.0_dp/8.0_dp*ROOTPI*(alfa(m)+alfa(n))**(-5.0_dp/2.0_dp)
               ELSE
                  S(m,n)=15.0_dp/16.0_dp*ROOTPI*(alfa(m)+alfa(n))**(-7.0_dp/2.0_dp)
               END IF
            END IF
         ENDDO
      ENDDO

      END
