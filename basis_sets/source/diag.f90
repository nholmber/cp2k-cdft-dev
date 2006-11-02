      SUBROUTINE diag(f,c,s,e,lda,n)
      USE basic_data_types, ONLY: dp
      USE atom
      USE mathlib
      USE upd
      IMPLICIT NONE
      INTEGER, PARAMETER :: naux = 1000
      INTEGER :: lda,n
      REAL(dp) f(lda,*),c(lda,*),s(lda,*),e(*),aux(naux)
      REAL(dp) w(lda),a(lda,lda),s12(lda,lda),b(lda,lda)
      INTEGER i,j,info,ndep

      if(naux.lt.3*n) stop 'diag'

      ndep = 0
      DO j=1,n
        DO i=1,n
          a(i,j) = s(i,j)
        END DO
      END DO
      CALL powmat(a,s12,lda,n,-0.5D0,1.0D-6,ndep)
      IF (ndep > 0) THEN
        print*,"ndep =",ndep
        print*,"Overlap matrix S"
        DO i=1,n
          WRITE (*,"(15F12.6)") (s(i,j),j=1,n)
          DO j=1,n
            a(j,i) = s(j,i)
          END DO
        END DO
        CALL diamat(a,w,lda,n)
        print*,"Eigenvalues and eigenvectors of the overlap matrix S"
        DO i=1,n
          WRITE (*,"(10F12.6)") w(i),(a(i,j),j=1,n)
        END DO
        WRITE (*,*) 'alpha:  '
        DO j=0,lmax
          WRITE (*,'(4F20.10)') (alpha(i,j),i=1,nalpha(j))
        END DO
        WRITE (*,*) 'galpha:  '
        DO j=0,lmax
          WRITE (*,'(4F20.10)') (galpha(i,j),i=1,nalpha(j))
        END DO
        STOP "diag: Overlap matrix has linear dependency"
      END IF
      CALL dgemm("N","N",n,n,n,1.0D0,f,lda,s12,lda,0.0D0,a,lda)
      CALL dgemm("N","N",n,n,n,1.0D0,s12,lda,a,lda,0.0D0,b,lda)
      CALL dsyev("V","U",n,b,lda,e,aux,naux,info)
      CALL dgemm("N","N",n,n,n,1.0D0,s12,lda,b,lda,0.0D0,c,lda)
      IF (ndep > 0) THEN
        print*,"Eigenvalues and Eigenvectors"
        DO i=1,n
          WRITE (*,"(10F12.6)") e(i),(c(i,j),j=1,n)
        END DO
      END IF
      END
