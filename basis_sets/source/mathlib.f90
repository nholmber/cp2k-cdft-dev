MODULE mathlib

! *** Needs LAPACK ***

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,200)

  PUBLIC :: diamat,invmat,powmat

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE diamat(a,eigval,lda,n)

!   Purpose: Diagonalize the symmetric n by n matrix a using the LAPACK
!            library. Only the upper triangle of matrix a is used.

!   History: - Creation (29.03.1999,MK)

!   ***************************************************************************

!   a       : Symmetric matrix to be diagonalized (input; upper triangle) ->
!             eigenvectors of the matrix a (output).
!   eigval  : Eigenvalues of the matrix a (output).

!   ***************************************************************************

    INTEGER, INTENT(IN)                         :: lda,n
    REAL(dp), DIMENSION(lda,lda), INTENT(INOUT) :: a
    REAL(dp), DIMENSION(lda), INTENT(OUT)       :: eigval

!   *** Local variables ***

    INTEGER :: info,istat,lwork,nb

    REAL(dp), DIMENSION(:), ALLOCATABLE :: work

!   *** Externals (LAPACK) ***

    INTEGER, EXTERNAL :: ilaenv

    EXTERNAL dsyev

!   ---------------------------------------------------------------------------

!   *** Get the optimal work storage size ***

    nb = ilaenv(1,"DSYTRD","U",n,-1,-1,-1)
    lwork = (nb + 2)*n

!   *** Allocate work storage ***

    ALLOCATE (work(lwork),STAT=istat)
    IF (istat /= 0) STOP "diamat: Allocate work"

!   *** Diagonalize the matrix a ***

    CALL dsyev("V","U",n,a(1,1),lda,eigval(1),work(1),lwork,info)

    IF (info /= 0) STOP "diamat: The matrix diagonalization with dsyev failed"

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) STOP "diamat: Deallocate work"

  END SUBROUTINE diamat

! *****************************************************************************

  SUBROUTINE invmat(a,a_inverse,error,option)

!   Purpose: Compute the inverse of the n by n matrix a using the LAPACK
!            library.

!   History: - Creation (23.03.1999,MK)

!   ***************************************************************************

!   a        : Matrix to be inverted (input).
!   a_inverse: Inverse of the matrix a (output).
!   a_lu     : LU factorization of matrix a.
!   a_norm   : Norm of matrix a.
!   error    : Estimated error of the inversion.
!   r_cond   : Reciprocal condition number of the matrix a.
!   trans    : "N" => invert a
!              "T" => invert transpose(a)

!   ***************************************************************************

    REAL(dp), DIMENSION(:,:), INTENT(IN)   :: a
    REAL(dp), DIMENSION(:,:), INTENT(OUT)  :: a_inverse
    REAL(dp), INTENT(OUT)                  :: error
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: option

!   *** Local variables ***

    CHARACTER(LEN=1)  :: norm,trans
    REAL(dp)          :: a_norm,old_error,r_cond
    INTEGER           :: info,istat,iter,n

    REAL(dp), DIMENSION(:), ALLOCATABLE :: berr,ferr,work
    INTEGER, DIMENSION(:), ALLOCATABLE  :: ipiv,iwork

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: a_lu,b

!   *** Externals (LAPACK) ***

    EXTERNAL dgecon,dgerfs,dgetrf,dgetrs

!   ---------------------------------------------------------------------------

!   *** Check for optional parameter ***

    IF (PRESENT(option)) THEN
      trans = option
    ELSE
      trans = "N"
    END IF

!   *** Get the dimension of matrix a ***

    n = SIZE(a,1)

!   *** Check array dimensions ***

    IF (n /= SIZE(a,2)) STOP "invmat: Check the array bounds of parameter #1"

    IF ((n /= SIZE(a_inverse,1)).OR.&
        (n /= SIZE(a_inverse,2))) THEN
      STOP "invmat: Check the array bounds of parameter #2"
    END IF

!   *** Allocate work storage ***

    ALLOCATE (a_lu(n,n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate a_lu"

    ALLOCATE (b(n,n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate b"

    ALLOCATE (berr(n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate berr"

    ALLOCATE (ferr(n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate ferr"

    ALLOCATE (ipiv(n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate ipiv"

    ALLOCATE (iwork(n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate iwork"

    ALLOCATE (work(4*n),STAT=istat)
    IF (istat /= 0) STOP "invmat: Allocate work"

    a_lu(1:n,1:n) = a(1:n,1:n)

!   *** Compute the LU factorization of the matrix a ***

    CALL dgetrf(n,n,a_lu(:,:),n,ipiv(:),info)

    IF (info /= 0) STOP "invmat: The LU factorization in dgetrf failed"

!   *** Compute the norm of the matrix a ***

    IF (trans == "N") THEN
      norm = '1'
    ELSE
      norm = 'I'
    END IF

!   *** Compute the reciprocal of the condition number of a ***

    CALL dgecon(norm,n,a_lu(:,:),n,a_norm,r_cond,work(:),iwork(:),info)

    IF (info /= 0) THEN
      STOP "invmat: The computation of the condition number in dgecon failed"
    END IF

    IF (r_cond < EPSILON(0.0_dp)) THEN
      WRITE (*,"(A,ES10.3,A)") "R_COND =",r_cond
      STOP "invmat: Bad condition number"
    END IF

!   *** Solve a system of linear equations using ***
!   *** the LU factorization computed by dgetrf  ***

    CALL unit_matrix(a_inverse(:,:))

    CALL dgetrs(trans,n,n,a_lu(:,:),n,ipiv(:),a_inverse(:,:),n,info)

    IF (info /= 0) THEN
      STOP "invmat: Solving the system of linear equations in dgetrs failed"
    END IF

!   *** Improve the computed solution iteratively ***

    CALL unit_matrix(b) ! Initialize right-hand sides

    error = 0.0_dp

    DO iter=1,10

      CALL dgerfs(trans,n,n,a(:,:),n,a_lu(:,:),n,ipiv(:),b(:,:),n,&
                  a_inverse(:,:),n,ferr(:),berr(:),work(:),iwork(:),info)

      IF (info /= 0) THEN
        WRITE (*,"(A,I2)") "iter =",iter
        STOP "invmat: Improving the computed solution in dgerfs failed"
      END IF

      old_error = error
      error = MAXVAL(ferr(:))

      IF (ABS(error - old_error) <= EPSILON(0.0_dp)) EXIT

    END DO

!   *** Release work storage ***

    DEALLOCATE (work,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate work"

    DEALLOCATE (iwork,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate iwork"

    DEALLOCATE (ipiv,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate ipiv"

    DEALLOCATE (ferr,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate ferr"

    DEALLOCATE (berr,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate berr"

    DEALLOCATE (b,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate b"

    DEALLOCATE (a_lu,STAT=istat)
    IF (istat /= 0) STOP "invmat: Deallocate a_lu"

  END SUBROUTINE invmat

! *****************************************************************************

  SUBROUTINE powmat(a,a_power,lda,n,exponent,threshold,n_dependent)

!   Purpose: Raise the real symmetric n by n matrix a to the power given by
!            exponent. All eigenvectors with a corresponding eigenvalue lower
!            than threshold are quenched. Only the upper triangle of matrix a
!            is used.

!   History: - Creation (29.03.1999,MK)

!   ***************************************************************************

!   a          : Symmetric matrix to be powered (input; upper triangle) ->
!                Destroyed on exit.
!   a_power    : Power of matrix a => a**exponent (output).
!   exponent   : Matrix exponent (input).
!   n_dependent: Number of the eigenvectors which are linear dependent due to
!                the defined eigval_eps (output).
!   threshold  : Threshold value for eigenvector quenching (input).

!   ***************************************************************************

    INTEGER, INTENT(IN)                         :: lda,n
    REAL(dp), DIMENSION(lda,lda), INTENT(INOUT) :: a
    REAL(dp), DIMENSION(lda,lda), INTENT(OUT)   :: a_power
    REAL(dp), INTENT(IN)                        :: exponent
    REAL(dp), INTENT(IN), OPTIONAL              :: threshold
    INTEGER, INTENT(OUT), OPTIONAL              :: n_dependent

!   *** Local variables ***

    REAL(dp) :: eps_eigval,expa
    INTEGER  :: i,n_dep

    REAL(dp), DIMENSION(lda) :: eigval

    EXTERNAL dsyrk

!   ---------------------------------------------------------------------------

!   *** Define the threshold for the eigenvalue quenching ***

    IF (PRESENT(threshold)) THEN
      eps_eigval = threshold
    ELSE
      eps_eigval = EPSILON(0.0_dp)
    END IF

!   *** Compute the eigenvectors and eigenvalues of the matrix a ***

    CALL diamat(a,eigval,lda,n)

!   *** Build a**exponent with eigenvector quenching ***

    expa = 0.5_dp*exponent

    n_dep = 0

    DO i=1,n
      IF (eigval(i) < eps_eigval) THEN
        a(1:n,i) = 0.0_dp
        n_dep = n_dep + 1
      ELSE
        eigval(i) = eigval(i)**expa
        a(1:n,i) = eigval(i)*a(1:n,i)
      END IF
    END DO

    IF (PRESENT(n_dependent)) THEN
      n_dependent = n_dep
    END IF

!   *** a_power <- a*Transpose(a) ***

    CALL dsyrk("U","N",n,n,1.0_dp,a(1,1),lda,0.0_dp,a_power(1,1),lda)

!   *** Copy upper triangle of matrix a_power to lower triangle ***

    CALL symmat(a_power,lda,n,"upper_to_lower")

  END SUBROUTINE powmat

! *****************************************************************************

  SUBROUTINE symmat(a,lda,n,option)

!   Purpose: Symmetrize the matrix a.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                         :: lda,n
    REAL(dp), DIMENSION(lda,lda), INTENT(INOUT) :: a
    CHARACTER(LEN=*), INTENT(IN)                :: option

!   *** Local variables ***

    INTEGER :: i

!   ---------------------------------------------------------------------------

    IF (option == "lower_to_upper") THEN
      DO i=1,n-1
        a(i,i+1:n) = a(i+1:n,i)
      END DO
    ELSE IF (option == "upper_to_lower") THEN
      DO i=1,n-1
        a(i+1:n,i) = a(i,i+1:n)
      END DO
    ELSE IF (option == "anti_lower_to_upper") THEN
      DO i=1,n-1
        a(i,i+1:n) = -a(i+1:n,i)
      END DO
    ELSE IF (option == "anti_upper_to_lower") THEN
      DO i=1,n-1
        a(i+1:n,i) = -a(i,i+1:n)
      END DO
    ELSE
      STOP "symmat: Invalid option"
    END IF

  END SUBROUTINE symmat

! *****************************************************************************

  SUBROUTINE unit_matrix(a)

!   Purpose: Set the matrix a to be a unit matrix.

!   History: - Creation (16.10.1998,MK)

!   ***************************************************************************

    REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: a

    INTEGER :: i

!   ---------------------------------------------------------------------------

    a = 0.0_dp

    DO i=1,SIZE(a)
      a(i,i) = 1.0_dp
    END DO

  END SUBROUTINE unit_matrix

! *****************************************************************************

END MODULE mathlib
