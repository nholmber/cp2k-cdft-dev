MODULE fft
  USE kinds, ONLY: dbl
  PRIVATE

  PUBLIC :: rdft

CONTAINS
! Fast Fourier/Cosine/Sine Transform
!     dimension   :one
!     data length :power of 2
!     decimation  :frequency
!     radix       :4, 2
!     data        :inplace
!     table       :use
! subroutines
!     cdft: Complex Discrete Fourier Transform
!     rdft: Real Discrete Fourier Transform
!     ddct: Discrete Cosine Transform
!     ddst: Discrete Sine Transform
!     dfct: Cosine Transform of RDFT (Real Symmetric DFT)
!     dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
!
!
! -------- Complex DFT (Discrete Fourier Transform) --------
!     [definition]
!         <case1>
!             X(k) = sum_j=0^n-1 x(j)*exp(2*pi*i*j*k/n), 0<=k<n
!         <case2>
!             X(k) = sum_j=0^n-1 x(j)*exp(-2*pi*i*j*k/n), 0<=k<n
!         (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call cdft(2*n, -1, a, ip, w)
!     [parameters]
!         2*n          :data length (integer)
!                       n >= 1, n = power of 2
!         a(0:2*n-1)   :input/output data (real(kind=dbl))
!                       input data
!                           a(2*j) = Re(x(j)), 
!                           a(2*j+1) = Im(x(j)), 0<=j<n
!                       output data
!                           a(2*k) = Re(X(k)), 
!                           a(2*k+1) = Im(X(k)), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n)  ; if mod(n,4) = 0
!                                       2+sqrt(n/2); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real(kind=dbl))
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call cdft(2*n, -1, a, ip, w)
!         is 
!             call cdft(2*n, 1, a, ip, w)
!             do j = 0, 2 * n - 1
!                 a(j) = a(j) / n
!             end do
!         .
!
!
! -------- Real DFT / Inverse of Real DFT --------
!     [definition]
!         <case1> RDFT
!             R(k) = sum_j=0^n-1 a(j)*cos(2*pi*j*k/n), 0<=k<=n/2
!             I(k) = sum_j=0^n-1 a(j)*sin(2*pi*j*k/n), 0<k<n/2
!         <case2> IRDFT (excluding scale)
!             a(k) = R(0)/2 + R(n/2)/2 + 
!                    sum_j=1^n/2-1 R(j)*cos(2*pi*j*k/n) + 
!                    sum_j=1^n/2-1 I(j)*sin(2*pi*j*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call rdft(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call rdft(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real(kind=dbl))
!                       <case1>
!                           output data
!                               a(2*k) = R(k), 0<=k<n/2
!                               a(2*k+1) = I(k), 0<k<n/2
!                               a(1) = R(n/2)
!                       <case2>
!                           input data
!                               a(2*j) = R(j), 0<=j<n/2
!                               a(2*j+1) = I(j), 0<j<n/2
!                               a(1) = R(n/2)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n/2-1)   :cos/sin table (real(kind=dbl))
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call rdft(n, 1, a, ip, w)
!         is 
!             call rdft(n, -1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
!     [definition]
!         <case1> IDCT (excluding scale)
!             C(k) = sum_j=0^n-1 a(j)*cos(pi*j*(k+1/2)/n), 0<=k<n
!         <case2> DCT
!             C(k) = sum_j=0^n-1 a(j)*cos(pi*(j+1/2)*k/n), 0<=k<n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call ddct(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call ddct(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real(kind=dbl))
!                       output data
!                           a(k) = C(k), 0<=k<n
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/4-1) :cos/sin table (real(kind=dbl))
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call ddct(n, -1, a, ip, w)
!         is 
!             a(0) = a(0) / 2
!             call ddct(n, 1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- DST (Discrete Sine Transform) / Inverse of DST --------
!     [definition]
!         <case1> IDST (excluding scale)
!             S(k) = sum_j=1^n A(j)*sin(pi*j*(k+1/2)/n), 0<=k<n
!         <case2> DST
!             S(k) = sum_j=0^n-1 a(j)*sin(pi*(j+1/2)*k/n), 0<k<=n
!     [usage]
!         <case1>
!             ip(0) = 0  ! first time only
!             call ddst(n, 1, a, ip, w)
!         <case2>
!             ip(0) = 0  ! first time only
!             call ddst(n, -1, a, ip, w)
!     [parameters]
!         n            :data length (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real(kind=dbl))
!                       <case1>
!                           input data
!                               a(j) = A(j), 0<j<n
!                               a(0) = A(n)
!                           output data
!                               a(k) = S(k), 0<=k<n
!                       <case2>
!                           output data
!                               a(k) = S(k), 0<k<n
!                               a(0) = S(n)
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/2); if mod(n,4) = 2
!                                       2+sqrt(n/4); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/4-1) :cos/sin table (real(kind=dbl))
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call ddst(n, -1, a, ip, w)
!         is 
!             a(0) = a(0) / 2
!             call ddst(n, 1, a, ip, w)
!             do j = 0, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- Cosine Transform of RDFT (Real Symmetric DFT) --------
!     [definition]
!         C(k) = sum_j=0^n a(j)*cos(pi*j*k/n), 0<=k<=n
!     [usage]
!         ip(0) = 0  ! first time only
!         call dfct(n, a, t, ip, w)
!     [parameters]
!         n            :data length - 1 (integer)
!                       n >= 2, n = power of 2
!         a(0:n)       :input/output data (real(kind=dbl))
!                       output data
!                           a(k) = C(k), 0<=k<=n
!         t(0:n/2)     :work area (real(kind=dbl))
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/4); if mod(n,4) = 0
!                                       2+sqrt(n/8); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/8-1) :cos/sin table (real(kind=dbl))
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             a(0) = a(0) / 2
!             a(n) = a(n) / 2
!             call dfct(n, a, t, ip, w)
!         is 
!             a(0) = a(0) / 2
!             a(n) = a(n) / 2
!             call dfct(n, a, t, ip, w)
!             do j = 0, n
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
! -------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
!     [definition]
!         S(k) = sum_j=1^n-1 a(j)*sin(pi*j*k/n), 0<k<n
!     [usage]
!         ip(0) = 0  ! first time only
!         call dfst(n, a, t, ip, w)
!     [parameters]
!         n            :data length + 1 (integer)
!                       n >= 2, n = power of 2
!         a(0:n-1)     :input/output data (real(kind=dbl))
!                       output data
!                           a(k) = S(k), 0<k<n
!                       (a(0) is used for work area)
!         t(0:n/2-1)   :work area (real(kind=dbl))
!         ip(0:*)      :work area for bit reversal (integer)
!                       length of ip >= 2+sqrt(n/4); if mod(n,4) = 0
!                                       2+sqrt(n/8); otherwise
!                       ip(0),ip(1) are pointers of the cos/sin table.
!         w(0:n*5/8-1) :cos/sin table (real(kind=dbl))
!                       w(),ip() are initialized if ip(0) = 0.
!     [remark]
!         Inverse of 
!             call dfst(n, a, t, ip, w)
!         is 
!             call dfst(n, a, t, ip, w)
!             do j = 1, n - 1
!                 a(j) = a(j) * 2 / n
!             end do
!         .
!
!
  SUBROUTINE cdft(n, isgn, a, ip, w)
    INTEGER n, isgn, ip(0 : *), j
    REAL(kind=dbl) a(0 : n - 1), w(0 : *)
    IF (n .GT. 4 * ip(0)) THEN
       CALL makewt(n / 4, ip, w)
    END IF
    IF (n .GT. 4) CALL bitrv2(n, ip(2), a)
    IF (n .GT. 4 .AND. isgn .LT. 0) THEN
       DO j = 1, n - 1, 2
          a(j) = -a(j)
       END DO
       CALL cftsub(n, a, w)
       DO j = 1, n - 1, 2
          a(j) = -a(j)
       END DO
    ELSE
       CALL cftsub(n, a, w)
    END IF
  END SUBROUTINE cdft
  !
  SUBROUTINE rdft(n, isgn, a, ip, w)
    INTEGER n, isgn, ip(0 : *), j, nw, nc
    REAL(kind=dbl) a(0 : n - 1), w(0 : *), xi
    nw = ip(0)
    IF (n .GT. 4 * nw) THEN
       nw = n / 4
       CALL makewt(nw, ip, w)
    END IF
    nc = ip(1)
    IF (n .GT. 4 * nc) THEN
       nc = n / 4
       CALL makect(nc, ip, w(nw))
    END IF
    IF (isgn .LT. 0) THEN
       a(1) = 0.5d0 * (a(1) - a(0))
       a(0) = a(0) + a(1)
       DO j = 3, n - 1, 2
          a(j) = -a(j)
       END DO
       IF (n .GT. 4) THEN
          CALL rftsub(n, a, nc, w(nw))
          CALL bitrv2(n, ip(2), a)
       END IF
       CALL cftsub(n, a, w)
       DO j = 1, n - 1, 2
          a(j) = -a(j)
       END DO
    ELSE
       IF (n .GT. 4) CALL bitrv2(n, ip(2), a)
       CALL cftsub(n, a, w)
       IF (n .GT. 4) CALL rftsub(n, a, nc, w(nw))
       xi = a(0) - a(1)
       a(0) = a(0) + a(1)
       a(1) = xi
    END IF
  END SUBROUTINE rdft
  !
  SUBROUTINE ddct(n, isgn, a, ip, w)
    INTEGER n, isgn, ip(0 : *), j, nw, nc
    REAL(kind=dbl) a(0 : n - 1), w(0 : *), xr
    nw = ip(0)
    IF (n .GT. 4 * nw) THEN
       nw = n / 4
       CALL makewt(nw, ip, w)
    END IF
    nc = ip(1)
    IF (n .GT. nc) THEN
       nc = n
       CALL makect(nc, ip, w(nw))
    END IF
    IF (isgn .LT. 0) THEN
       xr = a(n - 1)
       DO j = n - 2, 2, -2
          a(j + 1) = a(j - 1) - a(j)
          a(j) = a(j) + a(j - 1)
       END DO
       a(1) = xr - a(0)
       a(0) = a(0) + xr
       IF (n .GT. 4) THEN
          CALL rftsub(n, a, nc, w(nw))
          CALL bitrv2(n, ip(2), a)
       END IF
       CALL cftsub(n, a, w)
       DO j = 1, n - 1, 2
          a(j) = -a(j)
       END DO
    END IF
    CALL dctsub(n, a, nc, w(nw))
    IF (isgn .GE. 0) THEN
       IF (n .GT. 4) CALL bitrv2(n, ip(2), a)
       CALL cftsub(n, a, w)
       IF (n .GT. 4) CALL rftsub(n, a, nc, w(nw))
       xr = a(0) - a(1)
       a(0) = a(0) + a(1)
       DO j = 2, n - 2, 2
          a(j - 1) = a(j) - a(j + 1)
          a(j) = a(j) + a(j + 1)
       END DO
       a(n - 1) = xr
    END IF
  END SUBROUTINE ddct
  !
  SUBROUTINE ddst(n, isgn, a, ip, w)
    INTEGER n, isgn, ip(0 : *), j, nw, nc
    REAL(kind=dbl) a(0 : n - 1), w(0 : *), xr
    nw = ip(0)
    IF (n .GT. 4 * nw) THEN
       nw = n / 4
       CALL makewt(nw, ip, w)
    END IF
    nc = ip(1)
    IF (n .GT. nc) THEN
       nc = n
       CALL makect(nc, ip, w(nw))
    END IF
    IF (isgn .LT. 0) THEN
       xr = a(n - 1)
       DO j = n - 2, 2, -2
          a(j + 1) = a(j - 1) + a(j)
          a(j) = a(j) - a(j - 1)
       END DO
       a(1) = -xr - a(0)
       a(0) = a(0) - xr
       IF (n .GT. 4) THEN
          CALL rftsub(n, a, nc, w(nw))
          CALL bitrv2(n, ip(2), a)
       END IF
       CALL cftsub(n, a, w)
       DO j = 1, n - 1, 2
          a(j) = -a(j)
       END DO
    END IF
    CALL dstsub(n, a, nc, w(nw))
    IF (isgn .GE. 0) THEN
       IF (n .GT. 4) CALL bitrv2(n, ip(2), a)
       CALL cftsub(n, a, w)
       IF (n .GT. 4) CALL rftsub(n, a, nc, w(nw))
       xr = a(0) - a(1)
       a(0) = a(0) + a(1)
       DO j = 2, n - 2, 2
          a(j - 1) = -a(j) - a(j + 1)
          a(j) = a(j) - a(j + 1)
       END DO
       a(n - 1) = -xr
    END IF
  END SUBROUTINE ddst
  !
  SUBROUTINE dfct(n, a, t, ip, w)
    INTEGER n, ip(0 : *), j, k, l, m, mh, nw, nc
    REAL(kind=dbl) a(0 : n), t(0 : n / 2), w(0 : *), xr, xi
    nw = ip(0)
    IF (n .GT. 8 * nw) THEN
       nw = n / 8
       CALL makewt(nw, ip, w)
    END IF
    nc = ip(1)
    IF (n .GT. 2 * nc) THEN
       nc = n / 2
       CALL makect(nc, ip, w(nw))
    END IF
    m = n / 2
    xr = a(0) + a(n)
    a(0) = a(0) - a(n)
    t(0) = xr - a(m)
    t(m) = xr + a(m)
    IF (n .GT. 2) THEN
       mh = m / 2
       DO j = 1, mh - 1
          k = m - j
          xr = a(j) + a(n - j)
          a(j) = a(j) - a(n - j)
          xi = a(k) + a(n - k)
          a(k) = a(k) - a(n - k)
          t(j) = xr - xi
          t(k) = xr + xi
       END DO
       t(mh) = a(mh) + a(n - mh)
       a(mh) = a(mh) - a(n - mh)
       CALL dctsub(m, a, nc, w(nw))
       IF (m .GT. 4) CALL bitrv2(m, ip(2), a)
       CALL cftsub(m, a, w)
       IF (m .GT. 4) CALL rftsub(m, a, nc, w(nw))
       xr = a(0) + a(1)
       a(n - 1) = a(0) - a(1)
       DO j = m - 2, 2, -2
          a(2 * j + 1) = a(j) + a(j + 1)
          a(2 * j - 1) = a(j) - a(j + 1)
       END DO
       a(1) = xr
       l = 2
       m = mh
       DO WHILE (m .GE. 2)
          CALL dctsub(m, t, nc, w(nw))
          IF (m .GT. 4) CALL bitrv2(m, ip(2), t)
          CALL cftsub(m, t, w)
          IF (m .GT. 4) CALL rftsub(m, t, nc, w(nw))
          a(n - l) = t(0) - t(1)
          a(l) = t(0) + t(1)
          k = 0
          DO j = 2, m - 2, 2
             k = k + 4 * l
             a(k - l) = t(j) - t(j + 1)
             a(k + l) = t(j) + t(j + 1)
          END DO
          l = 2 * l
          mh = m / 2
          DO j = 0, mh - 1
             k = m - j
             t(j) = t(m + k) - t(m + j)
             t(k) = t(m + k) + t(m + j)
          END DO
          t(mh) = t(m + mh)
          m = mh
       END DO
       a(l) = t(0)
       a(n) = t(2) - t(1)
       a(0) = t(2) + t(1)
    ELSE
       a(1) = a(0)
       a(2) = t(0)
       a(0) = t(1)
    END IF
  END SUBROUTINE dfct
  !
  SUBROUTINE dfst(n, a, t, ip, w)
    INTEGER n, ip(0 : *), j, k, l, m, mh, nw, nc
    REAL(kind=dbl) a(0 : n - 1), t(0 : n / 2 - 1), w(0 : *), xr, xi
    nw = ip(0)
    IF (n .GT. 8 * nw) THEN
       nw = n / 8
       CALL makewt(nw, ip, w)
    END IF
    nc = ip(1)
    IF (n .GT. 2 * nc) THEN
       nc = n / 2
       CALL makect(nc, ip, w(nw))
    END IF
    IF (n .GT. 2) THEN
       m = n / 2
       mh = m / 2
       DO j = 1, mh - 1
          k = m - j
          xr = a(j) - a(n - j)
          a(j) = a(j) + a(n - j)
          xi = a(k) - a(n - k)
          a(k) = a(k) + a(n - k)
          t(j) = xr + xi
          t(k) = xr - xi
       END DO
       t(0) = a(mh) - a(n - mh)
       a(mh) = a(mh) + a(n - mh)
       a(0) = a(m)
       CALL dstsub(m, a, nc, w(nw))
       IF (m .GT. 4) CALL bitrv2(m, ip(2), a)
       CALL cftsub(m, a, w)
       IF (m .GT. 4) CALL rftsub(m, a, nc, w(nw))
       xr = a(0) + a(1)
       a(n - 1) = a(1) - a(0)
       DO j = m - 2, 2, -2
          a(2 * j + 1) = a(j) - a(j + 1)
          a(2 * j - 1) = -a(j) - a(j + 1)
       END DO
       a(1) = xr
       l = 2
       m = mh
       DO WHILE (m .GE. 2)
          CALL dstsub(m, t, nc, w(nw))
          IF (m .GT. 4) CALL bitrv2(m, ip(2), t)
          CALL cftsub(m, t, w)
          IF (m .GT. 4) CALL rftsub(m, t, nc, w(nw))
          a(n - l) = t(1) - t(0)
          a(l) = t(0) + t(1)
          k = 0
          DO j = 2, m - 2, 2
             k = k + 4 * l
             a(k - l) = -t(j) - t(j + 1)
             a(k + l) = t(j) - t(j + 1)
          END DO
          l = 2 * l
          mh = m / 2
          DO j = 1, mh - 1
             k = m - j
             t(j) = t(m + k) + t(m + j)
             t(k) = t(m + k) - t(m + j)
          END DO
          t(0) = t(m + mh)
          m = mh
       END DO
       a(l) = t(0)
    END IF
    a(0) = 0
  END SUBROUTINE dfst
  !
  ! -------- initializing routines --------
  !
  SUBROUTINE makewt(nw, ip, w)
    INTEGER nw, ip(0 : *), nwh, j
    REAL(kind=dbl) w(0 : nw - 1), delta, x, y
    ip(0) = nw
    ip(1) = 1
    IF (nw .GT. 2) THEN
       nwh = nw / 2
       delta = ATAN(1.0d0) / nwh
       w(0) = 1
       w(1) = 0
       w(nwh) = COS(delta * nwh)
       w(nwh + 1) = w(nwh)
       DO j = 2, nwh - 2, 2
          x = COS(delta * j)
          y = SIN(delta * j)
          w(j) = x
          w(j + 1) = y
          w(nw - j) = y
          w(nw - j + 1) = x
       END DO
       CALL bitrv2(nw, ip(2), w)
    END IF
  END SUBROUTINE makewt
  !
  SUBROUTINE makect(nc, ip, c)
    INTEGER nc, ip(0 : *), nch, j
    REAL(kind=dbl) c(0 : nc - 1), delta
    ip(1) = nc
    IF (nc .GT. 1) THEN
       nch = nc / 2
       delta = ATAN(1.0d0) / nch
       c(0) = 0.5d0
       c(nch) = 0.5d0 * COS(delta * nch)
       DO j = 1, nch - 1
          c(j) = 0.5d0 * COS(delta * j)
          c(nc - j) = 0.5d0 * SIN(delta * j)
       END DO
    END IF
  END SUBROUTINE makect
  !
  ! -------- child routines --------
  !
  SUBROUTINE bitrv2(n, ip, a)
    INTEGER n, ip(0 : *), j, j1, k, k1, l, m, m2
    REAL(kind=dbl) a(0 : n - 1), xr, xi
    ip(0) = 0
    l = n
    m = 1
    DO WHILE (4 * m .LT. l)
       l = l / 2
       DO j = 0, m - 1
          ip(m + j) = ip(j) + l
       END DO
       m = m * 2
    END DO
    IF (4 * m .GT. l) THEN
       DO k = 1, m - 1
          DO j = 0, k - 1
             j1 = 2 * j + ip(k)
             k1 = 2 * k + ip(j)
             xr = a(j1)
             xi = a(j1 + 1)
             a(j1) = a(k1)
             a(j1 + 1) = a(k1 + 1)
             a(k1) = xr
             a(k1 + 1) = xi
          END DO
       END DO
    ELSE
       m2 = 2 * m
       DO k = 1, m - 1
          DO j = 0, k - 1
             j1 = 2 * j + ip(k)
             k1 = 2 * k + ip(j)
             xr = a(j1)
             xi = a(j1 + 1)
             a(j1) = a(k1)
             a(j1 + 1) = a(k1 + 1)
             a(k1) = xr
             a(k1 + 1) = xi
             j1 = j1 + m2
             k1 = k1 + m2
             xr = a(j1)
             xi = a(j1 + 1)
             a(j1) = a(k1)
             a(j1 + 1) = a(k1 + 1)
             a(k1) = xr
             a(k1 + 1) = xi
          END DO
       END DO
    END IF
  END SUBROUTINE bitrv2
  !
  SUBROUTINE cftsub(n, a, w)
    INTEGER n, j, j1, j2, j3, k, k1, ks, l, m
    REAL(kind=dbl) a(0 : n - 1), w(0 : *)
    REAL(kind=dbl) wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
    REAL(kind=dbl) x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
    l = 2
    DO WHILE (2 * l .LT. n)
       m = 4 * l
       DO j = 0, l - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x0r - x2r
          a(j2 + 1) = x0i - x2i
          a(j1) = x1r - x3i
          a(j1 + 1) = x1i + x3r
          a(j3) = x1r + x3i
          a(j3 + 1) = x1i - x3r
       END DO
       IF (m .LT. n) THEN
          wk1r = w(2)
          DO j = m, l + m - 2, 2
             j1 = j + l
             j2 = j1 + l
             j3 = j2 + l
             x0r = a(j) + a(j1)
             x0i = a(j + 1) + a(j1 + 1)
             x1r = a(j) - a(j1)
             x1i = a(j + 1) - a(j1 + 1)
             x2r = a(j2) + a(j3)
             x2i = a(j2 + 1) + a(j3 + 1)
             x3r = a(j2) - a(j3)
             x3i = a(j2 + 1) - a(j3 + 1)
             a(j) = x0r + x2r
             a(j + 1) = x0i + x2i
             a(j2) = x2i - x0i
             a(j2 + 1) = x0r - x2r
             x0r = x1r - x3i
             x0i = x1i + x3r
             a(j1) = wk1r * (x0r - x0i)
             a(j1 + 1) = wk1r * (x0r + x0i)
             x0r = x3i + x1r
             x0i = x3r - x1i
             a(j3) = wk1r * (x0i - x0r)
             a(j3 + 1) = wk1r * (x0i + x0r)
          END DO
          k1 = 1
          ks = -1
          DO k = 2 * m, n - m, m
             k1 = k1 + 1
             ks = -ks
             wk1r = w(2 * k1)
             wk1i = w(2 * k1 + 1)
             wk2r = ks * w(k1)
             wk2i = w(k1 + ks)
             wk3r = wk1r - 2 * wk2i * wk1i
             wk3i = 2 * wk2i * wk1r - wk1i
             DO j = k, l + k - 2, 2
                j1 = j + l
                j2 = j1 + l
                j3 = j2 + l
                x0r = a(j) + a(j1)
                x0i = a(j + 1) + a(j1 + 1)
                x1r = a(j) - a(j1)
                x1i = a(j + 1) - a(j1 + 1)
                x2r = a(j2) + a(j3)
                x2i = a(j2 + 1) + a(j3 + 1)
                x3r = a(j2) - a(j3)
                x3i = a(j2 + 1) - a(j3 + 1)
                a(j) = x0r + x2r
                a(j + 1) = x0i + x2i
                x0r = x0r - x2r
                x0i = x0i - x2i
                a(j2) = wk2r * x0r - wk2i * x0i
                a(j2 + 1) = wk2r * x0i + wk2i * x0r
                x0r = x1r - x3i
                x0i = x1i + x3r
                a(j1) = wk1r * x0r - wk1i * x0i
                a(j1 + 1) = wk1r * x0i + wk1i * x0r
                x0r = x1r + x3i
                x0i = x1i - x3r
                a(j3) = wk3r * x0r - wk3i * x0i
                a(j3 + 1) = wk3r * x0i + wk3i * x0r
             END DO
          END DO
       END IF
       l = m
    END DO
    IF (l .LT. n) THEN
       DO j = 0, l - 2, 2
          j1 = j + l
          x0r = a(j) - a(j1)
          x0i = a(j + 1) - a(j1 + 1)
          a(j) = a(j) + a(j1)
          a(j + 1) = a(j + 1) + a(j1 + 1)
          a(j1) = x0r
          a(j1 + 1) = x0i
       END DO
    END IF
  END SUBROUTINE cftsub
  !
  SUBROUTINE rftsub(n, a, nc, c)
    INTEGER n, nc, j, k, kk, ks
    REAL(kind=dbl) a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
    ks = 4 * nc / n
    kk = 0
    DO k = n / 2 - 2, 2, -2
       j = n - k
       kk = kk + ks
       wkr = 0.5d0 - c(kk)
       wki = c(nc - kk)
       xr = a(k) - a(j)
       xi = a(k + 1) + a(j + 1)
       yr = wkr * xr - wki * xi
       yi = wkr * xi + wki * xr
       a(k) = a(k) - yr
       a(k + 1) = a(k + 1) - yi
       a(j) = a(j) + yr
       a(j + 1) = a(j + 1) - yi
    END DO
  END SUBROUTINE rftsub
  !
  SUBROUTINE dctsub(n, a, nc, c)
    INTEGER n, nc, j, k, kk, ks, m
    REAL(kind=dbl) a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
    ks = nc / n
    kk = ks
    m = n / 2
    DO k = 1, m - 1
       j = n - k
       wkr = c(kk) - c(nc - kk)
       wki = c(kk) + c(nc - kk)
       kk = kk + ks
       xr = wki * a(k) - wkr * a(j)
       a(k) = wkr * a(k) + wki * a(j)
       a(j) = xr
    END DO
    a(m) = 2 * c(kk) * a(m)
  END SUBROUTINE dctsub
  !
  SUBROUTINE dstsub(n, a, nc, c)
    INTEGER n, nc, j, k, kk, ks, m
    REAL(kind=dbl) a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
    ks = nc / n
    kk = ks
    m = n / 2
    DO k = 1, m - 1
       j = n - k
       wkr = c(kk) - c(nc - kk)
       wki = c(kk) + c(nc - kk)
       kk = kk + ks
       xr = wki * a(j) - wkr * a(k)
       a(j) = wkr * a(j) + wki * a(k)
       a(k) = xr
    END DO
    a(m) = 2 * c(kk) * a(m)
  END SUBROUTINE dstsub
  !
END MODULE fft
