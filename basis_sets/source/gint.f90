        FUNCTION gint(n,a)
        USE basic_data_types, ONLY: dp
        IMPLICIT NONE
        INTEGER n
        REAL(dp) :: a,gint
        REAL(dp) ::rootpi,dblfac,facul
        REAL(dp), PARAMETER :: pi=3.14159265358979323846264_dp
        external dblfac,facul

        rootpi=DSQRT(pi)
        if (2*(n/2).eq.n) then
          gint=rootpi*dblfac(n-1) / (2**(n/2+1) * a**(n/2+.5_dp))
        else
          gint=facul((n-1)/2) / (2 * a**((n+1)/2))
        endif

        end
