      SUBROUTINE calcnn(cnn)
      USE basic_data_types, ONLY: dp
      USE atom
      IMPLICIT NONE
      REAL(dp) cnn(namax,namax,0:lamax)
      REAL(dp) pi,rootpi,dblfac
      INTEGER l,kappa,lambda
      PARAMETER (pi=3.14159265358979323846264_dp)
      rootpi=dsqrt(pi)

      do l=0,lmax
        do kappa=1,nalpha(l)
          do lambda=kappa,nalpha(l)
            cnn(kappa,lambda,l)=2.D0**(2*l+7.D0/2.D0)/rootpi/dblfac(2*l+1)*&
                                (alpha(kappa,l)*alpha(lambda,l))**((2.D0*l+3.D0)/4.D0)
            cnn(lambda,kappa,l)=cnn(kappa,lambda,l)
          enddo
        enddo
      enddo

      return
      end
