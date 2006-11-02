      SUBROUTINE denmat(p,c)
      USE basic_data_types, ONLY: dp
      USE atom
      IMPLICIT NONE
      REAL(dp) :: p(namax,namax,0:lamax),c(namax,namax,0:lamax)
      INTEGER :: i,l,kappa,lambda
      do l=0,lmax
        do kappa=1,nalpha(l)
          do lambda=kappa,nalpha(l)
            p(kappa,lambda,l)=0.0_dp
            do i=1,nocc(l)
              p(kappa,lambda,l)=p(kappa,lambda,l) +&
                                occ(i,l)*c(kappa,i,l)*c(lambda,i,l)
            enddo
            p(lambda,kappa,l)=p(kappa,lambda,l)
         enddo
        enddo
      enddo
      END
