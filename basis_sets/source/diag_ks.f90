      SUBROUTINE diag_ks(ksmat,smat,cmat,ksener)
      USE basic_data_types, ONLY: dp
      USE atom
      IMPLICIT NONE
      REAL(dp) :: ksmat(namax,namax,0:lamax)
      REAL(dp) :: smat(namax,namax,0:lamax)
      REAL(dp) :: cmat(namax,namax,0:lamax)
      REAL(dp) :: x1mat(namax,namax,0:lamax)
      REAL(dp) :: x2mat(namax,namax,0:lamax)
      REAL(dp) :: ksener(namax,0:lamax)
      INTEGER :: l,n

      x1mat=ksmat
      x2mat=smat
      do l=0,lmax
        n=nalpha(l)
        call diag(ksmat(1,1,l),cmat(1,1,l),smat(1,1,l),ksener(1,l),namax,n)
      enddo
      ksmat=x1mat
      smat=x2mat
      end
