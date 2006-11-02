!------------------------------------------------------------------------------!
  subroutine v_add_d(s,sd,alpha,nalpha,lmax,lamax,namax,uadd_d)
!------------------------------------------------------------------------------!
  USE basic_data_types, ONLY: dp
  IMPLICIT NONE
!out:
  INTEGER    :: lmax,lamax,namax
  REAL(dp)   :: uadd_d(namax,namax,0:lamax)
!in:
  REAL(dp)   :: s(namax,namax,0:lamax)
  REAL(dp)   :: sd(namax,namax,0:lamax)
  REAL(dp)   :: alpha(namax,0:lamax)
  INTEGER    :: nalpha(0:lamax)
!local:
  REAL(dp)   :: raux,p,q
  INTEGER    :: i,j,l,k
!------------------------------------------------------------------------------!

  do l=0,lmax
    raux=l+1.5d0
    do k=1,nalpha(l)
      p=alpha(k,l)
      do j=1,nalpha(l)
        q=alpha(j,l)
        uadd_d(k,j,l)=raux/(p+q)*( sd(k,j,l) - s(k,j,l)/(p+q) )
      enddo
    enddo
  enddo

end subroutine v_add_d



