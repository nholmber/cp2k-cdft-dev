subroutine v_add(s,uadd)

  USE basic_data_types, ONLY: dp
  USE atom
  IMPLICIT NONE
  REAL(dp) :: uadd(namax,namax,0:lamax),&
              s(namax,namax,0:lamax),&
              p,q,raux
  integer :: i,j,l

  do l=0,lmax
    raux=l+1.5d0
    do i=1,nalpha(l)
      p=alpha(i,l)
      do j=i,nalpha(l)
        q=alpha(j,l)
        uadd(i,j,l)=raux/(p+q)*s(i,j,l)
        uadd(j,i,l)=uadd(i,j,l)
      enddo
    enddo
  enddo

end subroutine v_add

