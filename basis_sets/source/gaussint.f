!------------------------------------------------------------------------------!
  function gaussint(alpha,iexpor)
!------------------------------------------------------------------------------!
	!  calculates the 3-dimensional integral over a cartesian Gaussian
	!  of the form
	!     x**iexpor(1) * y**iexpor(2) * z**iexpor(3) * exp(alpha * r**2)
	!
	!----------------------------------------------------------------------!

  use constants
  implicit none

!in:
  integer(it6)		::  iexpor(3)
  real(rt15)		::  alpha

!out:
  real(rt15)		::  gaussint

!locals:
  integer(it6)		::  isum,i

!------------------------------------------------------------------------------!

  isum=-sum(iexpor)
  iexpor=iexpor-1_it6
  gaussint=PI32 * ROOT2**isum * idblfac(iexpor) * sqrt(alpha)**(isum-3_it6)


contains

  !----------------------------------------------------------------------------!
    function idblfac(ivector)
  !----------------------------------------------------------------------------!

    implicit none

  !in:
    integer(it6)		::  ivector(:)

  !out:
    integer(it6) 		::  idblfac

  !locals:
    integer(it6)		::  idfacl(-1:10) &
  	 			     = (/1,1,1,2,3,8,15,48,105,384,945,3840/)
    integer(it6)		::  idx,i
    integer(it6)		::  ifac(size(ivector))

  !----------------------------------------------------------------------------!

    do idx=1,size(ivector) 
      if(ivector(idx)<-1) then
        print*, ' double factorial of ',ivector(idx),' not defined'
        stop 'dblfac'
      elseif(ivector(idx)<=10) then
        ifac(idx)=idfacl(ivector(idx))
      elseif(modulo(ivector(idx),2_it6)==0) then
       ifac(idx)=3840
        do i=12,ivector(idx),2
          ifac(idx)=ifac(idx)*i
        enddo
      else
        ifac(idx)=945
        do i=11,ivector(idx),2
          ifac(idx)=ifac(idx)*i
        enddo
      endif
    enddo
    idblfac=product(ifac)

  end function idblfac

end function gaussint
