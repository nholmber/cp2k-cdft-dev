!------------------------------------------------------------------------------!
  module pspot
!------------------------------------------------------------------------------!
   USE basic_data_types, ONLY: dp
   USE atom, ONLY: namax,lamax
   IMPLICIT NONE
   PRIVATE
   SAVE
!..Parameters for the pseudopotential in the BHS form
!	Vbhs= SUM_i^ERFnr( -Zeff/r* PPLerf(i)*ERF(PPNerf(i)*r) )
!	 +SUM_i^EXPnr( PPLexp(i,l)*r**(PPRexp-2)*EXP(-PPNexp(i,l)*r**2 )
	
!..Effective core charge
    REAL(dp),public           :: Zeff
!..Highest angular momentum component in the pseudopotential
    integer,public           :: PPlmax
!..Limit for PPlmax
    integer,parameter,public :: PPlamax=3
!..Number of sets of coefficients for ERF and EXP part of Vbhs rsp. VKB_core
    integer,public           :: ERFnr,EXPnr(0:PPlamax)
!..Limits for ERFnr,EXPnr and PPlmax
    integer,parameter,public :: ERFmax=4,&
                                EXPmax=8
!..Linear coefficients for the BHS pseudopotential rsp. VKB_core
    REAL(dp),public           :: PPLerf(ERFmax),&
                                 PPLexp(EXPmax,0:PPlamax)
    REAL(dp),public           :: PPC(EXPmax)
!..Exponents of r in the general form of the pseudopotential rsp. VKB_core
    integer,public           :: PPRexp(EXPmax,0:PPlamax)
!..Nonlinear coefficients for the BHS pseudopotential rsp. VKB_core
    REAL(dp),public           :: PPNerf(ERFmax),PPNexp(EXPmax,0:PPlamax)
    REAL(dp),public           :: r_loc

!..Limits for KBPROJnr,KBEXPnr
    integer,parameter,public :: KBPROJmax=4,&
                                KBEXPmax=EXPmax
!..Number of Gaussians per projector
    integer,public           :: KBEXPnr(0:PPlamax)
!..Number of projectors in Vkb
    integer,public           :: KBPROJnr(0:PPlamax)
!..Coefficient matrix for the projectors
    REAL(dp),public           :: KBV(KBPROJmax,KBPROJmax,0:PPlamax)
!..Linear coefficients for the projectors in the KB-like pseudopotential
    REAL(dp),public           :: KBLexp(KBEXPmax,KBPROJmax,0:PPlamax)
!..Exponents of r for the projectors in the KB-like pseudopotential
    integer,public           :: KBRexp(KBEXPmax,KBPROJmax,0:PPlamax)
!..Gaussian exponents for the projectors in the KB-like pseudopotential
    REAL(dp),public           :: KBNexp(KBEXPmax,KBPROJmax,0:PPlamax)
    REAL(dp),public           :: r_proj(0:3)

    REAL(dp),public           :: GP(KBPROJmax,namax,0:lamax)


!..Pseudopotential type (1=BHSC,2=BHSA,3=CGEN,4=GRID,5=KBlike)
    integer,public           :: pptype
    character*50,public      :: ppstring

!------------------------------------------------------------------------------!
  end module pspot
!------------------------------------------------------------------------------!
