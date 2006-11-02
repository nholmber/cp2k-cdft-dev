!------------------------------------------------------------------------------!
  module gdiis_save
!------------------------------------------------------------------------------!

  USE basic_data_types, ONLY: dp

  real(dp),allocatable :: xbuf(:,:),gbuf(:,:),gmax(:),gnorm(:),bc(:,:),vc(:)

end module gdiis_save

!------------------------------------------------------------------------------!
  subroutine gdiis(x,g,h0,n,mdiis,init)
!------------------------------------------------------------------------------!
        !     GDIIS/LBFGS method                                               
        !     T.H.Fischer and J.Almlof J.Phys.Chem. 96, 9768, (1992)           
        !                                                                      
        !     On Input : x(*) is the current position vector                   
        !                g(*) is the current gradient vector                   
        !                h0(*) is the approximation to the inverse Hessian     
        !                n is the number of parameters                         
        !                mdiis the maximum number of DIIS vectors and BFGS     
        !                      history                                         
        !                init asks for a new start for the DIIS history        
        !                                                                      
        !     On Output : x(*) is the new parameter vector                     
        !                 g(*) is the negative displacement vector             
        !                                                                      
        !----------------------------------------------------------------------!

  USE basic_data_types, ONLY: dp
  use gdiis_save, only: xbuf,gbuf,gmax,gnorm,bc,vc
  implicit none
!in:
  integer, intent(in) :: n,mdiis
!inout:
  real(dp)              :: x(n),g(n),h0(n)
  logical                 :: init
!locals:
  integer,parameter  :: maxdis=20
  integer,save       :: now=0,ndiis=0,nbfgs=-1
  integer            :: i,ii,j,jj,nst
  real(dp)              :: r1,r2,rt,eold,enew
!------------------------------------------------------------------------------!

  if (.not.allocated( xbuf )) then
    allocate( xbuf(n,mdiis) )
    allocate( gbuf(n,mdiis) )
    allocate( gmax(mdiis) )
    allocate( gnorm(mdiis) )
    allocate( bc(mdiis+1,mdiis+1) )
    allocate( vc(mdiis+1) )
  endif

!reset DIIS history
  if (init) then
    ndiis=0
    eold=enew
    init=.false.
  endif
!update counters
  call ucount(now,ndiis,nbfgs,mdiis)
!update buffers
  xbuf(:,now)=x
  gbuf(:,now)=g
  gmax(now)=maxval(abs(g))
  gnorm(now)=sqrt(dot_product(g,g) / real(n,dp) )
!check for progress
  if (ndiis>=3) then
    rt=1._dp-ndiis/100._dp
    nst=now-ndiis+1
    if (nst<=0) nst=nst+mdiis
    r1=gmax(now)/gmax(nst)
    r2=gnorm(now)/gnorm(nst)
    if (r1>rt .and. r2>rt) ndiis=1
  endif
!setup DIIS matrix
  do i=1,ndiis
    ii=now-i+1
    if (ii<=0) ii=ii+mdiis
    x=gbuf(:,ii)
    call lbfgs(x,xbuf,gbuf,h0,n,n,nbfgs,now,mdiis)
    bc(i,i)=dot_product(x,x)
    do j=i+1,ndiis
      jj=now-j+1
      if(jj<=0) jj=jj+mdiis
      g=gbuf(:,jj)
      call lbfgs(g,xbuf,gbuf,h0,n,n,nbfgs,now,mdiis)
      bc(i,j)=dot_product(x,g)
      bc(j,i)=bc(i,j)
    enddo
    bc(ndiis+1,i)=-1._dp
    bc(i,ndiis+1)=-1._dp
    vc(i)=0._dp
  enddo
  bc(ndiis+1,ndiis+1)=0._dp
  vc(ndiis+1)=-1._dp
!solve the linear system:
  call ssolv(bc,mdiis+1,ndiis+1,vc)
!compute the interpolated coefficient and gradient vectors:
  x=0._dp
  g=0._dp
  do i=1,ndiis
    ii=now-i+1
    if (ii<=0) ii=ii+mdiis
    x=x+vc(i)*xbuf(:,ii)
    g=g+vc(i)*gbuf(:,ii)
  enddo
!estimate the new parameter vector:
  call lbfgs(g,xbuf,gbuf,h0,n,n,nbfgs,now,mdiis)
  x=x-g


contains

  !----------------------------------------------------------------------------!
    subroutine ucount(now,ndiis,nbfgs,mdiis)
  !----------------------------------------------------------------------------!
    implicit none
  !in:
    integer,intent(in):: mdiis
  !inout:
    integer           :: now,ndiis,nbfgs
  !----------------------------------------------------------------------------!

    now=now+1
    if (now>mdiis) now=now-mdiis
    ndiis=ndiis+1
    if (ndiis>mdiis) ndiis=mdiis
    nbfgs=nbfgs+1
    if (nbfgs>=mdiis) nbfgs=mdiis-1

  end subroutine ucount

    !--------------------------------------------------------------------------!
      subroutine ssolv(b,ldb,ndim,v)
    !--------------------------------------------------------------------------!
      USE basic_data_types, ONLY: dp
      implicit none
    !in:
      integer     ,intent(in):: ldb,ndim
    !inout:
      real(dp)             :: b(ldb,ldb)
      real(dp)             :: v(ldb)
    !locals:
      integer                :: rank,laux,info
      real(dp),parameter   :: toleig=0.2e-15_dp
!MK   real(dp)             :: sv(ndim),aux(4*ndim+4)
      real(dp)             :: sv(ndim),aux(5*ndim)
    !--------------------------------------------------------------------------!
!MK   laux=4*ndim+4
      laux = 5*ndim
      CALL dgelss(ndim,ndim,1,b,ldb,v,ldb,sv,toleig,rank,aux,laux,info)

    end subroutine ssolv

end subroutine gdiis
