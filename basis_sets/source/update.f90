      SUBROUTINE update(icall,ksmat,smat,cmat,pmat)
      USE basic_data_types, ONLY: dp
      USE atom
      IMPLICIT NONE
      INTEGER, PARAMETER :: nsdiis = 10000
      INTEGER :: icall
      REAL(dp) :: ksmat(namax,namax,0:lamax),&
                  smat(namax,namax,0:lamax),&
                  cmat(namax,namax,0:lamax),&
                  pmat(namax,namax,0:lamax)
!local
      REAL(dp) :: xmat(namax,namax)
      REAL(dp) :: zmat(namax,namax)
      REAL(dp) :: eemat(namax*namax*(lamax+1))
      REAL(dp) :: pvmat(namax*namax*(lamax+1))
      REAL(dp) :: eesto(nsdiis)
      REAL(dp) :: pvsto(nsdiis),eenow,eemax
      INTEGER :: ipar,i,j,l,n
      save eesto,pvsto

! error matrix
      ipar=0
      do l=0,lmax
        n=nalpha(l)
        call dgemm('N','N',n,n,n,1.d0,ksmat(1,1,l),namax,pmat(1,1,l),namax,0.0d0,xmat,namax)
        call dgemm('N','N',n,n,n,1.d0,xmat,namax,smat(1,1,l),namax,0.0d0,zmat,namax)
        call dgemm('N','N',n,n,n,1.d0,smat(1,1,l),namax,pmat(1,1,l),namax,0.0d0,xmat,namax)
        call dgemm('N','N',n,n,n,1.d0,xmat,namax,ksmat(1,1,l),namax,-1.0d0,zmat,namax)
        eemax=0.0d0
        do i=1,n
          do j=1,n
            ipar=ipar+1
            eemat(ipar)=zmat(i,j)
            eenow=abs(eemat(ipar))
            if(eenow.gt.eemax) eemax=eenow
            pvmat(ipar)=ksmat(i,j,l)
          enddo
        enddo
      enddo
      dmix1=0.0d0
      if(eemax.gt.0.1) dmix1=dmix
      call ksdiis(pvmat,eemat,eesto,pvsto,icall,nsdiis,ipar)
      ipar=0
      do l=0,lmax
        n=nalpha(l)
        do i=1,n
          do j=1,n
            ipar=ipar+1
            ksmat(i,j,l)=pvmat(ipar)
          enddo
        enddo
      enddo
      end

      SUBROUTINE ksdiis(pv,ee,eesto,pvsto,icall,nsdiis,n)
      USE basic_data_types, ONLY: dp
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxdiis = 10
      INTEGER :: icall,n,nsdiis
      REAL(dp) :: pv(*),ee(*),eesto(n,*),pvsto(n,*)
      REAL(dp) :: v(maxdiis+1),dm(maxdiis+1,maxdiis+1)
      REAL(dp) :: s1(maxdiis+1,maxdiis+1),s2(maxdiis+1,maxdiis+1)
      REAL(dp) :: toleig,ddot
      INTEGER :: i,j,m1,isto,mdiis,nx,nsize,nr,n1,lw,info
      m1=maxdiis+1
      nx=nsdiis/n
      mdiis=min(maxdiis,nx)
      isto=mod(icall-1,mdiis)+1
      call dcopy(n,pv,1,pvsto(1,isto),1)
      call dcopy(n,ee,1,eesto(1,isto),1)
      nsize=icall
      if(icall.gt.mdiis) nsize=mdiis
      do i=1,nsize
        do j=1,i
          dm(i,j)=ddot(n,eesto(1,i),1,eesto(1,j),1)
          dm(j,i)=dm(i,j)
        enddo
        v(i)=0.0d0
        dm(i,nsize+1)=-1.d0
        dm(nsize+1,i)=-1.d0
      enddo
      v(nsize+1)=-1.d0
      dm(nsize+1,nsize+1)=0.0d0
      toleig=0.2d-15
      lw=m1*m1
      n1=nsize+1
      call dgelss(n1,n1,1,dm,m1,v,m1,s1,toleig,nr,s2,lw,info)
      do i=1,n
        pv(i)=0.0d0
      enddo
      do i=1,nsize
        do j=1,n
          pv(j)=pv(j)+v(i)*pvsto(j,i)
        enddo
      enddo

      END
