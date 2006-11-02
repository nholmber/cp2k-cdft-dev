      SUBROUTINE expdiis(cmat,etot,iter)
      USE basic_data_types, ONLY: dp
      USE atom
      USE upd
      IMPLICIT NONE
      REAL(dp) etot
      REAL(dp) cmat(namax,namax,0:lamax)
      REAL(dp) ksener(namax,0:lamax)
      integer iter,i,l,mxdis,nh
      parameter(mxdis=10)
      REAL(dp) e2
      REAL(dp) cmat2(namax,namax,0:lamax)
      REAL(dp) alphaold(namax,0:lamax)
      REAL(dp) keepa(mxdis*namax*lamax),keepg(mxdis*namax*lamax)
      REAL(dp) error(mxdis*namax*lamax)
      save keepa,keepg

!     write(*,*) ' try DIIS step:'
      cmat2=cmat
      call wfn_ortho(cmat2)
      call lda_scf(cmat2,.false.,ksener,etot)
      call exp_grad (cmat2,ksener)
      call alpha2beta
      nh=namax*(lamax+1)
!..diis routine
      call ediis(iter,nbeta,beta,gbeta,error,hessinv,keepa,keepg,nh,mxdis)
!.backtransform
      alphaold=alpha
      do l=0,lmax
        do i=1,nalpha(l)
          if (alpp(i,l).ne.0) then
            alpha(i,l)=beta(alpp(i,l))
          endif
        enddo
      enddo
!..check for energy increase
      cmat2=cmat
      call wfn_ortho(cmat2)
      call lda_scf(cmat2,.false.,ksener,e2)
      if(e2.gt.etot) then
        alpha=alphaold
        write(*,*) ' DIIS UPDATE REJECTED ',e2-etot
      else
        write(*,*) ' DIIS UPDATE ACCEPTED ',e2-etot
        call exp_grad(cmat,ksener)
        etot=e2
        cmat=cmat2
      endif
      end

      subroutine ediis(icall,n,beta,gbeta,error,hessinv,keepa,keepg,nh,m)
      USE basic_data_types, ONLY: dp
      IMPLICIT NONE
      INTEGER icall,n,nh,m,n1
      REAL(dp) beta(n),gbeta(n),hessinv(nh,nh)
      REAL(dp) keepa(n,m),keepg(n,m)
      REAL(dp) error(n,m)
      !..local
      integer isto,nsize,i,j,info,nr
      REAL(dp) toleig,ddot
      REAL(dp) dm(20,20),v(20),s1(20),s2(400)

      if(m.gt.19) stop 'ediis'
      isto=mod(icall-1,m)+1
      call dcopy(n,beta,1,keepa(1,isto),1)
      call dcopy(n,gbeta,1,keepg(1,isto),1)
      nsize=icall
      if(icall.gt.m) nsize=m
      do i=1,nsize
        call dgemv('N',n,n,-1.0_dp,hessinv,nh,keepg(1,i),1,0.0_dp,error(1,i),1)
      enddo
      do i=1,nsize
        do j=1,i
          dm(i,j)=ddot(n,error(1,i),1,error(1,j),1)
          dm(j,i)=dm(i,j)
        enddo
        v(i)=0.0_dp
        dm(i,nsize+1)=-1.0_dp
        dm(nsize+1,i)=-1.0_dp
      enddo
      v(nsize+1)=-1.0_dp
      dm(nsize+1,nsize+1)=0.0_dp
      toleig=0.2d-15
      n1=nsize+1
      call dgelss(n1,n1,1,dm,20,v,20,s1,toleig,nr,s2,400,info)
      beta=0.0_dp
      gbeta=0.0_dp
      do i=1,nsize
        do j=1,n
          beta(j)=beta(j)+v(i)*keepa(j,i)
          gbeta(j)=gbeta(j)+v(i)*keepg(j,i)
        enddo
      enddo
      call dgemv('N',n,n,-1.0_dp,hessinv,nh,gbeta,1,0.0_dp,s1,1)
      do i=1,n
        if(abs(s1(i)).lt.(0.1*beta(i))) beta(i)=beta(i)+s1(i)
      enddo
      end
