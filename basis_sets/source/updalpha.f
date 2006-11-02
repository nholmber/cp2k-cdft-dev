      SUBROUTINE updalpha(dalpha)
      USE basic_data_types, ONLY: dp
      use atom
      use upd
      implicit none
      REAL(dp) :: dalpha(namax,0:lamax)
      logical firstopt,hout,check,updd
c..local
      integer i,j,k,l
c
      call alpha2beta
c
c..new search direction
      do i=1,nbeta
         xi(i)=0.D0
         do j=1,nbeta
            xi(i)=xi(i)-hessinv(i,j)*gbeta(j)
         enddo
      enddo
c
c..transform xi to dalpha
      do l=0,lmax
         do i=1,nalpha(l)
            if (alpp(i,l).eq.0) then
               dalpha(i,l)=0.D0
            else
               dalpha(i,l)=xi(alpp(i,l))
            endif
         enddo
      enddo
c
      return
      end
c
c     ------------------------------------------------------------------
      SUBROUTINE updalpha_eps(epslimit,converged,who,dalpha)
      USE basic_data_types, ONLY: dp
      USE atom
      USE upd
      IMPLICIT NONE
      REAL(dp) :: dalpha(namax,0:lamax),epslimit
      logical firstopt,hout,check,updd,converged
      integer who
c..   local
      integer i,j,k,l
c
      call alpha2beta
c
c..   new search direction
      converged=.true.
      who=800000000
      do i=1,nbeta
         xi(i)=0.D0
         do j=1,nbeta
            xi(i)=xi(i)-hessinv(i,j)*gbeta(j)
         enddo
         if (abs(gbeta(i))>epslimit) then
            who=who+10**(nbeta-i)
            converged=.false.
         else
            if (converged) then
               xi(i)=0.D0
            else
               who=who+10**(nbeta-i)
            endif
         endif
      enddo
c
c..   transform xi to dalpha
      do l=0,lmax
         do i=1,nalpha(l)
            if (alpp(i,l).eq.0) then
               dalpha(i,l)=0.D0
            else
               dalpha(i,l)=xi(alpp(i,l))
            endif
         enddo
      enddo
c
      return
      end
c
c     ------------------------------------------------------------------
      SUBROUTINE updalpha_allin1(converged)
      USE basic_data_types, ONLY: dp
      USE atom
      USE upd
      IMPLICIT NONE
      REAL(dp) :: dalpha(namax,0:lamax),epslimit,factor,dfactor
      logical firstopt,hout,check,updd,converged
      integer i,j,k,l
c
      call alpha2beta
c
      dfactor=0.d0
      do i=1,nbeta
         dfactor=dfactor+gbeta(i)*beta(i)
      enddo
      if (abs(dfactor)>0.1d0) dfactor=sign(0.1d0,dfactor)
      factor=1.d0-dfactor
      print*,'factor =',factor
      converged=abs(1.d0-factor)<1e-5
c
c..   transform xi to dalpha
      do l=0,lmax
         do i=1,nalpha(l)
            alpha(i,l)=alpha(i,l)*factor
         enddo
      enddo
c
      return
      end
c
c     ------------------------------------------------------------------
      subroutine updhess
      USE basic_data_types, ONLY: dp
      use atom
      use upd
      implicit none
c..   local
      REAL(dp) :: xii(namax*(lamax+1))
      REAL(dp) :: dg(namax*(lamax+1))
      REAL(dp) :: hdg(namax*(lamax+1))
      REAL(dp) :: sumdg,sumxi,fac,fae,fad,xtg
      integer i,j,k,l
c
c..   update hessian
      call alpha2beta
      xii=beta-betaold
      dg=gbeta-gbetaold
      do i=1,nbeta
         hdg(i)=0.D0
         do j=1,nbeta
            hdg(i)=hdg(i)+hessinv(i,j)*dg(j)
         enddo
      enddo
      fac=0.D0
      fae=0.D0
      sumdg=0.D0
      sumxi=0.D0
      do i=1,nbeta
         fac=fac+dg(i)*xii(i)
         fae=fae+dg(i)*hdg(i)
         sumdg=sumdg+dg(i)**2
         sumxi=sumxi+xii(i)**2
      enddo
      if (fac**2.gt.3.D-8*sumdg*sumxi) then
c     write (*,*) 'updating Hessian...'
         fac=1.D0/fac
         fad=1.D0/fae
         do i=1,nbeta
            dg(i)=fac*xii(i)-fad*hdg(i)
         enddo
         do i=1,nbeta
            do j=1,nbeta
               hessinv(i,j)=hessinv(i,j)+fac*xii(i)*xii(j)
     &              -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
            enddo
         enddo
      else
         write (*,*) 'cannot update Hessian.'
      endif
c
      betaold=beta
      gbetaold=gbeta
c
      return
      end
c
c     ------------------------------------------------------------------
      subroutine alpha2beta
      use atom
      use upd
      implicit none
      integer i,l
c..   project exponents and gradient on selected optimization subspace
      nbeta=0
      gbeta=0.D0
      do l=0,lmax
         do i=1,nalpha(l)
            nbeta=max(alpp(i,l),nbeta)
            if (alpp(i,l).ne.0) then
               beta(alpp(i,l))=alpha(i,l)
               gbeta(alpp(i,l))=gbeta(alpp(i,l))+galpha(i,l)
            endif
         enddo
      enddo
c
      return
      end
c
c     ------------------------------------------------------------------
      subroutine beta2dalpha(dalpha)
      USE basic_data_types, ONLY: dp
      use atom
      use upd
      implicit none
      integer i,l
      REAL(dp) :: dalpha(namax,0:lamax)
c
c..   project exponents and gradient
      dalpha=0.D0
      do l=0,lmax
         do i=1,nalpha(l)
            if (alpp(i,l).ne.0) then
               dalpha(i,l)=-gbeta(alpp(i,l))
            endif
         enddo
      enddo
c
      return
      end
c
c     ------------------------------------------------------------------
      subroutine hessunit
      use atom
      use upd
      implicit none
c..   local
      integer j
c
      call alpha2beta
      hessinv=0.D0
c
      write (*,*) 'initializing Hessian as diagonal matrix...'
c
      do j=1,nbeta
         hessinv(j,j)=1.D0*beta(j)**2
      enddo
c
      return
      end
c
c     ------------------------------------------------------------------
      subroutine hessinit(cmat)
      USE basic_data_types, ONLY: dp
      use atom
      use upd
      implicit none
      REAL(dp) :: dalpha(namax,0:lamax)
      REAL(dp) :: cmat(namax,namax,0:lamax)
      REAL(dp) :: ksener(namax,0:lamax)
      logical firstopt,posdef
c..   local
      REAL(dp) :: gbeta1(namax*(lamax+1))
      REAL(dp) :: gbeta2(namax*(lamax+1))
      REAL(dp) :: alphasave(namax,0:lamax)
      REAL(dp) :: galphasave(namax,0:lamax)
      REAL(dp) :: betasave(namax*(lamax+1))
      REAL(dp) :: gbetasave(namax*(lamax+1))
      REAL(dp) :: hessian(namax*(lamax+1),namax*(lamax+1))
      REAL(dp) :: w(namax*(lamax+1),namax*(lamax+1))
      REAL(dp) :: hess(namax*(lamax+1)*namax*(lamax+1))
      REAL(dp) :: hessev(namax*(lamax+1))
      REAL(dp) :: rcond,det(2),aux(200*nbeta),delta,etot,
     &       work(namax*(lamax+1))
      integer i,j,k,l,naux,info,ipiv(namax*(lamax+1))
      logical firstloop
c
      firstloop=.false.
      alphasave=alpha
      galphasave=galpha
      betasave=beta
      gbetasave=gbeta
      hessinv=0.D0
      hessian=0.D0
c
      write (*,*) 'calculating new Hessian...'
c
      do j=1,nbeta
         delta=1.D-3*beta(j)
         beta(j)=beta(j)+0.5D0*delta
c..   beta->alpha
         do l=0,lmax
            do i=1,nalpha(l)
               if (alpp(i,l).ne.0) alpha(i,l)=beta(alpp(i,l))
            enddo
         enddo
c..   gradient
         call wfn_ortho(cmat)
         call lda_scf(cmat,firstloop,ksener,etot)
         call exp_grad (cmat,ksener)
c..   galpha->gbeta1
         call alpha2beta
         gbeta1=gbeta
c
         beta(j)=beta(j)-delta
c..   beta->alpha
         do l=0,lmax
            do i=1,nalpha(l)
               if (alpp(i,l).ne.0) alpha(i,l)=beta(alpp(i,l))
            enddo
         enddo
c..   gradient
         call wfn_ortho(cmat)
         call lda_scf(cmat,firstloop,ksener,etot)
         call exp_grad (cmat,ksener)
c..   galpha->gbeta2
         call alpha2beta
         gbeta2=gbeta
c..   Hessian
         do i=1,nbeta
            hessian(i,j)=(gbeta1(i)-gbeta2(i))/delta
         enddo
c
         beta=betasave
         gbeta=gbetasave
      enddo
c..   restore alpha,galpha
      alpha=alphasave
      galpha=galphasave
c..   symmetrize Hessian
      do i=1,nbeta
         do j=i,nbeta
            hessian(i,j)=(hessian(i,j)+hessian(j,i))/2.D0
            hessian(j,i)=hessian(i,j)
         enddo
      enddo
c..   invert Hessian
      naux=200*nbeta
CMK   call dgeicd (hessian,namax*(lamax+1),nbeta,2,rcond,det,aux,naux)
      CALL dgetrf(nbeta,nbeta,hessian,namax*(lamax+1),ipiv,info)
      CALL dgetri(nbeta,hessian,namax*(lamax+1),ipiv,work,
     &            namax*(lamax+1),info)
      hessinv=hessian
c..   check positive definiteness
      call matcheck(hessinv,namax*(lamax+1),nbeta,posdef)
      if (.not.posdef) then
         write (*,*) 'Hessian has nonpositive eigenvalues !! ...',
     &        'corrected hessinv.'
      endif
c
      write (*,*) '...'
c
      return
      end
c
c     ------------------------------------------------------------------
      subroutine matcheck(mat,lda,dim,posdef)
      USE basic_data_types, ONLY: dp
      implicit none
      logical posdef
      integer lda,dim
      REAL(dp) :: mat(lda,*)
c..   local
      REAL(dp) :: matpack(lda*lda)
      REAL(dp) :: eval(dim)
      REAL(dp) :: evec(dim,dim)
      REAL(dp) :: matev(dim,dim)
      REAL(dp) :: mat1(dim,dim)
      REAL(dp) :: aux(3*lda)
      integer i,j,k,l,naux,info
      logical firstloop
c
      naux=3*lda
c
c..   pack mat
      do i=1,dim
         do j=i,dim
            matpack(i+((j*(j-1))/2))=mat(i,j)
         enddo
      enddo
c..   calculate eigenvectors and eigenvalues of mat
      call dspev("V","U",dim,matpack,eval,evec,dim,aux,info)
      if (info /= 0) STOP "dspev: info /= 0"
c..   check and correct eigenvalues
      posdef=.true.
      do i=1,dim
         if (eval(i).le.0) then
            posdef=.false.
            print*,"matcheck:",i,eval(i)
            eval(i)=abs(eval(i))
         endif
      enddo
      if (posdef) return
c
c..   construct eigenvalue matrix
      matev=0.D0
      do i=1,dim
         matev(i,i)=eval(i)
      enddo
c..   reconstruct corrected mat
c..   mat1 = matev evecT
      call dgemm('N','T',dim,dim,dim,1.D0,matev,dim,evec,dim,0.D0,
     &     mat1,dim)
c..   mat = evec mat1
      call dgemm('N','N',dim,dim,dim,1.D0,evec,dim,mat1,dim,0.D0,
     &     mat,lda)
c

      end
