      subroutine numd_wfn(cmat)
      use atom
      use rint
      use pspot
      implicit none
      logical firstloop
      real*8 cmat(namax,namax,0:lamax)
c..local arrays
      real*8 cdmat(namax,namax,0:lamax,-2:2)
      real*8 dermat(namax,0:lamax,3)
      real*8 smat(namax,namax,0:lamax)
      real*8 pmat(namax,namax,0:lamax)
      real*8 ksener(namax,0:lamax)
      real*8 eone,ekin,epot,etot,etot2,exc,edif,etoto,ecb,trace
      real*8 rmax,c,emax,emax1,dfn,dd,over,ff
      integer ldval,ndval
      integer k,l,m,i,j,n,ii,iter,lnn
c
      firstloop=.True.
      call lda_scf(cmat,firstloop,ksener,etot)
      firstloop=.False.
c..overlap matrix
      call overlap(smat)
c..HOMO
      ldval=0
      ndval=0
      emax=-9999.
      do l=0,lmax
        do n=1,nocc(l)
          emax1=ksener(n,l)
          if(emax1.gt.emax) then
            ldval=l
            ndval=n
            emax=emax1
          endif
        enddo
      enddo
      write(6,'(A,2i5,2F10.6)') ' derivative wrt occupation ',
     *           ldval,ndval,emax,occ(ndval,ldval)
c..finite differences
      lnn=namax*namax*(lamax+1)
      call dcopy(lnn,cmat,1,cdmat(1,1,0,0),1)
      dfn=0.005
      occ(ndval,ldval)=occ(ndval,ldval)-2*dfn
	write (*,*) '4'
      call lda_scf(cmat,firstloop,ksener,etot)
      call dcopy(lnn,cmat,1,cdmat(1,1,0,-2),1)
      occ(ndval,ldval)=occ(ndval,ldval)+dfn
	write (*,*) '3'
      call lda_scf(cmat,firstloop,ksener,etot)
      call dcopy(lnn,cmat,1,cdmat(1,1,0,-1),1)
      occ(ndval,ldval)=occ(ndval,ldval)+2*dfn
	write (*,*) '2'
      call lda_scf(cmat,firstloop,ksener,etot)
      call dcopy(lnn,cmat,1,cdmat(1,1,0,1),1)
      occ(ndval,ldval)=occ(ndval,ldval)+dfn
	write (*,*) '1'
      call lda_scf(cmat,firstloop,ksener,etot)
      call dcopy(lnn,cmat,1,cdmat(1,1,0,2),1)
      occ(ndval,ldval)=occ(ndval,ldval)-2*dfn
      call dcopy(lnn,cdmat(1,1,0,0),1,cmat,1)
c..derivatives
      dd=(2*ldval+1)*dfn
      do l=0,lmax
        n=nocc(l)
c..first
        do i=1,nalpha(l)
          dermat(i,l,1)=(cdmat(i,n,l,1)-cdmat(i,n,l,-1))/(2*dd)
        enddo
c..second
        do i=1,nalpha(l)
          dermat(i,l,2)=(cdmat(i,n,l,2)-2*cdmat(i,n,l,0)+
     *                   cdmat(i,n,l,-2))/(4*dd*dd)
        enddo
c..third
        do i=1,nalpha(l)
          dermat(i,l,3)=(cdmat(i,n,l,2)-2*cdmat(i,n,l,1)+
     *           2*cdmat(i,n,l,-1)-cdmat(i,n,l,-2))/(2*dd*dd*dd)
        enddo
c..normalization
        call gsortho(cmat(1,1,l),dermat(1,l,1),smat(1,1,l),
     *               nalpha(l),n,namax)
        call dcopy(nalpha(l),dermat(1,l,1),1,cmat(1,n+1,l),1)
        call gsortho(cmat(1,1,l),dermat(1,l,2),smat(1,1,l),
     *               nalpha(l),n+1,namax)
        call dcopy(nalpha(l),dermat(1,l,2),1,cmat(1,n+2,l),1)
        call gsortho(cmat(1,1,l),dermat(1,l,3),smat(1,1,l),
     *               nalpha(l),n+2,namax)
        call dcopy(nalpha(l),dermat(1,l,3),1,cmat(1,n+3,l),1)
cdeb        write(6,*) l,n

cdeb	write(6,*) ' exponent:   |       1s           2s           3s',
cdeb     &			'           4s           5s'
cdeb	write(6,*) ' exponent:   |       1p           2p           3p',
cdeb     &			'           4p'
cdeb	write(6,*) '--------------------------------------------------',
cdeb     &			'------------------------------'
cdeb        do i=1,nalpha(l)
cdeb          write(6,'(i5,10F14.8)') i,(cmat(i,j,l),j=1,n+3)
cdeb          write(6,'(F13.8,A,10F13.8)') alpha(i,l),' |',
cdeb     &						(cmat(i,j,l),j=1,n+3)
cdeb        enddo
        do i=n+4,nalpha(l)
          do j=1,nalpha(l)
            cmat(j,i,l)=0.0d0
          enddo
        enddo
      enddo
      return
      end
c
      subroutine gsortho(c,v,s,n,ns,namax)
      implicit none
      integer n,ns,namax,i
      real*8 c(namax,*),v(*),s(namax,*)
      real*8 ff,over
      do i=1,ns
        ff=over(c(1,i),v,s,n,namax)
        call daxpy(n,-ff,c(1,i),1,v,1)
      enddo
      ff=over(v,v,s,n,namax)
      ff=1./sqrt(ff)
      call dscal(n,ff,v,1)
      return
      end
      function over(c1,c2,s,n,namax)
      implicit none
      integer n,namax,i,j
      real*8 c1(*),c2(*),s(namax,*),over
      over=0.0
      do i=1,n
        do j=1,n
          over=over+c1(i)*s(i,j)*c2(j)
        enddo
      enddo
      return
      end
