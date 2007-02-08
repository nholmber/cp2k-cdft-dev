C     ==================================================================
      subroutine kbppfit(cmat)
C     ==================================================================
      use atom
      use pspot
      implicit none
      real*8 nn(namax,namax,0:lamax)
      real*8 cmat(namax,namax,0:lamax)
      real*8 umat(namax,namax,0:lamax)
      real*8 A(namax,0:lamax)
      integer nvkb(namax,0:lamax)
      real*8 avkb(namax,0:lamax)
      real*8 Vl(0:lamax)
      integer Vkbnr(0:lamax)
      common /Vkbpara/ A,nvkb,avkb,Vkbnr,Vl
C.....locals
      real*8 Cvpsi(2*namax*EXPmax,0:lamax)
      real*8 avpsi(2*namax*EXPmax,0:lamax)
      integer nvpsi(2*namax*EXPmax,0:lamax)
      integer vpsinr(0:lamax)
      common /vpsipara/ Cvpsi,avpsi,nvpsi,vpsinr
      real*8	G(namax,0:lamax),H(namax,namax,0:lamax)
      real*8	gint,amin,rmax,step,f(2*lamax)
      real*8	VPSIR,Vkbr,r,v1,Vl2(0:lamax)
      integer i,j,k,l,lambda,kappa,info,ir
      integer ipvt(namax)
C     ------------------------------------------------------------------
      call calcnn(nn)
C
C.....The exponents of the Gaussians and the exponents of r in the
C.....projectors are chosen to be the ones used for the pseudopotential,
C.....but with the Gaussian exponents multiplied by a factor.
C
      do l=0,lmax-1
         Vkbnr(l)=EXPnr(l)
         do i=1,EXPnr(l)
	    nvkb(i,l)=PPRexp(i,l)
	    avkb(i,l)=1.3*PPNexp(i,l)
         enddo
      enddo
      Vkbnr(lmax)=1
      nvkb(1,lmax)=0
      avkb(1,lmax)=0.d0
      A(1,lmax)=0.d0
C
C.....Initialize the coefficients for the projector to fit
C
      call vpsiinit(nn,cmat)
C
C.....Optimize the coefficients A for given X (least square fit)
C.....by computing a solution of  G = H A
C
C.....First compute G
C
      do l=0,lmax-1
         do k=1,Vkbnr(l)
	    g(k,l)=0.d0
	    do i=1,vpsinr(l)
               g(k,l)=g(k,l)+Cvpsi(i,l)*
     &              gint(2+nvkb(k,l)+nvpsi(i,l),avkb(k,l)+avpsi(i,l))
	    enddo
	    g(k,l)=2*g(k,l)
         enddo
      enddo
C
C.....Next compute H
C
      do l=0,lmax-1
         do k=1,Vkbnr(l)
	    do j=k,Vkbnr(l)
              h(k,j,l)=2*gint(2+nvkb(k,l)+nvkb(j,l),avkb(k,l)+avkb(j,l))        
              h(j,k,l)=h(k,j,l)
	    enddo
         enddo
      enddo
C     
C.....Now solve  H A = G  for A
C
      do l=0,lmax-1
!MK         call dgef(h(1,1,l),namax,Vkbnr(l),ipvt)
!MK         call dges(h(1,1,l),namax,Vkbnr(l),ipvt,g(1,l),0)
         CALL dgetrf(Vkbnr(l),Vkbnr(l),h(1,1,l),namax,ipvt,info)
         CALL dgetrs("N",Vkbnr(l),1,h(1,1,l),namax,ipvt,g(1,l),namax,
     &               info)
      enddo
      A=g
C
C.....Calculate the matrix elements <PSI(1,l)|Vkb(1,l)>
C
      do l=0,lmax-1
         Vl(l)=0.d0
         do kappa=1,nalpha(l)
	    v1=0.d0
	    do i=1,Vkbnr(l)
               v1=v1+A(i,l)*gint(l+2+nvkb(i,l),alpha(kappa,l)+avkb(i,l))
	    enddo
	    Vl(l)=Vl(l)+v1*sqrt(nn(kappa,kappa,l))*cmat(kappa,1,l)
         enddo
         Vl(l)=1.d0/Vl(l)
      enddo
c
C.....Calculate the matrix elements <PSI(1,l)|Vps(l)*PSI(1,l)>
C
      do l=0,lmax-1
         Vl2(l)=0.d0
         do kappa=1,nalpha(l)
	    v1=0.d0
	    do i=1,Vpsinr(l)
               v1=v1+Cvpsi(i,l)
     &              *gint(l+2+nvpsi(i,l),alpha(kappa,l)+avpsi(i,l))
	    enddo
	    Vl2(l)=Vl2(l)+v1*sqrt(nn(kappa,kappa,l))*cmat(kappa,1,l)
         enddo
         Vl2(l)=1.d0/Vl2(l)
      enddo
c
C.....Output
C
C.....Projector plots
c.....determine smallest exponent
      amin=avkb(1,0)
      do l=0,lmax-1
         do i=1,Vkbnr(l)
	    if (avkb(i,l).lt.amin) then
               amin=avkb(i,l)
	    endif
         enddo
      enddo
      rmax=sqrt(20.d0/amin)
      step=rmax/300.d0
c
      open(16,FILE='atom.proj')
      do ir=0,int(rmax/step)
         r=ir*step
         i=1
         do l=0,lmax-1
	    f(i)=VPSIR(l,r,cmat)
	    f(i+1)=Vkbr(l,r)
	    i=i+1
         enddo
         write (16,*) r,(f(i),i=1,2*lmax)
      enddo
      close(16)
      write (*,*) 'Created File "Atom.Proj".'
C
C.....Generate input file
C
      Open (17,FILE='atom.kbpp')
      write (17,*) '&PSEUDO'
      write (17,*) 'KBPP'
      write (17,*) 'ZEFF'
      write (17,'(1G24.17)') Zeff
      write (17,*) 'ERFC'
      write (17,*)  ERFnr
      write (17,'(3G24.17)') (PPNerf(j),j=1,ERFnr)
      write (17,'(3G24.17)') (PPLerf(j),j=1,ERFnr)
      write (17,*) 'GCFS'
      write (17,*)  EXPnr(lmax)
      write (17,'(3G24.17)') (PPNexp(j,lmax),j=1,EXPnr(lmax))
      write (17,'(3I16)') (PPRexp(j,lmax),j=1,EXPnr(lmax))
      write (17,'(3G24.17)') (PPLexp(j,lmax),j=1,EXPnr(lmax))
      write (17,*) 'PPLM'
      write (17,*) lmax
      write (17,*) 'NLCS'
      write (17,*) '1',Vkbnr(0)
      write (17,'(1G24.17)') Vl(0)
      write (17,'(3G24.17)') (avkb(j,0),j=1,Vkbnr(0))
      write (17,'(3I16)') (nvkb(j,0),j=1,Vkbnr(0))
      write (17,'(3G24.17)') (A(j,0),j=1,Vkbnr(0))
      if (lmax.ge.1) then
         write (17,*) 'NLCP'
         write (17,*) '1',Vkbnr(1)
         write (17,'(1G24.17)') Vl(1)
         write (17,'(3G24.17)') (avkb(j,1),j=1,Vkbnr(1))
         write (17,'(3I16)') (nvkb(j,1),j=1,Vkbnr(1))
         write (17,'(3G24.17)') (A(j,1),j=1,Vkbnr(1))
      endif
      if (lmax.ge.2) then
         write (17,*) 'NLCD'
         write (17,*) '1',Vkbnr(2)
         write (17,'(1G24.17)') Vl(2)
         write (17,'(3G24.17)') (avkb(j,2),j=1,Vkbnr(2))
         write (17,'(3I16)') (nvkb(j,2),j=1,Vkbnr(2))
         write (17,'(3G24.17)') (A(j,2),j=1,Vkbnr(2))
      endif
      if (lmax.ge.3) then
         write (17,*) 'NLCF'
         write (17,*) '1',Vkbnr(3)
         write (17,'(1G24.17)') Vl(3)
         write (17,'(3G24.17)') (avkb(j,3),j=1,Vkbnr(3))
         write (17,'(3I16)') (nvkb(j,3),j=1,Vkbnr(3))
         write (17,'(3G24.17)') (A(j,3),j=1,Vkbnr(3))
      endif
      write (17,*) '&END'
      write (17,*) ' '
      write (17,*) ' '
      write (17,*) ' '
C
C.....Generate alternative input file
C
      open (17,FILE='atom.kbpp')
      write (17,*) '&PSEUDO'
      write (17,*) 'KBPP'
      write (17,*) 'ZEFF'
      write (17,'(1G24.17)') Zeff
      write (17,*) 'ERFC'
      write (17,*)  ERFnr
      write (17,'(3G24.17)') (PPNerf(j),j=1,ERFnr)
      write (17,'(3G24.17)') (PPLerf(j),j=1,ERFnr)
      write (17,*) 'GCFS'
      write (17,*)  EXPnr(lmax)
      write (17,'(3G24.17)') (PPNexp(j,lmax),j=1,EXPnr(lmax))
      write (17,'(3I16)') (PPRexp(j,lmax),j=1,EXPnr(lmax))
      write (17,'(3G24.17)') (PPLexp(j,lmax),j=1,EXPnr(lmax))
      write (17,*) 'PPLM'
      write (17,*) lmax
      write (17,*) 'NLCS'
      write (17,*) '1',Vpsinr(0)
      write (17,'(1G24.17)') Vl2(0)
      write (17,'(3G24.17)') (avpsi(j,0),j=1,Vpsinr(0))
      write (17,'(3I16)') (nvpsi(j,0),j=1,Vpsinr(0))
      write (17,'(3G24.17)') (Cvpsi(j,0),j=1,Vpsinr(0))
      if (lmax.ge.1) then
         write (17,*) 'NLCP'
         write (17,*) '1',Vpsinr(1)
         write (17,'(1G24.17)') Vl2(1)
         write (17,'(3G24.17)') (avpsi(j,1),j=1,Vpsinr(1))
         write (17,'(3I16)') (nvpsi(j,1),j=1,Vpsinr(1))
         write (17,'(3G24.17)') (Cvpsi(j,1),j=1,Vpsinr(1))
      endif
      if (lmax.ge.2) then
         write (17,*) 'NLCD'
         write (17,*) '1',Vpsinr(2)
         write (17,'(1G24.17)') Vl2(2)
         write (17,'(3G24.17)') (avpsi(j,2),j=1,Vpsinr(2))
         write (17,'(3I16)') (nvpsi(j,2),j=1,Vpsinr(2))
         write (17,'(3G24.17)') (Cvpsi(j,2),j=1,Vpsinr(2))
      endif
      if (lmax.ge.3) then
         write (17,*) 'NLCF'
         write (17,*) '1',Vpsinr(3)
         write (17,'(1G24.17)') Vl2(3)
         write (17,'(3G24.17)') (avpsi(j,3),j=1,Vpsinr(3))
         write (17,'(3I16)') (nvpsi(j,3),j=1,Vpsinr(3))
         write (17,'(3G24.17)') (Cvpsi(j,3),j=1,Vpsinr(3))
      endif
      write (17,*) '&END'
c     
      close(17)
      write (*,*) 'created file "atom.kbpp".'
C
      return
      end
C
C     ==================================================================
      subroutine vpsiinit(nn,cmat)
C     ==================================================================
      use atom
      use pspot
      implicit none
      integer l
      real*8 r,cmat(namax,namax,0:lamax)
      real*8 nn(namax,namax,0:lamax)
      real*8 umat(namax,namax,0:lamax)
      integer i,kappa,lambda,n
      real*8 Cvpsi(2*namax*EXPmax,0:lamax)
      real*8 avpsi(2*namax*EXPmax,0:lamax)
      integer nvpsi(2*namax*EXPmax,0:lamax)
      integer vpsinr(0:lamax)
      common /vpsipara/ Cvpsi,avpsi,nvpsi,vpsinr
C     ------------------------------------------------------------------
      do l=0,lmax-1
         n=0
         do kappa=1,nalpha(l)
	    do i=1,EXPnr(l)
               n=n+1
               Cvpsi(n,l)=PPLexp(i,l)
     &              *cmat(kappa,1,l)*sqrt(nn(kappa,kappa,l))
               avpsi(n,l)=PPNexp(i,l)+alpha(kappa,l)
               nvpsi(n,l)=PPRexp(i,l)+l
	    enddo
	    do i=1,EXPnr(lmax)
               n=n+1
               Cvpsi(n,l)=-PPLexp(i,lmax)
     &              *cmat(kappa,1,l)*sqrt(nn(kappa,kappa,l))
               avpsi(n,l)=PPNexp(i,lmax)+alpha(kappa,l)
               nvpsi(n,l)=PPRexp(i,lmax)+l
	    enddo
         enddo
         vpsinr(l)=n
      enddo
c     
      Vpsinr(lmax)=1
      nvpsi(1,lmax)=0
      avpsi(1,lmax)=0.d0
      Cvpsi(1,lmax)=0.d0
c
      return
      end
c
C     ==================================================================
      REAL*8 FUNCTION VPSIR(l,r,cmat)
C     ==================================================================
      use atom
      use pspot
      implicit none
      integer l
      real*8 r,cmat(namax,namax,0:lamax)
      real*8 Cvpsi(2*namax*EXPmax,0:lamax)
      real*8 avpsi(2*namax*EXPmax,0:lamax)
      integer nvpsi(2*namax*EXPmax,0:lamax)
      integer vpsinr(0:lamax)
      common /vpsipara/ Cvpsi,avpsi,nvpsi,vpsinr
      integer i,kappa
C     ------------------------------------------------------------------
      vpsir=0.d0
      do i=1,vpsinr(l)
         vpsir=vpsir+Cvpsi(i,l)*r**nvpsi(i,l)*exp(-avpsi(i,l)*r**2)
      enddo
c
      return
      end
C     ==================================================================
      REAL*8 FUNCTION Vkbr(l,r)
C     ==================================================================
      use atom
      implicit none
      integer l
      integer i,j
      real*8 A(namax,0:lamax)
      integer nvkb(namax,0:lamax)
      real*8 avkb(namax,0:lamax)
      real*8 Vl(0:lamax)
      integer Vkbnr(0:lamax)
      common /Vkbpara/ A,nvkb,avkb,Vkbnr,Vl
      real*8 r
C     ------------------------------------------------------------------
      vkbr=0.d0
      do j=1,Vkbnr(l)
         vkbr=vkbr+A(j,l)*r**nvkb(j,l)*exp(-avkb(j,l)*r**2)
      enddo
c
      return
      end
