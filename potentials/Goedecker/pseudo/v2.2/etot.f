      subroutine etot(
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     aeval,
     :     rprb,zion,rloc,gpot,r_l,hsep,
     :     xp,ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,
     :     xcgrd,excgrd,rhogrd,occup,rr,rw,
     :     expxpr)

      implicit real*8 (a-h,o-z)
      dimension aeval(noccmx,lmx,nsmx),
     :     gpot(4),r_l(lmx),hsep(6,lpmx,nsmx),
     :     xp(0:ng),  ud(nint,((ng+1)*(ng+2))/2,lcx+1), 
     :     psi(0:ngmx,noccmx,lmx,nsmx),rho(((ng+1)*(ng+2))/2,lmax+1),
     :     pp1(0:ng,lmax+1),pp2(0:ng,lmax+1),pp3(0:ng,lmax+1),
     :     xcgrd(nint), rhogrd(nint),
     :     occup(noccmx,lmx,nsmx),expxpr(0:ng,nint)
      dimension scpr1(nsmx),scpr2(nsmx),scpr3(nsmx)

      dimension vexgrd(nint),excgrd(nint),vhgrd(nint)
      dimension rr(nint),rw(nint)

      external gamma

      fourpi = 16.d0*atan(1.d0)
      exc  = 0.0d0
      eext = 0.0d0
      ekin = 0.0d0
      vxc  = 0.0d0
      ehart= 0.0d0
      eigsum=0.0d0
      enl  = 0.0d0
c
c calc hartree potential
      do i=1,nint
         vhgrd(i)=0.0d0
      enddo
      call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint,
     &             rho,1,1.d0,vhgrd,1)
c   calc eext, ehart, exc
      do k=1,nint
         r=rr(k)
         vexgrd(k)=.5d0*(r/rprb**2)**2-zion*Derf(r/(sqrt(2.d0)*rloc))/r 
     1        + exp(-.5d0*(r/rloc)**2)*
     1        ( gpot(1) + gpot(2)*(r/rloc)**2 + gpot(3)*(r/rloc)**4 + 
     1        gpot(4)*(r/rloc)**6 )
         eext  = eext  +       vexgrd(k)*rhogrd(k)*rw(k)
         ehart = ehart + 0.5d0*vhgrd(k) *rhogrd(k)*rw(k)
         exc   = exc   +       excgrd(k)*rhogrd(k)*rw(k)
         vxc   = vxc+xcgrd(k)*rhogrd(k)
      enddo
      vxc=vxc*fourpi
      do ll=0,lmax
         do ispin=1,min(2*ll+1,nspin)
            if (ll.le.lpx) then
	rnrm1=1.d0/sqrt(.5d0*gamma(ll+1.5d0)*r_l(ll+1)**(2*ll+3))
	rnrm2=1.d0/sqrt(.5d0*gamma(ll+3.5d0)*r_l(ll+1)**(2*ll+7))
	rnrm3=1.d0/sqrt(.5d0*gamma(ll+5.5d0)*r_l(ll+1)**(2*ll+11))
            endif
            do iocc=1,noccmax
               zz = occup(iocc,ll+1,ispin)
               if (zz.gt.1.0d-8) then 
                  eigsum=eigsum + aeval(iocc,ll+1,ispin) *zz
c     separabel part
                  if (ll.le.lpx) then
                     scpr1(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin),
     :                    1,pp1(0,ll+1),1)
                     scpr2(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin),
     :                    1,pp2(0,ll+1),1)
                     scpr3(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin),
     :                    1,pp3(0,ll+1),1)
                  endif
                  do k=1,nint
                     r=rr(k)
c     wavefunction on grid
                     psigrd = wave2(ng,ll,psi(0,iocc,ll+1,ispin),
     :                    expxpr,r,k,nint)
c     kinetic energy	
                     rkin=0.d0
                     do i=0,ng
                        rkin=rkin+psi(i,iocc,ll+1,ispin)*(xp(i)*
     :                       (3.d0+2.d0*ll-2.d0*xp(i)*r**2)
     :                       *exp(-xp(i)*r**2))
                     enddo
                     rkin=rkin*r**ll
                     ekin = ekin+rkin*(psigrd*r**ll)*zz*rw(k)/fourpi
c     separabel part
                     if (ll.le.lpx) then
                        sep = (scpr1(ispin)*hsep(1,ll+1,ispin) 
     :                       + scpr2(ispin)*hsep(2,ll+1,ispin) 
     :                       + scpr3(ispin)*hsep(4,ll+1,ispin))
     :                       *rnrm1*r**ll*exp(-.5d0*(r/r_l(ll+1))**2)+
     :                       (scpr1(ispin)*hsep(2,ll+1,ispin) 
     :                       + scpr2(ispin)*hsep(3,ll+1,ispin) 
     :                       + scpr3(ispin)*hsep(5,ll+1,ispin))
     :                       *rnrm2*r**(ll+2)
     :                       *exp(-.5d0*(r/r_l(ll+1))**2)   
     :                       +(scpr1(ispin)*hsep(4,ll+1,ispin) 
     :                       + scpr2(ispin)*hsep(5,ll+1,ispin) 
     :                       + scpr3(ispin)*hsep(6,ll+1,ispin))
     :                       *rnrm3*r**(ll+4)
     :                       *exp(-.5d0*(r/r_l(ll+1))**2)
                        enl = enl+sep*(psigrd*r**ll)*zz*rw(k)/fourpi

                     else
                        sep=0.d0
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo
      etotal = eigsum - ehart - vxc + exc
      write(6,*) 
      write(6,*)' Pseudo atom energies'
      write(6,*)' --------------------'
      write(6,'(a,f16.10)')' kinetic energy            =',ekin
      write(6,'(a,f16.10)')' vxc    correction         =',vxc
      write(6,'(a,f16.10)')' sum of eigenvalues        =',eigsum
      write(6,'(a,f16.10)')' exchange + corr energy    =',exc
      write(6,'(a,f16.10)')' external energy           =',eext
      write(6,'(a,f16.10)')' non local energy          =',enl
      write(6,'(a,f16.10)')' el-el  interaction energy =',ehart
      write(6,'(a,f16.10)')' total energy              =',etotal
      denfull=0.0d0
      write(6,*)
      do k=1,nint
         denfull=denfull+rhogrd(k)*rw(k)
      enddo
      write(6,'(a,e10.4)')'atomic+electronic charge ',zion-denfull
      return
      end
