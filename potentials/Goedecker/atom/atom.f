      program main
c
c     run atomic program
c
      call atom
      end
c
       subroutine atom
c
       implicit double precision(a-h,o-z)
c      ******************************************************
c      *     program for atomic calculations                *
c      *     written by Sverre Froyen, February 1982        *
c      *     while at UC Berkeley, California               *
c      *                                                    *
c      *    this is a modified version from                 *
c      *    christian hartwigsen, february 1998             *
c      *    while at MPI Stuttgart, Germany                 *
c      *                                                    *
c      ******************************************************
c
c      some parameters are set inside the program
c      the most important ones are the tolerance for the self-
c      consistency in the potential (set in the main program)
c      and the accuracy in the eigenvalue for a given potential
c      (set in difrel,difnrl)
c
      parameter (nrmax=10000, maxorb=60, lmax=5, maxconf=19)
c
       dimension r(nrmax),rab(nrmax),
     2 no(maxorb),lo(maxorb),so(maxorb),zo(maxorb),
     3 cdd(nrmax),cdu(nrmax),cdc(nrmax),
     4 viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     5 vod(nrmax),vou(nrmax),
     6 etot(10),econf(maxconf),ev(maxorb),ek(maxorb),ep(maxorb)
c
       character*2 naold,itype,ityold,ispp*1,nameat,stop
       character*10 icorr,icold
C
c     c.hartwig: additioanl grids for modified integration
      dimension rw(10000),rd(10000)
      common /intgrd/ rw,rd
c----------------------------------------------------------------------
       tol = 1.0D-11
       naold = '  '
       icold = '  '
       ityold = '  '
       stop  = 'st'
       zsold = 0.D0
       nconf = 0
       dvold = 1.0D10
       nr    = 1
       norb  = 1
       open(35,file='atom.dat',status='unknown')
       open(unit=40,file='atom.ae',form='formatted')
       open(unit=50,file='psp.par',form='formatted')
       open(unit=60,file='weights.par',form='formatted')
c
c     begin main loop
c
 20    continue
c
c      read input data
c
       call input (itype,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,rprb,rcov,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep,nconf)

       if (itype .eq. stop) goto 140
       if (nconf .gt. maxconf) then
          write(6,*) 'too many configurations, max. is:',maxconf
          stop
       endif
c
c     r ...... radial mesh
c     nr ..... # mesh points
c     norb ... # orbitals
c     ncore .. # core orbitals (closed shells)
c     no ..... n quantum number
c     lo ..... l do.
c     so ..... spin (+/- 0.5, or 0 for unpolarized)
c     zo ..... # electrons
c     znuc ... atomic number
c

c       if (zsold .eq. zsh .and. naold .eq. nameat .and.
c     +    ityold .eq. itype) goto 45



c
c      set up initial charge density.
c      cdd and cdu  =  2 pi r**2 rho(r)
c
        aa = sqrt(sqrt(znuc))/2.0d0+1.0d0
        a2 = zel/4.0d0*aa**3
        do 30 i=1,nr
          cdd(i) = a2*exp(-aa*r(i))*r(i)**2
          cdu(i) = cdd(i)
 30     continue
c
c     cdd ..... charge density (spin down)
c     cdu ..... charge density (spin up)
c     cdc ..... core charge density (up to ncore orbitals)
c
c      set up ionic potentials
c
 40    call vionic(itype,icorr,ifcore,
     1 nrmax,nr,a,b,r,rab,rprb,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c     Potentials: always multiplied by r.
c     viod,u ..... ionic potential (down,up)
c     vid,u ...... input screening potential (down,up)
c     vod,u ...... output screening potential (down,up)
c
c      set up electronic potential
c
 45    call velect(0,0,icorr,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
       do 50 i=1,nr
       vid(i) = vod(i)
       viu(i) = vou(i)
 50    continue
c
c      start iteration loop
c
       iconv = 0
       icon2 = 0
       maxit = 5000
c     empirical function
       xmixo = 1.0d0/log(znuc+7.0d0)
c
       do 100 iter=1,maxit
c
       if (iter .eq. maxit) iconv=1
c
c      compute orbitals (solve Schrodinger equation)
c
       if (icon2 .eq. 0) then
c
c     finite difference solution (less accurate)
c
       call dsolv1(
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
      else
c
c     predictor - corrector method (more accurate)
c
       call dsolv2
     + (iter,iconv,icorr,ispp,ifcore,itype,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep,rcov,rprb,nconf)
c
      endif
c
c     etot ..... terms in Etotal
c     ev ....... eigenvalues
c     ek ....... kinetic energy for each orbital
c     ep ....... potential energy (Vionic*rho) for each orbital
c
c      set up output electronic potential from charge density
c
       call velect(iter,iconv,icorr,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
c      check for convergence (Vout - Vin)
c
       if (iconv .gt. 0) goto 120
       dvmax = 0.D0
       do 60 i=2,nr
          dv = (vod(i)-vid(i))/(1.D0+vod(i)+vou(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
          dv = (vou(i)-viu(i))/(1.D0+vou(i)+vod(i))
          if (abs(dv) .gt. dvmax) dvmax=abs(dv)
 60    continue
       iconv = 1
       icon2 = icon2+1
       if (dvmax .gt. tol) iconv=0
       if (dvmax .ge. dvold) xmixo=0.8D0*xmixo
c     diverging - reduce mixing coefficient
       if (xmixo .lt. 0.01D0) xmixo=0.01D0
       dvold = dvmax
      write(6,70) iter,dvmax,xmixo
 70    format(7h iter =,i5,9h dvmax = ,e9.3,8h xmixo =,e9.3)
c
c      mix input and output electronic potentials
c
       call mixer(iter,iconv,icon2,xmixo,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
c
 100   continue
c
       write(6,110) dvmax,xmixo
 110   format(/,34h potential not converged - dvmax =,e10.4,
     1 9h  xmixo =,f5.3)
       call ext(1)
c
c      find total energy
c
 120   call etotal(itype,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       if (naold .ne. nameat .or. icold .ne. icorr .or.
     +     ityold .ne. itype ) call prdiff(nconf,econf)
c       if (nconf .eq. 9) nconf=1
       nconf = nconf + 1
       econf(nconf) = etot(10)
       if (nconf .ne. 1) write(6,130) etot(10)-econf(1)
 130   format(//,28h excitation energy         =,f18.8,/,1x,45('-'))
       naold = nameat
       icold = icorr
       zsold = zsh
       ityold = itype
c
c     next configuration of the atom
c
       goto 20
c
 140    continue
c
c     write data to files psp.par/weights.par
c     for pseudopotential-fit
      if (ispp.ne.'r') ispp='n'
c     psp.par
      write(50,*) ' 10   2.0     ng, rij (initial guess) '
      write(50,'(2f15.10,a)') rcov, rprb, ' rcov, rprb '
      if (ispp.eq.'r') then
         write(50,*)'relativistic calculation'
      else
         write(50,*)'non relativistic calculation'
      endif
      write(50,'(t2,a10,t15,a)')icorr ,'XC-functional'
      write(50,'(3f7.3,2a)') znuc,zps,rcov/4.d0,
     :     '  0.0 0.0 0.0 0.0',
     :     '    znuc,zpseudo,rloc,gpot()'
      lpx=2
      write(50,*) lpx ,' lpx '
      do l=0,lpx
         write(50,'(f7.3,2a)')  rcov/4.0d0,
     :        '  0.0 0.0 0.0 0.0 0.0 0.0 ',
     :        'r_l(), hsep()'
         if (ispp.eq.'r' .and. l.ne.0 )
     :        write(50,'(tr7,2a)')
     :        '  0.0 0.0 0.0 0.0 0.0 0.0 ',
     :        ' hsep()'
      enddo
c     weights.par
      write(60,*) 'Initial weight-guess from atomic-program'
      write(60,*) '----------------------------------------'
      write(60,*)'   .1E+03     weight_psir0'
      if (nconf.gt.1) then
         write(60,*) ('1.00 ',i=1,nconf),' weights for configurations'
         write(60,*) '0.00 ',('1.00 ',i=2,nconf),
     :        ' weights for excitation-energies '
      endif
      write(60,*) '   n   l  so   eigval  chrg    dchrg  ddchrg',
     :     '  res    rnode dnode ddnode'
      do iorb=ncore+1,norb
         weight=0.0d0
         if (zo(iorb).gt.1.0d-4) then
            write(60,'(2i4,f5.2,tr3,a)') no(iorb),lo(iorb),so(iorb),
     :       '1.0e0   1.0e0   0.0e0  0.0e0   0.0e0  1.0e0 0.0e0 0.0e0'
         else
            write(60,'(2i4,f5.2,tr3,a)') no(iorb),lo(iorb),so(iorb),
     :       '0.0e0   0.0e0   0.0e0  0.0e0   0.0e0  0.0e0 0.0e0 0.0e0'
         endif
      enddo
c
c     append excitation energies (in hartree!) to file atom.ae
c     if more than one configuration
c
      if (nconf.gt.1) then
         write(40,*) 'EXCITATION ENERGIES:'
         write(40,*) ((econf(i)-econf(1))/2.d0,i=1,nconf)
      endif
c
      call prdiff(nconf,econf)
      call ext(0)
      end
c
c      *****************************************************************
c
       subroutine prdiff(nconf,econf)
       implicit double precision(a-h,o-z)
       dimension econf (*)
       if (nconf .le. 1) goto 40
       write(6,*)
       write(6,*)'---------------------------------------------'
       write(6,10) (i,i=1,nconf)
 10    format(25h Total energy differences,//,2x,19i9)
       do 30 i=1,nconf
       write(6,20) i,(econf(i)-econf(j),j=1,i)
 20    format(1x,i2,1x,19f9.5)
 30    continue
 40    nconf = 0
       return
       end
c
c      *****************************************************************
c
       subroutine mixer(iter,iconv,icon2,xmixo,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      subroutine computes the new exchange correlation potential
c      given the input and the output potential from the previous
c      iteration.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       character*2 ispp*1,nameat
       character*10 icorr
c
       xmixi = 1 - xmixo
       do 120 i=1,nr
       vid(i) = xmixo * vod(i) + xmixi * vid(i)
       viu(i) = xmixo * vou(i) + xmixi * viu(i)
 120   continue
       return
       end
c
c      *****************************************************************
c
       subroutine etotal(itype,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      etotal computes the total energy from the electron charge density.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       character*2 itype,nameat
c
c      etot(i)    i=1,10 contains various contributions to the total
c                 energy.
c                 (1)   sum of eigenvalues ev
c                 (2)   sum of orbital kinetic energies ek
c                 (3)   el-ion interaction from sum of orbital
c                       potential energies ep
c                 (4)   electrostatic el-el interaction  (from velect)
c                 (5)   vxc (exchange-correlation) correction to sum
c                       of eigenvalues                   (from velect)
c                 (6)   3 * vc - 4 * ec
c                       correction term for virial theorem
c                       when correlation is included     (from velect)
c                 (7)   exchange and correlation energy  (from velect)
c                 (8)   kinetic energy from eigenvalues  (1,3,4,5)
c                 (9)   potential energy
c                 (10)  total energy
c
       dimension il(5)
       character il*1
 1     format(/,1x,a10,30(/,1x,10e13.4))
c       pi = 4*atan(1.D0)
c
c      sum up eigenvalues ev, kinetic energies ek, and
c      el-ion interaction ep
c
       etot(1) = 0.D0
       etot(2) = 0.D0
       etot(3) = 0.D0
c     c.hartwig
c     subtract vshift
      vshift=-15.0d0
       do 10 i=1,norb
        etot(1) = etot(1) + zo(i)*(ev(i)-vshift)
        etot(2) = etot(2) + zo(i)*ek(i)
        etot(3) = etot(3) + zo(i)*(ep(i)-vshift)
c       etot(1) = etot(1) + zo(i)*ev(i)
c       etot(2) = etot(2) + zo(i)*ek(i)
c       etot(3) = etot(3) + zo(i)*ep(i)
 10    continue
c
c      compute interaction shell - (nucleus-core)
c
       esh = 0.D0
       if (zsh .ne. 0.D0) esh = 2*zsh*(znuc-zcore)/rsh
c
c      kinetic energy
c
       etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)
c
c      potential energy
c
       etot(9) = etot(3) + etot(4) + etot(7) + esh
c
c      total energy
c
       etot(10) = etot(1) - etot(4) - etot(5) + etot(7) + esh
c
c      printout
c
       il(1) = 's'
       il(2) = 'p'
       il(3) = 'd'
       il(4) = 'f'
       il(5) = 'g'
       write(6,*)
       write(6,20) nameat
 20    format(a3,25h output data for orbitals,/,1x,27('-'),//,
     1 17h nl    s      occ,9x,'eigenvalue',4x,14hkinetic energy,
     2 6x,'pot energy'/)
       do 40 i=1,norb
c     c.hartwig give energies in hartree
      ev(i) = ev(i) - vshift
      ep(i) = ep(i) - vshift
       write(6,30) no(i),il(lo(i)+1),so(i),zo(i),
     :         ev(i)/2,ek(i)/2,ep(i)/2
 30    format(1x,i1,a1,f6.1,f10.4,3f17.8)
 40    continue
c     c.hartwig give energies in hartree; no virial correction
      write(6,50) (etot(i)/2,i=1,5), (etot(i)/2,i=7,10)
 50    format(//,15h total energies,/,1x,14('-'),/,
     1 /,28h sum of eigenvalues        =,f18.8,
     2 /,28h kinetic energy from ek    =,f18.8,
     3 /,28h el-ion interaction energy =,f18.8,
     4 /,28h el-el  interaction energy =,f18.8,
     5 /,28h vxc    correction         =,f18.8,
     7 /,28h exchange + corr energy    =,f18.8,
     8 /,28h kinetic energy from ev    =,f18.8,
     9 /,28h potential energy          =,e18.8,/,1x,45('-'),
     X /,28h total energy              =,f18.8)
       return
       end
c
c      *****************************************************************
c
       subroutine ext(i)
c
c      i  is a stop parameter
c
c      000-099 main (0 is normal exit)
c      100-199 input
c      200-299 charge
c      300-399 vionic
c      400-499 velect
c      500-599 dsolv1
c      600-699 dsolv2 (including difnrl and difrel)
c      700-799 etotal
c      800-899 pseudo
c
       if (i .ne. 0) write(6,10) i
 10    format(17h1stop parameter =,i3)
       stop
       end
c
c      *****************************************************************
c
       subroutine vionic(itype,icorr,ifcore,
     1 nrmax,nr,a,b,r,rab,rprb,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      vionic sets up the ionic potential
c      note that vio is the ionic potential times r
c
       dimension r(*),rab(*),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(*),cdu(*),cdc(*),
     4 viod(lmax,*),viou(lmax,*),vid(*),viu(*),vod(*),vou(*),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       character*2 itype,nameat,icalc,cdtyp
       character*10 icorr
c
       dimension iray(6)
       character namef*6,iray*8,
     1 namet*2,icorrt*2,mcore*4,irel*3
c.....files
      common /files/iinput,iout,in290,in213,istore,iunit7,iunit8,istruc,
     +               ivnlkk,isumry,ikpts
c
c      2*znuc part (Rydberg units)
c
       ifcore = 0
       do 10 i=1,lmax
       do 12 j=1,nrmax
c  c.hartwig  add confining potential
          viod(i,j) = -2.0d0*(znuc -.5d0*(r(j)/rprb**2)**2*r(j))
          viou(i,j) = -2.0d0*(znuc -.5d0*(r(j)/rprb**2)**2*r(j))
c     viod(i,j) = -2.0*( -.5d0*(r(j)/rprb**2)**2*r(j))
c     viou(i,j) = -2.0*( -.5d0*(r(j)/rprb**2)**2*r(j))
c
c     c.hartwig  shift potential to avoid positive eigenvalues
c     and convergence problems
          vshift=-15.0d0*r(j)
          viod(i,j) = viod(i,j)+vshift
          viou(i,j) = viou(i,j)+vshift
 12    continue
 10   continue
c
c      add potential from shell charge
c
 105   if (zsh .eq. 0.D0) return
       do 110 i=1,lmax
       do 110 j=1,nr
       if (r(j) .ge. rsh) viod(i,j) = viod(i,j) - 2*zsh
       if (r(j) .ge. rsh) viou(i,j) = viou(i,j) - 2*zsh
       if (r(j) .lt. rsh) viod(i,j) = viod(i,j) - 2*zsh*r(i)/rsh
       if (r(j) .lt. rsh) viou(i,j) = viou(i,j) - 2*zsh*r(i)/rsh
 110   continue
       return
       end
c
c      *****************************************************************
c
       subroutine velect(iter,iconv,icorr,ispp,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      velect generates the electronic output potential from
c      the electron charge density.
c      the ionic part is added in dsolve.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       dimension vtemp(1000)
       character*2 ispp*1,nameat,itype
       character*10 icorr
c
      parameter ( mesh = 2000 )
       dimension y(mesh),yp(mesh),ypp(mesh),w(3*mesh),s1(mesh),s2(mesh)
       common  y,yp,ypp,w,s1,s2
c
c      for use in routine atomwr:
       parameter (ntitle = 60)
       character*40 text(ntitle)
       character irel*3, xccore*4, cdtyp*2

c     c.hartwig
      dimension rho(nr),excgrd(nr),vxcgrd(nr)
      dimension rw(10000),rd(10000)
      common /intgrd/ rw,rd
      INCLUDE 'func.inc'

c
c
       pi = 4*atan(1.D0)
c
c      fit cd/r by splines
c
       y(1) = 0.D0
       do 10 i=2,nr
       y(i) = (cdd(i)+cdu(i))/r(i)
       if (ifcore .eq. 2) y(i) = y(i) + cdc(i)/r(i)
 10    continue
       isx = 0
       a1 = 0.D0
       an = 0.D0
       b1 = 0.D0
       bn = 0.D0
       call splift(r,y,yp,ypp,nr,w,ierr,isx,a1,b1,an,bn)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
       xlo = 0.D0
       call spliq(r,y,yp,ypp,nr,xlo,r,nr,s2,ierr)
       do 20 i=1,nr
       ypp(i) = r(i)*ypp(i) + 2*yp(i)
       yp(i)  = r(i)*yp(i)  + y(i)
       y(i)   = r(i)*y(i)
 20    continue
       call spliq(r,y,yp,ypp,nr,xlo,r,nr,s1,ierr)
c
c      check normalization
c
       xnorm = 0.D0
       if (zel .ne. 0.D0) xnorm = zel/s1(nr)
       if (iter .gt. 0 .and. abs(zel-s1(nr)) .gt. 0.01D0)
     1 write(6,25) iter,xnorm
 25    format(/,46h warning *** charge density rescaled in velect,
     1 /,17h iteration number,i4,3x,16hscaling factor =,g10.3,/)
c
c      compute new hartree potential
c      renormalize the charge density
c
       do 30 i=2,nr
       vod(i) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
       vou(i) = vod(i)
       cdd(i) = xnorm*cdd(i)
       cdu(i) = xnorm*cdu(i)
 30    continue
c
       if (iconv .ne. 1) goto 50
c
c      compute hartree contribution to total energy
c
       ehart = 0.D0
       ll = 4
       do 40 i=2,nr
       ehart = ehart + ll * (cdd(i)+cdu(i)) * vod(i) * rab(i)
       ll = 6 - ll
 40    continue
       ehart = ehart / 6
c
c      find derivatives of the charge density
c
       do 45 i=2,nr
c
 45    continue
c
c      store the atomic Coulomb (ionic + Hartree) potential on file
c
c      first construct the total potential, store in array vtemp:
c
c       ifile = 2
c       irectp = 31
       do 300 l = 1, 3
c       do 310 i = 1, nr
c         vtemp(i) = viod(l,i) + vod(i) * r(i)
c310    continue
c       cdtyp = ' '
c       itype = ' '
c       call atomwr
c     +  (ifile,irectp,nameat,icorr,irel,xccore,zcore,norb,text,
c     +   nr,aa,bb,r,nql,delql,nqnl,delqnl,numnl,
c     +   itype,cdtyp,0,(l-1),mode,vtemp)
c
c       if (ispp .ne. ' ') then
c         do 220 i = 1, nr
c           vtemp(i) = viou(l,i) + vou(i) * r(i)
c220      continue
c         cdtyp = 'up'
c         call atomwr
c     +    (ifile,irectp,nameat,icorr,irel,xccore,zcore,norb,text,
c     +     nr,aa,bb,r,nql,delql,nqnl,delqnl,numnl,
c     +     itype,cdtyp,0,(l-1),mode,vtemp)
c        endif
300     continue

c        goto 50
c
c      add exchange and correlation
c
C..functionals
 50     continue
       salpha=2.D0/3.D0
       bbeta=0.0042D0
       betapp=0.0042D0
       IF (INDEX(icorr,"NONE").NE.0) THEN
         mfxcx = 0
         mfxcc = 0
         mgcx  = 0
         mgcc  = 0
       ELSE IF (INDEX(icorr,"SONLY").NE.0) THEN
         mfxcx = 1
         mfxcc = 0
         mgcx  = 0
         mgcc  = 0
       ELSE IF (INDEX(icorr,"VWN").NE.0) THEN
         mfxcx = 1
         mfxcc = 2
         mgcx  = 0
         mgcc  = 0
       ELSE IF (INDEX(icorr,"LDA").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 0
         mgcc  = 0
       ELSE IF (INDEX(icorr,"PADE").NE.0) THEN
         mfxcx = 0
         mfxcc = 9
         mgcx  = 0
         mgcc  = 0
       ELSE IF (INDEX(icorr,"BONLY").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 1
         mgcc  = 0
       ELSE IF (INDEX(icorr,"BP").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 1
         mgcc  = 1
       ELSE IF (INDEX(icorr,"BLYP").NE.0) THEN
         mfxcx = 1
         mfxcc = 3
         mgcx  = 1
         mgcc  = 2
       ELSE IF (INDEX(icorr,"XLYP").NE.0) THEN
         mfxcx = 1
         mfxcc = 3
         mgcx  = 7
         mgcc  = 2
       ELSE IF (INDEX(icorr,"PW91").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 2
         mgcc  = 3
       ELSE IF (INDEX(icorr,"PBE1W").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 3
         mgcc  = 6
       ELSE IF (INDEX(icorr,"REVPBE").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 4
         mgcc  = 4
       ELSE IF (INDEX(icorr,"PBES").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 9
         mgcc  = 7
       ELSE IF (INDEX(icorr,"PBE").NE.0) THEN
         mfxcx = 1
         mfxcc = 1
         mgcx  = 3
         mgcc  = 4
       ELSE IF (INDEX(icorr,"HCTH").NE.0) THEN
         mfxcx = 0
         mfxcc = 0
         IF (INDEX(icorr,"93").NE.0) THEN
           mgcx = 13
         ELSE IF (INDEX(icorr,"120").NE.0) THEN
           mgcx = 14
         ELSE IF (INDEX(icorr,"147").NE.0) THEN
           mgcx = 15
         ELSE IF (INDEX(icorr,"407").NE.0) THEN
           mgcx = 16
         ELSE
CMK        HCTH/120 implemented in CPMD
           icorr = "HCTH"
           mgcx = 5
         END IF
         mgcc = mgcx
       ELSE IF (INDEX(icorr,"OPTX").NE.0) THEN
         mfxcx = 0
         mfxcc = 0
         mgcx  = 6
         mgcc  = 6
       ELSE IF (INDEX(icorr,"OLYP").NE.0) THEN
         mfxcx = 0
         mfxcc = 3
         mgcx  = 6
         mgcc  = 2
       ELSE IF (INDEX(icorr,"B97").NE.0) THEN
         mfxcx = 0
         mfxcc = 0
         IF (INDEX(icorr,'GRIMME').NE.0) THEN
           mgcx = 12
         ELSE
           PRINT *,"B97 needs exact exchange"
           mgcx = 11
         END IF
         mgcc = mgcx
       ELSE
         WRITE (6,*) 'Unknown functional(s): ',icorr
         STOP
       END IF
       do i=2,nr
          rho(i)=(cdd(i)+cdu(i))/4.d0/pi/r(i)**2
       enddo
       rho(1)=rho(2)-(rho(3)-rho(2))*r(2)/(r(3)-r(2))

CMK   this should avoid problems with XC-functionals
CMK   with kinks (BLYP,....)
      if (iter.lt.30) then
        mfxcx=0
        mfxcc=9
        mgcc=0
        mgcx=0
      else if (iter.eq.30) then
        write(6,*) 'Switching from LDA to the requested functional'
      endif

c     hutter's routine
c       call evxc(nr,r,rho,vxcgrd,excgrd)
c     goedecker's routine
       call ggaenergy_15(nr,rw,rd,rho,enexc,vxcgrd,excgrd)
c
c     c.hartwig modified integration
       exc=0.d0
       vxc=0.d0
       do i=1,nr
c     need energy/potential in ryd
          exct =  2.d0*excgrd(i)
          vxcd =  2.d0*vxcgrd(i)
          vxcu = vxcd
          rhoup=rho(i)/2.d0
          rhodw=rhoup
          vod(i) = vod(i) + vxcd
          vou(i) = vou(i) + vxcu
          vxc = vxc + (vxcd*rhodw + vxcu*rhoup) * rw(i)
          exc = exc + exct * rw(i)
       enddo
       etot(4) = ehart
       etot(5) = vxc
       etot(7) = exc
       return
       end
c
c      *****************************************************************
c
       subroutine input (itype,icorr,ispp,
     1 nrmax,nr,a,b,r,rab,rprb,rcov,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep,nconf)
       implicit double precision(a-h,o-z)
c
c      subroutine to read input parameters
c
       dimension r(nrmax),rab(nrmax),
     2 no(*),lo(*),so(*),zo(*),
     3 cdd(nrmax),cdu(nrmax),cdc(nrmax),
     4 viod(lmax,nrmax),viou(lmax,nrmax),vid(nrmax),viu(nrmax),
     + vod(nrmax),vou(nrmax),
     5 etot(10),ev(*),ek(*),ep(*)
       character*2 itype,ispp*1,nameat,blank*1
       character*10 icorr

       dimension rw(10000),rd(10000)
       common /intgrd/ rw,rd
       save nvalo,ncoreo


c      for use in routine atomwr:
       parameter (ntitle = 60)
       character*40 text(ntitle)
       character*80 instrg
       character irel*3, icalc*2, cdtyp*2
       character spdf(5)
       dimension nc(15),lc(15),nomin(5),iray(2)
       character iray*8,name*3
c
       data nc /1,2,2,3,3,3,4,4,4,4,5,5,5,6,6/
       data lc /0,0,1,0,1,2,0,1,2,3,0,1,2,0,1/
       data nomin /5*10/
       data spdf /'s','p','d','f','g'/
      data blank /' '/
c------------------------------------------------------------------
      itype='ae'
 10   read(35,'(a)',err=998,end=999) instrg
      if (instrg.eq.' ') goto 10
      if (index(instrg,'NEXT CONFIGURATION').ne.0) goto 89
      if (nconf.ge.1) goto 10
      j1=1
      j2=2
      do i=len(instrg),1,-1
         if (instrg(i:i).ne.' ') j1=i
      enddo
      do i=len(instrg),j1,-1
         if (instrg(i:i).eq.' ') j2=i
      enddo
      j2=j2-1
      nameat=instrg(j1:j2)
      if (j2.eq.1) nameat(2:2)=' '
      read(35,'(a)',err=998,end=999) instrg
      j1=1
      j2=2
      do i=len(instrg),1,-1
         if (instrg(i:i).ne.' ') j1=i
      enddo
      do i=len(instrg),j1,-1
         if (instrg(i:i).eq.' ') j2=i
      enddo
      j2=j2-1
      icorr=instrg(j1:j2)
      read(35,'(a)',err=998,end=999) instrg
      do i=len(instrg),1,-1
         if (instrg(i:i).ne.' ') j1=i
      enddo
      ispp=instrg(j1:j1)
      if (ispp.eq.'R') ispp='r'
      if (ispp.ne.'r') ispp=' '

       if (ispp .ne. 's' .and. ispp  .ne. 'r')  ispp=blank
c      spin-polarization needs relativistic calculation
       znuc=0.d0
       read(35,*,err=998,end=999) rmax,aa,bb
       read(35,*,err=998,end=999) rcov,rprb
       znuc=charge(nameat)
c
c      set up grid
c
       if (abs(rmax) .lt. 0.00001) rmax=100.0d0
       if (abs(aa) .lt. 0.00001) aa = 3.0d0
       if (abs(bb) .lt. 0.00001) bb = 40.0d0
       if (znuc .eq. 0.0d0) then
          a = 10**(-aa)
          goto 29
       endif
       a=exp(-aa)/znuc
       b = 1/bb
c
c     modify grid-parameter, so that one grid-point matches
c     rcov exact
c
        do i=1,nrmax
           if (i .eq. nrmax) then
              write(6,50)
              stop 'input two'
           endif
           r(i) = a*(exp(b*(i-1))-1)
           if (r(i).ge.rcov) then
              a= rcov/(exp(b*(i-1))-1)
              aa=-log(a*znuc)
              goto 29
           endif
        enddo
 29     continue
        do 30 i=1,nrmax
           if (i .eq. nrmax) then
              write(6,50)
 50           format(/,' error in input - arraylimits',
     1             ' for radial array exceeded',/)
              call ext(100)
           endif
          r(i) = a*(exp(b*(i-1))-1)
          rab(i) = (r(i)+a)*b
c
c     c.hartwig: set up grids for modified integration
c
          rw(i) = b*(r(i)+a)
          rd(i) = 1.d0/rw(i)
          rw(i)=rw(i)*12.56637061435917d0*r(i)**2
          if (r(i) .gt. rmax) goto 60
 30     continue
 60     nr = i-1
c
c     modify weights at end point for improved accuracy
c
        rw(1)=rw(1)*17.d0/48.d0
        rw(2)=rw(2)*59.d0/48.d0
        rw(3)=rw(3)*43.d0/48.d0
        rw(4)=rw(4)*49.d0/48.d0
c
c      read the number of core and valence orbitals
c
 6011 read(35,*,err=998,end=999) ncore, nval
      nvalo=nval
      ncoreo=ncore
       if (ncore .gt. 15) then
          write (6,*) 'more than 15 core orbitals'
          call ext(101)
       endif

 89    continue
       ncore=ncoreo
       nval =nvalo
       if (ispp.eq.'R') ispp='r'
       if (ispp.ne.'r') ispp=' '
       if (ispp .ne. 's' .and. ispp  .ne. 'r')  ispp=blank
c
c      compute occupation numbers and orbital energies for the core
c
       zcore = 0.D0
       if (ispp .eq. blank) then
          jmax = 1
         sc = 0.0D0
       else
         jmax = 2
         sc = - 0.5 D0
         endif
      norb = 0
       if (ncore .eq. 0) goto 85
       do 80 i=1,ncore
       do 80 j=1,jmax
       norb = norb + 1
       no(norb) = nc(i)
       lo(norb) = lc(i)
       so(norb) = sc
       zo(norb) = 2*lo(norb)+1
       if (ispp .eq. blank) zo(norb) = 2*zo(norb)
       if (ispp .eq. 'r') zo(norb) = 2*(lo(norb)+sc)+1
       zcore = zcore + zo(norb)
       if (abs(zo(norb)) .lt. 0.1D0) norb=norb-1
       if (ispp .ne. blank) sc=-sc
 80    continue
       ncore = norb
 85    continue
       norb = ncore
       zval = 0.D0
       if (nval .eq. 0) goto 105
c
       do 90 i=1,nval
c
       read(35,*,err=998,end=999) ni,li,zd,zu
       si = 0.D0
       if (ispp .ne. blank) si=0.5D0
c
       do 90 j=1,jmax
c
       norb = norb + 1
       if (ispp .ne. blank) si=-si
       no(norb) = ni
       lo(norb) = li
       so(norb) = si
       zo(norb) = zd+zu
c     c.hartwig
          if (zo(norb) .eq. 0.0) then
             zo(norb)=1.0d-20
          endif
c
       if (ispp .eq. 's' .and. si .lt. 0.1D0) zo(norb) = zd
       if (ispp .eq. 's' .and. si .gt. 0.1D0) zo(norb) = zu
       if (ispp .eq. 'r') zo(norb)=zo(norb)*(2*(li+si)+1)/(4*li+2)
       zval = zval + zo(norb)
       if (ispp .eq. 'r' .and. li+si .lt. 0.D0) norb=norb-1
       if (norb .eq. 0) goto 90
       if (nomin(lo(norb)+1) .gt. no(norb)) nomin(lo(norb)+1)=no(norb)
 90    continue
       nval = norb - ncore
c
c      abort if two orbitals are equal
c
       if (norb .le. 0) call ext(110)
       do 100 i = 1, (norb - 1)
       do 100 j = (i + 1),norb
         if (no(i) .eq. no(j) .and. lo(i) .eq. lo(j)) then
            if (abs(so(i)-so(j)) .lt. 1.0D-3) then
               print*,'i,no(i),lo(i),so(i):',i,no(i),lo(i),so(i)
               print*,'j,no(j),lo(j),so(j):',j,no(j),lo(j),so(j)

               call ext(110+i)
            endif
          endif
 100     continue
 105   zion = znuc - zcore - zval
c       write(6,*)' zion = ',zion
c       write(6,*)' znuc = ',znuc
c       write(6,*)' zcore = ',zcore
c       write(6,*)' zval = ',zval
       zel = zval
       zel=zel+zcore
c
       write(6,120) nameat
 120   format(1x,a2,' all electron calculation  '/,1x,27('-'),/)
       if (ispp .eq. 'r') write(6,150)
 150   format(' r e l a t i v i s t i c ! !'/)
       name = '   '
       if (ispp .ne. 's') name = 'non'
       write(6,160) icorr,name
 160   format(' correlation = ',a10,3x,a3,' spin-polarized'/)
       write(6,170) znuc,ncore,nval,zel,zion
 170   format(' nuclear charge             =',f10.6,/,
     1        ' number of core orbitals    =',i3,/,
     2        ' number of valence orbitals =',i3,/,
     3        ' electronic charge          =',f10.6,/,
     4        ' ionic charge               =',f10.6,//)
       if (zsh .ne. 0.D0) write(6,175) zsh,rsh
 175   format(' shell charge =',f6.2,' at radius =',f6.2,//)
       write(6,180)
 180   format(' input data for orbitals'//,
     1        '  i    n    l    s     j     occ'/)
       xji = 0.D0
       do 200 i=1,norb
         if (ispp .eq. 'r') xji = lo(i) + so(i)
         write(6,190) i,no(i),lo(i),so(i),xji,zo(i)
 190     format(1x,i2,2i5,2f6.1,f10.4)
 200     continue
      write(6,210) r(2),nr,r(nr),aa,bb
 210  format(//,' radial grid parameters',//,
     1 ' r(1) = .0 , r(2) =',e12.6,' , ... , r(',i4,') =',f12.8,
     2 /,' a =',f12.8,'  b =',f12.8,/)
      irel   = 'nrl'
      if (ispp .eq. 'r') irel = 'rel'
      if (ispp .eq. 's') irel = 'spp'
      do 25 i = 1, norb
        write (text(i),24) no(i),spdf(lo(i)+1),so(i),zo(i),irel
24      format (1x,i1,a,' s=',f4.1,' (occ=',f6.3,') ',a)
25      continue
1000   return

 998   write(6,*) 'Error while reading atom.dat'
       stop
 999   write(6,*) 'Reached end of file atom.dat'
       itype='stop'
       return
       end
c
c      *****************************************************************
c
c      *****************************************************************
c
       double precision function charge(name)
c
c    function determines the nuclear charge of an element
c
       parameter ( nelem = 103 )
       character*2 name, elemnt, pertab(nelem)
       integer ic(2)
c      the periodic table
       data pertab /
     +  'H ','HE',
     +  'LI','BE','B ','C ','N ','O ','F ','NE',
     +  'NA','MG','AL','SI','P ','S ','CL','AR',
     +  'K ','CA',
     +       'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     +            'GA','GE','AS','SE','BR','KR',
     +  'RB','SR',
     +       'Y ','ZR','NB','MO','TC','RU','RH','PD','AG','CD',
     +            'IN','SN','SB','TE','I ','XE',
     +  'CS','BA',
     +       'LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     +                                'HO','ER','TM','YB','LU',
     +            'HF','TA','W ','RE','OS','IR','PT','AU','HG',
     +            'TL','PB','BI','PO','AT','RN',
     +  'FR','RA',
     +       'AC','TH','PA','U ','NP','PU','AM','CM','BK','CF',
     +                                'ES','FM','MD','NO','LR'/
c
c      convert the name to upper-case, and possibly left-justify
c
c      code 97-122: lower case
c      code 65-90:  upper case
c      code 32:     blank
c
       do 100 i = 1,2
c      get the ascii value
       ic(i) = ichar( name(i:i) )
       if (ic(i) .ge. 97 .and. ic(i) .le. 122) then
c        convert to upper case
         ic(i) = ic(i) - 32
       else if (ic(i) .ge. 65 .and. ic(i) .le. 90) then
c        upper-case - do nothing
       else if (ic(i) .eq. 32) then
c        'space' - do nothing
       else if (ic(i) .eq. 0) then
c        'nul' - replace by space
         ic(i) = 32
       else
         write (6,*) 'unrecognized element name:',name
         call ext(200)
         endif
100    continue
c
c      left justify
       if (ic(1) .eq. 32) then
         ic(1) = ic(2)
         ic(2) = 32
         endif
c      the standard name of the element:
       elemnt = char(ic(1))//char(ic(2))
c
c      find the element in the periodic table
c
       do 150 i = 1, nelem
         if (elemnt .eq. pertab(i)) then
           charge = i
           return
           endif
150      continue
       write (6,160) name,elemnt,ic
160    format (' could not locate name in list of elements:'/
     + ' name=',a,' converted to=',a,' ascii codes=',2i3)
       call ext (200)
       return
       end
c      *****************************************************************
c
       subroutine dsolv1(
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      dsolv1 finds the non relativistic wave function
c      using finite differences and matrix diagonalization
c      initial guess for the eigenvalues need not be supplied
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       character*2 nameat
c
      parameter ( mesh = 4000 , nvmax = 6*mesh )
       dimension nmax(2,5),dk(mesh),d(mesh),sd(mesh),sd2(mesh),e(10),
     1 ind(10),z(nvmax),
     2 rv1(mesh),rv2(mesh),rv3(mesh),rv4(mesh),rv5(mesh)
       common dk,d,sd,sd2,z,rv1,rv2,rv3,rv4,rv5
c.....files
      common /files/iinput,iout,in290,in213,istore,iunit7,iunit8,istruc,
     +               ivnlkk,isumry,ikpts
c
c
c      initialize charge density arrays
c
       do 10 i=1,nr
       cdd(i) = 0.D0
       cdu(i) = 0.D0
 10    continue
c
c      find max n given l and s
c      zero spin is treated as down
c
       do 20 i=1,2
       do 20 j=1,lmax
       nmax(i,j) = 0
       do 20 k=1,norb
       if (no(k) .le. 0) goto 20
       if (lo(k) .ne. j-1) goto 20
       if ((so(k)-0.1D0)*(i-1.5D0) .lt. 0.D0) goto 20
       nmax(i,j)=no(k)
       if (no(k)*(nr-1) .gt. nvmax) then
         print*,no(k),nr-1
         print*,no(k)*(nr-1)," > ",nvmax
         call ext(500)
       end if
 20    continue
c
c      set up hamiltonian matrix for kinetic energy
c      only the diagonal depends on the potential
c
       c2 = -1.D0/b**2
       c1 = -2.D0*c2 + 0.25D0
       dk(1)  = c1 / (r(2)+a)**2
       sd(1)  = 0.D0
       sd2(1) = 0.D0
       do 30 i=3,nr
       dk(i-1)  = c1 / (r(i)+a)**2
       sd(i-1)  = c2 / ((r(i)+a)*(r(i-1)+a))
       sd2(i-1) = sd(i-1)**2
 30    continue
c
c      start loop over spin down=1 and spin up=2
c
       nrm = nr - 1
       do 80 i=1,2
c
c      start loop over s p d... states
c
       do 80 j=1,lmax
       if (nmax(i,j) .eq. 0) goto 80
       llp = j*(j-1)
       do 40 k=2,nr
       if (i .eq. 1) d(k-1) = dk(k-1)
     1  + (viod(j,k) + llp/r(k))/r(k) + vid(k)
       if (i .eq. 2) d(k-1) = dk(k-1)
     1  + (viou(j,k) + llp/r(k))/r(k) + viu(k)
 40    continue
c
c      diagonalize
c
       eps = -1.D0
       call tridib(nrm,eps,d,sd,sd2,bl,bu,1,nmax(i,j),e,ind,ierr,
     1 rv4,rv5)
       if (ierr .ne. 0) write(6,50) ierr
 50    format(/,21h ****** error  ierr =,i3,/)
       call tinvit(nrm,nrm,d,sd,sd2,nmax(i,j),e,ind,z,ierr,
     1 rv1,rv2,rv3,rv4,rv5)
       if (ierr .ne. 0) write(6,50) ierr
c
c      save energy levels and add to charge density
c
       ki = 1
       kn = 0
       do 70 k=1,norb
       if (no(k) .le. 0) goto 70
       if (lo(k) .ne. j-1) goto 70
       if ((so(k)-0.1D0)*(i-1.5D0) .lt. 0.D0) goto 70
       ev(k) = e(ki)
       do 60 l=2,nr
       denr = zo(k) * z(kn+l-1)**2 / rab(l)
       if (i .eq. 1) cdd(l) = cdd(l) + denr
       if (i .eq. 2) cdu(l) = cdu(l) + denr
 60    continue
       ki = ki + 1
       kn = kn + nrm
 70    continue
 80    continue
c
c      end loop over s p and d states
c
       return
       end
c
c      *****************************************************************
c
       subroutine dsolv2
     + (iter,iconv,icorr,ispp,ifcore,itype,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep,rcov,rprb,nconf)
       implicit double precision(a-h,o-z)
c
c      dsolv2 finds the (non) relativistic wave function using
c      difnrl to intgrate the Scroedinger equation or
c      difrel to intgrate the Dirac equation
c      the energy level from the previous iteration is used
c      as initial guess, and it must therefore be reasonable
c      accurate.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       character*2 ispp*1,nameat,itype
       character*10 icorr
c
       parameter ( mesh = 2000 )
       dimension v(mesh),ar(mesh),br(mesh)
       common  v,ar,br
c.....files
      common /files/iinput,iout,in290,in213,istore,iunit7,iunit8,istruc,
     +               ivnlkk,isumry,ikpts
c
c
c      initialize arrays for charge density
c
       do 10 i=1,nr
       cdd(i) = 0.D0
       cdu(i) = 0.D0
       if (ifcore .ne. 1) cdc(i)=0.D0
 10    continue
c
c      start loop over orbitals
c      note that spin zero is treated as down
c
       do 50 i=1,norb
       if (no(i) .le. 0) goto 50
       if (zo(i) .eq. 0.D0 .and. iconv .eq. 0) goto 50
       if (ev(i) .ge. 0.D0) ev(i)=-1.D0
c
c      set up potential
c
       lp  = lo(i)+1
       llp = lo(i)*lp
       do 20 j=2,nr
       if (so(i) .lt. 0.1D0) v(j) = viod(lp,j)/r(j) + vid(j)
       if (so(i) .gt. 0.1D0) v(j) = viou(lp,j)/r(j) + viu(j)
       if (ispp .ne. 'r') v(j) = v(j) + llp/r(j)**2
 20    continue
c
c      call integration routine
c
       if (ispp .ne. 'r' ) call difnrl(iter,i,v,ar,br,
     1 lmax,nr,a,b,r,rab,
     2 norb,no,lo,so,
     3 znuc,
     4 viod,viou,vid,viu,ev)
       if (ispp .eq. 'r' ) call difrel(iter,i,v,ar,br,
     1 lmax,nr,a,b,r,rab,norb,
     2 no,lo,so,znuc,viod,viou,vid,viu,ev)
c
c      add to the charge density
c
       do 30 j=1,nr
       denr = zo(i) * ar(j) * ar(j)
       if (ispp .eq. 'r') denr = denr + zo(i) * br(j) * br(j)
       if (so(i) .lt. 0.1D0) cdd(j) = cdd(j) + denr
       if (so(i) .gt. 0.1D0) cdu(j) = cdu(j) + denr
       if (ifcore .ne. 1 .and. i .le. ncore) cdc(j)=cdc(j)+denr
 30    continue
c
c      compute various quantitities if last iteration
c
       if (iconv .eq. 1) call orban
     + (itype,icorr,ispp,i,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,v,ev,ek,ep,rcov,rprb,nconf)
 50    continue
c
c      end loop over orbitals
c
       return
       end
c
c      *****************************************************************
c
      subroutine difnrl(iter,iorb,v,ar,br,lmax,
     1 nr,a,b,r,rab,norb,no,lo,so,znuc,viod,viou,
     2 vid,viu,ev)
c
c    difnrl integrates the Schroedinger equation
c    if finds the eigenvalue ev, the wavefunction ar
c    and the derivative br = d(ar)/dr
c
      implicit real*8 (a-h,o-z)
c
c  Tolerence
c
      parameter(etol=-1.d-7)
      parameter(tol=1.0d-14)
c
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),no(norb),
     1 lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb)
c
c    Arrays added to gain speed.
c
      dimension rabrlo(5),rlp(5),rab2(nr),fa(nr),fb(nr)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c IBM
      expzer = 3.7D2
cIris     expzer = 3.7E2
cApollo   expzer = 3.7D2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer =  2.8E3
c
c     for numerical stability:
c
      expzer = expzer/2
c
c      integration coefficients
c
       abc1 = 1901.D0/720.D0
       abc2 = -1387.D0/360.D0
       abc3 = 109.D0/30.D0
       abc4 = -637.D0/360.D0
       abc5 = 251.D0/720.D0
       amc0 = 251.D0/720.D0
       amc1 = 323.D0/360.D0
       amc2 = -11.D0/30.D0
       amc3 = 53.D0/360.D0
       amc4 = -19.D0/720.D0
      itmax = 100
      lp = lo(iorb)+1
      ar(1) = 0.0d0
      if (lo(iorb) .eq. 0) then
        br(1) = b*a
      else
        br(1) = 0.0d0
      endif
      do 1 j=2,nr
        ar(j) = 0.0d0
 1    continue
      do 2 j=2,nr
        br(j) =0.0d0
 2    continue
      do 4 j=2,5
        rlp(j)=r(j)**lp
 4    continue
      do 5 j=2,5
        rabrlo(j)=rab(j)*r(j)**lo(iorb)
 5    continue
      do 6 j=1,nr
        rab2(j)=rab(j)*rab(j)
 6    continue
c
c   set underflow trap
c
      juflow=1
      do 42 j=2,nr
        if (lp*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
c
c   determine effective charge and vzero for startup of
c   outward integration
c   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
c   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
c
      zeff = 0.0d0
      if (so(iorb) .lt. 0.1 .and. viod(lp,2) .lt. -0.1) zeff=znuc
      if (so(iorb) .gt. 0.1 .and. viou(lp,2) .lt. -0.1) zeff=znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa
      if (zeff .eq. 0.0) then
        if (so(iorb) .lt. 0.1 ) then
          vzero=vzero+viod(lp,2)/r(2)
        else
          vzero=vzero+viou(lp,2)/r(2)
        endif
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vzero+vid(2)
      else
        vzero=vzero+viu(2)
      endif
      var0 = 0.0d0
      if (lo(iorb) .eq. 0) var0=-2*zeff
      if (lo(iorb) .eq. 1) var0=2.0d0
      emax = 0.0d0
      emin = -200000.0d0
      if (ev(iorb) .gt. emax) ev(iorb) = emax
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) return
      if (ev(iorb) .gt. 0.0) then
        write(6,1000)iorb
        stop 'difnrl one'
      endif
 1000 format(//,' error in difnrl - ev(',i2,
     1 ') greater then v(infinty)')
c
c   find practical infinity ninf and classical turning
c   point nctp for orbital
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) -ev(iorb)
        if (temp .lt. 0.0) temp = 0.0d0
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue
      if (ev(iorb) .ge. etol*10) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=0.0d0
      if (nctp .le. 6) then
        ev(iorb) = 0.9d0*ev(iorb)
        if (icount .gt. 100) then
          write(6,1010)iorb
          stop 'difnrl two'
        endif
        goto 20
      endif
 1010 format(//,'error in difnrl - cannot find the classical '
     1 ,/' turning point for orbital ',i2)
c
c   outward integration from 1 to nctp
c   startup
c
      bb = (vzero-ev(iorb))/(4*lp+2)
      do 35 j=2,5
        ar(j) = rlp(j) * (1+(aa+bb*r(j))*r(j))
        br(j) = rabrlo(j) * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
 35   continue
c
c    Predictor-corrector array added.
c
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-ev(iorb))*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-ev(iorb))*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-ev(iorb))*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-ev(iorb))*ar(5)
c
c   intergration loop
c
      nodes = 0
      do 40 j=6,nctp
c
c   predictor (Adams-Bashforth)
c
        j1=j-1
        j2=j-2
        j3=j-3
        j4=j-4
        j5=j-5
        vev=v(j)-ev(iorb)
        arp = ar(j1) + abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5)
        brp = br(j1) + abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5)
        fb1 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) + amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4)
        brc = br(j1) + amc0*fb1+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4)
        fb0 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc + amc0*(brc-brp)
        br(j) = brc + amc0*(fb0-fb1)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
c
c   count nodes - if no underflow
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.0.0)nodes=nodes+1
 40   continue
c
      arctp = ar(nctp)
      brctp = br(nctp)
c
c   end outward integration
c
c   if number of nodes correct, start inward integration
c   else modify energy stepwise and try again
c
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
c     c.hartwig
c         write(6,*) 'nodes,ev(iorb)',nodes,ev(iorb)
        if (nodes .lt. no(iorb)-lo(iorb)-1) then
c
c  too few nodes; increase ev
c
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        else
c
c  too many nodes; decrease ev
c
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c   inward integration from ninf to nctp
c   startup
c
      do 71 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. 0.0) alf = 0.0d0
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = -rab(j)*alf*ar(j)
 71   continue
c
c    Array for predictor-corrector added.
c
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*
     1 (v(ninf)-ev(iorb))*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*
     1       (v(ninf1)-ev(iorb))*ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*
     1       (v(ninf2)-ev(iorb))*ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*
     1       (v(ninf3)-ev(iorb))*ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*
     1       (v(ninf4)-ev(iorb))*ar(ninf4)
c
c   integration loop
c
      istop = ninf - nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c   predictor (Adams-Bashforth)
c
        j1 = j + 1
        j2 = j + 2
        j3 = j + 3
        j4 = j + 4
        j5 = j + 5
        vev = v(j)-ev(iorb)
        arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+
     1   abc4*fa(j4)+abc5*fa(j5))
        brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+
     1   abc4*fb(j4)+abc5*fb(j5))
        fb0 = b*brp + rab2(j)*vev*arp
c
c   corrector (Adams-Moulton)
c
        arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+
     1   amc3*fa(j3)+amc4*fa(j4))
        brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+
     1   amc3*fb(j3)+amc4*fb(j4))
c
        fb1 = b*brc + rab2(j)*vev*arc
c
c   error reduction step
c
        ar(j) = arc - amc0*(brc-brp)
        br(j) = brc - amc0*(fb1-fb0)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
 80   continue
c
c   end inward integration
c
c   rescale ar and br outside nctp to match ar(nctp) from
c   outward integration
c
  222 factor = arctp/ar(nctp)
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
c
c   find normalizing factor
c
      factor = 0.0d0
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*ar(j)*ar(j)*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c   modify eigenvalue ev
c
      dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=sign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) ev(iorb) = (evold + emax) / 2
      if (ev(iorb) .lt. emin) ev(iorb) = (evold + emin) / 2
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j) / rab(j)
 110  continue
 111  continue
      return
      end
c
c      *****************************************************************
c
      subroutine difrel(iter,iorb,v,ar,br,lmax,nr,a,b,r,rab,
     1 norb,no,lo,so,znuc,viod,viou,vid,viu,ev)
c
c  difrel integrates the relativistic Dirac equation
c  it finds the eigenvalue ev, the major and minor component
c  of the wavefunction, ar and br.  It uses an intial guess
c  for the eigenvalues from dsolv1
c
      implicit real*8 (a-h,o-z)
c
      parameter (ai=2*137.0360411d0)
c
c  Tolernce
c
      parameter (etol=-1.d-7)
      parameter (tol = 1.0d-14)
c
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr),
     1 no(norb),lo(norb),so(norb),viod(lmax,nr),viou(lmax,nr),
     2 vid(nr),viu(nr),ev(norb),rabkar(nr),rabai(nr),
     3 fa(nr),fb(nr)
c
      dimension rs(5)
c
c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c IBM
      expzer = 3.7D2
cIris     expzer =3.7E2
cApollo   expzer = 3.7E2
cSun      expzer = 3.7D2
cVax      expzer = 44.D0
Cray      expzer = 2.8E3
c
c     for numerical stability:
c
      expzer = expzer/2
c
c
c      integration coefficients
c
       abc1 = 1901.D0/720.D0
       abc2 = -1387.D0/360.D0
       abc3 = 109.D0/30.D0
       abc4 = -637.D0/360.D0
       abc5 = 251.D0/720.D0
       amc0 = 251.D0/720.D0
       amc1 = 323.D0/360.D0
       amc2 = -11.D0/30.D0
       amc3 = 53.D0/360.D0
       amc4 = -19.D0/720.D0
      itmax = 100
      ai2 = ai * ai
      az = znuc/(2*ai)
      ka = lo(iorb)+1
      if (so(iorb) .lt. 0.1 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c  determine effective charge and vzero for startup of
c  outward integration
c  ar = r**s * (1  + a1 r + a2 r**2 + ... )
c  br = r**s * (b0 + b1 r + b2 r**2 + ... )
c  s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
c  an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - ai**2) b(n-1))
c        / (n ai (2 s + n))
c  bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( ai (s + n + ka))
c
      s = sqrt(ka*ka-az*az)
      if (ka .gt. 0) then
        b0 = -az/(s+ka)
      else
        b0 = (s-ka)/az
      endif
      if (so(iorb) .lt. 0.1) then
        vzero=vid(2)
      else
        vzero=viu(2)
      endif
c
c    Loop data calculated only once.
c    Set ar() and br() to zero.
c
      do 1 j=1,nr
        ar(j) = 0.0d0
        br(j) = 0.0d0
 1    continue
      do 3 j=2,nr
        rabkar(j)=rab(j)*ka/r(j)
 3    continue
      do 4 j=2,nr
        rabai(j)=rab(j)/ai
 4    continue
      do 5 j=2,5
        rs(j)=r(j)**s
 5    continue
c
c  set the underflow trap
c
      juflow=1
      do 42 j=2,nr
        if (s*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
c

      emax = 0.0d0
      emin = -100000.0d0
      if (ev(iorb) .gt. emax) ev(iorb) = emax
 10   if (itmax .lt. 2) write(6,15) iorb,iter,ev(iorb),nodes
 15   format(' iorb =',i3,' iter =',i3,' ev =',e18.10,' nodes =',i2)
      if (itmax .eq. 0) return
      if (ev(iorb) .gt. 0.0) then
        write(6,1000)iorb
        stop 'difrel one'
      endif
 1000 format(//,' error in difrel - ev(',i2,
     1 ') greater then v(infinty)')
c
c  Find practical infinity ninf and classical turning
c  point nctp for orbital.
c
      icount=0
 20   icount=icount+1
      do 22 j=nr,2,-1
        temp = v(j) - ev(iorb)
        if (temp .lt. 0.0) temp = 0.0d0
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
 22   continue
 23   ninf=j
      nctp = ninf - 5
      do 25 j=2,ninf-5
        if (v(j) .lt. ev(iorb)) nctp = j
 25   continue
      if (ev(iorb) .ge. etol*100) nctp=ninf-5
      if (ev(iorb) .ge. etol) ev(iorb)=0.0d0

      if (nctp .le. 6) then
        ev(iorb) = 0.9d0*ev(iorb)
        if (icount .gt. 100) then
          write(6,1010)iorb
          stop 'difrel two'
        endif
        goto 20
      endif
 1010 format(//,'error in difrel - cannot find classical',
     1 /,'turning point in orbital ',i2)
c
c  Outward integration from 1 to nctp, startup.
c
      a1 = (az*(vzero-ev(iorb))-(s+1+ka)*(vzero-ev(iorb)-ai2)*b0)
     1   / (ai*(2*s+1))
      b1 = ((vzero-ev(iorb))-2*znuc*a1) / (ai*(s+1+ka))
      a2 = (az*(vzero-ev(iorb))*a1-(s+2+ka)*(vzero-ev(iorb)-ai2)*b1)
     1   / (2*ai*(2*s+2))
      b2 = ((vzero-ev(iorb))*a1-2*znuc*a2) / (ai*(s+2+ka))
      do 35 j=2,5
        ar(j) = rs(j) * (1 +(a1+a2*r(j))*r(j))
        br(j) = rs(j) * (b0+(b1+b2*r(j))*r(j))
 35   continue
      fa(1) = 0.0d0
      fb(1) = 0.0d0
      fa(2) = rabkar(2)*ar(2)+(ev(iorb)-v(2)+ai2)*br(2)*rabai(2)
      fb(2) = -rabkar(2)*br(2)-(ev(iorb)-v(2))*ar(2)*rabai(2)
      fa(3) = rabkar(3)*ar(3)+(ev(iorb)-v(3)+ai2)*br(3)*rabai(3)
      fb(3) = -rabkar(3)*br(3)-(ev(iorb)-v(3))*ar(3)*rabai(3)
      fa(4) = rabkar(4)*ar(4)+(ev(iorb)-v(4)+ai2)*br(4)*rabai(4)
      fb(4) = -rabkar(4)*br(4)-(ev(iorb)-v(4))*ar(4)*rabai(4)
      fa(5) = rabkar(5)*ar(5)+(ev(iorb)-v(5)+ai2)*br(5)*rabai(5)
      fb(5) = -rabkar(5)*br(5)-(ev(iorb)-v(5))*ar(5)*rabai(5)
c
c  Intergration loop.
c
      nodes = 0
      do 40 j=6,nctp
c
c  Predictor (Adams-Bashforth).
c
        evvai2=ev(iorb)-v(j)+ai2
        evv=ev(iorb)-v(j)
        arp = ar(j-1) + abc1*fa(j-1)+abc2*fa(j-2)+abc3*fa(j-3)
     1   +abc4*fa(j-4)+abc5*fa(j-5)
        brp = br(j-1) + abc1*fb(j-1)+abc2*fb(j-2)+abc3*fb(j-3)
     1   +abc4*fb(j-4)+abc5*fb(j-5)
        fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
        fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
        arc = ar(j-1) + amc0*fa(j)+amc1*fa(j-1)+amc2*fa(j-2)
     1   +amc3*fa(j-3)+amc4*fa(j-4)
        brc = br(j-1) + amc0*fb(j)+amc1*fb(j-1)+amc2*fb(j-2)
     1   +amc3*fb(j-3)+amc4*fb(j-4)
        faj = rabkar(j)*arc+evvai2*brc*rabai(j)
        fbj = -rabkar(j)*brc-evv*arc*rabai(j)
c
c  Error reduction step.
c
        ar(j) = arc + amc0*(faj-fa(j))
        br(j) = brc + amc0*(fbj-fb(j))
        fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
        fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
c
c  Count nodes - if no underflow.
c
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.0.0)nodes=nodes+1
 40   continue
       arout = ar(nctp)
       arpout = fa(nctp)
c
c  End outward integration.
c  If number of nodes correct, start inward integration
c  else modify energy stepwise and try again.
c
      if (nodes .ne. no(iorb)-lo(iorb)-1) then
c
c  too many nodes decrease ev
c
        if (nodes .gt. no(iorb)-lo(iorb)-1) then
          if (ev(iorb) .lt. emax) emax = ev(iorb)
          ev(iorb) = ev(iorb) + ev(iorb)/10
c
c  too few nodes increase ev
c
        else
          if (ev(iorb) .gt. emin) emin = ev(iorb)
          ev(iorb) = ev(iorb) - ev(iorb)/10
        endif
        itmax = itmax-1
        goto 10
      endif
c
c  Inward integration from ninf to nctp startup.
c
      do 70 j=ninf,ninf-4,-1
        alf = v(j) - ev(iorb)
        if (alf .lt. 0.0) alf = 0.0d0
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = ai*(alf+ka/r(j))*ar(j)/(v(j)-ev(iorb)-ai2)
 70   continue
      fa(ninf) = rabkar(ninf)*ar(ninf)+
     1    (ev(iorb)-v(ninf)+ai2)*br(ninf)*rabai(ninf)
      fb(ninf) = -rabkar(ninf)*br(ninf)
     1    -(ev(iorb)-v(ninf))*ar(ninf)*rabai(ninf)
      fa(ninf-1) = rabkar(ninf-1)*ar(ninf-1)+
     1    (ev(iorb)-v(ninf-1)+ai2)*br(ninf-1)*rabai(ninf-1)
      fb(ninf-1) = -rabkar(ninf-1)*br(ninf-1)
     1    -(ev(iorb)-v(ninf-1))*ar(ninf-1)*rabai(ninf-1)
      fa(ninf-2) = rabkar(ninf-2)*ar(ninf-2)
     1    +(ev(iorb)-v(ninf-2)+ai2)*br(ninf-2)*rabai(ninf-2)
      fb(ninf-2) = -rabkar(ninf-2)*br(ninf-2)
     1    -(ev(iorb)-v(ninf-2))*ar(ninf-2)*rabai(ninf-2)
      fa(ninf-3) = rabkar(ninf-3)*ar(ninf-3)
     1    +(ev(iorb)-v(ninf-3)+ai2)*br(ninf-3)*rabai(ninf-3)
      fb(ninf-3) = -rabkar(ninf-3)*br(ninf-3)
     1    -(ev(iorb)-v(ninf-3))*ar(ninf-3)*rabai(ninf-3)
      fa(ninf-4) = rabkar(ninf-4)*ar(ninf-4)
     1    +(ev(iorb)-v(ninf-4)+ai2)*br(ninf-4)*rabai(ninf-4)
      fb(ninf-4) = -rabkar(ninf-4)*br(ninf-4)
     1    -(ev(iorb)-v(ninf-4))*ar(ninf-4)*rabai(ninf-4)
c
c  Integration loop.
c
      istop = ninf-nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
c
c  Predictor (Adams-Bashforth).
c
        evvai2=ev(iorb)-v(j)+ai2
        evv=ev(iorb)-v(j)
        arp = ar(j+1)-(abc1*fa(j+1)+abc2*fa(j+2)+abc3*fa(j+3)
     1   +abc4*fa(j+4)+abc5*fa(j+5))
        brp = br(j+1)-(abc1*fb(j+1)+abc2*fb(j+2)+abc3*fb(j+3)
     1   +abc4*fb(j+4)+abc5*fb(j+5))
        fa(j) = rabkar(j)*arp+evvai2*brp*rabai(j)
        fb(j) = -rabkar(j)*brp-evv*arp*rabai(j)
c
c  Corrector (Adams-Moulton).
c
        arc = ar(j+1)-(amc0*fa(j)+amc1*fa(j+1)+amc2*fa(j+2)
     1   +amc3*fa(j+3)+amc4*fa(j+4))
        brc = br(j+1)-(amc0*fb(j)+amc1*fb(j+1)+amc2*fb(j+2)
     1   +amc3*fb(j+3)+amc4*fb(j+4))
        faj = rabkar(j)*arc+evvai2*brc*rabai(j)
        fbj = -rabkar(j)*brc-evv*arc*rabai(j)
c
c  Error reduction step.
c
        ar(j) = arc + amc0*(faj-fa(j))
        br(j) = brc + amc0*(fbj-fb(j))
        fa(j) = rabkar(j)*ar(j)+evvai2*br(j)*rabai(j)
        fb(j) = -rabkar(j)*br(j)-evv*ar(j)*rabai(j)
 80   continue
 222  arin = ar(nctp)
      arpin = fa(nctp)
c
c  End inward integration
c  Rescale ar and br outside nctp to match ar(nctp) from
c  outward integration.
c
      factor = arout/arin
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
      arpin = factor * arpin
c
c  Find the normalizing factor.
c
      factor = 0.0d0
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
c
c  Modify the eigenvalue ev.
c
      dev = arout * (arpout-arpin) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev(iorb)) dev=dsign(ev(iorb),dev)/5
      itmax = itmax-1
      evold = ev(iorb)
      ev(iorb) = ev(iorb) + dev
      if (ev(iorb) .gt. emax) then
        ev(iorb) = (evold + emax) / 2
      elseif (ev(iorb) .lt. emin) then
        ev(iorb) = (evold + emin) / 2
      endif
      if (abs(dev) .gt. tol*(1-ev(iorb))) goto 10
c
c  Normalize the wavefunction.
c
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j)
 110  continue
 111  continue
      return
      end
c
c      *****************************************************************
c
       subroutine orban
     + (itype,icorr,ispp,iorb,ar,br,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,v,ev,ek,ep,rcov,rprb,nconf)
       implicit double precision(a-h,o-z)
c
c      orban is used to analyze and printout data about the orbital
c
       dimension ar(nr),br(nr)
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 v(nr),etot(10),ev(norb),ek(norb),ep(norb)
       character*2 ispp*1,nameat,itype
c
       dimension rzero(10),rextr(10),aextr(10),bextr(10)
       dimension cg(100),gzero(10),gextr(10),cextr(10)
       character*10 name
       character icorr*10,irel*3,ifcore*4
c.....files
      common /files/iinput,iout,in290,in213,istore,iunit7,iunit8,istruc,
     +               ivnlkk,isumry,ikpts
c     c.hartwig
c     work-arrays for integration, and xc-potential
      dimension ttx(50000),tty(50000),ttyp(50000),ttypp(50000),
     :     ttw(150000),rho(nr),vxcgrd(nr),excgrd(nr)
      dimension rr(10000),rw(10000),rd(10000)
      common /intgrd/ rw,rd
      character*1 il(5)
      character*2 cnum
c
c
c       ai = 2*137.04D0
       ai=2*137.0360411d0
       pi = 4.D0 * atan(1.D0)
       ka = lo(iorb)+1
       lp = ka
       if (so(iorb) .lt. 0.1D0 .and. lo(iorb) .ne. 0) ka=-lo(iorb)
c
c      compute zeroes and extrema
c
       nzero = 0
       nextr = 0
       rzero(1) = 0.D0
       arp = br(2)
       if (ispp .eq. 'r' .and. so(iorb) .lt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(iorb) - viod(lp,2)/r(2) - vid(2) + ai*ai) * br(2) / ai
       if (ispp .eq. 'r' .and. so(iorb) .gt. 0.1D0) arp = ka*ar(2)/r(2)
     1  + (ev(iorb) - viou(lp,2)/r(2) - viu(2) + ai*ai) * br(2) / ai
       do 20 i=3,nr
       if (nextr .ge. no(iorb)-lo(iorb)) goto 30
       if (ar(i)*ar(i-1) .gt. 0.D0) goto 10
c
c      zero
c
       nzero = nzero + 1
       rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i)) / (ar(i)-ar(i-1))
 10    arpm = arp
       arp = br(i)
       if (ispp .eq. 'r' .and. so(iorb) .lt. 0.1D0) arp = ka*ar(i)/r(i)
     1  + (ev(iorb) - viod(lp,i)/r(i) - vid(i) + ai*ai) * br(i) / ai
       if (ispp .eq. 'r' .and. so(iorb) .gt. 0.1D0) arp = ka*ar(i)/r(i)
     1  + (ev(iorb) - viou(lp,i)/r(i) - viu(i) + ai*ai) * br(i) / ai
       if (arp*arpm .gt. 0.D0) goto 20
c
c      extremum
c
       nextr = nextr + 1
       rextr(nextr) = (arp*r(i-1)-arpm*r(i)) / (arp-arpm)
       aextr(nextr) = (ar(i)+ar(i-1))/2
     1 - (arp**2+arpm**2) * (r(i)-r(i-1)) / (4*(arp-arpm))
       bextr(nextr) = br(i)
 20    continue
c
c      find orbital kinetic and potential energy
c      the potential part includes only the interaction with
c      the nuclear part
c
 30    ek(iorb) = br(1)*br(1)*rab(1)
       ep(iorb) = 0.D0
       sa2 = 0.D0
       lp = lo(iorb)+1
       llp = lo(iorb)*lp
       ll = 2
       if (2*(nr/2) .eq. nr) ll=4
       do 40 ii=2,nr
       i = nr-ii+2
       ar2 = ar(i)*ar(i)
       br2 = br(i)*br(i)
       deni = ar2
       if (ispp .eq. 'r') deni=deni+br2
       ek(iorb) = ek(iorb) + ll * (br2 + ar2*llp/r(i)**2)*rab(i)
       if (so(iorb) .lt. 0.1D0) ep(iorb) = ep(iorb)
     1    + ll * deni*viod(lp,i)*rab(i)/r(i)
       if (so(iorb) .gt. 0.1D0) ep(iorb) = ep(iorb)
     1    + ll * deni*viou(lp,i)*rab(i)/r(i)
       ll = 6 - ll
       if (sa2 .gt. 0.10D0) goto 40
       sa2 = sa2 + deni*rab(i)
       if (sa2 .le. 0.01D0) i99 = i
       i90 = i
 40    continue
       ek(iorb) = ek(iorb) / 3
       ep(iorb) = ep(iorb) / 3
       if (ispp .eq. 'r') ek(iorb) = 0.D0
c
c      fourier analyze orbital
c
c       if (iorb .lt. ncore) return
c       kzero = 0
c       kextr = 0
c       iextr = 1
c       delg = 0.2D0*pi/r(i90)
c       do 60 i=1,100
c       g = delg * (i-1)
c       cg(i) = 0.D0
c       if (i .eq. 1 .and. lp .ne. 1) goto 60
c       ll = 4
c       do 50 j=2,nr
c       rg = r(j) * g
c       bsl = 1.D0
c       if (i  .ne. 1) bsl = sin(rg) / rg
c       if (lp .eq. 2) bsl = (bsl - cos(rg)) / rg
c       if (lp .eq. 3) bsl = 3.D0 * (bsl - cos(rg)) / rg**2 -  bsl
c       cg(i) = cg(i) + ll * r(j) * ar(j) * bsl * rab(j)
c       ll = 6 - ll
c 50    continue
c       cg(i) = cg(i) / (6.D0*pi**2)
cc      write(6,'(2i3,3f13.6)') lo(iorb),i,g,cg(i),cg(i)*g**2
c       if (i .eq. 1) goto 60
c
c      find extremum
c
c       if (abs(cg(i)) .gt. abs(cg(iextr))) iextr = i
c       if (i .eq. 2) goto 60
c
c      zero
c
c       if (cg(i)*cg(i-1) .gt. 0.D0) goto 60
c
c      zero found - update arrays
c
c       if (i-iextr .lt. 4) goto 70
c       kzero = kzero + 1
c       gzero(kzero) = delg*(cg(i)*(i-2)-cg(i-1)*(i-1))/(cg(i)-cg(i-1))
c       kextr = kextr + 1
c       cextr(kextr) = Dlog10(abs(cg(iextr)))
c       gextr(kextr) = delg * (iextr-1)
c       if (kextr .eq. 5) goto 70
c       iextr = i
c 60    continue
c       kextr = kextr + 1
c       cextr(kextr) = Dlog10(abs(cg(iextr)))
c       gextr(kextr) = delg * iextr
c
c      printout
c
c      vshift=-15.d0
c 70    if (iorb .lt. ncore) return
c       write(6,80) no(iorb),lo(iorb),so(iorb)
c 80    format(/' n =',i2,'  l =',i2,'  s =',f4.1)
c       write(6,90)  (ev(iorb)-vshift)/2.,ek(iorb)/2.,ep(iorb)/2.
c 90    format(8x,'ev =',e15.8,'  ek =',e14.8,'  ep =',e15.8)
c       name = 'a extr    '
c       write(6,100) name,(aextr(i),i=1,nextr)
c       name = 'b extr    '
c       if (ispp .eq. 'r') write(6,100) name,(bextr(i),i=1,nextr)
c       name = 'r extr    '
c       write(6,100) name,(rextr(i),i=1,nextr)
c       name = 'r zero    '
c       write(6,100) name,(rzero(i),i=1,nzero)
c       name = 'r 90/99 % '
c       write(6,100) name,r(i90),r(i99)
c       name = 'c extr lg '
c       write(6,100) name,(cextr(i),i=1,kextr)
c       name = 'g extr    '
c       write(6,100) name,(gextr(i),i=1,kextr)
c       name = 'g zero    '
c       write(6,100) name,(gzero(i),i=1,kzero)
c 100   format(8x,a10,8f8.3)

c------Machine dependent parameter-
c------Require exp(-2*expzer) to be within the range of the machine
c IBM
      expzer = 3.7D2
c     c.hartwig for numerical stability:
      expzer = expzer/2
c
c  Find practical infinity ninf and classical turning
c  point nctp for orbital.
      do  j=nr,2,-1
         temp = v(j) - ev(iorb)
         if (temp .lt. 0.0) temp = 0.0
         if (r(j)*sqrt(temp) .lt. expzer) goto 23
      enddo
 23   ninf=j
c
c     compute charge at rcov + higher moments
c     spline interpolation/integration

c     some additional points for the spline
      npoint=min(ninf+5,nr)
c     charge(rcov)= int_0^rcov g^2 r^2 dr + int_0^infinity f^2 r^2 dr
c
      a1=0
      an=0
      b1=0
      bn=0
      isx=0
c     int_0^rcov g^2 r^2 dr
      do i=1,npoint
         ttx(i)=r(i)
         tty(i)=ar(i)*ar(i)
         if (r(i).le.rcov) ircov=i
      enddo
      if (ircov.gt.ninf) then
         ircov=ninf
         write(6,*) 'warning: ircov > ninf ! (ircov set to ninf)'
      endif
      call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,a1,b1,an,bn)
      if(ierr.ne.1) stop 'splift'
      isx=1
      ttxup=ttx(ircov)
      call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,crcov,ierr)
      if(ierr.ne.1) stop 'spliq'
      if (ispp .eq. 'r') then
c     int_0^infinity f^2 r^2 dr
         cmin=0.
         do i=1,npoint
            tty(i) = br(i)*br(i)
         enddo
         call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,
     :        a1,b1,an,bn)
         if(ierr.ne.1) stop 'splift'
c         ttxup=ttx(ircov)
c         call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,cmin,ierr)
c         if(ierr.ne.1) stop 'spliq'
c         print*,'crcov+cmin:',crcov+cmin
         ttxup=ttx(ninf)
         call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,cmin,ierr)
         if(ierr.ne.1) stop 'spliq'
c         print*,'crcov+cmin:',crcov+cmin
         crcov=crcov+cmin
      endif
c
c     dcharge      = int_0^infinity (f^2+g^2) r^4 dr
c
      ttxup=ttx(ninf)
      do i=1,npoint
         tty(i)=ar(i)*ar(i)
         if (ispp.eq.'r')tty(i)=tty(i)+br(i)*br(i)
         tty(i)=tty(i)*r(i)*r(i)
      enddo
      call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,a1,b1,an,bn)
      if(ierr.ne.1) stop 'splift'
c     ddd =  = int_0^rcov (f^2+g^2) r^4 dr
      call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttx(ircov),
     :     1,ddd,ierr)
      if(ierr.ne.1) stop 'spliq'
      call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,dcrcov,ierr)
      if(ierr.ne.1) stop 'spliq'
c
c     int_0^infinity (f^2+g^2) r^6 dr
c
      do i=1,npoint
         tty(i)=tty(i)*r(i)*r(i)
      enddo
      call splift(ttx,tty,ttyp,ttypp,npoint,ttw,ierr,isx,a1,b1,an,bn)
      if(ierr.ne.1) stop 'splift'
      call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttxup,1,ddcrcov,ierr)
      if(ierr.ne.1) stop 'spliq'
c     dddd =  = int_0^rcov (f^2+g^2) r^6 dr
      call spliq(ttx,tty,ttyp,ttypp,npoint,ttxlo,ttx(ircov),
     :     1,dddd,ierr)
      if(ierr.ne.1) stop 'spliq'
c
c   printout
c
      vshift=-15d0
      il(1) = 's'
      il(2) = 'p'
      il(3) = 'd'
      il(4) = 'f'
      il(5) = 'g'

      if (iorb.eq.1)then
         write(6,*)
         write(6,*) 'rcov         = ',rcov
         if (ispp .ne. 'r' ) then
            write(6,*) 'charge(rcov) = int_0^rcov psi^2 r^2 dr'
            write(6,*) 'dcharge      = int_0^infinity psi^2 r^4 dr'
            write(6,*) 'ddcharge     = int_0^infinity psi^2 r^6 dr'
         else
            write(6,*) 'charge(rcov) = int_0^rcov g^2 r^2 dr ',
     1           '  +  int_0^infinity f^2 r^2 dr'
            write(6,*) 'dcharge      = int_0^infinity (f^2+g^2) r^4 dr '
            write(6,*) 'ddcharge     = int_0^infinity (f^2+g^2) r^6 dr '
         endif
         write(6,21)
 21      format(/,' nl   s    occ',4x,'eigenvalue',3x,'charge(rcov)',
     1        4 x,'dcharge',4x,'ddcharge')
      endif
      write(6,31) no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),
     :     (ev(iorb)-vshift)/2.,crcov,dcrcov,ddcrcov
 31   format(1x,i1,a1,f4.1,f8.3,2e14.7,2e12.5)
c      write(6,*) 'drcov at rcov :',ddd
c      write(6,*) 'ddrcov at rcov:',dddd
       name = 'r extr    '
       write(6,100) name,(rextr(i),i=1,nextr)
 100   format(5x,a10,9f7.2)
c
c     write data to files atom.ae for pseudopotential-fit
c     only valence electrons
      if (ispp.ne.'r') ispp='n'
      if (iorb.gt.ncore) then
         if (iorb.eq.ncore+1)then
            zps=0
            do jj=iorb,norb
               zps=zps+zo(jj)
            enddo
            if (nconf .eq. 0 ) then
               write(40,*) norb-ncore,'Number of orbitals'
               write(40,*) znuc,zps,rcov,rprb,
     :              'znuc, zpseudo, rcov, rprb'
               if (ispp.eq.'r') then
                  write(40,*)'relativistic calculation'
               else
                  write(40,*)'non relativistic calculation'
               endif
               write(40,'(a10,a)') icorr, '   XC-functional'
               write(40,*) nr,        'number of gridpoints'
            else
               write(40,'(a,i2,a)') ' NEXT CONFIGURATION (',nconf,')'
            endif
         endif
c
         write(40,*) no(iorb),lo(iorb),so(iorb),zo(iorb),
     :        (ev(iorb)-vshift)/2.,crcov,dcrcov,ddcrcov,
     :        ' n l s z eval, charge, dcharge, ddcharge residue'
         do i=1,nr
            if (ispp.eq.'r') then
               write(40,*) r(i),ar(i),br(i)
            else
               write(40,*) r(i),ar(i)
            endif
         enddo
      endif

c     c.chartwig:
c     save data for plots
      write(cnum,'(i2)') nconf
      inum=2
      if (nconf.gt.9) inum=1
      if (ispp.eq.'r') then
         open(unit=33,file='ae.'//
     :        char(ichar('0')+no(iorb))//
     :        il(lo(iorb)+1)//
     :        char(ichar('0')+int(2*(lo(iorb)+so(iorb))))//'by2'//
     :        '.conf'//cnum(inum:2)//'.dat',
     :        form='formatted',status='unknown')
c     :        char(ichar('A')-ichar('a')+ichar(il(lo(iorb)+1)))//
         else
         open(unit=33,file='ae.'//
     :        char(ichar('0')+no(iorb))//il(lo(iorb)+1)//
     :        '.conf'//cnum(inum:2)//'.dat',
     :        form='formatted',status='unknown')
      endif
      dena=0
      denb=0
      if (ispp.eq.'r') then
         write(33,*)'# r , major , minor , den(major) , den(minor)'
      else
         write(33,*)'# r , psi , den(major)'
      endif
      do i=1,npoint
         dena = dena +  ar(i)*ar(i)*rab(i)
         denb = denb +  br(i)*br(i)*rab(i)
         if (ispp.eq.'r') then
            write(33,'(20e20.10)') r(i),ar(i),br(i),dena,denb
         else
            write(33,'(20e20.10)') r(i),ar(i),dena
         endif
      enddo
      close(33)




       return
       end
c
c      *****************************************************************
c
      subroutine splift(x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)
      implicit double precision(a-h,o-z)
c
c     sandia mathematical program library
c     applied mathematics division 2613
c     sandia laboratories
c     albuquerque, new mexico  87185
c     control data 6600/7600  version 7.2  may 1978
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    issued by sandia laboratories                     *
c  *                   a prime contractor to the                       *
c  *                united states department of energy                 *
c  * * * * * * * * * * * * * * * notice  * * * * * * * * * * * * * * * *
c  * this report was prepared as an account of work sponsored by the   *
c  * united states government.  neither the united states nor the      *
c  * united states department of energy nor any of their employees,    *
c  * nor any of their contractors, subcontractors, or their employees  *
c  * makes any warranty, express or implied, or assumes any legal      *
c  * liability or responsibility for the accuracy, completeness or     *
c  * usefulness of any information, apparatus, product or process      *
c  * disclosed, or represents that its use would not infringe          *
c  * owned rights.                                                     *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * the primary document for the library of which this routine is     *
c  * part is sand77-1441.                                              *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     written by rondall e. jones
c
c     abstract
c         splift fits an interpolating cubic spline to the n data points
c         given in x and y and returns the first and second derivatives
c         in yp and ypp.  the resulting spline (defined by x, y, and
c         ypp) and its first and second derivatives may then be
c         evaluated using splint.  the spline may be integrated using
c         spliq.  for a smoothing spline fit see subroutine smoo.
c
c     description of arguments
c         the user must dimension all arrays appearing in the call list,
c         e.g.   x(n), y(n), yp(n), ypp(n), w(3n)
c
c       --input--
c
c         x    - array of abscissas of data (in increasing order)
c         y    - array of ordinates of data
c         n    - the number of data points.  the arrays x, y, yp, and
c                ypp must be dimensioned at least n.  (n .ge. 4)
c         isx  - must be zero on the initial call to splift.
c                if a spline is to be fitted to a second set of data
c                that has the same set of abscissas as a previous set,
c                and if the contents of w have not been changed since
c                that previous fit was computed, then isx may be
c                set to one for faster execution.
c         a1,b1,an,bn - specify the end conditions for the spline which
c                are expressed as constraints on the second derivative
c                of the spline at the end points (see ypp).
c                the end condition constraints are
c                        ypp(1) = a1*ypp(2) + b1
c                and
c                        ypp(n) = an*ypp(n-1) + bn
c                where
c                        abs(a1).lt. 1.0  and  abs(an).lt. 1.0.
c
c                the smoothest spline (i.e., least integral of square
c                of second derivative) is obtained by a1=b1=an=bn=0.
c                in this case there is an inflection at x(1) and x(n).
c                if the data is to be extrapolated (say, by using splint
c                to evaluate the spline outside the range x(1) to x(n)),
c                then taking a1=an=0.5 and b1=bn=0 may yield better
c                results.  in this case there is an inflection
c                at x(1) - (x(2)-x(1)) and at x(n) + (x(n)-x(n-1)).
c                in the more general case of a1=an=a  and b1=bn=0,
c                there is an inflection at x(1) - (x(2)-x(1))*a/(1.0-a)
c                and at x(n) + (x(n)-x(n-1))*a/(1.0-a).
c
c                a spline that has a given first derivative yp1 at x(1)
c                and ypn at y(n) may be defined by using the
c                following conditions.
c
c                a1=-0.5
c
c                b1= 3.0*((y(2)-y(1))/(x(2)-x(1))-yp1)/(x(2)-x(1))
c
c                an=-0.5
c
c                bn=-3.0*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)/(x(n)-x(n-1))
c
c       --output--
c
c         yp   - array of first derivatives of spline (at the x(i))
c         ypp  - array of second derivatives of spline (at the x(i))
c         ierr - a status code
c              --normal code
c                 1 means that the requested spline was computed.
c              --abnormal codes
c                 2 means that n, the number of points, was .lt. 4.
c                 3 means the abscissas were not strictly increasing.
c
c       --work--
c
c         w    - array of working storage dimensioned at least 3n.
      dimension x(n),y(n),yp(n),ypp(n),w(n,3)
c
      if (n.lt.4) go to 200
      nm1  = n-1
      nm2  = n-2
      if (isx.gt.0) go to 40
      do 5 i=2,n
      if (x(i)-x(i-1)) 300,300,5
    5 continue
c
c     define the tridiagonal matrix
c
      w(1,3) = x(2)-x(1)
      do 10 i=2,nm1
      w(i,2) = w(i-1,3)
      w(i,3) = x(i+1)-x(i)
   10 w(i,1) = 2.D0*(w(i,2)+w(i,3))
      w(1,1) = 4.D0
      w(1,3) =-4.D0*a1
      w(n,1) = 4.D0
      w(n,2) =-4.D0*an
c
c     l u decomposition
c
      do 30 i=2,n
      w(i-1,3) = w(i-1,3)/w(i-1,1)
   30 w(i,1)   = w(i,1) - w(i,2)*w(i-1,3)
c
c     define *constant* vector
c
   40 ypp(1) = 4.D0*b1
      dold   = (y(2)-y(1))/w(2,2)
      do 50 i=2,nm2
      dnew   = (y(i+1) - y(i))/w(i+1,2)
      ypp(i) = 6.D0*(dnew - dold)
      yp(i)  = dold
   50 dold   = dnew
      dnew   = (y(n)-y(n-1))/(x(n)-x(n-1))
      ypp(nm1) = 6.D0*(dnew - dold)
      ypp(n) = 4.D0*bn
      yp(nm1)= dold
      yp(n)  = dnew
c
c     forward substitution
c
      ypp(1) = ypp(1)/w(1,1)
      do 60 i=2,n
   60 ypp(i) = (ypp(i) - w(i,2)*ypp(i-1))/w(i,1)
c
c     backward substitution
c
      do 70 j=1,nm1
      i = n-j
   70 ypp(i) = ypp(i) - w(i,3)*ypp(i+1)
c
c     compute first derivatives
c
      yp(1)  = (y(2)-y(1))/(x(2)-x(1)) - (x(2)-x(1))*(2.D0*ypp(1)
     1         + ypp(2))/6.D0
      do 80 i=2,nm1
   80 yp(i)  = yp(i) + w(i,2)*(ypp(i-1) + 2.D0*ypp(i))/6.D0
      yp(n)  = yp(n) + (x(n)-x(nm1))*(ypp(nm1) + 2.D0*ypp(n))/6.D0
c
      ierr = 1
      return
  200 ierr = 2
      write(6,210)
  210 format(47h in splift, there were less than 4 data values.)
      return
  300 ierr = 3
      write(6,310)
  310 format(11h in splift,,
     144h the abscissas were not strictly increasing.)
      return
      end
c
c      *****************************************************************
c
      subroutine spliq(x,y,yp,ypp,n,xlo,xup,nup,ans,ierr)
      implicit double precision(a-h,o-z)
      dimension x(n),y(n),yp(n),ypp(n),xup(nup),ans(nup)
c
c     sandia mathematical program library
c     applied mathematics division 2613
c     sandia laboratories
c     albuquerque, new mexico  87185
c     control data 6600/7600  version 7.2  may 1978
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    issued by sandia laboratories                     *
c  *                   a prime contractor to the                       *
c  *                united states department of energy                 *
c  * * * * * * * * * * * * * * * notice  * * * * * * * * * * * * * * * *
c  * this report was prepared as an account of work sponsored by the   *
c  * united states government.  neither the united states nor the      *
c  * united states department of energy nor any of their employees,    *
c  * nor any of their contractors, subcontractors, or their employees  *
c  * makes any warranty, express or implied, or assumes any legal      *
c  * liability or responsibility for the accuracy, completeness or     *
c  * usefulness of any information, apparatus, product or process      *
c  * disclosed, or represents that its use would not infringe          *
c  * owned rights.                                                     *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * the primary document for the library of which this routine is     *
c  * part is sand77-1441.                                              *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     this routine was written by m. k. gordon
c
c     abstract
c
c     subroutine spliq integrates a cubic spline (generated by
c     splift, smoo, etc.) on the intervals (xlo,xup(i)), where xup
c     is a sequence of upper limits on the intervals of integration.
c     the only restrictions on xlo and xup(*) are
c                xlo .lt. xup(1),
c                xup(i) .le. xup(i+1)   for each i .
c     endpoints beyond the span of abscissas are allowed.
c     the spline over the interval (x(i),x(i+1)) is regarded
c     as a cubic polynomial expanded about x(i) and is integrated
c     analytically.
c
c     description of arguments
c         the user must dimension all arrays appearing in the call list,
c         e.g.  x(n), y(n), yp(n), ypp(n), xup(nup), ans(nup)
c
c      --input--
c
c        x    - array of abscissas (in increasing order) that define the
c               spline.  usually x is the same as x in splift or smoo.
c        y    - array of ordinates that define the spline.  usually y is
c               the same as y in splift or as r in smoo.
c        yp   - array of first derivatives of the spline at abscissas.
c               usually yp is the same as yp in splift or r1 in smoo.
c        ypp  - array of second derivatives that define the spline.
c               usually ypp is the same as ypp in splift or r2 in smoo.
c        n    - the number of data points that define the spline.
c        xlo  - left endpoint of integration intervals.
c        xup  - right endpoint or array of right endpoints of
c               integration intervals in ascending order.
c        nup  - the number of right endpoints.  if nup is greater than
c               1, then xup and ans must be dimensioned at least nup.
c
c      --output--
c
c        ans -- array of integral values, that is,
c               ans(i) = integral from xlo to xup(i)
c        ierr -- error status
c                = 1 integration successful
c                = 2 improper input - n.lt.4 or nup.lt.1
c                = 3 improper input - abscissas not in
c                        strictly ascending order
c                = 4 improper input - right endpoints xup not
c                        in ascending order
c                = 5 improper input - xlo.gt.xup(1)
c                = 6 integration successful but at least one endpoint
c                        not within span of abscissas
c
c   check for improper input
c
      ierr = 2
      if(n .ge. 4  .and.  nup .ge. 1) go to 1
      write(6,110)
 110  format(36h in spliq, either n.lt.4 or nup.lt.1)
      return
 1    nm1 = n-1
      nm2 = n-2
      ierr = 3
      do 2 i = 1,nm1
        if(x(i) .lt. x(i+1)) go to 2
        write(6,120)
 120    format(43h in spliq, abscissas not in ascending order)
        return
 2      continue
      if(nup .eq. 1) go to 4
      ierr = 4
      do 3 i = 2,nup
        if(xup(i-1) .le. xup(i)) go to 3
        write(6,130)
 130    format(49h in spliq, right endpoints not in ascending order)
        return
 3      continue
 4    ierr = 5
      if(xlo .le. xup(1)) go to 5
      write(6,140)
 140  format(26h in spliq, xlo .gt. xup(1))
      return
    5 ierr = 1
      if(xlo .lt. x(1)  .or.  xup(nup) .gt. x(n)) ierr = 6
c
c   locate xlo in interval (x(i),x(i+1))
c
      do 10 i = 1,nm2
        if(xlo .lt. x(i+1)) go to 20
 10     continue
      i = nm1
 20   hlo = xlo-x(i)
      hlo2 = hlo*hlo
      hi = x(i+1)-x(i)
      hi2 = hi*hi
      do 30 j = 1,nup
        if(xup(j) .gt. x(i+1)  .and.  xlo .lt. x(nm1)) go to 40
c
c   compute special cases of xup in interval with xlo
c
        hup = xup(j)-x(i)
        hsum = hup+hlo
        hdiff = hup-hlo
        hup2 = hup*hup
        sum = (ypp(i+1)-ypp(i))*hsum*hdiff*(hup2+hlo2)/(24.D0*hi)
        sum = sum + ypp(i)*hdiff*(hup2+hlo*hup+hlo2)/6.D0
        sum = sum + yp(i)*hdiff*hsum/2.D0
        sum = sum + y(i)*hdiff
 30     ans(j) = sum
      return
c
c   compute integral between xlo and x(i+1) as four terms in taylor
c   polynomial and advance i to i+1
c
 40   hdiff = hi-hlo
      hsum = hi+hlo
      sum0 = y(i)*hdiff
      sum1 = yp(i)*hdiff*hsum
      sum2 = ypp(i)*hdiff*(hi2+hi*hlo+hlo2)
      sum3 = (ypp(i+1)-ypp(i))*hdiff*hsum*(hi2+hlo2)/hi
      i = i+1
c
c   locate each xup(m) in interval (x(i),x(i+1))
c
      do 80 m = j,nup
 50     if(xup(m) .lt. x(i+1)  .or.  i .eq. nm1) go to 60
c
c   augment integral between abscissas to include interval
c   (x(i),x(i+1)) and advance i to i+1
c
        hi = x(i+1)-x(i)
        hi2 = hi*hi
        hi3 = hi2*hi
        sum0 = sum0 + y(i)*hi
        sum1 = sum1 + yp(i)*hi2
        sum2 = sum2 + ypp(i)*hi3
        sum3 = sum3 + (ypp(i+1)-ypp(i))*hi3
        i = i+1
        go to 50
c
c   integral between x(i) and xup(m) is zero
c
 60     if(xup(m) .ne. x(i)) go to 70
        sum = ((sum3/24.D0 + sum2/6.D0) + sum1/2.D0) + sum0
        go to 80
c
c   compute integral between x(i) and xup(m) and evaluate
c   taylor polynomial in reverse order
c
 70     hup = xup(m)-x(i)
        hup2 = hup*hup
        hup3 = hup2*hup
        hup4 = hup3*hup
        hi = x(i+1)-x(i)
        psum0 = y(i)*hup
        psum1 = yp(i)*hup2
        psum2 = ypp(i)*hup3
        psum3 = (ypp(i+1)-ypp(i))*hup4/hi
        sum = (sum3+psum3)/24.D0 + (sum2+psum2)/6.D0
        sum = sum + (sum1+psum1)/2.D0
        sum = sum + (sum0+psum0)
 80     ans(m) = sum
      return
      end
c
c      *****************************************************************
c
      subroutine tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5)
c
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm
      double precision d(n),e(n),e2(n),w(m),rv4(n),rv5(n)
      double precision u,v,lb,t1,t2,ub,xu,x0,x1,eps1,machep
c     real abs,max,min,DBLE
      integer ind(m)
c
c     this subroutine is a translation of the algol procedure bisect,
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971).
c
c     this subroutine finds those eigenvalues of a tridiagonal
c     symmetric matrix between specified boundary indices,
c     using bisection.
c
c     on input-
c
c        n is the order of the matrix,
c
c        eps1 is an absolute error tolerance for the computed
c          eigenvalues.  if the input eps1 is non-positive,
c          it is reset for each submatrix to a default value,
c          namely, minus the product of the relative machine
c          precision and the 1-norm of the submatrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary,
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary,
c
c        m11 specifies the lower boundary index for the desired
c          eigenvalues,
c
c        m specifies the number of eigenvalues desired.  the upper
c          boundary index m22 is then obtained as m22=m11+m-1.
c
c     on output-
c
c        eps1 is unaltered unless it has been reset to its
c          (last) default value,
c
c        d and e are unaltered,
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero,
c
c        lb and ub define an interval containing exactly the desired
c          eigenvalues,
c
c        w contains, in its first m positions, the eigenvalues
c          between indices m11 and m22 in ascending order,
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.,
c
c        ierr is set to
c          zero       for normal return,
c          3*n+1      if multiple eigenvalues at index m11 make
c                     unique selection impossible,
c          3*n+2      if multiple eigenvalues at index m22 make
c                     unique selection impossible,
c
c        rv4 and rv5 are temporary storage arrays.
c
c     note that subroutine tql1, imtql1, or tqlrat is generally faster
c     than tridib, if more than n/4 eigenvalues are to be found.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
      machep = 2.D0**(-47)
c
      ierr = 0
      tag = 0
      xu = d(1)
      x0 = d(1)
      u = 0.D0
c     ********** look for small sub-diagonal entries and determine an
c                interval containing all the eigenvalues **********
      do 40 i = 1, n
         x1 = u
         u = 0.D0
         if (i .ne. n) u = abs(e(i+1))
         xu = min(d(i)-(x1+u),xu)
         x0 = max(d(i)+(x1+u),x0)
         if (i .eq. 1) go to 20
         if (abs(e(i)) .gt. machep * (abs(d(i)) + abs(d(i-1))))
     x      go to 40
   20    e2(i) = 0.D0
   40 continue
c
      x1 = max(abs(xu),abs(x0)) * machep * DBLE(n)
      xu = xu - x1
      t1 = xu
      x0 = x0 + x1
      t2 = x0
c     ********** determine an interval containing exactly
c                the desired eigenvalues **********
      p = 1
      q = n
      m1 = m11 - 1
      if (m1 .eq. 0) go to 75
      isturm = 1
   50 v = x1
      x1 = xu + (x0 - xu) * 0.5D0
      if (x1 .eq. v) go to 980
      go to 320
   60 if (s - m1) 65, 73, 70
   65 xu = x1
      go to 50
   70 x0 = x1
      go to 50
   73 xu = x1
      t1 = x1
   75 m22 = m1 + m
      if (m22 .eq. n) go to 90
      x0 = t2
      isturm = 2
      go to 50
   80 if (s - m22) 65, 85, 70
   85 t2 = x1
   90 q = 0
      r = 0
c     ********** establish and process next submatrix, refining
c                interval by the gerschgorin bounds **********
  100 if (r .eq. m) go to 1001
      tag = tag + 1
      p = q + 1
      xu = d(p)
      x0 = d(p)
      u = 0.D0
c
      do 120 q = p, n
         x1 = u
         u = 0.D0
         v = 0.D0
         if (q .eq. n) go to 110
         u = abs(e(q+1))
         v = e2(q+1)
  110    xu = min(d(q)-(x1+u),xu)
         x0 = max(d(q)+(x1+u),x0)
         if (v .eq. 0.D0) go to 140
  120 continue
c
  140 x1 = max(abs(xu),abs(x0)) * machep
      if (eps1 .le. 0.D0) eps1 = -x1
      if (p .ne. q) go to 180
c     ********** check for isolated root within interval **********
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940
      m1 = p
      m2 = p
      rv5(p) = d(p)
      go to 900
  180 x1 = x1 * DBLE(q-p+1)
      lb = max(t1,xu-x1)
      ub = min(t2,x0+x1)
      x1 = lb
      isturm = 3
      go to 320
  200 m1 = s + 1
      x1 = ub
      isturm = 4
      go to 320
  220 m2 = s
      if (m1 .gt. m2) go to 940
c     ********** find roots by bisection **********
      x0 = ub
      isturm = 5
c
      do 240 i = m1, m2
         rv5(i) = ub
         rv4(i) = lb
  240 continue
c     ********** loop for k-th eigenvalue
c                for k=m2 step -1 until m1 do --
c                (-do- not used to legalize -computed go to-) **********
      k = m2
  250    xu = lb
c     ********** for i=k step -1 until m1 do -- **********
         do 260 ii = m1, k
            i = m1 + k - ii
            if (xu .ge. rv4(i)) go to 260
            xu = rv4(i)
            go to 280
  260    continue
c
  280    if (x0 .gt. rv5(k)) x0 = rv5(k)
c     ********** next bisection step **********
  300    x1 = (xu + x0) * 0.5D0
         if ((x0 - xu) .le. (2.D0 * machep *
     x      (abs(xu) + abs(x0)) + abs(eps1))) go to 420
c     ********** in-line procedure for sturm sequence **********
  320    s = p - 1
         u = 1.D0
c
         do 340 i = p, q
            if (u .ne. 0.D0) go to 325
            v = abs(e(i)) / machep
            if (e2(i) .eq. 0.D0) v = 0.D0
            go to 330
  325       v = e2(i) / u
  330       u = d(i) - x1 - v
            if (u .lt. 0.D0) s = s + 1
  340    continue
c
         go to (60,80,200,220,360), isturm
c     ********** refine intervals **********
  360    if (s .ge. k) go to 400
         xu = x1
         if (s .ge. m1) go to 380
         rv4(m1) = x1
         go to 300
  380    rv4(s+1) = x1
         if (rv5(s) .gt. x1) rv5(s) = x1
         go to 300
  400    x0 = x1
         go to 300
c     ********** k-th eigenvalue found **********
  420    rv5(k) = x1
      k = k - 1
      if (k .ge. m1) go to 250
c     ********** order eigenvalues tagged with their
c                submatrix associations **********
  900 s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
c
      do 920 l = 1, r
         if (j .gt. s) go to 910
         if (k .gt. m2) go to 940
         if (rv5(k) .ge. w(l)) go to 915
c
         do 905 ii = j, s
            i = l + s - ii
            w(i+1) = w(i)
            ind(i+1) = ind(i)
  905    continue
c
  910    w(l) = rv5(k)
         ind(l) = tag
         k = k + 1
         go to 920
  915    j = j + 1
  920 continue
c
  940 if (q .lt. n) go to 100
      go to 1001
c     ********** set error -- interval cannot be found containing
c                exactly the desired eigenvalues **********
  980 ierr = 3 * n + isturm
 1001 lb = t1
      ub = t2
      return
c     ********** last card of tridib **********
      end
c
c      *****************************************************************
c
      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      double precision d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      double precision u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,machep
c     real sqrt,abs,DBLE
      integer ind(m)
c     level 2, z
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input-
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement,
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary,
c
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0 if the eigenvalues are in ascending order, or 2.0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here,
c
c        m is the number of specified eigenvalues,
c
c        w contains the m eigenvalues in ascending or descending order,
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c
c     on output-
c
c        all input arrays are unaltered,
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero,
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations,
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c     ------------------------------------------------------------------
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
      machep = 2.D0**(-47)
c
      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.D0 - e2(1)
      q = 0
c     ********** establish and process next submatrix **********
  100 p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.D0) go to 140
  120 continue
c     ********** find vectors by inverse iteration **********
  140 tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
c     ********** check for isolated root **********
         xu = 1.D0
         if (p .ne. q) go to 490
         rv6(p) = 1.D0
         go to 870
  490    norm = abs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = norm + abs(d(i)) + abs(e(i))
c     ********** eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow **********
         eps2 = 1.0D-3 * norm
         eps3 = machep * norm
         uk = DBLE(q-p+1)
         eps4 = uk * eps3
         uk = eps4 / sqrt(uk)
         s = p
  505    group = 0
         go to 520
c     ********** look for close or coincident roots **********
  510    if (abs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.D0) x1 = x0 + order * eps3
c     ********** elimination with interchanges and
c                initialization of vector **********
  520    v = 0.D0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (abs(e(i)) .lt. abs(u)) go to 540
c     ********** warning -- a divide check may occur here if
c                e2 array has not been specified correctly **********
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.D0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.D0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.D0) u = eps3
         rv1(q) = u
         rv2(q) = 0.D0
         rv3(q) = 0.D0
c     ********** back substitution
c                for i=q step -1 until p do -- **********
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     ********** orthogonalize with respect to previous
c                members of group **********
         if (group .eq. 0) go to 700
         j = r
c
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.D0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.D0
c
         do 720 i = p, q
  720    norm = norm + abs(rv6(i))
c
         if (norm .ge. 1.D0) go to 840
c     ********** forward substitution **********
         if (its .eq. 5) go to 830
         if (norm .ne. 0.D0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     ********** elimination operations on next vector
c                iterate **********
  780    do 820 i = ip, q
            u = rv6(i)
c     ********** if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process **********
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     ********** set error -- non-converged eigenvector **********
  830    ierr = -r
         xu = 0.D0
         go to 870
c     ********** normalize so that sum of squares is
c                1 and expand to full order **********
  840    u = 0.D0
c
         do 860 i = p, q
  860    u = u + rv6(i)**2
c
         xu = 1.D0 / sqrt(u)
c
  870    do 880 i = 1, n
  880    z(i,r) = 0.D0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (q .lt. n) go to 100
 1001 return
c     ********** last card of tinvit **********
      end
CMK begin include functionals.f from CPMD
C     ==================================================================
      SUBROUTINE XC(RHO,EX,EC,VX,VC)
C     ==--------------------------------------------------------------==
C     ==  LDA EXCHANGE AND CORRELATION FUNCTIONALS                    ==
C     ==                                                              ==
C     ==  EXCHANGE  :  SLATER alpha                                   ==
C     ==  CORRELATION : CEPERLEY & ALDER (PERDEW-ZUNGER PARAMETERS)   ==
C     ==                VOSKO, WILK & NUSSAIR                         ==
C     ==                LEE, YANG & PARR                              ==
C     ==                PERDEW & WANG                                 ==
C     ==                WIGNER                                        ==
C     ==                HEDIN & LUNDQVIST                             ==
C     ==                ORTIZ & BALLONE (PERDEW-ZUNGER FORMULA)       ==
C     ==                ORTIZ & BALLONE (PERDEW-WANG FORMULA)         ==
C     ==                HCTH/120                                      ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'func.inc'
      PARAMETER (SMALL=1.D-10)
      PARAMETER (PI34= 0.75D0 / 3.141592653589793D+00
     +         ,THIRD=1.D0/3.D0)
C     ==--------------------------------------------------------------==
C..Exchange
      IF(MFXCX.EQ.1) THEN
        CALL SLATERX(RHO,EX,VX,SALPHA)
      ELSE
        EX=0.0D0
        VX=0.0D0
      ENDIF
      IF(RHO.LE.SMALL) THEN
        EC = 0.0D0
        VC = 0.0D0
        EX = 0.0D0
        VX = 0.0D0
      ELSE IF(MFXCC.EQ.1) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL PZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.2) THEN
        RS = (PI34/RHO)**THIRD
        CALL VWN(RS,EC,VC)
      ELSEIF(MFXCC.EQ.3) THEN
        CALL LYP(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.4) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL PW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.5) THEN
        CALL WIGNER(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.6) THEN
        CALL HEDIN(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.7) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL OBPZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.8) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL OBPW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.9) THEN
        RS=(PI34/RHO)**THIRD
        CALL PADE(RS,EC,VC)
      ELSE
        EC=0.0D0
        VC=0.0D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GCXC(RHO,GRHO,SX,SC,V1X,V2X,V1C,V2C)
      use xc_b97, only: eval_b97
C     ==--------------------------------------------------------------==
C     ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
C     ==                                                              ==
C     ==  EXCHANGE  :  BECKE88                                        ==
C     ==               GGAX                                           ==
C     ==               PBEX                                           ==
C     ==               PBESX                                           ==
C     ==               revPBEX                                        ==
C     ==               HCTH/120                                       ==
C     ==               OPTX                                           ==
C     ==  CORRELATION : PERDEW86                                      ==
C     ==                LEE, YANG & PARR                              ==
C     ==                GGAC                                          ==
C     ==                PBEC                                          ==
C     ==                PBESC                                          ==
C     ==                HCTH/120                                      ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SMALL=1.D-10)
      INCLUDE 'func.inc'
C     ==--------------------------------------------------------------==
C..Exchange
      IF(RHO.LE.SMALL) THEN
        SX  = 0.0D0
        V1X = 0.0D0
        V2X = 0.0D0
      ELSEIF(MGCX.EQ.1) THEN
        CALL BECKE88(bbeta,RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.2) THEN
        CALL GGAX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.3) THEN
        CALL PBEX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.4) THEN
        CALL revPBEX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.5.AND.MGCC.EQ.5) THEN
        CALL HCTH120(RHO,GRHO,SX,V1X,V2X) ! x&c
        SC=0.0D0
        V1C=0.0D0 
        V2C=0.0D0
      ELSEIF(MGCX.EQ.6) THEN
        CALL OPTX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.7) THEN
        CALL BECKE88(BBETA,RHO,GRHO,SXA,V1XA,V2XA)
        CALL GGAX(RHO,GRHO,SXB,V1XB,V2XB)
        SX=0.722D0*SXA+0.347D0*SXB
        V1X=0.722D0*V1XA+0.347D0*V1XB
        V2X=0.722D0*V2XA+0.347D0*V2XB
      ELSEIF(MGCX.EQ.8) THEN
        CALL BECKE88(BBETA,RHO,GRHO,SXA,V1XA,V2XA)
        CALL GGAX(RHO,GRHO,SXB,V1XB,V2XB)
        SX=0.542D0*SXA+0.167D0*SXB
        V1X=0.542D0*V1XA+0.167D0*V1XB
        V2X=0.542D0*V2XA+0.167D0*V2XB
      ELSEIF(MGCX.EQ.9) THEN
        CALL PBESX(RHO,GRHO,SX,V1X,V2X)
CMK Extra functionals not yet available with CPMD
      ELSEIF(MGCX.EQ.13.AND.MGCC.EQ.13) THEN
        CALL HCTH(93,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.14.AND.MGCC.EQ.14) THEN
        CALL HCTH(120,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.15.AND.MGCC.EQ.15) THEN
        CALL HCTH(147,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.16.AND.MGCC.EQ.16) THEN
        CALL HCTH(407,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.11) THEN
        CALL eval_b97(1,RHO,GRHO,SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.12) THEN
        CALL eval_b97(2,RHO,GRHO,SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSE
        SX=0.0D0
        V1X=0.0D0
        V2X=0.0D0
      ENDIF
C..Correlation
      IF(RHO.LE.SMALL) THEN
        SC  = 0.0D0
        V1C = 0.0D0
        V2C = 0.0D0
      ELSEIF(MGCC.EQ.1) THEN
        CALL PERDEW86(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.2) THEN
        CALL GLYP(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.3) THEN
        CALL GGAC(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.4) THEN
        W1=1.D0
        CALL PBEC(RHO,GRHO,W1,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.6) THEN
        W1=0.74D0
        CALL PBEC(RHO,GRHO,W1,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.7) THEN
        CALL PBESC(RHO,GRHO,W1,SC,V1C,V2C)
      ELSE
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SLATERX(RHO,EX,VX,ALPHA)
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SMALL=1.D-10)
      PARAMETER (F1 = -1.10783814957303361D0)
      PARAMETER (THIRD=1.D0/3.D0,F43=4.D0/3.D0)
C     ==--------------------------------------------------------------==
      IF(RHO.LE.SMALL) THEN
        EX = 0.0D0
        VX = 0.0D0
      ELSE
        RS = RHO**THIRD
        EX = F1*ALPHA*RS
        VX = F43*F1*ALPHA*RS
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PZ(RS,EPZ,VPZ,IFLG)
C     ==--------------------------------------------------------------==
C     ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.0311D0,B=-0.048D0,C=0.0020D0,D=-0.0116D0,
     *           GC=-0.1423D0,B1=1.0529D0,B2=0.3334D0)
C     ==--------------------------------------------------------------==
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+
     *              (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7.D0/6.D0*B1*RS1+4.D0/3.D0*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE VWN(RS,EVWN,VVWN)
C     ==--------------------------------------------------------------==
C     ==  S.H VOSKO, L.WILK, AND M. NUSAIR,                           ==
C     ==                 CAN. J. PHYS. 58 1200  (1980)                ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.0310907D0,B=3.72744D0,C=12.9352D0,X0=-0.10498D0)
      PARAMETER (TWO=2.0D0)
C     ==--------------------------------------------------------------==
      Q  = SQRT(4.D0*C-B*B)
      F1 = TWO*B/Q
      F2 = B*X0/(X0*X0+B*X0+C)
      F3 = TWO*(TWO*X0+B)/Q
      X  = SQRT(RS)
      FX = X*X+B*X+C
      QX = ATAN(Q/(TWO*X+B))
      EVWN=A*(LOG(RS/FX)+F1*QX-F2*(LOG((X-X0)**2/FX)+F3*QX))
      TXPB=TWO*X+B
      TTQQ=TXPB*TXPB+Q*Q
      VVWN=EVWN - X*A/6.D0*(TWO/X-TXPB/FX-4.D0*B/TTQQ-F2*(TWO/(X-X0)
     *          -TXPB/FX-4.D0*(TWO*X0+B)/TTQQ))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE LYP(RHO,ELYP,VLYP)
C     ==--------------------------------------------------------------==
C     ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
C     ==  THIS IS ONLY THE LDA PART                                   ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0)
      PARAMETER (CF=2.87123400018819108D0)
C     ==--------------------------------------------------------------==
      RS=RHO**(-1.D0/3.D0)
      ECRS=B*CF*EXP(-C*RS)
      OX=1.D0/(1.D0+D*RS)
      ELYP=-A*OX*(1.D0+ECRS)
      VLYP=ELYP-RS/3.D0*A*OX*(D*OX+ECRS*(D*OX+C))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PW(RS,EPWC,VPWC,IFLG)
C     ==--------------------------------------------------------------==
C     ==  J.P. PERDEW AND YUE WANG PRB 45, 13244 (1992)               ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,A1=0.21370D0,B1=7.5957D0,B2=3.5876D0,
     *           B3=1.6382D0,B4=0.49294D0,C0=A,C1=0.046644D0,
     *           C2=0.00664D0,C3=0.01043D0,D0=0.4335D0,D1=1.4408D0)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-
     *              (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5D0*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1.D0+A1*RS)*OLOG
        VPWC=-2.D0*A*(1.D0+2.D0/3.D0*A1*RS)*OLOG
     *       -2.D0/3.D0*A*(1.D0+A1*RS)*DOM/(OM*(OM+1.D0))
      ELSEIF(IFLG.EQ.3) THEN
C..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WIGNER(RHO,EXC,FXC)
      IMPLICIT REAL*8 (A-H,O-Z)
      RH=RHO
      X=RH**0.33333333333333333D0
      FXC=-X*((0.943656D0+8.8963D0*X)/(1.0D0+12.57D0*X)**2)
      EXC=-0.738D0*X*(0.959D0/(1.0D0+12.57D0*X))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE HEDIN(RHO,ECP,FCP)
      IMPLICIT REAL*8 (A-H,O-Z)
cmb-ike   VARIABLES with SAVE attribute hinder the vectorization
cmb-ike      SAVE RH
cmb-ike      IF(RH .EQ. 0.0D0) RETURN
      RH=RHO
      RSM1=0.62035049D0*RH**(0.3333333333333333D0)
      ALN=DLOG(1.0D0 + 21.0D0*RSM1)
      X=21.0D0/RSM1
      ECP = ALN+(X**3*ALN-X*X)+X/2.0D0-1.0D0/3.0D0
      ECP = -0.0225D0*ECP
      FCP = -0.0225D0*ALN
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE OBPZ(RS,EPZ,VPZ,IFLG)
C     ==--------------------------------------------------------------==
C     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
C     ==  PERDEW-ZUNGER FORMULA                                       ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,B=-0.046644D0,C=0.00419D0,D=-0.00983D0,
     *           GC=-0.103756D0,B1=0.56371D0,B2=0.27358D0)
C     ==--------------------------------------------------------------==
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+
     *              (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7.D0/6.D0*B1*RS1+4.D0/3.D0*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE OBPW(RS,EPWC,VPWC,IFLG)
C     ==--------------------------------------------------------------==
C     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
C     ==  PERDEW-WANG FORMULA                                         ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,A1=0.026481D0,B1=7.5957D0,B2=3.5876D0,
     *           B3=-0.46647D0,B4=0.13354D0,C0=A,C1=0.046644D0,
     *           C2=0.00664D0,C3=0.01043D0,D0=0.4335D0,D1=1.4408D0)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-
     *              (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5D0*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1.0D0+A1*RS)*OLOG
        VPWC=-2.D0*A*(1.D0+2.D0/3.D0*A1*RS)*OLOG
     *       -2.D0/3.D0*A*(1.D0+A1*RS)*DOM/(OM*(OM+1.D0))
      ELSEIF(IFLG.EQ.3) THEN
C..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PADE(RS,EC,VC)
C     ==--------------------------------------------------------------==
C     ==  PADE APPROXIMATION                                          ==
C     ==  S. GOEDECKER, M. TETER, J. HUTTER, PRB in press             ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A0=0.4581652932831429D0,A1=2.217058676663745D0,
     *           A2=0.7405551735357053D0,A3=0.01968227878617998D0)
      PARAMETER (B1=1.0000000000000000D0,B2=4.504130959426697D0,
     *           B3=1.110667363742916D0,B4=0.02359291751427506D0)
      PARAMETER (O3=1.D0/3.D0)
C     ==--------------------------------------------------------------==
      TOP=A0+RS*(A1+RS*(A2+RS*A3))
      DTOP=A1+RS*(2.D0*A2+3.D0*A3*RS)
      BOT=RS*(B1+RS*(B2+RS*(B3+RS*B4)))
      DBOT=B1+RS*(2.D0*B2+RS*(3.D0*B3+RS*4.D0*B4))
      EC=-TOP/BOT
      VC=EC+RS*O3*(DTOP/BOT-TOP*DBOT/(BOT*BOT))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE BECKE88(B1,RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C BECKE EXCHANGE: PRA 38, 3098 (1988)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(OB3=1.D0/3.D0)
C     ==--------------------------------------------------------------==
      TWO13 = 2.0D0**(1.D0/3.D0)
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      XS    = TWO13*A/BR4
      XS2   = XS*XS
      SA2B8 = SQRT(1.0D0+XS2)
      SHM1  = LOG(XS+SA2B8)
      DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
      DD2   = DD*DD
      EE    = 6.0D0*B1*XS2/SA2B8 - 1.D0
      SX    = TWO13*AA/BR4*(-B1/DD)
      V1X   = -(4.D0/3.D0)/TWO13*XS2*B1*BR1*EE/DD2
      V2X   = TWO13*B1*(EE-DD)/(BR4*DD2)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GGAX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C J.P.PERDEW ET AL. PRB 46 6671 (1992)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(F1=0.19645D0,F2=7.7956D0,F3=0.2743D0,F4=0.1508D0,
     *          F5=0.004D9,PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      FP1   = -3.D0/(16.D0*PI)*(3.D0*PI*PI)**(-1.D0/3.D0)
      FP2   = 0.5D0*(3.D0*PI*PI)**(-1.D0/3.D0)
      AA    = GRHO
      A     = SQRT(AA)
      RR    = RHO**(-4.D0/3.D0)
      S     = FP2*A*RR
      S2    = S*S
      S3    = S2*S
      S4    = S2*S2
      EXPS  = F4*EXP(-100.D0*S2)
      AS    = F3-EXPS-F5*S2
      SA2B8 = SQRT(1.0D0+F2*F2*S2)
      SHM1  = LOG(F2*S+SA2B8)
      BS    = 1.D0+F1*S*SHM1+F5*S4
      DAS   = 200.D0*S*EXPS-2.D0*S*F5
      DBS   = F1*(SHM1+F2*S/SA2B8)+4.D0*F5*S3
      DLS   = (DAS/AS-DBS/BS)
      SX    = FP1*AA*RR*AS/BS
      V1X   = -4.D0/3.D0*SX/RHO*(1.D0+S*DLS)
      V2X   = FP1*RR*AS/BS*(2.D0+S*DLS)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBESX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C PBESol functional
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0,
     *          UM=0.123456790123456789D0,UK=0.8040D0,UL=UM/UK)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4.D0/3.D0)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBEX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C J.P.PERDEW ET AL. PRL 77 3865 (1996)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0,
     *          UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4.D0/3.D0)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE revPBEX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C Y. ZHANG ET AL. PRL 80 890 (1998)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0,
     *          UM=0.2195149727645171D0,UK=1.2450D0,UL=UM/UK)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4.D0/3.D0)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PERDEW86(RHO,GRHO,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PERDEW CORRELATION: PRB 33, 8822 (1986)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(P1=0.023266D0,P2=7.389D-6,P3=8.723D0,P4=0.472D0)
      PARAMETER(PC1=0.001667D0,PC2=0.002568D0,PCI=PC1+PC2)
      PARAMETER(OB3=1.D0/3.D0, FPI=4.0D0*3.141592653589793D0)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      RS    = (3.D0/(FPI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      CNA   = PC2+P1*RS+P2*RS2
      CNB   = 1.D0+P3*RS+P4*RS2+1.D4*P2*RS3
      CN    = PC1 + CNA/CNB
      DRS   = -OB3*(3.D0/FPI)**OB3 / BR4
      DCNA  = (P1+2.D0*P2*RS)*DRS
      DCNB  = (P3+2.D0*P4*RS+3.D4*P2*RS2)*DRS
      DCN   = DCNA/CNB - CNA/(CNB*CNB)*DCNB
      PHI   = 0.192D0*PCI/CN*A*RHO**(-7.D0/6.D0)
      EPHI  = EXP(-PHI)
      SC    = AA/BR4*CN*EPHI
      V1C   = SC*((1.D0+PHI)*DCN/CN -((4.D0/3.D0)-(7.D0/6.D0)*PHI)/RHO)
      V2C   = CN*EPHI/BR4*(2.D0-PHI)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GLYP(RHO,GRHO,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C LEE, YANG PARR: GRADIENT CORRECTION PART
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      R     = RHO**(-1.d0/3.d0)
      OM    = EXP(-C*R)/(1.d0+D*R)
      R5    = R**5
      XL    = 1.d0+(7.d0/3.d0)*(C*R + D*R/(1.d0+D*R))
      FF    = A*B*AA/24.d0
      SC    = FF*R5*OM*XL
      DR5   = 5.d0*R*R*R*R
      DOM   = -OM*(C+D+C*D*R)/(1.d0+D*R)
      DXL   = (7.d0/3.d0)*(C+D+2.d0*C*D*R+C*D*D*R*R)/(1.d0+D*R)**2
      V1C   = -FF*(R*R*R*R)/3.d0*( DR5*OM*XL + R5*DOM*XL + R5*OM*DXL)
      V2C   = A*B*R5*OM*XL/12.d0
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GGAC(RHO,GRHO,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PERDEW & WANG GGA CORRELATION PART      
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(AL=0.09D0,PA=0.023266D0,PB=7.389D-6,PC=8.723D0,
     *          PD=0.472D0,CX=-0.001667D0,CXC0=0.002568D0,CC0=-CX+CXC0)
      PARAMETER(OB3=1.D0/3.D0, PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      XNU   = 16.D0/PI*(3.D0*PI*PI)**OB3
      BE    = XNU*CC0
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3.D0/(4.D0*PI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      XKF   = (9.D0*PI/4.D0)**OB3/RS
      XKS   = SQRT(4.D0*XKF/PI)
      T     = A/(2.D0*XKS*RHO)
      EXPE  = EXP(-2.D0*AL*EC/(BE*BE))
      AF    = 2.D0*AL/BE * (1.D0/(EXPE-1.D0))
      BF    = EXPE*(VC-EC)
      Y     = AF*T*T
      XY    = (1.D0+Y)/(1.D0+Y+Y*Y)
      QY    = Y*Y*(2.D0+Y)/(1.D0+Y+Y*Y)**2
      S1    = 1.D0+2.D0*AL/BE*T*T*XY
      H0    = BE*BE/(2.D0*AL) * LOG(S1)
      DH0   = BE*T*T/S1*(-7.D0/3.D0*XY-QY*(AF*BF/BE-7.D0/3.D0))
      DDH0  = BE/(2.D0*XKS*XKS*RHO)*(XY-QY)/S1
      EE    = -100.D0*(XKS/XKF*T)**2
      CNA   = CXC0+PA*RS+PB*RS2
      DCNA  = -(PA*RS+2.D0*PB*RS2)/3.D0
      CNB   = 1.D0+PC*RS+PD*RS2+1.D4*PB*RS3
      DCNB  = -(PC*RS+2.D0*PD*RS2+3.D4*PB*RS3)/3.D0
      CN    = CNA/CNB - CX
      DCN   = DCNA/CNB - CNA*DCNB/(CNB*CNB)
      H1    = XNU*(CN-CC0-3.D0/7.D0*CX)*T*T*EXP(EE)
      DH1   = -OB3*(H1*(7.D0+8.D0*EE)+XNU*T*T*EXP(EE)*DCN)
      DDH1  = 2.D0*H1*(1.D0+EE)*RHO/AA
      SC    = RHO*(H0+H1)
      V1C   = H0+H1+DH0+DH1
      V2C   = DDH0+DDH1
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBESC(RHO,GRHO,W1,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PBESol Correlation functional
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(BE=0.046D0,GA=0.031090690869654895D0)
      PARAMETER(OB3=1.D0/3.D0,PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3.D0/(4.D0*PI*RHO))**OB3
      XKF   = (9.D0*PI/4.D0)**OB3/RS
      XKS   = SQRT(4.D0*XKF/PI)
      T     = A/(2.D0*XKS*RHO)
      EXPE  = EXP(-EC/GA)
      AF    = BE/GA * (1.D0/(EXPE-1.D0))
      Y     = AF*T*T
      XY    = (1.D0+Y)/(1.D0+Y+Y*Y)
      S1    = 1.D0+BE/GA*T*T*XY
      H0    = GA * LOG(S1)
      DTDR  = -T*7.D0/(6.D0*RHO)
      DADR  = AF*AF*EXPE/BE*(VC-EC)/RHO
      DSDA  = -BE/GA * AF * T**6 * (2.D0+Y) / (1.D0+Y+Y*Y)**2
      DSDT  = 2.D0*BE/GA * T * (1.D0+2.D0*Y) / (1.D0+Y+Y*Y)**2
      DSDR  = DSDA*DADR + DSDT*DTDR
      DHDT  = GA/S1*DSDT
      DHDR  = GA/S1*DSDR
      SC    = W1*RHO*H0
      V1C   = W1*H0+W1*DHDR*RHO
      V2C   = W1*RHO*DHDT*T/AA
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBEC(RHO,GRHO,W1,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PBE Correlation functional
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(BE=0.06672455060314922D0,GA=0.031090690869654895D0)
      PARAMETER(OB3=1.D0/3.D0,PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3.D0/(4.D0*PI*RHO))**OB3
      XKF   = (9.D0*PI/4.D0)**OB3/RS
      XKS   = SQRT(4.D0*XKF/PI)
      T     = A/(2.D0*XKS*RHO)
      EXPE  = EXP(-EC/GA)
      AF    = BE/GA * (1.D0/(EXPE-1.D0))
      Y     = AF*T*T
      XY    = (1.D0+Y)/(1.D0+Y+Y*Y)
      S1    = 1.D0+BE/GA*T*T*XY
      H0    = GA * LOG(S1)
      DTDR  = -T*7.D0/(6.D0*RHO)
      DADR  = AF*AF*EXPE/BE*(VC-EC)/RHO
      DSDA  = -BE/GA * AF * T**6 * (2.D0+Y) / (1.D0+Y+Y*Y)**2
      DSDT  = 2.D0*BE/GA * T * (1.D0+2.D0*Y) / (1.D0+Y+Y*Y)**2
      DSDR  = DSDA*DADR + DSDT*DTDR
      DHDT  = GA/S1*DSDT
      DHDR  = GA/S1*DSDR
      SC    = W1*RHO*H0
      V1C   = W1*H0+W1*DHDR*RHO
      V2C   = W1*RHO*DHDT*T/AA
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE hcth120(rho,grho,sx,v1x,v2x)
C     HCTH, JCP 109, 6264 (1998)
C     Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
C     Present release: Tsukuba, 09/02/2005
c--------------------------------------------------------------------------
c     rhoa = rhob = 0.5 * rho
c     grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
c     sx  : total exchange correlation energy at point r 
c     v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
c     v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
c--------------------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER(o3=1.0d0/3.0d0,fr83=8.d0/3.d0)
      DIMENSION cg0(6),cg1(6),caa(6),cab(6),cx(6)
      r3q2=DEXP(-o3*0.69314718055994531d0)
      r3pi=DEXP(-o3*0.04611759718129048d0)
c.....coefficients for PW correlation......................................
      cg0(1)= 0.031091d0
      cg0(2)= 0.213700d0
      cg0(3)= 7.595700d0
      cg0(4)= 3.587600d0
      cg0(5)= 1.638200d0
      cg0(6)= 0.492940d0
      cg1(1)= 0.015545d0
      cg1(2)= 0.205480d0
      cg1(3)=14.118900d0
      cg1(4)= 6.197700d0
      cg1(5)= 3.366200d0
      cg1(6)= 0.625170d0
C......HCTH-19-4.....................................
      caa(1)=  0.489508D+00
      caa(2)= -0.260699D+00
      caa(3)=  0.432917D+00
      caa(4)= -0.199247D+01
      caa(5)=  0.248531D+01
      caa(6)=  0.200000D+00
      cab(1)=  0.514730D+00
      cab(2)=  0.692982D+01
      cab(3)= -0.247073D+02
      cab(4)=  0.231098D+02
      cab(5)= -0.113234D+02
      cab(6)=  0.006000D+00
      cx(1) =  0.109163D+01
      cx(2) = -0.747215D+00
      cx(3) =  0.507833D+01
      cx(4) = -0.410746D+01
      cx(5) =  0.117173D+01
      cx(6) =  0.004000D+00
c...........................................................................
      gr=DSQRT(grho)
      rho_o3=rho**(o3) 
      rho_o34=rho*rho_o3
      xa=1.25992105d0*gr/rho_o34
      xa2=xa*xa
      ra=0.781592642d0/rho_o3
      rab=r3q2*ra
      dra_drho=-0.260530881d0/rho_o34
      drab_drho=r3q2*dra_drho
      CALL pwcorr(ra,cg1,g,dg)
      era1=g
      dera1_dra=dg
      CALL pwcorr(rab,cg0,g,dg)
      erab0=g
      derab0_drab=dg
      ex=-0.75d0*r3pi*rho_o34
      dex_drho=-r3pi*rho_o3
      uaa=caa(6)*xa2
      uaa=uaa/(1.0d0+uaa)
      uab=cab(6)*xa2
      uab=uab/(1.0d0+uab)
      ux=cx(6)*xa2
      ux=ux/(1.0d0+ux)
      ffaa=rho*era1
      ffab=rho*erab0-ffaa
      dffaa_drho=era1+rho*dera1_dra*dra_drho
      dffab_drho=erab0+rho*derab0_drab*drab_drho-dffaa_drho
cmb-> i-loop removed
      denaa=1.d0/(1.0d0+caa(6)*xa2)
      denab=1.d0/(1.0d0+cab(6)*xa2)
      denx =1.d0/(1.0d0+cx(6)*xa2)
      f83rho=fr83/rho
      bygr=2.0d0/gr
      gaa=caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))
      gab=cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))
      gx=cx(1)+ux*(cx(2)+ux*(cx(3)+ux*(cx(4)+ux*cx(5))))
      taa=denaa*uaa*(caa(2)+uaa*(2.d0*caa(3)+uaa 
     &    *(3.d0*caa(4)+uaa*4.d0*caa(5))))
      tab=denab*uab*(cab(2)+uab*(2.d0*cab(3)+uab
     &    *(3.d0*cab(4)+uab*4.d0*cab(5))))
      txx=denx*ux*(cx(2)+ux*(2.d0*cx(3)+ux
     &    *(3.d0*cx(4)+ux*4.d0*cx(5))))
      dgaa_drho=-f83rho*taa
      dgab_drho=-f83rho*tab
      dgx_drho=-f83rho*txx
      dgaa_dgr=bygr*taa
      dgab_dgr=bygr*tab
      dgx_dgr=bygr*txx
cmb
      sx=ex*gx+ffaa*gaa+ffab*gab
      v1x=dex_drho*gx+ex*dgx_drho
     .   +dffaa_drho*gaa+ffaa*dgaa_drho
     .   +dffab_drho*gab+ffab*dgab_drho
      v2x=(ex*dgx_dgr+ffaa*dgaa_dgr+ffab*dgab_dgr)/gr
      RETURN
      END
C =-------------------------------------------------------------------=
      SUBROUTINE pwcorr(r,c,g,dg)
      IMPLICIT real*8 (a-h,o-z)
      DIMENSION c(6)
      r12=DSQRT(r)
      r32=r*r12
      r2=r*r
      rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
      sb=1.0d0+1.0d0/(2.0d0*c(1)*rb)
      g=-2.0d0*c(1)*(1.0d0+c(2)*r)*DLOG(sb)
      drb=c(3)/(2.0d0*r12)+c(4)+1.5d0*c(5)*r12+2.0d0*c(6)*r
      dg=(1.0d0+c(2)*r)*drb/(rb*rb*sb)-2.0d0*c(1)*c(2)*DLOG(sb)
      RETURN
      END 
C     ==================================================================
      SUBROUTINE OPTX(rho,grho,sx,v1x,v2x)
C     OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
C     Present release: Tsukuba, 20/6/2002
c--------------------------------------------------------------------------
c     rhoa = rhob = 0.5 * rho in LDA implementation
c     grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
c     sx  : total exchange correlation energy at point r
c     v1x : d(sx)/drho
c     v2x : 1/gr*d(sx)/d(gr)
c--------------------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER(SMALL=1.D-20,SMAL2=1.D-08)
C.......coefficients and exponents....................
      PARAMETER(o43=4.0d0/3.0d0,two13=1.259921049894873D0
     .         ,two53=3.174802103936399D0,gam=0.006D0
     .         ,a1cx=0.9784571170284421D0,a2=1.43169D0)
C.......OPTX in compact form..........................
      IF(RHO.LE.SMALL) THEN
       sx=0.0D0
       v1x=0.0D0
       v2x=0.0D0
      ELSE
       gr=DMAX1(grho,SMAL2)
       rho43=rho**o43
       xa=two13*DSQRT(gr)/rho43
       gamx2=gam*xa*xa
       uden=1.d+00/(1.d+00+gamx2)
       uu=a2*gamx2*gamx2*uden*uden
       uden=rho43*uu*uden
       sx=-rho43*(a1cx+uu)/two13
       v1x=o43*(sx+two53*uden)/rho
       v2x=-two53*uden/gr
      ENDIF
C
      RETURN
      END
C =-------------------------------------------------------------------=
CMK end include functionals.f from CPMD
        subroutine ggaenergy_15(nrad,rw,rd,rho,enexc,pot,eps)
        implicit real*8 (a-h,o-z)
c calculate exc energy enexc
        dimension rho(nrad),rw(nrad),rd(nrad),pot(nrad),eps(nrad)
        dimension c(-8:8)

        enexc=0.d0
        call zero(nrad,pot)
        call zero(nrad,eps)

        j=1
         c(0)=-2.717857142857143d0
         c(1)=8.d0
         c(2)=-14.d0
         c(3)=18.66666666666667d0
         c(4)=-17.5d0
         c(5)=11.2d0
         c(6)=-4.666666666666667d0
         c(7)=1.142857142857143d0
         c(8)=-0.125d0
        rder=0.d0
        do i=-0,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-0,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=2
         c(-1)=-0.1111111111111111d0
         c(0)=-1.717857142857143d0
         c(1)=4.d0
         c(2)=-4.666666666666667d0
         c(3)=4.666666666666667d0
         c(4)=-3.5d0
         c(5)=1.866666666666667d0
         c(6)=-0.6666666666666666d0
         c(7)=0.1428571428571428d0
         c(8)=-0.01388888888888889d0
        rder=0.d0
        do i=-1,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-1,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=3
         c(-2)=0.01111111111111111d0
         c(-1)=-0.2222222222222222d0
         c(0)=-1.217857142857143d0
         c(1)=2.666666666666666d0
         c(2)=-2.333333333333333d0
         c(3)=1.866666666666667d0
         c(4)=-1.166666666666667d0
         c(5)=0.5333333333333333d0
         c(6)=-0.1666666666666666d0
         c(7)=0.03174603174603174d0
         c(8)=-0.2777777777777778d-2
        rder=0.d0
        do i=-2,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-2,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=4
         c(-3)=-0.202020202020202d-2
         c(-2)=0.03333333333333333d0
         c(-1)=-0.3333333333333333d0
         c(0)=-0.88452380952381d0
         c(1)=2.d0
         c(2)=-1.4d0
         c(3)=0.933333333333333d0
         c(4)=-0.5d0
         c(5)=0.2d0
         c(6)=-0.05555555555555556d0
         c(7)=0.952380952380952d-2
         c(8)=-0.7575757575757577d-3
        rder=0.d0
        do i=-3,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-3,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=5
         c(-4)=0.5050505050505051d-3
         c(-3)=-0.808080808080808d-2
         c(-2)=0.06666666666666666d0
         c(-1)=-0.4444444444444445d0
         c(0)=-0.6345238095238095d0
         c(1)=1.6d0
         c(2)=-0.933333333333333d0
         c(3)=0.5333333333333333d0
         c(4)=-0.25d0
         c(5)=0.0888888888888889d0
         c(6)=-0.02222222222222222d0
         c(7)=0.3463203463203463d-2
         c(8)=-0.2525252525252525d-3
        rder=0.d0
        do i=-4,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-4,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=6
         c(-5)=-0.1554001554001554d-3
         c(-4)=0.2525252525252525d-2
         c(-3)=-0.0202020202020202d0
         c(-2)=0.1111111111111111d0
         c(-1)=-0.5555555555555556d0
         c(0)=-0.4345238095238095d0
         c(1)=1.333333333333333d0
         c(2)=-0.6666666666666666d0
         c(3)=0.3333333333333333d0
         c(4)=-0.1388888888888889d0
         c(5)=0.04444444444444445d0
         c(6)=-0.0101010101010101d0
         c(7)=0.1443001443001443d-2
         c(8)=-0.971250971250971d-4
        rder=0.d0
        do i=-5,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-5,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=7
         c(-6)=0.555000555000555d-4
         c(-5)=-0.932400932400932d-3
         c(-4)=0.7575757575757577d-2
         c(-3)=-0.04040404040404041d0
         c(-2)=0.1666666666666666d0
         c(-1)=-0.6666666666666666d0
         c(0)=-0.2678571428571428d0
         c(1)=1.142857142857143d0
         c(2)=-0.5d0
         c(3)=0.2222222222222222d0
         c(4)=-0.0833333333333333d0
         c(5)=0.02424242424242424d0
         c(6)=-0.5050505050505051d-2
         c(7)=0.6660006660006659d-3
         c(8)=-0.4162504162504162d-4
        rder=0.d0
        do i=-6,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-6,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=8
         c(-7)=-0.222000222000222d-4
         c(-6)=0.3885003885003884d-3
         c(-5)=-0.3263403263403263d-2
         c(-4)=0.01767676767676768d0
         c(-3)=-0.07070707070707071d0
         c(-2)=0.2333333333333333d0
         c(-1)=-0.7777777777777778d0
         c(0)=-0.125d0
         c(1)=1.d0
         c(2)=-0.3888888888888889d0
         c(3)=0.1555555555555556d0
         c(4)=-0.05303030303030303d0
         c(5)=0.01414141414141414d0
         c(6)=-0.2719502719502719d-2
         c(7)=0.3330003330003329d-3
         c(8)=-0.1942501942501942d-4
        rder=0.d0
        do i=-7,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-7,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo


         c(-8)=9.71250971250971d-6
         c(-7)=-0.1776001776001776d-3
         c(-6)=0.1554001554001554d-2
         c(-5)=-0.87024087024087d-2
         c(-4)=0.3535353535353535d-1
         c(-3)=-0.1131313131313131d0
         c(-2)=0.3111111111111111d0
         c(-1)=-0.888888888888889d0
         c(0)=0.d0
         c(1)=0.888888888888889d0
         c(2)=-0.3111111111111111d0
         c(3)=0.1131313131313131d0
         c(4)=-0.3535353535353535d-1
         c(5)=0.87024087024087d-2
         c(6)=-0.1554001554001554d-2
         c(7)=0.1776001776001776d-3
         c(8)=-9.71250971250971d-6
        do 100,j=9,nrad-8
        rder=0.d0
        do i=-8,8
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,8
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo
100     continue

        j=nrad-7
        rder=0.d0
        do i=-8,7
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,7
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-6
        rder=0.d0
        do i=-8,6
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,6
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-5
        rder=0.d0
        do i=-8,5
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,5
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-4
        rder=0.d0
        do i=-8,4
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,4
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-3
        rder=0.d0
        do i=-8,3
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,3
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-2
        rder=0.d0
        do i=-8,2
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,2
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-1
        rder=0.d0
        do i=-8,1
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,1
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        j=nrad-0
        rder=0.d0
        do i=-8,0
        rder=rder+c(i)*rho(j+i)
        enddo
        if (rder.ge.0.d0) then
        sign=rd(j)
        else
        sign=-rd(j)
        endif
        rder=sign*rder
        call XCFUNCTION(epsxc,rho(j),der1,der2,rder)
        enexc=enexc+epsxc*rw(j)
        eps(j)=eps(j)+epsxc
        pot(j)=pot(j)+der1*rw(j)
        do i=-8,0
        pot(j+i)=pot(j+i)+(sign*c(i)*der2)*rw(j)
        enddo

        do j=2,nrad
           pot(j)=pot(j)/rw(j)
        enddo
        pot(1)=pot(2)

        return
        end

      SUBROUTINE XCFUNCTION(EXC,RHOE,V,VTMP,GRAD)
c     input:  rho=RHO, GRAD=grad(rho)
c     output: exc=energy density: energy = \int EXC dr
c             v=dEXC/dRHO, vtmp=dEXC/dGRAD
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
C     ==--------------------------------------------------------------==

c     using hutter's routine
      exc=0.0d0
      v=0.0d0
      rhoe=rhoe
      call xc(rhoe,ex,ec,vx,vc)
      exc=(ex+ec)*rhoe
      v=vx+vc
      Vtmp=0.d0
c     c.hartwig: don't calc. gradients for rho>1.0d-18
c     because of numerical noise in the density!!!!
CMK      if (rhoe.gt.1.0d-18) then
      IF (rhoe.gt.1.0d-10) then
         dr2=grad*grad
         call gcxc(rhoe,dr2,sx,sc,v1x,v2x,v1c,v2c)
         exc=exc+(sx+sc)
         v=v+v1x+v1c
         VTMP= (V2X+V2C)*abs(grad)
      endif
      return
      end
      subroutine zero(n,x)
      implicit real*8 (a-h,o-z)
      dimension x(n)
      do 10,i=1,n-1,2
         x(i)=0.d0
         x(i+1)=0.d0
 10   continue
      istart=i
      do 11,i=istart,n
         x(i)=0.d0
 11   continue
      return
      end
