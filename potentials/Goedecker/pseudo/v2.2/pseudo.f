C finds approximate pseudopotential parameters by simplex downhill method

      implicit real*8 (a-h,o-z)
      logical fullac, avgl1,avgl2,avgl3,plotwf,denbas,ignore,
     :     ortprj, litprj, energ,igrad,info

      parameter ( norbmx=40, nrmax=10000, maxdim=30 )
      parameter ( lmx=5, lpmx= 4, noccmx=7, nsmx=2 )
      parameter ( ngmx=27, nintmx=5*(ngmx+14) )

      dimension aeval(noccmx,lmx,nsmx),chrg(noccmx,lmx,nsmx),
     1     dhrg(noccmx,lmx,nsmx),ehrg(noccmx,lmx,nsmx),
     1     res(noccmx,lmx,nsmx),
     1     wfnode(noccmx,lmx,nsmx,3),
     1     psi(0:ngmx,noccmx,lmx,nsmx),wght(noccmx,lmx,nsmx,8),
     1     gpot(4),hsep(6,lpmx,nsmx),r_l(lpmx),
     1     occup(noccmx,lmx,nsmx),
     1     vh((lmx+1)*((ngmx+1)*(ngmx+2))/2,
     1     (lmx+1)*((ngmx+1)*(ngmx+2))/2),xp(0:ngmx),
     1     rmt(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1),
     1     rmtg(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1),
     1     ud(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1)

      dimension pp(maxdim*(maxdim+1)),yp(maxdim+1),rae(nrmax)
      dimension psiold(nrmax,noccmx,lmx+1,nsmx)
      dimension rr(nintmx),rw(nintmx),rd(nintmx)

      dimension no(norbmx),noae(norbmx),lo(norbmx),so(norbmx),
     :     zo(norbmx),ev(norbmx),crcov(norbmx),dcrcov(norbmx),
     :     ddcrcov(norbmx),
     :     gf(nrmax,norbmx,nsmx),
     :     nomax(0:lmx),nomin(0:lmx),hso(6),havg(6),hhsep(6)

      character*1 il(5),ispp,isppp
      character*10 icorr,icorrp
      character*7 is(2)
      character fname*30,tname*10,string*25,strcyc*10,form*25
      character*80 label
      integer namoeb,nsmplx


      common /ecalc/ energ
      INCLUDE 'func.inc'
c     initialize som counters
      ntime=0
      itertot=0
c
      write(6,*) '*********************************************'
      write(6,*) '***           pseudo_2.2                  ***'
      write(6,*) '***           fitting of                  ***'
      write(6,*) '***   goedecker type pseudopotentials     ***'
      write(6,*) '***   last changes: 14.7.98               ***'
      write(6,*) '*********************************************'
      write(6,*)
c
c     default values:
c
      NAMOEB=0
      nsmplx=2000
      ng=0
      rij=0
      FULLAC=.false.
      avgl1=.false.
      avgl2=.false.
      avgl3=.false.
      plotwf=.false.
      ortprj=.false.
      litprj=.false.
      energ=.false.
      igrad=.false.
      info=.false.
c
c     read command line options
c
      ICARG=IARGC()
      IF(ICARG.LT.1) THEN
        WRITE(6,*) 'No value for NAMOEB specified.'
      else
         do iicarg=1,icarg
c     IBM/DEC (1 line)
            CALL GETARG(iicarg,STRING)
c     CRAYT3D (2 lines)
c         CALL PXFGETARG(iicarg,STRING,L,IERROR)
c         IF(IERROR.NE.0) STRING=' '
            ii=index(string,'-c')
            if (ii.ne.0) then
               ndigit=index(string,' ') - 3
               jmax=min(9999,10**ndigit-1)
               do j=0,jmax
                  write(strcyc,'(i'//char(ichar('0')+ndigit)//')') j
                  if (string(ii+2:ii+1+ndigit).eq.strcyc) namoeb=j
               enddo
               write(6,*) namoeb, 'fit cycles'
               goto 11
            endif
            ii=index(string,'-orth')
            if (ii.ne.0) then
               ortprj=.true.
               write(6,*) 'orthogonalize the projectors'
               goto 11
            endif
            ii=index(string,'-lith')
            if (ii.ne.0) then
               litprj=.true.
               write(6,*) 'transform the projectors as in ',
     :              'literature'
               goto 11
            endif
            ii=index(string,'-n')
            if (ii.ne.0) then
               ndigit=index(string,' ') - 3
               jmax=min(99999,10**ndigit-1)
               do j=0,jmax
                  write(strcyc,'(i'//char(ichar('0')+ndigit)//')') j
                  if (string(ii+2:ii+1+ndigit).eq.strcyc) nsmplx=j
               enddo
               write(6,*) nsmplx, 'max. simplex iterations'
               goto 11
            endif
            ii=index(string,'-g')
            if (ii.ne.0) then
               ndigit=index(string,' ') - 3
               jmax=min(99,10**ndigit-1)
               do j=0,jmax
                  write(strcyc,'(i'//char(ichar('0')+ndigit)//')') j
                  if (string(ii+2:ii+1+ndigit).eq.strcyc) ng=j
               enddo
               write(6,*)ng,'gaussians (don''t use value from psp.par)'
               goto 11
            endif
            ii=index(string,'-r')
            if (ii.ne.0) then
               ndigit=index(string,' ') - 3
               jmax=min(99,10**ndigit-1)
               do j=0,jmax
                  write(strcyc,'(i'//char(ichar('0')+ndigit)//')') j
                  if (string(ii+2:ii+1+ndigit).eq.strcyc) rij=j
               enddo
               rij=rij/10.d0
               write(6,*)rij,' rij (don''t use value from psp.par)'
               goto 11
            endif
            ii=index(string,'-fullacc')
            if (ii.ne.0) then
               fullac=.true.
               write(6,*) 'Use max. number of gaussians'
               goto 11
            endif
            ii=index(string,'-plot')
            if (ii.ne.0) then
               plotwf=.true.
               write(6,*) 'plot wfs after each iteration'
               goto 11
            endif
            ii=index(string,'-denbas')
            if (ii.ne.0) then
               denbas=.true.
               write(6,*) 'use dense gaussian basis'
               goto 11
            endif
            ii=index(string,'-ignore')
            if (ii.ne.0) then
               ignore=.true.
               write(6,*) 'override several warnings'
               goto 11
            endif
            ii=index(string,'-info')
            if (ii.ne.0) then
               info=.true.
               write(6,*) 'append gaussian coefficients of final'
               write(6,*) 'wavefunctions to psp.par'
              goto 11
            endif
            ii=index(string,'-l1so')
            if (ii.ne.0) then
               avgl1=.true.
               write(6,*) 'average nonlocal potential zero for l=1'
               write(6,*) '(only for the highest projector )'
               goto 11
            endif
            ii=index(string,'-l2so')
            if (ii.ne.0) then
               avgl2=.true.
               write(6,*) 'average nonlocal potential zero for l=2'
               write(6,*) '(only for the highest projector)'
               goto 11
            endif
            ii=index(string,'-l3so')
            if (ii.ne.0) then
               avgl3=.true.
               write(6,*) 'average nonlocal potential zero for l=3'
               write(6,*) '(only for the highest projector)'
               goto 11
            endif
            write(6,*) 'unknown commandline option: '//string
            write(6,*) 'possible options: -cN       ',
     :           'number of fitting cycles (-1< N< 10000)'
            write(6,*) '                  -nN       ',
     :           'number of iterations per simplex'
            write(6,*) '                  -gN       ',
     :           'number gaussians'
            write(6,*) '                  -rN       ',
     :           'set rij=N/10 '
            write(6,*) '                  -fullacc  ',
     :           'use max. number of gaussians'
            write(6,*) '                  -orth     ',
     :           'orthogonalisation of the projectors'
            write(6,*) '                  -lith     ',
     :           'transformation of the projectors as in lit.'
            write(6,*) '                  -denbas   ',
     :           'use dense gaussian basis    '
            write(6,*) '                  -ignore   ',
     :           'override several warnings   '
            write(6,*) '                  -info     ',
     :           'append gaussian coefficients of final'
            write(6,*) '                            ',
     :           'wavefunctions to psp.par'
            write(6,*) '                  -plot     ',
     :           'plot wfs after each iteration'
            write(6,*) '                  -lNso     ',
     :           'average nonlocal potential = 0 for l=N'
            write(6,*) '                            ',
     :           '(only for the highest projector)'
            stop
 11         continue
         enddo
      endif
      if (namoeb.eq.0) then
         write(6,*) 'Do one pseudopotential calculation.'
         write(6,*) 'No fitting.'
      endif

      if (ortprj .and. litprj ) then
         write(6,*) 'use only one option -orth or -lith!'
         stop
      endif
c
c     ----------------- read data from AE calculation -----------------
c
      write(6,*) '***        Reading data from atom.ae      *** '
      open(unit=40,file='atom.ae',form='formatted',status='unknown')
      read(40,*,err=456) norb
      read(40,*,err=456) znucp,zionp,rcovp,rprbp
      read(40,'(a)',err=456) label
      j=1
      do i=len(label),1,-1
         if (label(i:i).ne.' ') j=i
      enddo
      isppp=label(j:j)
      read(40,'(a)',err=456) label
      j1=1
      j2=2
      do i=len(label),1,-1
         if (label(i:i).ne.' ') j1=i
      enddo
      do i=len(label),j1,-1
         if (label(i:i).eq.' ') j2=i
      enddo
      j2=j2-1
      icorrp=label(j1:j2)
c      read(40,'(t2,a)',err=456) isppp
c      read(40,'(t2,a)',err=456) icorrp
      read(40,*,err=456) ngrid
      if (ngrid .gt. nrmax ) then
         write(6,*) 'array dimension problem: ngrid,nrmax:',
     :        ngrid,nrmax
         stop
      endif
      write(6,*)'pseudo states = ', norb
      write(6,*)'znuc          = ', znucp
      write(6,*)'zpseudo       = ', zionp
      write(6,*)'r_covalent    = ', rcovp
      write(6,*)'r_confining   = ', rprbp
      write(6,*)'ispp          = ', isppp
      write(6,*)'icorr         = ', icorrp
      write(6,*)'gridpoints    = ', ngrid
      il(1) = 's'
      il(2) = 'p'
      il(3) = 'd'
      il(4) = 'f'
      il(5) = 'g'
      nspin=1
      is(1) = '  so=0'
      if (isppp.eq.'r') then
         nspin=2
         is(1)= 'so=+0.5'
         is(2)= 'so=-0.5'
      endif
      write(6,*)' nl    s      occ        ',
     :     'eigenvalue     charge(rcov)    '
c      write(6,*)' nl    s      occ        ',
c     :     'eigenvalue     charge(rcov)    ',
c     :     'dcharge         ddcharge'
      do iorb=1,norb
         read(40,*,err=456) no(iorb),lo(iorb),so(iorb),zo(iorb),
     :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)
         write(6,30) no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),
     :        ev(iorb),crcov(iorb)
c      write(6,30) no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),
c     :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)
 30      format(1x,i1,a1,f6.1,f10.4,2f16.10,2f16.7)
         read(40,*,err=456) (rae(igrid),(gf(igrid,iorb,igf),
     :        igf=1,nspin),igrid=1,ngrid)
      enddo
      goto 457
 456  write(6,*) 'error during reading atom.ae'
      stop
 457  continue
      lmax=0
      lcx=0
      do iorb=1,norb
         lmax=max(lo(iorb),lmax)
         if (zo(iorb).gt.1.d-10)lcx=max(lo(iorb),lcx)
      enddo
c      print*,'lmax=',lmax
c      print*,'lcx=',lcx, '( charge > 1.0d-10)'
      if (lmax.gt.lmx+1) then
         write(6,*) 'array dimension problem:lmax,;lmx+1:',
     1        lmax,lmx+1
         stop
      endif
c     compute corresponding n-quantum numbers of the pseudostates
c     no()   will contain n quantum numbers of the pseudostates afterwards
c     noae() will contain n quantum numbers of the ae-states afterwards
c     no() starts from n=1 for each(!) l
      noccmax=0
      do l=0,lmax
         nomin(l)=100
         nomax(l)=0
      enddo
      do iorb=1,norb
         nomin(lo(iorb))=min(no(iorb),nomin(lo(iorb)))
         nomax(lo(iorb))=max(no(iorb),nomax(lo(iorb)))
      enddo
      do iorb=1,norb
         noae(iorb)=no(iorb)
         no(iorb)=no(iorb)-nomin(lo(iorb))+ 1
      enddo
      do l=0,lmax
         noccmax= max(noccmax,nomax(l)-nomin(l)+1)
      enddo
      if (noccmax.gt.noccmx) then
         write(6,*) 'array dimension problem:noccmax,noccmx:',
     1        noccmax,noccmx
      endif
c      print*,'noccmax=',noccmax
      do nocc=1,noccmax
         do l=0,lmax
            do ispin=1,nspin
               occup(nocc,l+1,ispin)=0.0d0
               aeval(nocc,l+1,ispin)=0.0d0
               chrg (nocc,l+1,ispin)=0.0d0
               dhrg (nocc,l+1,ispin)=0.0d0
               ehrg (nocc,l+1,ispin)=0.0d0
            enddo
         enddo
      enddo
      do iorb=1,norb
         nocc=no(iorb)
         l=lo(iorb)
         ispin=1
         if (so(iorb).lt.0) ispin=2
         occup(nocc,l+1,ispin)=zo(iorb)
         do j=1,ngrid
            psiold(j,nocc,l+1,ispin) = 0.0d0
c     use major comp. as reference
            if (rae(j).ne.0.0)
     :           psiold(j,nocc,l+1,ispin)=gf(j,iorb,1)/rae(j)
         enddo
      enddo

      write(6,*) '***        Quantum numbers        ***'
      write(6,*) 'AE             <->  PSEUDO'
      do iorb=1,norb
         write(6,*) 'n,l(AE): ', noae(iorb),lo(iorb),
     :        '      n,l(PS): ', no(iorb),lo(iorb)
      enddo


c     ---------------------------------------------------------------
c     main loop begins here
c     ---------------------------------------------------------------

      do iiter=1,max(namoeb,1)
c
c     read initial pseudopotential parameter from psp.par
c     test if parameter from atom.ae and psp.par are consistent
c
         if (iiter.eq.1) then
            write(6,*) '***        Reading data from psp.par      ***'
            open(unit=23,file='psp.par',form='formatted',
     1           status='unknown')
            read(23,*) nng,tmprij
            if (ng.eq.0) ng=nng
            if (rij.eq.0) rij=tmprij
            if (fullac) then
               ng=ngmx
            elseif (ng .gt. ngmx ) then
               write(6,*) 'gaussians:',ng
               write(6,*) 'max is   :',ngmx
               stop
            endif
            if (noccmax.gt.ng+1) stop 'noccmax>ng+1'
            read(23,*) rcov,rprb
            if (rcov.ne.rcovp) then
               write(6,*) 'rcov from atom.ae and psp.par not identical'
               write(6,*) 'atom.ae    rcov=',rcovp
               write(6,*) 'psp.par    rcov=',rcov
               if (ignore) then
                  write(6,*) 'Warning! continue program using rcov',
     :                 ' from psp.par'
               else
                  write(6,*) 'option ''-ignore'' ignores this warning'
                  stop
               endif
            endif
            if (rprb.ne.rprbp) then
               write(6,*) 'rbrb from atom.ae and psp.par not identical'
               write(6,*) 'atom.ae    rprb=',rprbp
               write(6,*) 'psp.par    rprb=',rprb
               if (ignore) then
                  write(6,*) 'Warning! continue program using rprb',
     :                 ' from psp.par'
               else
                  write(6,*) 'option ''-ignore'' ignores this warning'
                  stop
               endif
            endif
c            read(23,'(t2,a,t10,a)') ispp,icorr
            read(23,'(a)') label
            j=1
            do i=len(label),1,-1
               if (label(i:i).ne.' ') j=i
            enddo
            ispp=label(j:j)
            if (isppp.ne.ispp) then
               write(6,*) 'ispp from atom.ae and psp.par not identical'
               write(6,*) 'atom.ae    ispp=',isppp
               write(6,*) 'psp.par    ispp=',ispp
               if (ignore) then
                  write(6,*) 'Warning! continue program using ispp',
     :                 ' from psp.par'
               else
                  write(6,*) 'option ''-ignore'' ignores this warning'
                  stop
               endif
            endif


            is(1) = '  so=0'
            if (ispp.eq.'r') then
               nspin=2
               is(1)= 'so=+0.5'
               is(2)= 'so=-0.5'
            endif
            read(23,'(a)') label
            j1=1
            j2=2
            do i=len(label),1,-1
               if (label(i:i).ne.' ') j1=i
            enddo
            do i=len(label),j1,-1
               if (label(i:i).eq.' ') j2=i
            enddo
            j2=j2-1
            icorr=label(j1:j2)
            if (icorr.ne.icorrp) then
               write(6,*) 'xc-functional from  atom.ae and psp.par ',
     :              'not identical'
               write(6,*) 'atom.ae    icorr=',icorrp
               write(6,*) 'psp.par    icorr=',icorr
               if (ignore) then
                  write(6,*) 'Warning! continue program using icorr',
     :                 ' from psp.par'
               else
                  write(6,*) 'option ''-ignore'' ignores this warning'
                  stop
               endif
            endif
c     set the parameter/variables for the xc-functional(s)
C..   functionals
            salpha=2.D0/3.D0
            bbeta=0.0042D0
            betapp=0.0042D0
            if(index(icorr,'BP/PADE').ne.0) then
               mfxcx=0
               mfxcc=9
               mgcx=1
               mgcc=1
               igrad=.true.
            elseif(index(icorr,'PBE/PADE').ne.0) then
               mfxcx=0
               mfxcc=9
               mgcx=3
               mgcc=4
               igrad=.true.
            elseif(index(icorr,'LDA').ne.0) then
               mfxcx=1
               mfxcc=1
               mgcx=0
               mgcc=0
            elseif(index(icorr,'PADE').ne.0) then
               mfxcx=0
               mfxcc=9
               mgcx=0
               mgcc=0
            elseif(index(icorr,'BONL').ne.0) then
               mfxcx=1
               mfxcc=1
               mgcx=1
               mgcc=0
               igrad=.true.
            elseif(index(icorr,'BP').ne.0) then
               mfxcx=1
               mfxcc=1
               mgcx=1
               mgcc=1
               igrad=.true.
            elseif(index(icorr,'PW').ne.0) then
               mfxcx=1
               mfxcc=1
               mgcx=2
               mgcc=3
               igrad=.true.
            elseif(index(icorr,'PBE').ne.0) then
               mfxcx=1
               mfxcc=1
               mgcx=3
               mgcc=4
               igrad=.true.
            elseif(index(icorr,'BLYP').ne.0) then
               mfxcx=1
               mfxcc=3
               mgcx=1
               mgcc=2
               igrad=.true.
            elseif(index(icorr,'HCTH').ne.0) then
               mfxcx=0
               mfxcc=0
               IF (INDEX(icorr,'93').NE.0) THEN
                 mgcx=6
               ELSE IF (INDEX(icorr,'120').NE.0) THEN
                 mgcx=7
               ELSE IF (INDEX(icorr,'147').NE.0) THEN
                 mgcx=8
               ELSE IF (INDEX(icorr,'407').NE.0) THEN
                 mgcx=9
               ELSE
                 icorr='HCTH'
                 mgcx=5
               END IF
               mgcc=mgcx
               igrad=.true.
            else
               write(6,*) 'Unknown functional(s): ',icorr
               stop
            endif
            read(23,*) znuc, zion, rloc, gpot(1),gpot(2),gpot(3),gpot(4)
            if (znucp.ne.znuc) then
               write(6,*) 'znuc from atom.ae and psp.par not identical'
               write(6,*) 'atom.ae    znuc=',znucp
               write(6,*) 'psp.par    znuc=',znuc
               if (ignore) then
                  write(6,*) 'Warning! continue program using znuc',
     :                 ' from psp.par'
               else
                  write(6,*) 'option ''-ignore'' ignores this warning'
                  stop
               endif
            endif
            if (zionp.ne.zion) then
               write(6,*) 'zion from atom.ae and psp.par not identical'
               write(6,*) 'atom.ae    zion=',zionp
               write(6,*) 'psp.par    zion=',zion
               if (ignore) then
                  write(6,*) 'Warning! continue program using zion',
     :                 ' from psp.par'
               else
                  write(6,*) 'option ''-ignore'' ignores this warning'
                  stop
               endif
            endif
            read(23,*) lpx
            if (lpx .gt. lmx ) then
               write(6,*) 'array dimension problem: lpx,lpmx',lpx,lpmx
               stop
            endif
            do l=0,lpx
               read(23,*) r_l(l+1),(hsep(i,l+1,1),i=1,6)
               if (l.gt.0.and.nspin.eq.2)read(23,*)(hsep(i,l+1,2),i=1,6)
            enddo
            a0in=rloc
            write(6,*)
            write(6,*) ' Initial Pseudpotential Parameter from psp.par'
            write(6,*) ' ---------------------------------------------'
            write(6,'(2f10.3,t30,a)') rcov,rprb,
     :           'rcov,rprb'
            if (ispp.eq.'r') then
               write(6,'(t30,a)')'relativistic calculation'
            else
               write(6,'(t30,a)')'non relativistic calculation'
            endif
            write(6,'(t2,a10,t30,a)')icorr ,'XC-functional'
            write(6,*) 'local part'
            write(6,'(f5.0,f7.2,f7.3,4e11.3,t65,a)')
     :           znuc,zion,rloc,gpot(1),gpot(2),gpot(3),gpot(4),
     :           'znuc,zion, rloc, gpot() '
            if (lpx.ge.0) then
               write(6,*) 'nonlocal part'
               write(6,'(i4,t60,a)') lpx ,
     :              'lpx, (Projectors for l=0..lpx)'
               do l=0,lpx
                  write(6,*) 'l=',l
                  write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
     :                 (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
                  if (l.gt.0 .and. nspin.eq.2)
     :                 write(6,'(t8,6e11.3,t76,a)')
     :                 (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
               enddo
            endif
            close(23)
         endif
c
c     weights will be read from weights.par
c
         write(6,*)
         write(6,*) 'Reading actual weights from file weights.par'
         open(unit=24,file='weights.par',form='formatted')
         read(24,*) label
         read(24,*) label
         read(24,*) whgtp0
         read(24,*) label
         do iorb=1,norb
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            ss = so(iorb)
            if (ss.lt.0) ispin=2
            read(24,*) nw,lw,sw,(wght(nocc,l+1,ispin,i),i=1,8)
            if (noae(iorb).ne.nw .or. l.ne.lw .or. ss.ne.sw) then
               write(6,*) 'error in file ''weights.par'' '
               write(6,*) 'need weights for n,l,s:',
     :              noae(iorb),l,so(iorb)
               write(6,*) 'found            n,l,s:',nw,lw,sw
               stop
            endif
         enddo
         close(24)
c
c     calc. exponents of gaussians
         a0=a0in/rij
c     take this for an crude initial fit:
c     tt=2.d0**.4d0
         if (denbas) then
c     fine fit:
            tt=sqrt(sqrt(2.d0))
         else
c     normal fit:
            tt=2.d0**.3d0
         endif
         do i=0,ng
            a=a0*tt**i
            xp(i)=.5d0/a**2
         enddo
         write(6,*)
         write(6,*)'Basis:'
         write(6,*)'------'
         write(6,'(a,4e11.4)') ' amin,a0in,amax',a0,a0in,a
         write(6,'(a,t10,3(e11.4),a,2(e11.4))') ' gaussians ',
     &        xp(1),xp(2),xp(3),' .... ',xp(ng-1),xp(ng)
         write(6,*)'gaussians:',ng
c     set up radial grid
         nint=5*(ng+14)
         rmax=min(15.d0*rprb,120.d0)
         a_grd=a0/400.d0
         b_grd=log(rmax/a_grd)/(nint-2)
         call radgrid(nint,rr,rw,rd,a_grd,b_grd,rmax)
         write(6,'(a,t10,3(e11.4),a,2(e11.4))') ' r-grid: ',
     &        rr(1),rr(2),rr(3),' .... ',rr(nint-1),rr(nint)
         write(6,*)'gridpoints:',nint

         call crtvh(ng,lcx,lmax,xp,vh,nint,rmt,rmtg,ud,rr)
         if (namoeb.gt.0) then

c     refine simplex only every 10.th step
c     start of if-block
c          if (mod(iter,10).eq.0) then



c
c     pack initial guess
c
            call  ppack (rloc,gpot,hsep,r_l,pp(1),
     :           lpx,lpmx,nspin,nsmx,maxdim,nfit,'init',
     :           avgl1,avgl2,avgl3,ortprj,litprj)
c
c     initial simplex
c
            do i=1,nfit
               do j=1,nfit
                  pp(j+i*nfit)=pp(j)
               enddo
            enddo
            dh=0.2d0
            do i=1,nfit
CMK
              CALL random_number(rrand)
              pp(i+i*nfit)=pp(i)+dh*1.0d0*(rrand-.5d0)
c     IBM/DEC
c            pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand())-.5d0)
c     CRAY
c            pp(i+i*nfit)=pp(i)+dh*1.0d0*(ranf()-.5d0)
            enddo
            write(6,*) 'penalties for start simplex:'
            do i=1,nfit+1
               call penalty(nfit,pp(1+(i-1)*nfit),yp(i),
     :              noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :              no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :              occup,aeval,chrg,dhrg,ehrg,res,wght,
     :              wfnode,psir0,whgtp0,
     :              rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :              vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :              avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
     :              ntime,itertot)
               write(6,*)'penalty:',yp(i)
CMK            CALL FLUSH(6)
            enddo
c     refine simplex only every 10.th step
c     end of if-block
c         endif
c
c     starting amoeba
c
            ftol=1.d-7
            write(6,*) 'starting amoeba, iter=',iiter
            call AMOEBA(pp,yp,nfit,FTOL,ITER,nsmplx,namoeb,
     :           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :           no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :           occup,aeval,chrg,dhrg,ehrg,res,wght,
     :           wfnode,psir0,whgtp0,
     :           rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :           vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :           avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
     :           ntime,itertot)
            write(6,*) 'Finished amoeba with ',iter,'iterations'
         else
            energ=.true.
            call gatom(
     :           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :           occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
     :           rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :           vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
     :           rr,rw,rd,ntime,itertot)
            energ=.false.
         endif
c
c     print results
c
         write(6,*)
         write(6,*) ' Pseudpotential Parameter '
         write(6,*) ' -------------------------'
         write(6,'(2f10.3,t60,a)') rcov,rprb,'rcov,rprb'
         write(6,*) 'local part'
            write(6,'(f5.0,f7.2,f7.3,4e11.3,t65,a)')
     :        znuc,zion,rloc,gpot(1),gpot(2),gpot(3),gpot(4),
     :        'znuc, zion, rloc, gpot() '
         if (lpx.ge.0) then
            write(6,*) 'nonlocal part'
            write(6,'(i4,t60,a)') lpx ,
     :           'lpx, (Projectors for l=0..lpx)'
            do l=0,lpx
               write(6,*) 'l=',l
               write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
     :              (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
               if (l.gt.0 .and. nspin.eq.2)
     :              write(6,'(t8,6e11.3,t76,a)')
     :              (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
            enddo
            if (nspin.eq.2 .and. lpx.gt.0 ) then
               write(6,*) 'nonlocal part as V_average + V_so '
               do l=0,lpx
                  write(6,*) 'l=',l
                  do i=1,6
                     if (l.gt.0) then
                        havg(i)=((l+1)*hsep(i,l+1,1)+l*hsep(i,l+1,2))
     :                       /(2*l+1)
                        hso(i)=2*(hsep(i,l+1,1)-hsep(i,l+1,2))
     :                       /(2*l+1)
                     else
                        havg(i)=hsep(i,l+1,1)
                     endif
                  enddo
                  write(6,'(f7.3,(t8,6e11.3,t76,a))')
     :                 r_l(l+1),(havg(i),i=1,6),'r_l(),hsep_avg()'
                  if (l.gt.0)    write(6,'(t8,6e11.3,t76,a)')
     :                 (hso(i),i=1,6),'       hsep_so()'
               enddo
            endif
         endif
         write(6,*)
         write(6,'(2(tr10,a,e12.4))')
     :        'psir0 =',psir0,'weight_psir0 =',abs(psir0*whgtp0)
         write(6,*)
         write(6,'(a,t32,a,t42,a,t55,a,t64,a)')
     :        ' nl    s      occ','ae','pseudo','diff','diff*weight'

         do iorb=1,norb
            write(6,31) noae(iorb),il(lo(iorb)+1),so(iorb),zo(iorb)
 31         format(1x,i1,a1,f6.1,f10.4)
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            if (so(iorb).lt.0) ispin=2
            write(6,32) 'eigenvalue ',
     :           ev(iorb),aeval(nocc,l+1,ispin),
     :           aeval(nocc,l+1,ispin)-ev(iorb),
     :           abs(wght(nocc,l+1,ispin,1)*
     :           (aeval(nocc,l+1,ispin)-ev(iorb)))
            write(6,32) 'charge    ',
     :           crcov(iorb),chrg(nocc,l+1,ispin),
     :           chrg(nocc,l+1,ispin)-crcov(iorb),
     :           abs(wght(nocc,l+1,ispin,2)*
     :           (chrg(nocc,l+1,ispin)-crcov(iorb)))
            if (wght(nocc,l+1,ispin,3).ne.0.0d0)
     :           write(6,32) 'dcharge   ',
     :           dcrcov(iorb),dhrg(nocc,l+1,ispin),
     :           100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb)),
     :           abs(wght(nocc,l+1,ispin,3))*
     :           100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))
            if (wght(nocc,l+1,ispin,4).ne.0.0d0)
     :           write(6,32) 'echarge   ',
     :           ddcrcov(iorb),ehrg(nocc,l+1,ispin),
     :           100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb)),
     :           abs(wght(nocc,l+1,ispin,4))*
     :           100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))
            write(6,33) 'residue   ',
     :           res(nocc,l+1,ispin),
     :           abs(wght(nocc,l+1,ispin,5)*res(nocc,l+1,ispin))
            if (wght(nocc,l+1,ispin,6).ne.0.0d0)
     :           write(6,33) 'rnode     ',
     :           wfnode(nocc,l+1,ispin,1),
     :           abs(wght(nocc,l+1,ispin,6)*wfnode(nocc,l+1,ispin,1))
            if (wght(nocc,l+1,ispin,7).ne.0.0d0)
     :           write(6,33) 'dnode     ',
     :           wfnode(nocc,l+1,ispin,2),
     :           abs(wght(nocc,l+1,ispin,7)*wfnode(nocc,l+1,ispin,2))
            if (wght(nocc,l+1,ispin,8).ne.0.0d0)
     :           write(6,33) 'ddnode    ',
     :           wfnode(nocc,l+1,ispin,3),
     :           abs(wght(nocc,l+1,ispin,8)*wfnode(nocc,l+1,ispin,3))
            write(6,*)
         enddo
 32      format (t10,a,t25,4e12.4)
 33      format (t10,a,t25,2e24.4)
         write(6,*) 'diff for dcharg and echarge is given in (%)'
c
c     overwrite old values of 'psp.par' with the current ones
c
         if (namoeb.gt.0) then
            open(unit=23,file='psp.par',form='formatted',
     1           status='unknown')
            write(23,'(i5,f20.5,a)') ng,rij,   ' ng, rij '
            write(23,'(2f10.3,t60,a)') rcov,rprb ,'rcov,rprb '
            write(23,'(t2,a,t30,a)') ispp,
     :           '(non)relativistic calculation'
            write(23,'(t2,a,t30,a)') icorr,'XC-functional'
            write(23,'(2f7.3,f16.10,4e20.10,tr2,a)')
     :           znuc,zion,rloc, gpot, 'znuc,zion,rloc,gpot()'
            write(23,'(i4,t60,a)') lpx ,'lpx, (Projectors for l=0..lpx)'
            do l=0,lpx
               write(23,'(f16.10,6e20.10,tr2,a)') r_l(l+1),
     :              (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
               if (l.gt.0 .and. nspin.eq.2)
     :              write(23,'(6e20.10,tr2,a)')
     :             (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
            enddo
            write(23,*) '---------------------------------------------'
            if (info) then
               write(23,*) 'Additional information (last calculation):'
               write(23,*) 'gaussian exponents:'
               write(23,*) 'xp():',(xp(i),i=1,ng)
               write(23,*) 'orbital coefficients ',
     :              '(one orbital per column)'
               do l=0,lmax
                  do ispin=1,min(2*l+1,nspin)
                     write(23,*) 'l,ispin:',l,ispin
                     write(23,*) 'psi n=1,noccmax(l):'
                     do i=0,ng
                        write(23,'(10e20.10)')
     :                       (psi(i,nocc,l+1,ispin),
     :                       nocc=1,nomax(l))
                     enddo
                  enddo
               enddo
            endif
         endif
         close(unit=23)


c
c     PLOT WAVEFUNCTIONS (up to 5*rcov)   c.hartwigsen: version for gnuplot
c
         if (plotwf) then
            call detnp(ngrid,rae,5*rcov,np)
            open(32,file='pswf.gnu',form='formatted',status='unknown')
            write (32,*) 'set data style lines'
            do iorb=1,norb
               nocc=no(iorb)
               l=lo(iorb)
               ispin=1
               if (so(iorb).lt.0) ispin=2
               if (ispp.eq.'r') then
                  open(unit=2,file='ps.'//
     :                 char(ichar('0')+noae(iorb))//
     :                 il(lo(iorb)+1)//
     :                 char(ichar('0')+int(2*(lo(iorb)+so(iorb))))//
     :                 'by2.dat',form='formatted',status='unknown')
c     :        char(ichar('A')-ichar('a')+ichar(il(lo(iorb)+1)))//
               else
                  open(unit=2,file='ps.'//
     :                 char(ichar('0')+noae(iorb))//
     :                 il(lo(iorb)+1)//'.dat',
     :                 form='formatted',status='unknown')
               endif
c     find outer max of psi (approx), search from 10 bohr down
                  ra=10.d0
                  ttrmax=ra
                  ttmax= dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
                  do i=100,0, -1
                     ra= 0.1d0 * i
                     ttpsi=dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
c                     print*,ra,ttpsi
                     if ( ttpsi .gt. ttmax
     :                    .and. ttpsi .gt. 1.0d-4 ) then
                        ttmax=ttpsi
                        ttrmax=ra
                     endif
                  if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) goto 3456
                  enddo
 3456             continue
c     ae/pseudo wfs should have the same sign for large r when plotted
               call detnp(ngrid,rae,ttrmax,nsign)
c     never use first gridpoint! (only relevant for H and He)
               if (nsign.eq.1) nsign=nsign+1
               tt=psiold(nsign,nocc,l+1,ispin)
               sign1=tt/abs(tt)
               tt= wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(nsign))
               sign2=tt/abs(tt)
               do i=2,np
                  ttold=psiold(i,nocc,l+1,ispin)*sign1*rae(i)
                  ttold=max(min(3.d0,ttold),-3.d0)
                  tt=wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))
     :                 *sign2*rae(i)
                  tt=max(min(3.d0,tt),-3.d0)
                  ttdiff=psiold(i,nocc,l+1,ispin)*sign1-
     :                 wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))*sign2
                  ttdiff= ttdiff*rae(i)
                  ttdiff=log(max(abs(ttdiff),1.d-8))/log(10.d0)
                  write(2,'(7e20.10)') rae(i),ttold,tt,ttdiff
c     plot of the wavefunction and the higher derivatives
c                  write(2,'(7e20.10)') rae(i),
c     :                 psiold(i,nocc,l+1,ispin),
c     :                 wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)),
c     :                 dwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)),
c     :                 ddwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))
               enddo
               close(2)
               if (ispp.eq.'r')then
                  fname= 'ps.'//char(ichar('0')+noae(iorb))
     :                 //il(lo(iorb)+1)
     :                 //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))
     :                 //'by2.dat'
                  tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)
     :                 //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))
     :                 //'/2'
               else
                  fname= 'ps.'//char(ichar('0')+noae(iorb))
     :                 //il(lo(iorb)+1)//'.dat'
                  tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)
               endif
               write (32,*)'  plot "'//fname(1:index(fname,' ')-1)
     :              //'"     title "'//tname(1:index(tname,' ')-1)
     :              //'"'
               write (32,*)'replot "'//fname(1:index(fname,' ')-1)
     :              //'" using 1:3 title "pseudo"'
               write(32,*) 'pause -10 "Hit return to continue"'
               write (32,*)'replot "'//fname(1:index(fname,' ')-1)
     :              //'" using 1:4 title "diff"'
               write(32,*) 'pause -10 "Hit return to continue"'
            enddo
            write (32,*) 'set nokey'
            close(unit=32)
            write(6,*) 'to plot wfs type ''gnuplot pswf.gnu'' '
         endif
c    -----------------------------------------------
c                     MAIN LOOP END
      if (namoeb.eq.0) goto 1000
      enddo
 1000 continue


c
c     write psp-parameter in cpmd-form
c
        ngpot=0
        do j=1,4
           if (gpot(j).ne.0.d0) ngpot=j
        enddo
        open(unit=3,file='XX',form='formatted',status='unknown')
        write(3,*) '&ATOM'
        write(3,*) ' Z  =  ',znuc
        write(3,*) ' ZV =  ' ,zion
        write(3,'(a,4i1,f15.10)')
     :       '  XC = ',mfxcx,mfxcc,mgcx,mgcc,salpha
        write(3,*) ' TYPE = NORMCONSERVING GOEDECKER'
        write(3,*) '&END'
        write(3,*) '&INFO'
        write(3,*) '  Goedecker/Hartwigsen s ? PP'
        write(3,*) '&END'
        write(3,*) '&POTENTIAL'
        write(3,*) '    GOEDECKER'
        write(3,*) lpx+1  ,'                                   LMAX'
8       format(1x,f16.9,a)
        write(3,8) rloc,'                                 RC'
 9      format(1x,i3,a)
19      format(1x,i3,f16.9,a)
29      format(1x,i3,2(f16.9),a)
39      format(1x,i3,3(f16.9),a)
49      format(1x,i3,4(f16.9),a)
        if (ngpot.eq.0) then
           write(3,9) ngpot,'   #C '
        elseif (ngpot.eq.1) then
           write(3,19)   ngpot,(gpot(j),j=1,ngpot), '   #C  C1'
        else if (ngpot.eq.2) then
           write(3,29)   ngpot,(gpot(j),j=1,ngpot), '   #C  C1 C2'
        else if (ngpot.eq.3) then
           write(3,39)   ngpot,(gpot(j),j=1,ngpot), '   #C  C1 C2 C3'
        else
           write(3,49)   ngpot,(gpot(j),j=1,ngpot), '   #C  C1 C2 C3 C4'
        endif
        if (lpx.ge.0) then
           do l=0,lpx
              do i=1,6
                 if (nspin.eq.1 .or. l.eq.0 ) then
                    havg(i)=hsep(i,l+1,1)
                 else
                    havg(i)=((l+1)*hsep(i,l+1,1)+l*hsep(i,l+1,2))
     :                   /(2*l+1)
                 endif
              enddo
              npj=0
              if ( abs(havg(1)).gt.1.0d-8) npj=1
              if ( abs(havg(3)).gt.1.0d-8) npj=2
              if ( abs(havg(6)).gt.1.0d-8) npj=3
              if (npj .eq. 1 ) then
                 string=  ' H('//il(l+1)//') 11'
                 form  = '(1x,f16.9,i3,1(f14.9),a)'
                 nhsep = 1
              elseif (npj .eq. 2 ) then
                 string = ' H('//il(l+1)//') 11 12 22'
                 form  = '(1x,f16.9,i3,3(f14.9),a)'
                 nhsep = 3
              elseif (npj .eq. 3 ) then
                 string = ' H('//il(l+1)//') 11 12 13 22 23 33'
                 form  = '(1x,f16.9,i3,6(f14.9),a)'
                 nhsep = 6
              else
                 string = ' H('//il(l+1)//') no projector'
                 form  = '(1x,f16.9,i3,a)'
                 nhsep = 0
              endif
c
c     cpmd uses different ordering of the H-matrix
c     only relevant for npj > 2
              do i=1,6
                 hhsep(i) = havg(i)
              enddo
              if (npj.eq.3) then
                 hhsep(1) =havg(1)
                 hhsep(2) =havg(2)
                 hhsep(3) =havg(4)
                 hhsep(4) =havg(3)
                 hhsep(5) =havg(5)
                 hhsep(6) =havg(6)
              endif
              write(3,form) r_l(l+1),npj,(hhsep(i),i=1,nhsep),string
           enddo
        endif
        write(3,*) '&END'
c
c     test orthogonality of the projectors
c
      if (lpx.ge.0.and.ortprj)
     :     call pj2test(hsep,lpx,lpmx,lmx,nspin,nsmx,r_l,is)
c
      write(6,*) 'Total SCF-cycles:',itertot
      write(6,*) 'Pseudoatom calculations:',ntime
      write(6,*) '*************finished*****************'
      end


c
c     CRAY: no derf() -> user erf()
c
c      real*8 function derf(x)
c      REAL*8 X
c      DERF=ERF(X)
c      RETURN
c      END
