      subroutine input
      use atom
      use rint
      use pspot
      use xcfcn
      implicit real*8(a-h,o-z)
C..local
      parameter(pi=3.14159265358979323846264D0)
      character*8 string,estr
      character*80 line
      logical err,ppsg,pp_read
C..defaults
      optexp=.false.
      allelectron=.true.
      add_pot_factor=0.d0
      zval=0.d0
      zeff=0.d0
      ippn=100
      dmix=0.25
      epsscf=1.d-6
      maxscf=200
      lmax=-1
      do l=0,lamax
        nocc(l)=0
        nalpha(l)=0
        do i=1,namax
          occ(i,l)=0.0d0
          alpha(i,l)=0.0d0
        enddo
      enddo

      pptype=4
      ERFnr=0
      EXPnr=0
      PPlmax=0
      KBPROJnr=0
      KBV=0.d0
      pp_read=.false.
      ppsg=.false.
      r_proj=0.d0
      r_loc=0.d0

      xcstring=' '
      xcfstring=' '
      gcxstring=' '
      gccstring=' '
      mfxcx=1
      mfxcc=1
      mgcx=0
      mgcc=0
      bbeta=0.0042D0
      betapp=0.0042D0
      salpha=2.D0/3.D0





C..look for section ATOM
      ierr=inscan(5,'&ATOM')
      if(ierr.ne.0) then
         write(*,*) ' could not find section &ATOM '
         stop
      endif
C..input section &ATOM
 100  continue
      call cfield(string,4)
      if(index(string,'&END').ne.0) then
         goto 101
      elseif(index(string,'NAME').ne.0) then
         call cfield(atomname,2)
         goto 100
      elseif(index(string,'ZVAL').ne.0) then
         call realfl(rnumber,err,estr)
         if(err) call rerror(string,estr)
         zval=rnumber
         goto 100
      elseif(index(string,'ZEFF').ne.0) then
         call realfl(rnumber,err,estr)
         if(err) call rerror(string,estr)
         zeff=rnumber
         goto 100
      elseif(index(string,'RCOV').ne.0) then
         call realfl(rnumber,err,estr)
         print*,' RCOV =',rnumber
         add_pot_factor=0.5d0*4.d0/(2.d0*rnumber)**4
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'RPRB').ne.0) then
         call realfl(rnumber,err,estr)
         print*,' RPRB =',rnumber
         add_pot_factor=0.5D0/rnumber**4
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'ADDP').ne.0) then
         call realfl(rnumber,err,estr)
         print*,' ADDP =',rnumber
         add_pot_factor=rnumber*0.5d0
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'PPOT').ne.0) then
         allelectron=.false.
         goto 100
      elseif(index(string,'OPTI').ne.0) then
         optexp=.true.
         call realfl(epsgrad,err,estr)
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'LMAX').ne.0) then
         call intgfl(lmax,err,estr)
         print*,' LMAX =',lmax
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'OCCU').ne.0) then
         if(lmax.lt.0) then
            write(*,*) ' Specify LMAX before OCCUPATION'
            call rerror(string,string)
         endif
         do l=0,lmax
            call intgfl(nocc(l),err,estr)
            if(err) call rerror(string,estr)
            do i=1,nocc(l)
               call realfl(oo,err,estr)
               print*,' OCCU(',i,',',l,') = ',oo
               if(err) call rerror(string,estr)
               occ(i,l)=oo/(2.D0*l+1.D0)
            enddo
         enddo
         goto 100
      elseif(index(string,'XCFN').ne.0) then
         string=' '
         call cfield(string,10)
         xcstring=string
         xcinfo = ""
         if (index(string,'NONE').ne.0) then
            mfxcx=0
            mfxcc=0
            mgcx=0
            mgcc=0
            xcinfo = "No exchange-correlation functional"
         elseif (index(string,'PZ').ne.0) then
            mfxcx=1
            mfxcc=1
            mgcx=0
            mgcc=0
            xcinfo = "Perdew, Zunger, PRB 23, 5048 (1981)"
         elseif (index(string,'VWN').ne.0) then
            mfxcx=1
            mfxcc=2
            mgcx=0
            mgcc=0
            xcinfo ="Vosko, Wilk, Nusair, Can. J. Phys. 58, 1200 (1980)"
         elseif ((index(string,'PADE').ne.0).or.
     &           (index(string,'LDA').ne.0)) then
            mfxcx=0
            mfxcc=9
            mgcx=0
            mgcc=0
            xcstring = "LDA"
            xcinfo = "Goedecker, Teter, Hutter, PRB 54, 1703 (1996)"
         elseif (index(string,'BONL').ne.0) then
            mfxcx=1
            mfxcc=1
            mgcx=1
            mgcc=0
            xcinfo = "Becke (1988); Perdew-Zunger (1981)"
         elseif (index(string,'BP').ne.0) then
            mfxcx=1
            mfxcc=1
            mgcx=1
            mgcc=1
            xcinfo = "Becke (1988); Perdew (1986)"
         elseif (index(string,'PW').ne.0) then
            mfxcx=1
            mfxcc=1
            mgcx=2
            mgcc=3
            xcinfo = "Perdew, Wang, PRB 46, 6671 (1992)"
         elseif (index(string,'PBE').ne.0) then
            mfxcx=1
            mfxcc=1
            mgcx=3
            mgcc=4
            xcinfo = "Perdew, Burke, Ernzerhof (1996)"
         elseif (index(string,'BLYP').ne.0) then
            mfxcx=1
            mfxcc=3
            mgcx=1
            mgcc=2
            xcinfo ="Becke (1988); Lee, Yang, Parr (1988)"
         elseif (index(string,'HCTH').ne.0) then
            mfxcx=0
            mfxcc=0
            IF (INDEX(string,'93').NE.0) THEN
              mgcx=6
            ELSE IF (INDEX(string,'120').NE.0) THEN
              mgcx=7
            ELSE IF (INDEX(string,'147').NE.0) THEN
              mgcx=8
            ELSE IF (INDEX(string,'407').NE.0) THEN
              mgcx=9
            ELSE
              mgcx=5
            END IF
            mgcc=mgcx
            xcinfo ="Hamprecht, Cohen, Tozer, Handy"
         else
            call rerror(string,string)
         endif
         goto 100
      elseif(index(string,'EXCF').ne.0) then
         call cfield(string,4)
         xcfstring=string
         if (index(string,'NONE').ne.0) then
            mfxcc=0
            mfxcx=0
         elseif (index(string,'CA').ne.0) then
            mfxcc=1
         elseif (index(string,'VWN').ne.0) then
            mfxcc=2
         elseif (index(string,'LYP').ne.0) then
            mfxcc=3
         elseif (index(string,'PW').ne.0) then
            mfxcc=4
         elseif (index(string,'WIGN').ne.0) then
            mfxcc=5
         elseif (index(string,'HELU').ne.0) then
            mfxcc=6
         elseif (index(string,'OBPZ').ne.0) then
            mfxcc=7
         elseif (index(string,'OBPW').ne.0) then
            mfxcc=8
         else
            write(*,*) ' Unknown keyword :',string
            call rerror(string,string)
         endif
         goto 100
      elseif(index(string,'EXGC').ne.0) then
         call cfield(string,4)
         gcxstring=string
         if (index(string,'NONE').ne.0) then
            mgcx=0
         elseif (index(string,'BECK').ne.0) then
            mgcx=1
         elseif (index(string,'GGA').ne.0) then
            mgcx=2
         else
            write(*,*) ' Unknown keyword :',string
         endif
         goto 100
      elseif(index(string,'ECGC').ne.0) then
         call cfield(string,4)
         gccstring=string
         if (index(string,'NONE').ne.0) then
            mgcc=0
         elseif (index(string,'PERD').ne.0) then
            mgcc=1
         elseif (index(string,'LYP').ne.0) then
            mgcc=2
         elseif (index(string,'GGA').ne.0) then
            mgcc=3
         else
            call rerror(string,string)
         endif
         goto 100
      elseif(index(string,'MIXI').ne.0) then
         call realfl(dmix,err,estr)
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'CONV').ne.0) then
         call realfl(epsscf,err,estr)
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'ITER').ne.0) then
         call intgfl(maxscf,err,estr)
         if(err) call rerror(string,estr)
         goto 100
      elseif(index(string,'IPPN').ne.0) then
         call intgfl(ippn,err,estr)
         if(err) call rerror(string,estr)
         if(ippn.gt.ippmax) then
            write(*,*) ' IPPN must be less than',ippmax+1,'!'
            call rerror(string,string)
         endif
         goto 100
      else
         call rerror(string,string)
         goto 100
      endif
 101  continue
      if (mfxcc.eq.3) then
         if (.not.(mgcc.eq.2.or.mgcc.eq.0)) then
            write (*,*) ' LYP must be chosen for EXCF and for ECGC !'
	    call rerror(string,string)
         endif
      endif
      if (xcstring(1:1)==' ') then
         xcstring=trim(xcfstring)//trim(gcxstring)//trim(gccstring)
      endif


C..look for section BASIS
      ierr=inscan(5,'&BASIS')
      if(ierr.ne.0) then
         write(*,*) ' could not find section &BASIS '
         stop
      endif
C..input section &BASIS
 200  continue
      call cfield(string,4)
      if(index(string,'&END').ne.0) then
         goto 201
      elseif(index(string,'GAUS').ne.0) then
         do l=0,lmax
            call intgfl(nalpha(l),err,estr)
            print*,' NALPHA(',l,') = ',nalpha(l)
            if(err) call rerror(string,estr)
            do i=1,nalpha(l)
               call realfl(aa,err,estr)
               print*,' GAUS(',i,',',l,') = ',aa
               if(err) call rerror(string,estr)
               alpha(i,l)=aa
            enddo
            do i=1,nalpha(l)
               call intgfl(alpp(i,l),err,estr)
               print*,' ALPP(',i,',',l,') = ',alpp(i,l)
               if(err) call rerror(string,estr)
            enddo
         enddo
         goto 200
      else
         write(*,*) ' Unknown keyword :',string
         goto 200
      endif
 201  continue




C..look for section PSEUDO    ATTENTION: There is another part for keyword
c                                        POTENTIAL !!!
      if (.not.allelectron) then
         ierr=inscan(5,'&PSEUDO')
         if (ierr.eq.0) then
            pp_read=.true.
C..input section &PSEUDO
 400        continue
            call cfield(string,4)
            if(index(string,'&END').ne.0) then
               goto 401
            elseif(index(string,'KBPP').ne.0) then
               pptype=5
               ppstring='Kleinman-Bylander type'
               goto 400
            elseif(index(string,'PGTH').ne.0) then
               ppsg=.true.
               pptype=5
               ppstring='Goedecker'
               ERFnr=1
               PPLerf(1)=1.d0
               PPlmax=-1
               KBEXPnr(0)=1
               KBEXPnr(1)=1
               KBEXPnr(2)=1
               KBEXPnr(3)=1
               KBREXP(1,1,0)=0
               KBREXP(1,2,0)=2
               KBREXP(1,3,0)=4
               KBREXP(1,1,1)=1
               KBREXP(1,2,1)=3
               KBREXP(1,3,1)=5
               KBREXP(1,1,2)=2
               KBREXP(1,2,2)=4
               goto 400
            elseif(index(string,'BHSA').ne.0) then
               pptype=2
               ppstring='Bachelet-Hamann-Schl"uter type'
               goto 400
            elseif(index(string,'CGEN').ne.0) then
               pptype=3
               ppstring='Bachelet-Hamann-Schl"uter type'
               goto 400
            elseif(index(string,'BHSC').ne.0) then
               pptype=1
               ppstring='Bachelet-Hamann-Schl"uter type'
               goto 400
            elseif(index(string,'ZEFF').ne.0) then
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               zeff=rnumber
               goto 400
            elseif(index(string,'RLOC').ne.0) then
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               r_loc=rnumber
               PPNerf(1)=1.d0/(sqrt(2.d0)*r_loc)
               PPNexp(1,0)=PPNerf(1)**2
               goto 400
            elseif(index(string,'C1').ne.0) then
               EXPnr=max(EXPnr,1)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               PPC(1)=rnumber
               PPLexp(1,0)=rnumber
               goto 400
            elseif(index(string,'C2').ne.0) then
               EXPnr=max(EXPnr,2)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               PPC(2)=rnumber
               PPLexp(2,0)=rnumber*2*PPNexp(1,0)
               goto 400
            elseif(index(string,'C3').ne.0) then
               EXPnr=max(EXPnr,3)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               PPC(3)=rnumber
               PPLexp(3,0)=rnumber*4*PPNexp(1,0)**2
               goto 400
            elseif(index(string,'C4').ne.0) then
               EXPnr=max(EXPnr,4)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               PPC(4)=rnumber
               PPLexp(4,0)=rnumber*8*PPNexp(1,0)**3
               goto 400
            elseif(index(string,'HS1').ne.0) then
               KBPROJnr(0)=max(KBPROJnr(0),1)
               PPlmax=max(PPlmax,0)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               KBV(1,1,0)=rnumber
               goto 400
            elseif(index(string,'HS2').ne.0) then
               KBPROJnr(0)=max(KBPROJnr(0),2)
               PPlmax=max(PPlmax,0)
               call realfl(rnumber,err,estr)
               if (err) call rerror(string,estr)
               KBV(2,2,0)=rnumber
               goto 400
            elseif(index(string,'HP1').ne.0) then
               KBPROJnr(1)=max(KBPROJnr(1),1)
               PPlmax=max(PPlmax,1)
               call realfl(rnumber,err,estr)
               if (err) call rerror(string,estr)
               KBV(1,1,1)=rnumber
               goto 400
            elseif(index(string,'HS').ne.0) then
               PPlmax=max(PPlmax,0)
               call intgfl(inumber,err,estr)
               KBPROJnr(0)=inumber
               do ih1=1,KBPROJnr(0)
                  do ih2=ih1,KBPROJnr(0)
                     call realfl(rnumber,err,estr)
                     if(err) call rerror(string,estr)
                     KBV(ih1,ih2,0)=rnumber
                     KBV(ih2,ih1,0)=rnumber
                  enddo
               enddo
               goto 400
            elseif(index(string,'HP').ne.0) then
               PPlmax=max(PPlmax,1)
               call intgfl(inumber,err,estr)
               KBPROJnr(1)=inumber
               do ih1=1,KBPROJnr(1)
                  do ih2=ih1,KBPROJnr(1)
                     call realfl(rnumber,err,estr)
                     if(err) call rerror(string,estr)
                     KBV(ih1,ih2,1)=rnumber
                     KBV(ih2,ih1,1)=rnumber
                  enddo
               enddo
               goto 400
            elseif(index(string,'HD').ne.0) then
               PPlmax=max(PPlmax,2)
               call intgfl(inumber,err,estr)
               KBPROJnr(2)=inumber
               do ih1=1,KBPROJnr(2)
                  do ih2=ih1,KBPROJnr(2)
                     call realfl(rnumber,err,estr)
                     if(err) call rerror(string,estr)
                     KBV(ih1,ih2,2)=rnumber
                     KBV(ih2,ih1,2)=rnumber
                  enddo
               enddo
               goto 400
            elseif(index(string,'RS').ne.0) then
               PPlmax=max(PPlmax,0)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               r_proj(0)=rnumber
               goto 400
            elseif(index(string,'RP').ne.0) then
               PPlmax=max(PPlmax,1)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               r_proj(1)=rnumber
               goto 400
            elseif(index(string,'RD').ne.0) then
               PPlmax=max(PPlmax,2)
               call realfl(rnumber,err,estr)
               if(err) call rerror(string,estr)
               r_proj(2)=rnumber
               goto 400
            elseif(index(string,'ERFC').ne.0) then
               call intgfl(ERFnr,err,estr)
               if(err) call rerror(string,estr)
               do i=1,ERFnr
                  call realfl(PPNerf(i),err,estr)
                  if(err) call rerror(string,estr)
               enddo
               do i=1,ERFnr
                  call realfl(PPLerf(i),err,estr)
                  if(err) call rerror(string,estr)
               enddo
               goto 400
            elseif(index(string,'PPLM').ne.0) then
               call intgfl(PPlmax,err,estr)
               if(err) call rerror(string,estr)
               goto 400
            elseif(index(string,'NLCS').ne.0) then
               call intgfl(KBPROJnr(0),err,estr)
               if(err) call rerror(string,estr)
               call intgfl(KBEXPnr(0),err,estr)
               if(err) call rerror(string,estr)
               do i=1,KBPROJnr(0)
                  do j=i,KBPROJnr(0)
                     call realfl(KBV(i,j,0),err,estr)
                     if(err) call rerror(string,estr)
                     KBV(j,i,0)=KBV(i,j,0)
                  enddo
               enddo
               do i=1,KBPROJnr(0)
                  do j=1,KBEXPnr(0)
                     call realfl(KBNexp(j,i,0),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               do i=1,KBPROJnr(0)
                  do j=1,KBEXPnr(0)
                     call intgfl(KBRexp(j,i,0),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               do i=1,KBPROJnr(0)
                  do j=1,KBEXPnr(0)
                     call realfl(KBLexp(j,i,0),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               goto 400
            elseif(index(string,'NLCP').ne.0) then
               call intgfl(KBPROJnr(1),err,estr)
               if(err) call rerror(string,estr)
               call intgfl(KBEXPnr(1),err,estr)
               if(err) call rerror(string,estr)
               do i=1,KBPROJnr(1)
                  do j=i,KBPROJnr(1)
                     call realfl(KBV(i,j,1),err,estr)
                     if(err) call rerror(string,estr)
                     KBV(j,i,1)=KBV(i,j,1)
                  enddo
               enddo
               do i=1,KBPROJnr(1)
                  do j=1,KBEXPnr(1)
                     call realfl(KBNexp(j,i,1),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               do i=1,KBPROJnr(1)
                  do j=1,KBEXPnr(1)
                     call intgfl(KBRexp(j,i,1),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               do i=1,KBPROJnr(1)
                  do j=1,KBEXPnr(1)
                     call realfl(KBLexp(j,i,1),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               goto 400
            elseif(index(string,'NLCD').ne.0) then
               call intgfl(KBPROJnr(2),err,estr)
               if(err) call rerror(string,estr)
               call intgfl(KBEXPnr(2),err,estr)
               if(err) call rerror(string,estr)
               do i=1,KBPROJnr(2)
                  do j=i,KBPROJnr(2)
                     call realfl(KBV(i,j,2),err,estr)
                     if(err) call rerror(string,estr)
                     KBV(j,i,2)=KBV(i,j,2)
                  enddo
               enddo
               do i=1,KBPROJnr(2)
                  do j=1,KBEXPnr(2)
                     call realfl(KBNexp(j,i,2),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               do i=1,KBPROJnr(2)
                  do j=1,KBEXPnr(2)
                     call intgfl(KBRexp(j,i,2),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               do i=1,KBPROJnr(2)
                  do j=1,KBEXPnr(2)
                     call realfl(KBLexp(j,i,2),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               enddo
               goto 400
            elseif(index(string,'GCFS').ne.0) then
               call intgfl(EXPnr(0),err,estr)
               if(err) call rerror(string,estr)
               do i=1,EXPnr(0)
                  call realfl(PPNexp(i,0),err,estr)
                  if(err) call rerror(string,estr)
               enddo
               if (pptype.ge.3) then
                  do i=1,(EXPnr(0))
                     call intgfl(PPRexp(i,0),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
                  do i=1,(EXPnr(0))
                     call realfl(PPLexp(i,0),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               else
                  do i=1,(2*EXPnr(0))
                     call realfl(PPLexp(i,0),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               endif
               goto 400
            elseif(index(string,'GCFP').ne.0) then
               call intgfl(EXPnr(1),err,estr)
               if(err) call rerror(string,estr)
               do i=1,EXPnr(1)
                  call realfl(PPNexp(i,1),err,estr)
                  if(err) call rerror(string,estr)
               enddo
               if (pptype.ge.3) then
                  do i=1,(EXPnr(1))
                     call intgfl(PPRexp(i,1),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
                  do i=1,(EXPnr(1))
                     call realfl(PPLexp(i,1),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               else
                  do i=1,(2*EXPnr(1))
                     call realfl(PPLexp(i,1),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               endif
               goto 400
            elseif(index(string,'GCFD').ne.0) then
               call intgfl(EXPnr(2),err,estr)
               if(err) call rerror(string,estr)
               do i=1,EXPnr(2)
                  call realfl(PPNexp(i,2),err,estr)
                  if(err) call rerror(string,estr)
               enddo
               if (pptype.ge.3) then
                  do i=1,(EXPnr(2))
                     call intgfl(PPRexp(i,2),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
                  do i=1,(EXPnr(2))
                     call realfl(PPLexp(i,2),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               else
                  do i=1,(2*EXPnr(2))
                     call realfl(PPLexp(i,2),err,estr)
                     if(err) call rerror(string,estr)
                  enddo
               endif
               goto 400
            else
               write(*,*) ' Unknown keyword :',string
               goto 400
            endif
 401        continue
         endif
      endif




C..look for section POTENTIAL
CMK    if (.not.allelectron) then
         ierr=inscan(5,'&POTENTIAL')
         if (ierr.eq.0) then
            string='&POTENTI'
            pp_read=.true.
            ppsg=.true.
            pptype=5
            ppstring='Goedecker'
            ERFnr=1
            EXPnr=0
            PPLerf(1)=1.d0
            KBEXPnr(0)=1
            KBEXPnr(1)=1
            KBEXPnr(2)=1
            KBREXP(1,1,0)=0
            KBREXP(1,2,0)=2
            KBREXP(1,3,0)=4
            KBREXP(1,1,1)=1
            KBREXP(1,2,1)=3
            KBREXP(1,3,1)=5
            KBREXP(1,1,2)=2
            KBREXP(1,2,2)=4

            read(5,*) line

            call intgfl(inumber,err,estr)
            if (err) call rerror(string,estr)
            pplmax=inumber-1
            call newline

            call realfl(rnumber,err,estr)
            if (err) call rerror(string,estr)
            r_loc=rnumber
            PPNerf(1)=1.d0/(sqrt(2.d0)*r_loc)
            PPNexp(1,0)=PPNerf(1)**2
            call newline

            call intgfl(inumber,err,estr)
            if (err) call rerror(string,estr)
            EXPnr=inumber
            do i=1,EXPnr(0)
               call realfl(rnumber,err,estr)
               if (err) call rerror(string,estr)
               PPC(i)=rnumber
               PPLexp(i,0)=rnumber*(2*PPNexp(1,0))**(i-1)
            enddo
            call newline

            do l=0,pplmax
               call realfl(rnumber,err,estr)
               if (err) call rerror(string,estr)
               r_proj(l)=rnumber

               call intgfl(inumber,err,estr)
               if (err) call rerror(string,estr)
               KBPROJnr(l)=inumber

               do i=1,KBPROJnr(l)
                  do j=i,KBPROJnr(l)
                     call realfl(rnumber,err,estr)
                     if (err) call rerror(string,estr)
                     KBV(i,j,l)=rnumber
                     KBV(j,i,l)=rnumber
                  enddo
               enddo
               call newline

            enddo

         endif
CMK    endif


      if ((.not.allelectron).and.(.not.pp_read)) then
CMK   if (.not.pp_read) then
         write(*,*) ' could not find a pseudopotential section!'
         stop
      endif
      if (allelectron) Zeff=Zval
      if (ppsg) then
         do l=0,PPlmax
            KBNexp(1,1,l)=1.d0/(2*r_proj(l)**2)
            do ip=2,KBPROJnr(l)
               KBNexp(1,ip,l)=KBNexp(1,1,l)
            enddo
            do ip=1,KBPROJnr(l)
               iexp=l+2*(ip-1)
               KBLexp(1,ip,l)
     &              =1.d0/sqrt(gint(2*iexp+2,2*KBNexp(1,ip,l)))        
            enddo
         enddo
      endif
cdeb:
c      print*,Zeff,PPlmax,ERFnr,EXPnr
c      print*,PPNerf(1:ERFnr)
c      print*,PPLerf(1:ERFnr)
c      do l=0,PPlmax
c         do i=1,EXPnr(l)
c            print*,PPLexp(i,l)
c            print*,PPRexp(i,l)
c         enddo
c      enddo
c      print*,KBV
c      print*,KBEXPnr
c      print*,KBPROJnr
c      print*,KBLexp
c      print*,KBRexp
c      print*,KBNexp
c      print*,pptype,ppstring
c
       print *, 'Zval = ',zval
      return
      end
c
c
      subroutine rerror(s1,s2)
      character*8 s1,s2
      write(*,*) ' Error while reading option for keyword ',s1
      write(*,*) ' Last token read was ',s2
      stop
c
CMK   return
      end
