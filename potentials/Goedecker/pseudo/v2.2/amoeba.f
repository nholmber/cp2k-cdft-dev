c     cwh: if amoeba doesn't improve within trymax iterations
c     amoeba also returns's

      SUBROUTINE AMOEBA(P,Y,ndim,FTOL,ITER,itmax,namoeb,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,whgtp0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :     ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,ntime,itertot)

      implicit real*8 (a-h,o-z)
      logical   avgl1,avgl2,avgl3,ortprj,litprj,igrad,lexit,lnext
      PARAMETER (NMAX=50,ALPHA=1.0d0,BETA=0.5d0,GAMMA=2.0d0)
      parameter ( trymax = 200 )

      DIMENSION P(ndim,ndim+1),Y(ndim+1),PR(NMAX),PRR(NMAX),PBAR(NMAX)
      dimension no(*),lo(*),so(*),ev(*),crcov(*),dcrcov(*),ddcrcov(*),
     :     occup(*),aeval(*),chrg(*),dhrg(*),ehrg(*),res(*),wght(*),
     :     gpot(*),r_l(*),hsep(*),vh(*),xp(*),rmt(*),rmtg(*),
     :     ud(*),psi(*),wfnode(*),rr(*),rw(*),rd(*)
c      print*,'entered amoeba with nfit=',ndim
      if (ndim.gt.nmax-1) stop 'nmax'
      if (ndim.lt.1) then
         write(6,*)'entered amoeba with nfit=',ndim
         write(6,*)'no fit!'
         return
      endif
      MPTS=NDIM+1
      ITER=0
      trycount = 0
      iloold = 1
 1    ILO=1
      IF(Y(1).GT.Y(2))THEN
         IHI=1
         INHI=2
      ELSE
         IHI=2
         INHI=1
      ENDIF
      DO 11 I=1,MPTS
         IF(Y(I).LT.Y(ILO)) ILO=I
         IF(Y(I).GT.Y(IHI))THEN
            INHI=IHI
            IHI=I
         ELSE IF(Y(I).GT.Y(INHI))THEN
            IF(I.NE.IHI) INHI=I
         ENDIF
 11   CONTINUE

c     cwh
      if (ilo .ne. iloold) trycount = 0
      iloold = ilo 

c     RTOL=min(Y(IHI),(y(ihi)-y(ilo))/y(ilo)**2)
      RTOL=min(Y(IHI),(y(ihi)-y(ilo))/y(ilo))
      IF(RTOL.LT.FTOL) then
c     Check   
         write(6,*) 'values at the edges of the simplex:'
         write(6,'(40e15.8)') y
         do i=1,ndim
            if (y(i).lt.y(ilo)) write(6,*) 'WARNING ilo not lowest'
         enddo
c 	unpack variables
c         print*,'call with ndim:',ndim
         call  ppack (rloc,gpot,hsep,r_l,p(1,ilo),
     :        lpx,lpmx,nspin,nsmx,ndim,ndim,'unpack',avgl1,avgl2,avgl3,
     :        ortprj,litprj)
c     cwh
c        write(6,*) 'best values in amoeba:'
c        write(6,*) rloc,r_l(1),gpot(1),gpot(2),gpot(3)
         RETURN
      endif
c      if (mod(iter,100).eq.0) write(6,*) 'iter=',iter
c
      IF(ITER.gt.ITMAX) then
         write(6,*) 'WARNING: NO CONVERGENCE IN AMOEBA'
c         ftol=10.d0*ftol
c         write(6,*) 'FTOL SET TO ',FTOL
         write(6,*) 'values at the edges of the simplex:'
         write(6,'(40e15.8)') y
c         write(6,*)'ilo:',ilo
         do i=1,ndim
            if (y(i).lt.y(ilo)) write(6,*) 'WARNING ilo not lowest'
         enddo

c     unpack variables
         call  ppack (rloc,gpot,hsep,r_l,p(1,ilo),
     :        lpx,lpmx,nspin,nsmx,ndim,ndim,'unpack',avgl1,avgl2,avgl3,
     :        ortprj,litprj)
         RETURN
      endif
C     CWH
      IF(trycount.ge.trymax) then
         write(6,*) 'WARNING: NO IMPROVEMENT IN AMOEBA FOR THE LAST',
     :        TRYCOUNT,'STEPS'

c     unpack variables
         call  ppack (rloc,gpot,hsep,r_l,p(1,ilo),
     :        lpx,lpmx,nspin,nsmx,ndim,ndim,'unpack',avgl1,avgl2,avgl3,
     :        ortprj,litprj)
         RETURN
      ENDIF
      trycount = trycount +1
      if (mod(iter,10).eq.0) then
         write(6,*) 'iter (amoeba,gatom):',iter,ntime,
     :        ' y(ilo):',y(ilo)
CMK      CALL FLUSH(6)
         INQUIRE ( FILE = 'NEXT', EXIST = lnext )
         INQUIRE ( FILE = 'EXIT', EXIST = lexit )
         if (lnext .or. lexit) then
            write(6,*) 'leaving amoeba'
            call  ppack (rloc,gpot,hsep,r_l,p(1,ilo),
     :        lpx,lpmx,nspin,nsmx,ndim,ndim,'unpack',avgl1,avgl2,avgl3,
     :        ortprj,litprj)
         endif
         if (lnext) then
            open(99,file='NEXT')
            close(99, status='DELETE')
            return
         endif
         if (lexit) then
c     reset counter for minimsation-cycles
c     so that main programm finishes
            namoeb=0
            open(99,file='EXIT')
            close(99, status='DELETE')
            return
         endif
      endif

      ITER=ITER+1
      DO 12 J=1,NDIM
         PBAR(J)=0.d0
 12   CONTINUE
      DO 14 I=1,MPTS
         IF(I.NE.IHI)THEN
            DO 13 J=1,NDIM
               PBAR(J)=PBAR(J)+P(j,i)
 13         CONTINUE
         ENDIF
 14   CONTINUE
      DO 15 J=1,NDIM
         PBAR(J)=PBAR(J)/NDIM
         PR(J)=(1.d0+ALPHA)*PBAR(J)-ALPHA*P(j,IHI)
 15   CONTINUE
      call penalty(ndim,pr,ypr,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,whgtp0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :     ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,ntime,itertot)
      IF(YPR.LE.Y(ILO))THEN
        DO 16 J=1,NDIM
          PRR(J)=GAMMA*PR(J)+(1.d0-GAMMA)*PBAR(J)
16      CONTINUE
        call penalty(ndim,prr,yprr,
     :       noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :       no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :       occup,aeval,chrg,dhrg,ehrg,res,wght,
     :       wfnode,psir0,whgtp0,
     :       rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :       vh,xp,rmt,rmtg,
     :       ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :       ortprj,litprj,igrad,rr,rw,rd,ntime,itertot)
        IF(YPRR.LT.Y(ILO))THEN
          DO 17 J=1,NDIM
            P(j,IHI)=PRR(J)
17        CONTINUE
c          if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,a,i2,e15.7)')
c     :         'iter',iter,' found',YPRR,' reject',ihi,Y(IHI)
          Y(IHI)=YPRR
        ELSE
          DO 18 J=1,NDIM
            P(j,IHI)=PR(J)
18        CONTINUE
c          if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,a,i2,e15.7)')
c     :         'iter',iter,' found',YPR,' reject',ihi,Y(IHI)
          Y(IHI)=YPR
        ENDIF
      ELSE IF(YPR.GE.Y(INHI))THEN
        IF(YPR.LT.Y(IHI))THEN
          DO 19 J=1,NDIM
            P(j,IHI)=PR(J)
19        CONTINUE
c          if (mod(iter,10).eq.0) 
c     :         write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))')
c     :         'iter',iter,' found',YPR,' reject',ihi,Y(IHI),
c     :         ' best:',ilo,Y(Ilo)
          Y(IHI)=YPR
        ENDIF
        DO 21 J=1,NDIM
          PRR(J)=BETA*P(j,IHI)+(1.d0-BETA)*PBAR(J)
21      CONTINUE
      call penalty(ndim,prr,yprr,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :           wfnode,psir0,whgtp0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :       ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,ntime,itertot)
        IF(YPRR.LT.Y(IHI))THEN
          DO 22 J=1,NDIM
            P(j,IHI)=PRR(J)
22        CONTINUE
c          if (mod(iter,10).eq.0)  
c     :         write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))')
c     :         'iter',iter,' found',YPRR,' reject',ihi,Y(IHI),
c     :         ' best:',ilo,Y(Ilo)
          Y(IHI)=YPRR
        ELSE
          DO 24 I=1,MPTS
            IF(I.NE.ILO)THEN
              DO 23 J=1,NDIM
                PR(J)=0.5d0*(P(j,I)+P(j,ILO))
                P(j,I)=PR(J)
23            CONTINUE
      call penalty(ndim,pr,ypr,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :             wfnode,psir0,whgtp0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :             ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,ntime,itertot)
            ENDIF
24        CONTINUE
        ENDIF
      ELSE
        DO 25 J=1,NDIM
          P(j,IHI)=PR(J)
25      CONTINUE
c          if (mod(iter,10).eq.0) 
c     :         write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))') 
c     :         'iter',iter,' found',YPR,' reject',ihi,Y(IHI),
c     :         ' best:',ilo,Y(Ilo)
        Y(IHI)=YPR
      ENDIF
      GO TO 1
      END


