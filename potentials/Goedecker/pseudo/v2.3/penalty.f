      subroutine penalty(maxdim,pp,penal,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,norbmx,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,wghtp0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :     avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,nconf,
     :     eexcit,wghtcf,wghtex,ntime,itertot)

      implicit real*8 (a-h,o-z)
      logical avgl1,avgl2,avgl3,ortprj,litprj,igrad
      dimension pp(maxdim),no(norb),lo(norb),so(norb),
     :     ev(norbmx,nconf),crcov(norbmx,nconf),dcrcov(norbmx,nconf),
     :     ddcrcov(norbmx,nconf),
     :     occup(noccmx,lmx,nsmx,nconf),aeval(noccmx,lmx,nsmx,nconf),
     :     chrg(noccmx,lmx,nsmx,nconf),dhrg(noccmx,lmx,nsmx,nconf),
     :     ehrg(noccmx,lmx,nsmx,nconf),res(noccmx,lmx,nsmx,nconf),
     :     wght(noccmx,lmx,nsmx,8),
     :     wfnode(noccmx,lmx,nsmx,3,nconf),
     :     gpot(*),r_l(*),hsep(*),
     :     vh(*),xp(*),rmt(*),rmtg(*),ud(*),psi(*),
     :     rr(*),rw(*),rd(*),
     :     eexcit(2,nconf),wghtcf(nconf),wghtex(nconf),psir0(nconf)

c     set res() =-1 for orbitals with zero charge and
c      wght(nocc,l+1,ispin,5) set to zero
c     this avoids unneccessay computations of res() in the
c     routine resid()
      do iconf=1,nconf
         do iorb=1,norb
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            if (so(iorb).lt.0) ispin=2
            if ( (wght(nocc,l+1,ispin,5).eq.0.0d0) .and.
     :           (occup(nocc,l+1,ispin,iconf).lt.1.0d-8)) then
               res(nocc,l+1,ispin,iconf) = -1.d0
            else
               res(nocc,l+1,ispin,iconf) = 0.d0
            endif
         enddo
      enddo

c  unpack variables
      call  ppack (rloc,gpot,hsep,r_l,pp(1),
     :     lpx,lpmx,nspin,nsmx,maxdim,maxdim,'unpack',
     :     avgl1,avgl2,avgl3,ortprj,litprj)

      call gatom(
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
     :     rr,rw,rd,nconf,eexcit,ntime,itertot)
      penal=0

c     loop over all states
      do iconf=1,nconf
         ppenal=0.d0
c     diff for dcharg and echarge is calc. in (%)
         do iorb=1,norb
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            if (so(iorb).lt.0) ispin=2
            ppenal=ppenal + 
     :           ((aeval(nocc,l+1,ispin,iconf)-ev(iorb,iconf))
     :           *wght(nocc,l+1,ispin,1))**2 +
     :           ((chrg(nocc,l+1,ispin,iconf)-crcov(iorb,iconf))
     :           *wght(nocc,l+1,ispin,2))**2 +
     :           (100.d0*(1.d0-dhrg(nocc,l+1,ispin,iconf)
     :           /dcrcov(iorb,iconf))
     :           *wght(nocc,l+1,ispin,3))**2 +
     :           (100.d0*(1.d0-ehrg(nocc,l+1,ispin,iconf)
     :           /ddcrcov(iorb,iconf))
     :           *wght(nocc,l+1,ispin,4))**2 +
     :           (res(nocc,l+1,ispin,iconf)*wght(nocc,l+1,ispin,5))**2 +
     :           (wfnode(nocc,l+1,ispin,1,iconf)
     :           *wght(nocc,l+1,ispin,6))**2 +
     :           (wfnode(nocc,l+1,ispin,2,iconf)
     :           *wght(nocc,l+1,ispin,7))**2 +
     :           (wfnode(nocc,l+1,ispin,3,iconf)
     :           *wght(nocc,l+1,ispin,8))**2 +
     :           (psir0(iconf)*wghtp0)**2
         enddo
         if (ppenal.ne.0.0d0)ppenal = wghtcf(iconf)*sqrt(ppenal)
         penal = penal+ppenal
         penal = penal 
     :        + abs(eexcit(2,iconf)-eexcit(1,iconf))*wghtex(iconf)
      enddo
c      print*,'penal:',penal
      return
      end



