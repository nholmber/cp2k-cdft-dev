



      subroutine penalty(maxdim,pp,penal,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,whgtp0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :     avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
     :     ntime,itertot)

      implicit real*8 (a-h,o-z)
      logical avgl1,avgl2,avgl3,ortprj,litprj,igrad
      dimension pp(maxdim),no(norb),lo(norb),so(norb),
     :     ev(norb),crcov(norb),dcrcov(norb),ddcrcov(norb),
     :     occup(noccmx,lmx,nsmx),aeval(noccmx,lmx,nsmx),
     :     chrg(noccmx,lmx,nsmx),dhrg(noccmx,lmx,nsmx),
     :     ehrg(noccmx,lmx,nsmx),res(noccmx,lmx,nsmx),
     :     wght(noccmx,lmx,nsmx,8),
     :     wfnode(noccmx,lmx,nsmx,3),
     :     gpot(*),r_l(*),hsep(*),
     :     vh(*),xp(*),rmt(*),rmtg(*),ud(*),psi(*),
     :     rr(*),rw(*),rd(*)

c     set res() =-1 for orbitals with zero charge and
c      wght(nocc,l+1,ispin,5) set to zero
c     this avoids unneccessay computations of res() in the
c     routine resid()
      do iorb=1,norb
         nocc=no(iorb)
         l=lo(iorb)
         ispin=1
         if (so(iorb).lt.0) ispin=2
         if ( (wght(nocc,l+1,ispin,5).eq.0.0d0) .and.
     :        (occup(nocc,l+1,ispin).lt.1.0d-8)) then
            res(nocc,l+1,ispin) = -1.d0
         else
            res(nocc,l+1,ispin) = 0.d0
         endif
      enddo
c  unpack variables
c      print*,'penalty: maxdim=',maxdim
      call  ppack (rloc,gpot,hsep,r_l,pp(1),
     :     lpx,lpmx,nspin,nsmx,maxdim,maxdim,'unpack',
     :     avgl1,avgl2,avgl3,ortprj,litprj)

      call gatom(
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
     :     rcov,rprb,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
     :     rr,rw,rd,ntime,itertot)
      penal=0
c     diff for dcharg and echarge is calc. in (%)
      do iorb=1,norb
         nocc=no(iorb)
         l=lo(iorb)
         ispin=1
         if (so(iorb).lt.0) ispin=2
         penal=penal + 
     :        ((aeval(nocc,l+1,ispin)-ev(iorb))
     :        *wght(nocc,l+1,ispin,1))**2 +
     :        ((chrg(nocc,l+1,ispin)-crcov(iorb))
     :        *wght(nocc,l+1,ispin,2))**2 +
     :        (100.d0*(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))
     :        *wght(nocc,l+1,ispin,3))**2 +
     :        (100.d0*(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))
     :        *wght(nocc,l+1,ispin,4))**2 +
     :        (res(nocc,l+1,ispin)*wght(nocc,l+1,ispin,5))**2 +
     :        (wfnode(nocc,l+1,ispin,1)*wght(nocc,l+1,ispin,6))**2 +
     :        (wfnode(nocc,l+1,ispin,2)*wght(nocc,l+1,ispin,7))**2 +
     :        (wfnode(nocc,l+1,ispin,3)*wght(nocc,l+1,ispin,8))**2 +
     :        (psir0*whgtp0)**2
      enddo
      penal=sqrt(penal)
c      print*,'penal:',penal
      return
      end
