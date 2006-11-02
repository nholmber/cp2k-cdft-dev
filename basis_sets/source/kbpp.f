C     ==================================================================
      SUBROUTINE KBPP(cmat)
C     ==================================================================
c..This routine calculates the product of the pseudowavefunctions and the
c..pseudopotential basis functions to give the coefficients and exponents
c..for the Kleinman-Bylander form of the PP.
c..V_KBlocal = V_BHSlocal + V_BHSnonlocal^lmax 
c         = SUM_i^ERFnr( -Zeff/r* PPLerf(i)*ERF(PPNerf(i)*r) )
c  +SUM_i^EXPnr( PPLexp(i,lmax)*r**(PPRexp(i,lmax)-2)*EXP(-PPNexp(i,lmax)*r**2 )
c
c..the form of the products is
c	KB_coeff(i,n,l) * r**KB_rexp(i,l) * exp(-KB_exp(i,l)*r**2) * Y_lm
c
      use atom
      use pspot
      IMPLICIT none
c..Matrix of coefficients for the wavefunctions
      REAL*8	cmat(namax,namax,0:lamax)
      INTEGER KB_rexp(4*EXPmax*Namax,0:Lamax)
      REAL*8  KB_coeff(4*EXPmax*Namax,Namax,0:Lamax)
      REAL*8  KB_exp(4*EXPmax*Namax,0:Lamax)
      INTEGER KB_nr(0:Lamax)
      INTEGER i,j,kappa,l,n
c     ------------------------------------------------------------------
      do l=0,lmax
         j=0
         do kappa=1,Nalpha(l)
	    do i=1,Expnr(l)
               KB_exp(j,l)=alpha(kappa,l)+PPNexp(i,l)
               KB_exp(j+1,l)=alpha(kappa,l)+PPNexp(i,lmax)
               do n=1,Nocc(l)+3
                  KB_coeff(j,l,n)=cmat(kappa,n,l)*PPLexp(i,l)
                  KB_coeff(j+1,l,n)=-cmat(kappa,n,l)*PPLexp(i,lmax)
               enddo
               KB_rexp(j,l)=PPRexp(i,l)
               KB_rexp(j+1,l)=PPRexp(i,l)
               j=j+2
               KB_nr(l)=j
	    enddo
         enddo
      enddo
c
c..wenn dies nicht soviele Terme erg"abe, w"urde man jetzt noch die Ausgabe-
c..Prozedur anf"ugen
c
      return
      end


