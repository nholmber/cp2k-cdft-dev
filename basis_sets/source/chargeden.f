      subroutine charge_density(P,cnn,nr,r,rho,drho,ddrho)
      use atom
      implicit real*8(a-h,o-z)
      integer nr
      dimension P(namax,namax,0:lamax),cnn(namax,namax,0:lamax)
      dimension r(nr),rho(nr),drho(nr),ddrho(nr)
      parameter (PI=3.14159265358979323846264D0)
      do k=1,nr
	rho(k)=0
	drho(k)=0
	ddrho(k)=0
        do l=0,lmax
          do kap=1,nalpha(l)
            do lam=1,nalpha(l)
              h1=alpha(kap,l)
              h2=alpha(lam,l)
              h=(2.D0*l+1.D0)/(4.D0*PI)*exp(-(h1+h2)*r(k)*r(k))
     &		*cnn(kap,lam,l)*P(kap,lam,l)
C
              rho(k)=rho(k)+h*r(k)**(2*l)
C
              drho(k)=drho(k)+h*2.D0
     &		   *(l*r(k)**(2*l-1)-(h1+h2)*r(k)**(2*l+1))
C
              ddrho(k)=ddrho(k)+h*( 2*l*(2*l-1)*r(k)**(2*l-2)
     &			   -2.D0*(h1+h2)*(4*l+1)*r(k)**(2*l)
     &			   +4.D0*(h1+h2)**2*r(k)**(2*l+2) )
            enddo
          enddo
        enddo
      enddo
c
      return
      end
