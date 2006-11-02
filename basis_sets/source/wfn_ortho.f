      subroutine wfn_ortho(cmat)
      use atom
      implicit none
      real*8 cmat(namax,namax,0:lamax)
c..local arrays
      real*8 smat(namax,namax,0:lamax)
      real*8 h2pack(namax*namax,0:lamax)
      real*8 evec(namax,namax,0:lamax)
      real*8 h1mat(namax,namax,0:lamax)
      real*8 h2mat(namax,namax,0:lamax)
      real*8 w(namax,0:lamax),aux(3*namax)
      integer naux,l,i,j,k,info
c
c..Overlap matrix
      call overlap(smat) 
c
      do l=0,lmax
c
c..h1 = CT S
        call dgemm('T','N',nalpha(l),nalpha(l),nalpha(l),1.D0,
     &	  cmat(1,1,l),namax,smat(1,1,l),namax,0.D0,h1mat(1,1,l),namax) 
c..h2 = h1 C = CT S C
        call dgemm('N','N',nalpha(l),nalpha(l),nalpha(l),1.D0,
     &	  h1mat(1,1,l),namax,cmat(1,1,l),namax,0.D0,h2mat(1,1,l),namax)
c..pack CT S C = h2
	do i=1,nalpha(l)
	  do j=i,nalpha(l)
	    h2pack(i+((j*(j-1))/2),l)=h2mat(i,j,l)
	  enddo
	enddo
c..diagonalize h2
        naux=3*namax
        call dspev("V","U",nalpha(l),h2pack(1,l),w(1,l),evec(1,1,l),
     &             namax,aux,info)
        if (info /= 0) STOP "dspev: info /= 0"
c..h2 = E^(-1/2) evec
        do i=1,nalpha(l)
          do j=1,nalpha(l)
	    h2mat(i,j,l)=evec(i,j,l)/sqrt(w(i,l))
	  enddo
        enddo
c..h1 = evecT E^(-1/2) evec = (CT S C)^(-1/2)
        call dgemm('T','N',nalpha(l),nalpha(l),nalpha(l),1.D0,
     &	  evec(1,1,l),namax,h2mat(1,1,l),namax,0.D0,h1mat(1,1,l),namax) 
c..h2 = C h1 = C (CT S C)^(-1/2)
        call dgemm('N','N',nalpha(l),nalpha(l),nalpha(l),1.D0,
     &	  cmat(1,1,l),namax,h1mat(1,1,l),namax,0.D0,h2mat(1,1,l),namax)
c..C = h2 = C (CT S C)^(-1/2)
	cmat(:,:,l)=h2mat(:,:,l)
c
c.. CT S C
c..h1 = CT S
c        call dgemm('T','N',nalpha(l),nalpha(l),nalpha(l),1.D0,
c     &	  cmat(1,1,l),namax,smat(1,1,l),namax,0.D0,h1mat(1,1,l),namax) 
c..h2 = h1 C = CT S C
c        call dgemm('N','N',nalpha(l),nalpha(l),nalpha(l),1.D0,
c     &	  h1mat(1,1,l),namax,cmat(1,1,l),namax,0.D0,h2mat(1,1,l),namax)
c
      enddo
c
      return
      end
