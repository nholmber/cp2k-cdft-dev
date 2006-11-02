      subroutine diag_oneeh(h1mat,smat,cmat)
      use atom
      implicit real*8(a-h,o-z)
      dimension h1mat(namax,namax,0:lamax)
      dimension smat(namax,namax,0:lamax)
      dimension cmat(namax,namax,0:lamax)
      dimension x1mat(namax,namax,0:lamax)
      dimension x2mat(namax,namax,0:lamax)
      dimension ksener(namax,0:lamax)
      x1mat=h1mat
      x2mat=smat
      do l=0,lmax
        n=nalpha(l)
        call diag(h1mat(1,1,l),cmat(1,1,l),smat(1,1,l),ksener(1,l),
     &		  namax,n)
      enddo
      h1mat=x1mat
      smat=x2mat
c
      return
      end
