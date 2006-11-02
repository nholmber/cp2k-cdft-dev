      subroutine ks_matrix(ksmat,smat,pmat,h1mat,h2mat,vxcmat,excmat,
     *                     cnn)
      use atom
      use rint
      implicit none
      real*8 ksmat(namax,namax,0:lamax)
      real*8 smat(namax,namax,0:lamax)
      real*8 pmat(namax,namax,0:lamax)
      real*8 h1mat(namax,namax,0:lamax)
      real*8 h2mat(namax,namax,0:lamax)
      real*8 vxcmat(namax,namax,0:lamax)
      real*8 excmat(namax,namax,0:lamax)
      real*8 cnn(namax,namax,0:lamax)
      real*8 ksener(namax,0:lamax)
      integer i,j,k
c
      call hermite_mat(smat,pmat,cnn,h2mat,excmat,vxcmat)
c
c..Kohn-Sham-Matrix
      do k=0,lmax
        do i=1,nalpha(k)
          do j=1,nalpha(k)
            ksmat(i,j,k)=h1mat(i,j,k)+h2mat(i,j,k)+vxcmat(i,j,k)
          enddo
        enddo
      enddo
c
      return
      end
