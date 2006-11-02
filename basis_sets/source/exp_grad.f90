      SUBROUTINE exp_grad (cmat,ksener)
      USE basic_data_types, ONLY: dp
      USE atom
      USE pspot
      IMPLICIT NONE
      REAL(dp) cmat(namax,namax,0:lamax)
      REAL(dp) pmat(namax,namax,0:lamax)
      REAL(dp) nn(namax,namax,0:lamax)
      REAL(dp) smat(namax,namax,0:lamax)
      REAL(dp) umat(namax,namax,0:lamax)
      REAL(dp) vxcmat(namax,namax,0:lamax)
      REAL(dp) excmat(namax,namax,0:lamax)
      REAL(dp) h2mat(namax,namax,0:lamax)
      REAL(dp) ksener(namax,0:lamax)
      REAL(dp) sdmat(namax,namax,0:lamax)
      REAL(dp) tdmat(namax,namax,0:lamax)
      REAL(dp) tmat(namax,namax,0:lamax)
      REAL(dp) udmat(namax,namax,0:lamax)
      REAL(dp) uadddmat(namax,namax,0:lamax)
      REAL(dp) h2dmat(namax,namax,0:lamax)
      REAL(dp) vxcdmat(namax,namax,0:lamax)
      REAL(dp) g2,f1,f2
      INTEGER k,l,i,j

!..Normalization factors
      call calcnn(nn)
!..Density matrix
      call denmat(pmat,cmat)
!..overlap matrix
      call overlap(smat)
!..differentiated overlap matrix
      call overlapd(smat,sdmat)
!..Kinetic energy matrix
      call kinetic_int(tmat,smat)
!..differentiated kinetic energy matrix
      call kinetic_intd(smat,sdmat,tdmat)
!..Potential energy matrix
      if (allelectron) then
        call nuc_att(umat,smat)
        do k=0,lmax
          do i=1,nalpha(k)
            do j=1,nalpha(k)
                umat(i,j,k)=-zval*umat(i,j,k)
             enddo
          enddo
        enddo
      else
        call pseudopot(umat,nn)
      endif
!..differentiated potential energy matrix
      if (allelectron) then
        call nuc_attd(smat,sdmat,udmat)
      else
        call pseudopotd(udmat,nn)
      endif
!..additional parabolic potential derivative
      call v_add_d(smat,sdmat,alpha,nalpha,lmax,lamax,namax,uadddmat)
!..two electron energy matrices
      call hermite_mat(smat,pmat,nn,h2mat,excmat,vxcmat)
!..r^2 multiplied coulomb and exchange-correlation energy matrices
      call hermite_matd(smat,pmat,nn,h2dmat,vxcdmat)

!     ------------------------------------------------------------------
!..calculate gradient
      do l=0,lmax

        f1=2.D0*(2*l+1)

        do j=1,nalpha(l)

          f2=(2*l+3)/(4.D0*alpha(j,l))

          galpha(j,l)=0.D0
          g_constr(j,l)=0.D0
          g_ekin(j,l)=0.D0
          g_epot(j,l)=0.D0
          g_epotadd(j,l)=0.D0
          g_exc(j,l)=0.D0
          g_ecoul(j,l)=0.D0

          do k=1,nalpha(l)

!.........part of the constraint:
            g2=0.D0
            do i=1,nocc(l)
              g2=g2-occ(i,l)*ksener(i,l)*cmat(j,i,l)*cmat(k,i,l)
            enddo
            g_constr(j,l)=g_constr(j,l)+f1*g2*sdmat(k,j,l)

!.........part of kinetic energy:
            g_ekin(j,l)=g_ekin(j,l)+f1*pmat(k,j,l)*tdmat(k,j,l)

!.........part of additional potential energy:
            g_epotadd(j,l)=g_epotadd(j,l)+f1*pmat(k,j,l)*add_pot_factor*uadddmat(k,j,l)

!.........part of potential energy:
            if (allelectron) then
               g_epot(j,l)=g_epot(j,l)-f1*Zval*pmat(k,j,l)*udmat(k,j,l)
            else
               g_epot(j,l)=g_epot(j,l)+f1*pmat(k,j,l)*(f2*umat(k,j,l)-udmat(k,j,l))
            endif

!.........part of coulomb energy:
            g_ecoul(j,l)=g_ecoul(j,l)+f1*pmat(k,j,l)*(f2*h2mat(k,j,l)-h2dmat(k,j,l))

!.........part of XC energy:
            g_exc(j,l)=g_exc(j,l)+f1*pmat(k,j,l)*(f2*vxcmat(k,j,l)-vxcdmat(k,j,l))

!.........total energy gradient:
            galpha(j,l)=g_constr(j,l)+g_ekin(j,l)+g_epot(j,l)+g_ecoul(j,l)+g_exc(j,l)+g_epotadd(j,l)

          enddo
        enddo
      enddo
      end
