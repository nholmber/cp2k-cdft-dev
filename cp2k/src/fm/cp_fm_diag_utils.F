!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary tools to redistribute cp_fm_type matrices before and after diagonalization.
!>        Heuristics are used to determine the optimal number of CPUs for diagonalization and the
!>        input matrices are redistributed if necessary
!> \par History
!>      - [01.2018] moved redistribution related code from cp_fm_syevd here
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
MODULE cp_fm_diag
   USE cp_blacs_calls,                  ONLY: cp_blacs_gridexit,&
                                              cp_blacs_gridinit
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_gemm,&
                                              cp_fm_scale,&
                                              cp_fm_syrk,&
                                              cp_fm_triangular_invert,&
                                              cp_fm_triangular_multiply,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_elpa,                      ONLY: cp_fm_diag_elpa,&
                                              set_elpa_kernel,&
                                              set_elpa_qr,&
                                              set_elpa_print
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_get_unit_nr,&
                                              cp_logger_type
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_memory
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_comm_free,&
                                              mp_comm_split,&
                                              mp_sync
#if defined (__HAS_IEEE_EXCEPTIONS)
  USE ieee_exceptions,                 ONLY: ieee_get_halting_mode,&
                                             ieee_set_halting_mode,&
                                             ieee_all
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag_utils'

   ! Container for redistribution settings and temporary work structs
   TYPE cp_fm_redistribute_type
      INTEGER                                  :: subgroup
      INTEGER, DIMENSION(:), POINTER           :: group_distribution => null()
      INTEGER, DIMENSION(:), POINTER           :: group_partition => null()
      TYPE(cp_blacs_env_type), POINTER         :: blacs_env_new => null()
      TYPE(cp_para_env_type), POINTER          :: para_env_new => null()
   END TYPE cp_fm_redistribute_type

   ! Permanent instance of redistribute type
   TYPE(cp_fm_redistribute_type), PRIVATE, &
      SAVE                                     :: work_redistribute

   ! Public subroutines

   PUBLIC :: cp_fm_redistribute_start, &
             cp_fm_redistribute_end

CONTAINS

! **************************************************************************************************
!> \brief  Initializes temporary storage needed when redistributing arrays
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_work_init()

      work_redistribute%subgroup = -1
      NULLIFY(work_redistribute%group_distribution)
      NULLIFY(work_redistribute%group_partition)
      NULLIFY(work_redistribute%blacs_env_new)
      NULLIFY(work_redistribute%para_env_new)

   END SUBROUTINE cp_fm_redistribute_work_init

! **************************************************************************************************
!> \brief  Releases the temporary storage needed when redistributing arrays
!> \param  has_redistributed flag that determines if the processors holds a part of the
!>                           redistributed array
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_work_finalize(has_redistributed)
      LOGICAL                                  :: has_work

      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         IF (has_work) THEN
            CALL cp_blacs_env_release(work_redistribute%blacs_env_new)
            CALL cp_para_env_release(work_redistribute%para_env_new)
         ELSE
            CALL mp_comm_free(work_redistribute%subgroup)
         ENDIF
         DEALLOCATE(work_redistribute%group_distribution)
         DEALLOCATE(work_redistribute%group_partition)
      END IF
      ! Return work to its initial state
      CALL cp_fm_redistribute_work_init()

   END SUBROUTINE cp_fm_redistribute_work_finalize

! **************************************************************************************************
!> \brief   Determines the optimal number of CPUs for matrix diagonalization and redistributes
!>          the input matrices if necessary
!> \param matrix           the input cp_fm_type matrix to be diagonalized
!> \param eigenvectors     the cp_fm_type matrix that will hold the eigenvectors of the input matrix
!> \param matrix_new       the redistributed input matrix which will subsequently be diagonalized,
!>                         or a pointer to the original matrix if no redistribution is required
!> \param eigenvector_new  the redistributed eigenvectors matrix, or a pointer to the original
!>                         matrix if no redistribution is required
!> \par History
!>      - [01.2018] created by moving redistribution related code from cp_fm_syevd here
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new)

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      TYPE(cp_fm_type), POINTER, INTENT(OUT)   :: matrix_new, eigenvectors_new

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_start', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle
#if defined(__SCALAPACK)
      REAL(KIND=dp)                            :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, &
                                                  n, ngroups, num_pe_new, &
                                                  num_pe_old
      TYPE(cp_para_env_type), POINTER          :: para_env
#endif

      CALL timeset(routineN, handle)

#if defined(__SCALAPACK)

      NULLIFY(matrix_new)
      NULLIFY(eigenvectors_new)
      n = matrix%matrix_struct%nrow_global

      ! first figure out the optimal number of cpus
      ! this is pure heuristics, based on rosa timings
      ! that demonstrate that timings go up sharply if too many tasks are used
      ! we take a multiple of 4, and approximately n/60
      num_pe_old = matrix%matrix_struct%para_env%num_pe
      num_pe_new = ((n+60*4-1)/(60*4))*4
      para_env => matrix%matrix_struct%para_env
      mepos_old = para_env%mepos

      ! if the optimal is smaller than num_pe, we will redistribute the input matrix
      IF (num_pe_new < num_pe_old) THEN

         ! split comm, the first num_pe_new tasks will do the work
         ALLOCATE (work_redistribute%group_distribution(0:num_pe_old-1))
         ALLOCATE (work_redistribute%group_partition(0:1))
         work_redistribute%group_partition = (/num_pe_new, num_pe_old-num_pe_new/)
         CALL mp_comm_split(comm=para_env%group, sub_comm=work_redistribute%subgroup, &
                            ngroups=ngroups, group_distribution=work_redistribute%group_distribution, &
                            n_subgroups=2, group_partition=work_redistribute%group_partition)

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN

            ! create para_env, might need a proper bound to this para_env
            NULLIFY (work_redistribute%para_env_new)
            CALL cp_para_env_create(work_redistribute%para_env_new, work_redistribute%subgroup)
            ! test a sync
            CALL mp_sync(work_redistribute%para_env_new%group)

            ! create blacs, should inherit the preferences for the layout and so on, from the higher level
            NULLIFY (work_redistribute%blacs_env_new)
            CALL cp_blacs_env_create(blacs_env=work_redistribute%blacs_env_new, para_env=work_redistribute%para_env_new)

            ! create new matrix
            NULLIFY (fm_struct_new)
            CALL cp_fm_struct_create(fmstruct=fm_struct_new, &
                                     para_env=work_redistribute%para_env_new, &
                                     context=work_redistribute%blacs_env_new, &
                                     nrow_global=n, ncol_global=n)
            CALL cp_fm_create(matrix_new, matrix_struct=fm_struct_new, name="yevd_new_mat")
            CALL cp_fm_create(eigenvectors_new, matrix_struct=fm_struct_new, name="yevd_new_vec")
            CALL cp_fm_struct_release(fm_struct_new)

            ! redistribute old
            CALL pdgemr2d(n, n, matrix%local_data(1, 1), 1, 1, matrix%matrix_struct%descriptor, &
                          matrix_new%local_data(1, 1), 1, 1, matrix_new%matrix_struct%descriptor, &
                          matrix%matrix_struct%context%group)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pdgemr2d
            fake_descriptor = -1
            CALL pdgemr2d(n, n, matrix%local_data(1, 1), 1, 1, matrix%matrix_struct%descriptor, &
                          fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          matrix%matrix_struct%context%group)
         ENDIF
      ELSE
         ! No need to redistribute, just return pointers to the original arrays
         matrix_new => matrix
         eigenvectors_new => eigenvectors

      ENDIF

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(matrix_new)
      MARK_USED(eigenvectors_new)
      CPABORT("Routine called in non-parallel case.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_redistribute_start

! **************************************************************************************************
!> \brief Redistributes eigenvectors and eigenvalues  back to the original communicator group
!> \param matrix           the input cp_fm_type matrix to be diagonalized
!> \param eigenvectors     the cp_fm_type matrix that will hold the eigenvectors of the input matrix
!> \param eig              global array holding the eigenvalues of the input matrixmatrix
!> \param matrix_new       the redistributed input matrix which will subsequently be diagonalized,
!>                         or a pointer to the original matrix if no redistribution is required
!> \param eigenvector_new  the redistributed eigenvectors matrix, or a pointer to the original
!>                         matrix if no redistribution is required
!> \par History
!>      - [01.2018] created by moving redistribution related code from cp_fm_syevd here
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_end(matrix, eigenvectors, eig, matrix_new, eigenvectors_new)

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eig
      TYPE(cp_fm_type), POINTER, INTENT(OUT)   :: matrix_new, eigenvectors_new

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_end', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle
#if defined(__SCALAPACK)
      REAL(KIND=dp)                            :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, n
      TYPE(cp_para_env_type), POINTER          :: para_env
#endif

      CALL timeset(routineN, handle)

#if defined(__SCALAPACK)

      ! Check if matrix was redistributed
      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         n = matrix%matrix_struct%nrow_global
         para_env => matrix%matrix_struct%para_env
         mepos_old = para_env%mepos

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN
            ! redistribute results on CPUs that hold the redistributed matrix
            CALL pdgemr2d(n, n, eigenvectors_new%local_data(1, 1), 1, 1, eigenvectors_new%matrix_struct%descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context%group)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pdgemr2d
            fake_descriptor = -1
            CALL pdgemr2d(n, n, fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context%group)
         ENDIF
         ! free work
         CALL cp_fm_redistribute_work_finalize(work_redistribute%group_distribution(mepos_old) == 0)

         ! finally, also the eigenvalues need to end up on the non-group member tasks
         CALL mp_bcast(eig, 0, para_env%group)

      ELSE
         ! Just nullify pointers to the original matrices
         NULLIFY (matrix_new)
         NULLIFY (eigenvectors_new)

      ENDIF

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eig)
      MARK_USED(matrix_new)
      MARK_USED(eigenvectors_new)
      CPABORT("Routine called in non-parallel case.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_redistribute_end

END MODULE cp_fm_diag_utils
