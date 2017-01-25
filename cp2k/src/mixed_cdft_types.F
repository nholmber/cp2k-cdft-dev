!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for mixed CDFT calculations
!> \par   History
!>                 Separated CDFT routines from mixed_environment_types
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
MODULE mixed_cdft_types
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: becke_control_release,&
                                              becke_restraint_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_release,&
                                              pw_env_type
   USE qs_kind_types,                   ONLY: deallocate_qs_kind_set,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

! **************************************************************************************************
!> \brief Buffers for load balancing
! **************************************************************************************************
   TYPE buffers
      INTEGER                                          :: rank(2), tag(2)
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                              :: cavity, weight
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)                            :: gradients
   END TYPE buffers
! **************************************************************************************************
!> \brief To build array of buffers
! **************************************************************************************************
   TYPE p_buffers
      TYPE(buffers), DIMENSION(:), POINTER             :: buffs
   END TYPE p_buffers
! **************************************************************************************************
!> \brief Information about load balancing
! **************************************************************************************************
   TYPE repl_info
      INTEGER, DIMENSION(:), POINTER                   :: matrix_info
      INTEGER, DIMENSION(:, :), POINTER                 :: target_list
   END TYPE repl_info
! **************************************************************************************************
!> \brief Load balancing control for mixed CDFT calculation
! **************************************************************************************************
   TYPE mixed_cdft_dlb_type
      INTEGER                                          :: my_source, distributed(2), &
                                                          my_dest_repl(2), dest_tags_repl(2), &
                                                          more_work
      INTEGER, DIMENSION(:), POINTER                   :: bo, expected_work, &
                                                          prediction_error
      INTEGER, DIMENSION(:, :), POINTER                 :: target_list
      LOGICAL                                          :: recv_work, send_work
      LOGICAL, DIMENSION(:), POINTER                   :: recv_work_repl
      REAL(KIND=dp)                                    :: load_scale, very_overloaded
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                              :: cavity, weight
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :, :)                            :: gradients
      ! Should convert to TYPE(p_buffers), POINTER
      TYPE(buffers), DIMENSION(:), POINTER             :: sendbuff
      TYPE(p_buffers), DIMENSION(:), POINTER           :: recvbuff
      TYPE(repl_info), DIMENSION(:), POINTER           :: recv_info
   END TYPE mixed_cdft_dlb_type
! **************************************************************************************************
!> \brief Main mixed CDFT control type
! **************************************************************************************************
   TYPE mixed_cdft_type
      INTEGER                                          :: sim_step, multiplicity, &
                                                          constraint_type, combined_type
      INTEGER, POINTER, DIMENSION(:)                   :: source_list, dest_list, &
                                                          recv_bo, source_list_save, &
                                                          dest_list_save
      INTEGER, POINTER, DIMENSION(:, :)                 :: source_list_bo, dest_list_bo, &
                                                           source_bo_save, dest_bo_save
      LOGICAL                                          :: is_pencil, dlb, &
                                                          is_special, first_iteration, &
                                                          calculate_metric, &
                                                          wfn_overlap_method, &
                                                          has_unit_metric, &
                                                          is_parallel
      REAL(KIND=dp)                                    :: eps_rho_rspace, sim_dt, &
                                                          eps_svd
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :)         :: weight, cavity
      TYPE(becke_restraint_type), POINTER              :: becke_control
      TYPE(buffers), DIMENSION(:), POINTER             :: sendbuff
      TYPE(cp_blacs_env_type), POINTER                 :: blacs_env
      TYPE(mixed_cdft_dlb_type), POINTER               :: dlb_control
      TYPE(qs_kind_type), DIMENSION(:), &
         POINTER                                       :: qs_kind_set
      TYPE(pw_env_type), POINTER                      :: pw_env
   END TYPE mixed_cdft_type

! **************************************************************************************************
!> \brief Container for constraint settings to check consistency of force_evals
! **************************************************************************************************
   TYPE mixed_cdft_settings_type
      LOGICAL                                            :: is_spherical, &
                                                            is_odd
      INTEGER                                            :: nbecke
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, PARAMETER                                 :: max_nkinds = 30
      INTEGER, DIMENSION(:), POINTER                     :: npts, &
                                                            grid_span, &
                                                            spherical, &
                                                            odd
      INTEGER, DIMENSION(:, :), POINTER                  :: si, &
                                                            rs_dims, &
                                                            atoms
      REAL(KIND=dp)                                      :: radius
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoff, &
                                                            rel_cutoff
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sr, &
                                                            coeffs, &
                                                            confine_bounds, &
                                                            cutoffs, &
                                                            radii
   END TYPE mixed_cdft_settings_type

! *** Public data types ***

   PUBLIC :: mixed_cdft_type, &
             mixed_cdft_settings_type

! *** Public subroutines ***

   PUBLIC :: mixed_cdft_type_create, &
             mixed_cdft_type_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_cdft_types'

CONTAINS

! **************************************************************************************************
!> \brief inits the given mixed_cdft_type
!> \param cdft_control the object to init
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_type_create(cdft_control)
      TYPE(mixed_cdft_type), POINTER                     :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_type_create', &
         routineP = moduleN//':'//routineN

      NULLIFY (cdft_control%pw_env, cdft_control%blacs_env, cdft_control%qs_kind_set)
      NULLIFY (cdft_control%dlb_control, cdft_control%dest_list_bo, cdft_control%dest_list)
      NULLIFY (cdft_control%dest_bo_save, cdft_control%dest_list_save, cdft_control%source_list)
      NULLIFY (cdft_control%source_list_save, cdft_control%source_bo_save, cdft_control%source_list_bo)
      NULLIFY (cdft_control%cavity, cdft_control%weight, cdft_control%sendbuff)
      NULLIFY (cdft_control%becke_control, cdft_control%recv_bo)

   END SUBROUTINE mixed_cdft_type_create

! **************************************************************************************************
!> \brief releases the given mixed_cdft_type
!> \param cdft_control the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_type_release(cdft_control)
      TYPE(mixed_cdft_type), POINTER                     :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_type_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      CALL pw_env_release(cdft_control%pw_env)
      IF (ASSOCIATED(cdft_control%dest_list)) &
         DEALLOCATE (cdft_control%dest_list)
      IF (ASSOCIATED(cdft_control%dest_list_save)) &
         DEALLOCATE (cdft_control%dest_list_save)
      IF (ASSOCIATED(cdft_control%dest_list_bo)) &
         DEALLOCATE (cdft_control%dest_list_bo)
      IF (ASSOCIATED(cdft_control%dest_bo_save)) &
         DEALLOCATE (cdft_control%dest_bo_save)
      IF (ASSOCIATED(cdft_control%source_list)) &
         DEALLOCATE (cdft_control%source_list)
      IF (ASSOCIATED(cdft_control%source_list_save)) &
         DEALLOCATE (cdft_control%source_list_save)
      IF (ASSOCIATED(cdft_control%source_list_bo)) &
         DEALLOCATE (cdft_control%source_list_bo)
      IF (ASSOCIATED(cdft_control%source_bo_save)) &
         DEALLOCATE (cdft_control%source_bo_save)
      IF (ASSOCIATED(cdft_control%recv_bo)) &
         DEALLOCATE (cdft_control%recv_bo)
      IF (ASSOCIATED(cdft_control%weight)) &
         DEALLOCATE (cdft_control%weight)
      IF (ASSOCIATED(cdft_control%cavity)) &
         DEALLOCATE (cdft_control%cavity)
      IF (ASSOCIATED(cdft_control%dlb_control)) &
         CALL mixed_cdft_dlb_release(cdft_control%dlb_control)
      IF (ASSOCIATED(cdft_control%sendbuff)) THEN
         DO i = 1, SIZE(cdft_control%sendbuff)
            CALL mixed_cdft_buffers_release(cdft_control%sendbuff(i))
         END DO
         DEALLOCATE (cdft_control%sendbuff)
      END IF
      CALL becke_control_release(cdft_control%becke_control)
      IF (ASSOCIATED(cdft_control%blacs_env)) &
         CALL cp_blacs_env_release(cdft_control%blacs_env)
      IF (ASSOCIATED(cdft_control%qs_kind_set)) &
         CALL deallocate_qs_kind_set(cdft_control%qs_kind_set)
      DEALLOCATE (cdft_control)

   END SUBROUTINE mixed_cdft_type_release

! **************************************************************************************************
!> \brief releases the given load balancing control
!> \param dlb_control the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_dlb_release(dlb_control)
      TYPE(mixed_cdft_dlb_type), POINTER                 :: dlb_control

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_dlb_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(dlb_control%recv_work_repl)) &
         DEALLOCATE (dlb_control%recv_work_repl)
      IF (ASSOCIATED(dlb_control%sendbuff)) THEN
         DO i = 1, SIZE(dlb_control%sendbuff)
            CALL mixed_cdft_buffers_release(dlb_control%sendbuff(i))
         END DO
         DEALLOCATE (dlb_control%sendbuff)
      END IF
      IF (ASSOCIATED(dlb_control%recvbuff)) THEN
         DO i = 1, SIZE(dlb_control%recvbuff)
            CALL mixed_cdft_p_buffers_release(dlb_control%recvbuff(i))
         END DO
         DEALLOCATE (dlb_control%recvbuff)
      END IF
      IF (ASSOCIATED(dlb_control%recv_info)) THEN
         DO i = 1, SIZE(dlb_control%recv_info)
            IF (ASSOCIATED(dlb_control%recv_info(i)%matrix_info)) &
               DEALLOCATE (dlb_control%recv_info(i)%matrix_info)
            IF (ASSOCIATED(dlb_control%recv_info(i)%target_list)) &
               DEALLOCATE (dlb_control%recv_info(i)%target_list)
         END DO
         DEALLOCATE (dlb_control%recv_info)
      END IF
      IF (ASSOCIATED(dlb_control%bo)) &
         DEALLOCATE (dlb_control%bo)
      IF (ASSOCIATED(dlb_control%expected_work)) &
         DEALLOCATE (dlb_control%expected_work)
      IF (ASSOCIATED(dlb_control%prediction_error)) &
         DEALLOCATE (dlb_control%prediction_error)
      IF (ASSOCIATED(dlb_control%target_list)) &
         DEALLOCATE (dlb_control%target_list)
      IF (ASSOCIATED(dlb_control%cavity)) &
         DEALLOCATE (dlb_control%cavity)
      IF (ASSOCIATED(dlb_control%weight)) &
         DEALLOCATE (dlb_control%weight)
      IF (ASSOCIATED(dlb_control%gradients)) &
         DEALLOCATE (dlb_control%gradients)
      DEALLOCATE (dlb_control)

   END SUBROUTINE mixed_cdft_dlb_release

! **************************************************************************************************
!> \brief releases the given buffers
!> \param buffer the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_buffers_release(buffer)
      TYPE(buffers)                                      :: buffer

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_buffers_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(buffer%cavity)) &
         DEALLOCATE (buffer%cavity)
      IF (ASSOCIATED(buffer%weight)) &
         DEALLOCATE (buffer%weight)
      IF (ASSOCIATED(buffer%gradients)) &
         DEALLOCATE (buffer%gradients)

   END SUBROUTINE mixed_cdft_buffers_release

! **************************************************************************************************
!> \brief releases the given pointer of buffers
!> \param p_buffer the object to release
!> \author Nico Holmberg [01.2017]
! **************************************************************************************************
   SUBROUTINE mixed_cdft_p_buffers_release(p_buffer)
      TYPE(p_buffers)                                    :: p_buffer

      CHARACTER(len=*), PARAMETER :: routineN = 'mixed_cdft_p_buffers_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      IF (ASSOCIATED(p_buffer%buffs)) THEN
         DO i = 1, SIZE(p_buffer%buffs)
            CALL mixed_cdft_buffers_release(p_buffer%buffs(i))
         END DO
         DEALLOCATE (p_buffer%buffs)
      END IF

   END SUBROUTINE mixed_cdft_p_buffers_release

END MODULE mixed_cdft_types
