!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines for performing a constrained DFT SCF run with Quickstep.
!> \par History
!>      - Separated some routines from qs_scf (03.2018) [Nico Holmberg]
!> \author Nico Holmberg (03.2018)
! **************************************************************************************************
MODULE qs_scf
   USE cp_log_handling,                 ONLY: cp_to_string
   USE input_constants,                 ONLY: jacobian_fd1,&
                                              jacobian_fd1_backward,&
                                              jacobian_fd1_central,&
                                              jacobian_fd2,&
                                              jacobian_fd2_backward
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invert_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf'

   PUBLIC :: prepare_jacobian_stencil, build_diagonal_jacobian

CONTAINS

! **************************************************************************************************
!> \brief Prepares the finite difference stencil for computing the Jacobian. The constraints
!>        are re-evaluated by perturbing each constraint.
!> \param jacobian_type the type of stencil to use
!> \param coeff list of coefficients that determine how to sum up the various perturbations
!> \param step_multiplier list of values that determine how large steps to take for each perturbatio
!> \param dh finite difference step size
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE prepare_jacobian_stencil(jacobian_type, coeff, step_multiplier, dh)
      INTEGER                                            :: jacobian_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff, step_multiplier
      REAL(KIND=dp)                                      :: dh

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_jacobian_stencil', &
         routineP = moduleN//':'//routineN

      SELECT CASE (scf_control%outer_scf%jacobian_type)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown Jacobian type: "// &
                       cp_to_string(scf_control%outer_scf%jacobian_type))
      CASE (jacobian_fd1)
         ! f'(x0) = [ f(x0+h) - f(x0) ] / h
         nwork = 0
         pwork = 1
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(nwork) = -1.0_dp
         coeff(pwork) = 1.0_dp
         step_multiplier = 1.0_dp
         dh = scf_control%outer_scf%jacobian_step
      CASE (jacobian_fd1_backward)
         ! f'(x0) = [ f(x0) - f(x0-h) ] / h
         nwork = -1
         pwork = 0
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(nwork) = -1.0_dp
         coeff(pwork) = 1.0_dp
         step_multiplier = -1.0_dp
         dh = scf_control%outer_scf%jacobian_step
      CASE (jacobian_fd2)
         ! f'(x0) = [ -f(x0+2h) + 4f(x0+h) - 3f(x0) ] / 2h
         nwork = 0
         pwork = 2
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(0) = -3.0_dp
         coeff(1) = 4.0_dp
         coeff(2) = -1.0_dp
         step_multiplier(0) = 0.0_dp
         step_multiplier(1) = 1.0_dp
         step_multiplier(2) = 2.0_dp
         dh = 2.0_dp*scf_control%outer_scf%jacobian_step
      CASE (jacobian_fd2_backward)
         ! f'(x0) = [ 3f(x0) - 4f(x0-h) + f(x0-2h) ] / 2h
         nwork = -2
         pwork = 0
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(0) = 3.0_dp
         coeff(-1) = -4.0_dp
         coeff(-2) = 1.0_dp
         step_multiplier(0) = 0.0_dp
         step_multiplier(-1) = -1.0_dp
         step_multiplier(-2) = -2.0_dp
         dh = 2.0_dp*scf_control%outer_scf%jacobian_step
      CASE (jacobian_fd1_central)
         ! f'(x0) = [ f(x0+h) - f(x0-h) ] / 2h
         nwork = -1
         pwork = 1
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(0) = 0.0_dp
         coeff(nwork) = -1.0_dp
         coeff(pwork) = 1.0_dp
         step_multiplier(0) = 0.0_dp
         step_multiplier(nwork) = -1.0_dp
         step_multiplier(pwork) = 1.0_dp
         dh = 2.0_dp*scf_control%outer_scf%jacobian_step
      END SELECT

   END SUBROUTINE prepare_jacobian_stencil
! **************************************************************************************************
!> \brief Builds a strictly diagonal inverse Jacobian from MD/SCF history.
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \param used_history flag that determines if history was actually used to prepare the Jacobian
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE build_diagonal_jacobian(qs_env, used_history)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: used_history

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_diagonal_jacobian', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ihistory, nvar, outer_scf_ihistory
      LOGICAL                                            :: use_md_history
      REAL(KIND=dp)                                      :: inv_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: jacobian
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient_history, inv_jacobian, &
                                                            variable_history
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (scf_control, scf_env)

      CALL get_qs_env(qs_env, scf_env=scf_env, &
                      scf_control=scf_control, &
                      outer_scf_ihistory=outer_scf_ihistory)
      ihistory = scf_env%outer_scf%iter_count

      IF (outer_scf_ihistory .GE. 3 .AND. .NOT. used_history) THEN
         ! First, lets try using the history from previous energy evaluations
         CALL get_qs_env(qs_env, gradient_history=gradient_history, &
                         variable_history=variable_history)
         nvar = SIZE(scf_env%outer_scf%variables, 1)
         use_md_history = .TRUE.
         ! Check that none of the history values are identical in which case we should try something different
         DO i = 1, nvar
            IF (ABS(variable_history(i, 2)-variable_history(i, 1)) .LT. 1.0E-12_dp) &
               use_md_history = .FALSE.
         END DO
         IF (use_md_history) THEN
            ALLOCATE (jacobian(nvar, nvar))
            DO i = 1, nvar
               jacobian(i, i) = (gradient_history(i, 2)-gradient_history(i, 1))/ &
                                (variable_history(i, 2)-variable_history(i, 1))
            END DO
            IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
               ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
            inv_jacobian => scf_env%outer_scf%inv_jacobian
            CALL invert_matrix(jacobian, inv_jacobian, inv_error)
            DEALLOCATE (jacobian)
            ! Mark that an inverse Jacobian was just built and the next outer_loop_optimize should not perform
            ! a Broyden update of it
            scf_control%outer_scf%broyden_update = .FALSE.
            ! Mark that the MD history has been used and should not be reused anymore on this energy evaluation
            used_history = .TRUE.
         END IF
      END IF
      IF (ihistory .GE. 2 .AND. .NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
         ! Next, try history from current SCF procedure
         nvar = SIZE(scf_env%outer_scf%variables, 1)
         IF (SIZE(scf_env%outer_scf%gradient, 2) .LT. 3) &
            CALL cp_abort(__LOCATION__, &
                          "Keyword EXTRAPOLATION_ORDER in section OUTER_SCF must be greater than or equal "// &
                          "to 3 for optimizers that build the Jacobian from SCF history.")
         ALLOCATE (jacobian(nvar, nvar))
         DO i = 1, nvar
            jacobian(i, i) = (scf_env%outer_scf%gradient(i, ihistory)-scf_env%outer_scf%gradient(i, ihistory-1))/ &
                             (scf_env%outer_scf%variables(i, ihistory)-scf_env%outer_scf%variables(i, ihistory-1))
         END DO
         IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
            ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
         inv_jacobian => scf_env%outer_scf%inv_jacobian
         CALL invert_matrix(jacobian, inv_jacobian, inv_error)
         DEALLOCATE (jacobian)
         scf_control%outer_scf%broyden_update = .FALSE.
      ELSE
         ! No history => will fall back to SD optimizer in outer_loop_optimize
      END IF

   END SUBROUTINE build_diagonal_jacobian

END MODULE qs_cdft_scf_utils
