!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutines for CDFT calculations
!> \par   History
!>                 separated from et_coupling [03.2017]
!> \author Nico Holmberg [03.2017]
! **************************************************************************************************
MODULE qs_cdft_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: cdft_control_type,&
                                              dft_control_type, &
                                              becke_group_type, &
                                              becke_constraint_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_transpose
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_get_el_struct
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_deallocate_matrix_set,&
                                              dbcsr_p_type
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: get_hirshfeld_info,&
                                              hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: &
        becke_cutoff_element, becke_cutoff_global, cdft_combined_acceptor, cdft_combined_all, &
        cdft_combined_constraint, cdft_combined_donor, cdft_density_constraint, &
        cdft_magnetization_constraint, ref_charge_atomic, ref_charge_mulliken, cdft_atomic_magnetization_constraint
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_sum,&
                                              accurate_sum_arrays_product
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_set
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_modify_pab_block,             ONLY: FUNC_AB
   USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                              mpole_rho_atom,&
                                              rho0_mpole_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE realspace_grid_cube,             ONLY: cube_to_pw
   USE realspace_grid_types,            ONLY: &
        pw2rs, realspace_grid_desc_type, realspace_grid_type, rs2pw, rs_grid_create, &
        rs_grid_release, rs_grid_retain, rs_grid_zero, rs_pw_transfer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_utils'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

! *** Public subroutines ***
   PUBLIC :: becke_constraint_init, hfun_scale

! **************************************************************************************************
!> \brief calculates a Becke constraint
!> \param qs_env the qs_env where to build the constraint
!> \param calc_pot if the potential needs to be recalculated or just integrated
!> \param calculate_forces logical if potential has to be calculated or only_energy
!> \par   History
!>        Created 01.2007 [fschiff]
!>        Extended functionality 12/15-12/16 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_constraint_init(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: calc_pot, calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER :: atom_a, bounds(2), dir, handle, i, iatom, iex, ikind, ind(3), ip, ithread, ivar, &
         iw, j, jatom, k, katom, lb_index, natom, nelectron_total, nkind, np(3), npme, &
         nskipped, nthread, numexp, nvar, tmp_index(2), ub_index, unit_nr, igroup, jp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: catom
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, stride
      LOGICAL                                            :: build, in_memory, paw_atom, store_vectors
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_constraint, skip_me
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: atom_in_group
      REAL(kind=dp) :: alpha, chi, coef, dist1, dist2, dmyexp, dvol, dynamic_radius, eps_cavity, &
         ircov, jrcov, lb, my1, my1_homo, myexp, nelectron_frag, sign, sum_cell_f_all, &
         sum_cell_f_combined, sum_cell_f_constr, tc, th, tmp_const, ub, uij, zeff, nelectrons
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: cell_functions, coefficients, dE, &
                                                            distances, ds_dR_i, ds_dR_j, &
                                                            strength, target_val, &
                                                            sum_cell_f_group
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_sum_Pm_dR, &
                                                            distance_vecs, dP_i_dRi, &
                                                            electronic_charge, position_vecs, R12, gapw_offset
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dP_i_dRj, pair_dist_vecs
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec, dmy_dR_i, dmy_dR_j, &
                                                            dr, dr1_r2, dr_i_dR, dr_ij_dR, &
                                                            dr_j_dR, grid_p, r, r1, ra, shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoffs
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: gradients
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(becke_constraint_type), POINTER                :: becke_control
      TYPE(becke_group_type), POINTER, DIMENSION(:)      :: group
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_type), POINTER                      :: cavity_env
      TYPE(mpole_rho_atom), DIMENSION(:), POINTER        :: mp_rho
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_frag
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, charge
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(realspace_grid_type), POINTER                 :: rs_cavity
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole
      TYPE(section_vals_type), POINTER                   :: becke_constraint_section

      NULLIFY (cores, stride, atom_list, cell, para_env, dft_control, particle_set, rho_r, energy, rho, logger, &
               becke_constraint_section, qs_kind_set, particles, subsys, pab, pw_env, rs_cavity, cavity_env, &
               auxbas_pw_pool, atomic_kind_set, group)
      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, 
                      cell=cell, &
                      particle_set=particle_set, &
                      natom=natom, &
                      dft_control=dft_control, &
                      para_env=para_env)
      becke_constraint_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS%BECKE_CONSTRAINT")
      iw = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", extension=".cdftLog")
      becke_control => dft_control%qs_control%becke_control
      group => becke_control%group
      ! Sanity checks
      DO igroup = 1, SIZE(group)
         IF (.NOT. group(igroup)%is_charge_constraint .AND. dft_control%nspins == 1) &
            CALL cp_abort(__LOCATION__, &
                          "Becke magnetization density constraint requires UKS calculation.")
      END DO
      IF (natom < becke_control%natoms) &
         CALL cp_abort(__LOCATION__, &
                       "The number of constraint atoms exceeds the total number of atoms.")
      in_memory = .FALSE.
      IF (becke_control%save_pot) THEN
         in_memory = becke_control%in_memory
      END IF
      ! Setup atomic radii for adjusting cell boundaries
      IF (becke_control%adjust) THEN
         IF (.NOT. ASSOCIATED(becke_control%radii)) THEN
            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
            IF (.NOT. SIZE(atomic_kind_set) == SIZE(becke_control%radii_tmp)) &
               CALL cp_abort(__LOCATION__, &
                             "Length of keyword BECKE_CONSTRAINT\ATOMIC_RADII does not "// &
                             "match number of atomic kinds in the input coordinate file.")
            ALLOCATE (becke_control%radii(SIZE(atomic_kind_set)))
            becke_control%radii(:) = becke_control%radii_tmp(:)
            DEALLOCATE (becke_control%radii_tmp)
         END IF
      END IF
      ! Setup cutoff scheme
      IF (.NOT. ASSOCIATED(becke_control%cutoffs)) THEN
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
         ALLOCATE (becke_control%cutoffs(natom))
         SELECT CASE (becke_control%cutoff_type)
         CASE (becke_cutoff_global)
            becke_control%cutoffs(:) = becke_control%rglobal
         CASE (becke_cutoff_element)
            IF (.NOT. SIZE(atomic_kind_set) == SIZE(becke_control%cutoffs_tmp)) &
               CALL cp_abort(__LOCATION__, &
                             "Length of keyword BECKE_CONSTRAINT\ELEMENT_CUTOFFS does not "// &
                             "match number of atomic kinds in the input coordinate file.")
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
               DO iatom = 1, katom
                  atom_a = atom_list(iatom)
                  becke_control%cutoffs(atom_a) = becke_control%cutoffs_tmp(ikind)
               END DO
            END DO
            DEALLOCATE (becke_control%cutoffs_tmp)
         END SELECT
      END IF
      ! Zero weight functions
      DO igroup = 1, SIZE(group)
         group(igroup)%weight%pw%cr3d = 0.0_dp
      END DO
      IF (becke_control%atomic_charges) THEN
         DO iatom = 1, becke_control%natoms
            becke%charge(iatom)%pw%cr3d = 0.0_dp
         END DO
      END IF
      ! Setup static/dynamic confinement
      IF (becke_control%confine) THEN
         dir = becke_control%confine_dir
         lb = becke_control%confine_bounds(1)
         ub = becke_control%confine_bounds(2)
         ! Convert [0,L] to [-L/2, L/2]
         lb = lb-cell%hmat(dir, dir)/2._dp
         ub = ub-cell%hmat(dir, dir)/2._dp
      ELSE IF (becke_control%dynamic_confine) THEN
         dir = becke_control%confine_dir
         lb = HUGE(0.0_dp)
         ub = -HUGE(0.0_dp)
         dynamic_radius = becke_control%dynamic_radius
      END IF
      ! Allocate storage for cell adjustment coefficients and needed distance vectors
      build = .FALSE.
      IF (becke_control%adjust .AND. .NOT. ASSOCIATED(becke_control%aij)) THEN
         ALLOCATE (becke_control%aij(natom, natom))
         build = .TRUE.
      END IF 
      IF (becke_control%vector_buffer%store_vectors) THEN
         ALLOCATE (becke_control%vector_buffer%distances(natom))
         ALLOCATE (becke_control%vector_buffer%distance_vecs(3, natom))
         IF (in_memory) ALLOCATE (becke_control%vector_buffer%pair_dist_vecs(3, natom, natom))
         ALLOCATE (becke_control%vector_buffer%position_vecs(3, natom))
      END IF
      ALLOCATE (becke_control%vector_buffer%R12(natom, natom))
      ! Calculate pairwise distances between each atom pair
      DO i = 1, 3
         cell_v(i) = cell%hmat(i, i)
      END DO
      DO iatom = 1, natom-1
         DO jatom = iatom+1, natom
            r = particle_set(iatom)%r
            r1 = particle_set(jatom)%r
            DO i = 1, 3
               r(i) = MODULO(r(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
               r1(i) = MODULO(r1(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
            END DO
            dist_vec = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
            ! Store pbc corrected position and pairwise distance vectors for later reuse
            IF (becke_control%vector_buffer%store_vectors) THEN
               becke_control%vector_buffer%position_vecs(:, iatom) = r(:)
               IF (iatom == 1 .AND. jatom == natom) becke_control%vector_buffer%position_vecs(:, jatom) = r1(:)
               IF (in_memory) THEN
                  becke_control%vector_buffer%pair_dist_vecs(:, iatom, jatom) = dist_vec(:)
                  becke_control%vector_buffer%pair_dist_vecs(:, jatom, iatom) = -dist_vec(:)
               END IF
            END IF
            becke_control%vector_buffer%R12(iatom, jatom) = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
            becke_control%vector_buffer%R12(jatom, iatom) = becke_control%vector_buffer%R12(iatom, jatom)
            IF (becke_control%dynamic_confine) THEN
               IF (is_constraint(iatom) .AND. is_constraint(jatom)) THEN
                  ! ub = max(r(dir)+dynamic_radius), lb = min(r(dir)-dynamic_radius)
                  IF ((r(dir)+dynamic_radius) .GT. ub) ub = r(dir)+dynamic_radius
                  IF ((r1(dir)+dynamic_radius) .GT. ub) ub = r1(dir)+dynamic_radius
                  IF ((r(dir)-dynamic_radius) .LT. lb) lb = r(dir)-dynamic_radius
                  IF ((r1(dir)-dynamic_radius) .LT. lb) lb = r1(dir)-dynamic_radius
               END IF
            END IF
            ! Set up heteronuclear cell partitioning using user defined radii
            IF (build) THEN
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, kind_number=ikind)
               ircov = becke_control%radii(ikind)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, kind_number=ikind)
               jrcov = becke_control%radii(ikind)
               IF (ircov .NE. jrcov) THEN
                  chi = ircov/jrcov
                  uij = (chi-1.0_dp)/(chi+1.0_dp)
                  becke_control%aij(iatom, jatom) = uij/(uij**2-1.0_dp)
                  IF (becke_control%aij(iatom, jatom) .GT. 0.5_dp) THEN
                     becke_control%aij(iatom, jatom) = 0.5_dp
                  ELSE IF (becke_control%aij(iatom, jatom) .LT. -0.5_dp) THEN
                     becke_control%aij(iatom, jatom) = -0.5_dp
                  END IF
               ELSE
                  becke_control%aij(iatom, jatom) = 0.0_dp
               END IF
                ! Note change of sign
               becke_control%aij(jatom, iatom) = -becke_control%aij(iatom, jatom)
            END IF
         END DO
      END DO
      ! Dump some additional information about the calculation
      IF (becke_control%first_iteration) THEN
         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A,T60)') &
               '----------------------------- Becke atomic parameters -----------------------------'
            IF (becke_control%adjust) THEN
               WRITE (iw, '(T3,A)') &
                  'Atom  Element     Cutoff (angstrom)    CDFT Radius (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, element_symbol=element_symbol, &
                                       kind_number=ikind)
                  ircov = cp_unit_from_cp2k(becke_control%radii(ikind), "angstrom")
                  WRITE (iw, "(i7,T15,A2,T25,F8.3,T47,F8.3,T61,F8.3)") &
                     iatom, element_symbol, cp_unit_from_cp2k(becke_control%cutoffs(iatom), "angstrom"), &
                     ircov
               END DO
            ELSE
               WRITE (iw, '(T3,A)') &
                  'Atom  Element     Cutoff (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, element_symbol=element_symbol)
                  WRITE (iw, "(i7,T15,A2,T25,F8.3,T47,F8.3)") &
                     iatom, element_symbol, cp_unit_from_cp2k(becke_control%cutoffs(iatom), "angstrom")
               END DO
            END IF
            WRITE (iw, '(T3,A)') &
               '-----------------------------------------------------------------------------------'
            WRITE (iw, '(/,T3,A,T60)') &
               '----------------------------- Becke group definitions -----------------------------'
            DO igroup = 1, SIZE(group)
               IF (igroup > 1) WRITE (iw, '(T3,A)') ' '
               WRITE (iw, '(T3,A,I4,A,I4)') &
                  'Atomic group', igroup, ' of ', SIZE(group)
               WRITE (iw, '(T13,A)') ' Atom  Element'
               DO ip = 1, SIZE(group(igroup)%atoms)
                  iatom = group(igroup)%atoms(ip)
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, element_symbol=element_symbol)
                  WRITE (iw, '(T11,i7,T25,A2)') iatom, element_symbol
               END DO
            END DO
            WRITE (iw, '(T3,A)') &
               '-----------------------------------------------------------------------------------'
         END IF
         becke_control%first_iteration = .FALSE.
      END IF
      ! Convert confinement bounds to RS grid indices
      IF (becke_control%dynamic_confine .OR. &
          becke_control%confine) THEN
         tmp_const = (ub-shift(dir))/dr(dir)
         IF (tmp_const .GE. 0.0_dp) THEN
            ub_index = CEILING(tmp_const)
         ELSE
            ub_index = FLOOR(tmp_const)
         END IF
         tmp_const = (lb-shift(dir))/dr(dir)
         IF (tmp_const .GE. 0.0_dp) THEN
            lb_index = CEILING(tmp_const)
         ELSE
            lb_index = FLOOR(tmp_const)
         END IF
         ! Store confinement bound indices for use in force calculation
         becke_control%confine_bounds_int(1) = lb_index
         becke_control%confine_bounds_int(2) = ub_index
      END IF
      ! Setup cavity confinement using spherical Gaussians
      IF (becke_control%cavity_confine) THEN
         cavity_env => becke_control%cavity_env
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pw_env=pw_env, qs_kind_set=qs_kind_set)
         CPASSERT(ASSOCIATED(qs_kind_set))
         nkind = SIZE(qs_kind_set)
         ! Setup the Gaussian shape function
         IF (.NOT. ASSOCIATED(cavity_env%kind_shape_fn)) &
            CALL create_shape_function(cavity_env, qs_kind_set, atomic_kind_set, &
                                       radius=becke_control%rcavity, &
                                       radii_list=becke_control%radii)
         ! Form cavity by summing isolated Gaussian densities over constraint atoms
         NULLIFY (rs_cavity)
         CALL pw_env_get(pw_env, auxbas_rs_grid=rs_cavity, auxbas_pw_pool=auxbas_pw_pool)
         CALL rs_grid_retain(rs_cavity)
         CALL rs_grid_zero(rs_cavity)
         ALLOCATE (pab(1, 1))
         nthread = 1
         ithread = 0
         DO ikind = 1, SIZE(atomic_kind_set)
            numexp = cavity_env%kind_shape_fn(ikind)%numexp
            IF (numexp <= 0) CYCLE
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
            ALLOCATE (cores(katom))
            DO iex = 1, numexp
               alpha = cavity_env%kind_shape_fn(ikind)%zet(iex)
               coef = cavity_env%kind_shape_fn(ikind)%coef(iex)
               npme = 0
               cores = 0
               DO iatom = 1, katom
                  IF (rs_cavity%desc%parallel .AND. .NOT. rs_cavity%desc%distributed) THEN
                     ! replicated realspace grid, split the atoms up between procs
                     IF (MODULO(iatom, rs_cavity%desc%group_size) == rs_cavity%desc%my_pos) THEN
                        npme = npme+1
                        cores(npme) = iatom
                     ENDIF
                  ELSE
                     npme = npme+1
                     cores(npme) = iatom
                  ENDIF
               END DO
               DO j = 1, npme
                  iatom = cores(j)
                  atom_a = atom_list(iatom)
                  pab(1, 1) = coef
                  IF (store_vectors) THEN
                     ra(:) = position_vecs(:, atom_a)+cell_v(:)/2._dp
                  ELSE
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                  END IF
                  IF (is_constraint(atom_a)) &
                     CALL collocate_pgf_product_rspace(0, alpha, 0, 0, 0.0_dp, 0, ra, &
                                                       (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, 1.0_dp, &
                                                       pab, 0, 0, rs_cavity, cell, pw_env%cube_info(1), &
                                                       dft_control%qs_control%eps_rho_rspace, &
                                                       ga_gb_function=FUNC_AB, ithread=ithread, &
                                                       use_subpatch=.TRUE., subpatch_pattern=0_int_8, &
                                                       lmax_global=0)
               END DO
            END DO
            DEALLOCATE (cores)
         END DO
         DEALLOCATE (pab)
         CALL pw_pool_create_pw(auxbas_pw_pool, becke_control%cavity%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL rs_pw_transfer(rs_cavity, becke_control%cavity%pw, rs2pw)
         CALL rs_grid_release(rs_cavity)
         ! Ignore grid points where the Gaussian density falls below eps_cavity
         IF ((in_memory .OR. becke_control%save_pot) .AND. .NOT. becke_control%dynamic_confine) THEN
            CALL hfun_zero(becke_control%cavity%pw%cr3d, eps_cavity, just_bounds=.TRUE., bounds=bounds)
            becke_control%confine_bounds_int(1) = bounds(1)
            becke_control%confine_bounds_int(2) = bounds(2)
         END IF
         ! Optional printing of cavity (meant for testing, so options currently hardcoded...)
         IF (becke_control%print_cavity) THEN
            CALL hfun_zero(becke_control%cavity%pw%cr3d, eps_cavity, just_bounds=.FALSE.)
            ALLOCATE (stride(3))
            stride = (/2, 2, 2/)
            ! Note PROGRAM_RUN_INFO section neeeds to be active!
            unit_nr = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", &
                                           middle_name="BECKE_CAVITY", &
                                           extension=".cube", file_position="REWIND", &
                                           log_filename=.FALSE.)
            IF (para_env%mepos == para_env%source .AND. unit_nr .LT. 1) &
               CALL cp_abort(__LOCATION__, &
                             "Please turn on PROGRAM_RUN_INFO to print cavity")
            CALL get_qs_env(qs_env, subsys=subsys)
            CALL qs_subsys_get(subsys, particles=particles)
            CALL cp_pw_to_cube(becke_control%cavity%pw, unit_nr, "CAVITY", particles=particles, stride=stride)
            DEALLOCATE (stride)
         END IF
      END IF

   END SUBROUTINE becke_constraint_init

! **************************************************************************************************
!> \brief Calculate fout = fun1/fun2 or fout = fun1*fun2
!> \param fout the output 3D potential
!> \param fun1 the first input 3D potential
!> \param fun2 the second input 3D potential
!> \param divide logical that decides whether to divide or multiply the input potentials
! **************************************************************************************************
   SUBROUTINE hfun_scale(fout, fun1, fun2, divide)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: fout
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: fun1, fun2
      LOGICAL, INTENT(IN)                                :: divide

      CHARACTER(len=*), PARAMETER :: routineN = 'hfun_scale', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: small = 1.0e-12_dp

      INTEGER                                            :: i1, i2, i3, n1, n2, n3

      n1 = SIZE(fout, 1)
      n2 = SIZE(fout, 2)
      n3 = SIZE(fout, 3)
      CPASSERT(n1 == SIZE(fun1, 1))
      CPASSERT(n2 == SIZE(fun1, 2))
      CPASSERT(n3 == SIZE(fun1, 3))
      CPASSERT(n1 == SIZE(fun2, 1))
      CPASSERT(n2 == SIZE(fun2, 2))
      CPASSERT(n3 == SIZE(fun2, 3))

      IF (divide) THEN
         DO i3 = 1, n3
            DO i2 = 1, n2
               DO i1 = 1, n1
                  IF (fun2(i1, i2, i3) > small) THEN
                     fout(i1, i2, i3) = fun1(i1, i2, i3)/fun2(i1, i2, i3)
                  ELSE
                     fout(i1, i2, i3) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO
      ELSE
         DO i3 = 1, n3
            DO i2 = 1, n2
               DO i1 = 1, n1
                  fout(i1, i2, i3) = fun1(i1, i2, i3)*fun2(i1, i2, i3)
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE hfun_scale

! **************************************************************************************************
!> \brief Determine confinement bounds along confinement dir (hardcoded to be z)
!>        and optionally zero entries below a given threshold
!> \param fun input 3D potential (real space)
!> \param th threshold for screening values
!> \param just_bounds if the bounds should be computed without zeroing values
!> \param bounds the confinement bounds: fun is nonzero only between these values along 3rd dimension
! **************************************************************************************************
   SUBROUTINE hfun_zero(fun, th, just_bounds, bounds)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: fun
      REAL(KIND=dp), INTENT(IN)                          :: th
      LOGICAL                                            :: just_bounds
      INTEGER, OPTIONAL                                  :: bounds(2)

      CHARACTER(len=*), PARAMETER :: routineN = 'hfun_zero', routineP = moduleN//':'//routineN

      INTEGER                                            :: i1, i2, i3, lb, n1, n2, n3, nzeroed, ub
      LOGICAL                                            :: lb_final, ub_final

      n1 = SIZE(fun, 1)
      n2 = SIZE(fun, 2)
      n3 = SIZE(fun, 3)
      IF (just_bounds) THEN
         CPASSERT(PRESENT(bounds))
         lb = 1
         lb_final = .FALSE.
         ub_final = .FALSE.
      END IF

      DO i3 = 1, n3
         IF (just_bounds) nzeroed = 0
         DO i2 = 1, n2
            DO i1 = 1, n1
               IF (fun(i1, i2, i3) < th) THEN
                  IF (just_bounds) THEN
                     nzeroed = nzeroed+1
                  ELSE
                     fun(i1, i2, i3) = 0.0_dp
                  END IF
               ELSE
                  IF (just_bounds) EXIT
               END IF
            END DO
            IF (just_bounds .AND. nzeroed < n1) EXIT
         END DO
         IF (just_bounds) THEN
            IF (nzeroed == (n2*n1)) THEN
               IF (.NOT. lb_final) THEN
                  lb = i3
               ELSE IF (.NOT. ub_final) THEN
                  ub = i3
                  ub_final = .TRUE.
               END IF
            ELSE
               IF (.NOT. lb_final) lb_final = .TRUE.
               IF (ub_final) ub_final = .FALSE. ! Safeguard against "holes"
            END IF
         END IF
      END DO
      IF (just_bounds) THEN
         IF (.NOT. ub_final) ub = n3
         bounds(1) = lb
         bounds(2) = ub
         bounds = bounds-(n3/2)-1
      END IF

   END SUBROUTINE hfun_zero
END MODULE qs_cdft_utils