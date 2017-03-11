!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Converts mutable data to linear (array) type.
!>
!> \param[in,out] wm      work matrix to convert
! **************************************************************************************************
  SUBROUTINE tree_to_linear_d(wm)
    USE btree,&
        ONLY: btree_2d_data_d => btree_data_dp2d,&
              btree_destroy_d => btree_delete,&
              btree_size_d => btree_get_entries
    TYPE(dbcsr_work_type), INTENT(INOUT)     :: wm

    CHARACTER(len=*), PARAMETER :: routineN = 'tree_to_linear_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, blk_p, treesize, &
                                                error_handler, needed_size
    INTEGER(KIND=int_8), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    REAL(kind=real_8), DIMENSION(:), POINTER           :: target_data
    REAL(kind=real_8), DIMENSION(:, :), POINTER        :: block_2d
    TYPE(btree_2d_data_d), ALLOCATABLE, &
      DIMENSION(:)                           :: values

!   ---------------------------------------------------------------------------

    CALL timeset(routineN, error_handler)
    ! srt = .TRUE. ! Not needed because of the copy
    treesize = btree_size_d(wm%mutable%m%btree_d)
    IF(wm%lastblk .NE. treesize)&
       CPABORT("Mismatch in number of blocks")
    ALLOCATE (keys(treesize), values(treesize))
    CALL btree_destroy_d (wm%mutable%m%btree_d, keys, values)
    CALL ensure_array_size (wm%row_i, ub=treesize)
    CALL ensure_array_size (wm%col_i, ub=treesize)
    CALL dbcsr_unpack_i8_2i4 (keys, wm%row_i,&
         wm%col_i)
    ! For now we also fill the data, sloooowly, but this should
    ! be avoided and the data should be copied directly from the
    ! source in the subroutine's main loop.
    CALL ensure_array_size (wm%blk_p, ub=treesize)
    needed_size=0
    DO blk= 1, treesize
       block_2d => values(blk)%p
       needed_size=needed_size+SIZE(block_2d)
    ENDDO
    wm%datasize=needed_size
    CALL dbcsr_data_ensure_size (wm%data_area,&
         wm%datasize)
    target_data => dbcsr_get_data_p_d (wm%data_area)
    blk_p = 1
    DO blk = 1, treesize
       block_2d => values(blk)%p
       IF (.NOT. values(blk)%tr) THEN
          wm%blk_p(blk) = blk_p
       ELSE
          wm%blk_p(blk) = -blk_p
       ENDIF
       CALL block_copy_d (target_data, block_2d,&
            SIZE (block_2d), blk_p, 1)
       blk_p = blk_p + SIZE(block_2d)
       DEALLOCATE (block_2d)
    ENDDO
    DEALLOCATE (keys, values)
    CALL dbcsr_mutable_release (wm%mutable)
    CALL timestop(error_handler)
  END SUBROUTINE tree_to_linear_d

