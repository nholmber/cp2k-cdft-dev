!!****** gfit/global_variables [1.0] *
!!
!!   NAME
!!     kinds
!!
!!   FUNCTION
!!     Defines the global variables used through the GFIT code
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************
MODULE global_variables
  USE kinds,                           ONLY: default_string_length
  IMPLICIT NONE
  PRIVATE
  
  CHARACTER (LEN=default_string_length), PARAMETER,  PUBLIC :: IONAME="input.dat"
  
END MODULE global_variables
