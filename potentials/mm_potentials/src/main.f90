!!****** gfit [1.0] *
!!
!!   NAME
!!     gfit, The main program
!!
!!   COPYRIGHT
!!I---------------------------------------------------------------------------I
!!I                                                                           I
!!I GFIT: is part of the cp2k distribution.                                   I
!!I A general program to perform gaussian fit of MM potential energy functionsI
!!I Copyright (C) 2000,2001,2002,2003,2004  CP2K developers group             I
!!I                                                                           I
!!I This program is free software; you can redistribute it and/or modify      I
!!I it under the terms of the GNU General Public License as published by      I
!!I the Free Software Foundation; either version 2 of the License, or         I
!!I (at your option) any later version.                                       I
!!I                                                                           I
!!I This program is distributed in the hope that it will be useful,           I
!!I but WITHOUT ANY WARRANTY; without even the implied warranty of            I
!!I MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             I
!!I GNU General Public License for more details.                              I
!!I                                                                           I
!!I You should have received a copy of the GNU General Public License         I
!!I along with this program; if not, write to the Free Software               I
!!I Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   FUNCTION
!!      starts the program
!!
!!   AUTHORL
!!      Teodoro Laino 
!!
!!   NOTES
!!
!------------------------------------------------------------------------------

PROGRAM GFIT
  USE mathconstants, ONLY: init_mathcon
  USE gaussian_fit_types,  ONLY : gaussian_fit_p_type  
  USE input, ONLY: read_input_file
  USE gaussian_fit, ONLY: fit 
  USE output,       ONLY: punch_results
  IMPLICIT NONE
  TYPE(gaussian_fit_p_type), DIMENSION(:), POINTER :: pgfs
  INTEGER ::  stat

  NULLIFY(pgfs)

  CALL init_mathcon
  CALL read_input_file(pgfs)
  CALL fit(pgfs)
  CALL punch_results(pgfs)

  IF (ASSOCIATED(pgfs)) DEALLOCATE(pgfs,stat=stat)

END PROGRAM GFIT
