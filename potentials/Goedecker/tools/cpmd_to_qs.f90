PROGRAM cpmd_to_qs

! Purpose: Convert the output file XX of the program pseudo.x for the
!          generation of Goedecker-Teter-Hutter (GTH) pseudo potentials which
!          is written in CPMD-format to the Quickstep potential database
!          format.

! Input file:  XX
! Output file: QS

! History: - Creation (12.12.2003,MK)

! *****************************************************************************

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(14,200),&
                        max_ppl = 4,&  ! Max. number of local projectors
                        max_ppnl = 4   ! Max. number of non-local projectors

  CHARACTER(LEN=2), DIMENSION(103) :: elesym =&
    (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg","Al","Si",&
      "P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni",&
      "Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr","Nb","Mo",&
      "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe","Cs","Ba",&
      "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",&
      "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",&
      "At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf",&
      "Es","Fm","Md","No","Lr"/)

  CHARACTER(LEN=200) :: line
  CHARACTER(LEN=12)  :: fmtstr,xc_string
  REAL(KIND=wp)      :: rloc,z,zeff
  INTEGER            :: i,ippnl,istat,j,nppl,nppnl_max,xc_code

  REAL(KIND=wp), DIMENSION(max_ppl)      :: cppl
  REAL(KIND=wp), DIMENSION(max_ppnl)     :: rppnl
  REAL(KIND=wp), DIMENSION(max_ppnl,3,3) :: cppnl

  INTEGER, DIMENSION(max_ppnl) :: nppnl

! -----------------------------------------------------------------------------

  z = 0.0_wp
  zeff = 0.0_wp

  nppnl(:) = 0

  cppl(:) = 0.0_wp
  rppnl(:) = 0.0_wp
  cppnl(:,:,:) = 0.0_wp

! *** Open I/O units ***

  OPEN (UNIT=5,&
        FILE="XX",&
        STATUS="OLD",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="READ",&
        IOSTAT=istat)

  IF (istat /= 0) STOP "*** Could not open input file XX ***"

  OPEN (UNIT=6,&
        FILE="QS",&
        STATUS="REPLACE",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="WRITE",&
        IOSTAT=istat)

  IF (istat /= 0) STOP "*** Could not open output file QS ***"

! *** Read XX file (CPMD format) ***

  DO
    READ (UNIT=5,FMT="(A)",IOSTAT=istat) line
    IF (istat /= 0) STOP "*** Error reading input file XX ***"
    IF (INDEX(line,"Z  =") > 0) THEN
      READ (UNIT=line(7:),FMT=*) z
    ELSE IF (INDEX(line,"ZV =") > 0) THEN
      READ (UNIT=line(7:),FMT=*) zeff
    ELSE IF (INDEX(line,"XC =") > 0) THEN
      READ (UNIT=line(7:),FMT=*) xc_code
    ELSE IF (INDEX(line,"&POTENTIAL") > 0) THEN
      EXIT
    END IF
  END DO

  IF (z == 0.0_wp) STOP "*** Z is zero ***"
  IF (zeff == 0.0_wp) STOP "*** Z(eff) is zero ***"

  SELECT CASE (xc_code)
  CASE (1312)
    xc_string = "BLYP"
  CASE (1111)
    xc_string = "BP"
  CASE (0900)
    xc_string = "PADE"
  CASE (1134)
    xc_string = "PBE"
  CASE DEFAULT
    STOP "*** Invalid XC code found in input file XX ***"
  END SELECT

  READ (UNIT=5,FMT="(A)") line
  IF (INDEX(line,"GOEDECKER") == 0) STOP "*** Keyword GOEDECKER not found ***"
  READ (UNIT=5,FMT=*) nppnl_max
  READ (UNIT=5,FMT=*) rloc
  READ (UNIT=5,FMT=*) nppl,(cppl(i),i=1,nppl)
  DO ippnl=1,nppnl_max
    READ (UNIT=5,FMT=*) rppnl(ippnl),nppnl(ippnl),&
                        ((cppnl(ippnl,i,j),j=i,nppnl(ippnl)),i=1,nppnl(ippnl))
  END DO

! *** Quickstep database format ***

  WRITE (UNIT=6,FMT="(A,1X,A)") elesym(NINT(z)),"GTH-"//TRIM(xc_string)
  WRITE (UNIT=6,FMT="(I5,A)")&
    NINT(z)," -> split into total number of s p d ... electrons (all-elec. atom)"
  WRITE (UNIT=6,FMT="(I5,A)")&
    NINT(zeff)," -> split into total number of s p d ... electrons (pseudo atom)"
  WRITE (UNIT=6,FMT="(F15.8,I5,4F15.8)") rloc,nppl,(cppl(i),i=1,nppl)
  WRITE (UNIT=6,FMT="(I5)") nppnl_max
  DO ippnl=1,nppnl_max
    WRITE (UNIT=6,FMT="(F15.8,I5,4F15.8)")&
      rppnl(ippnl),nppnl(ippnl),(cppnl(ippnl,1,j),j=1,nppnl(ippnl))
    fmtstr = "(T  ,4F15.8)"
    DO i=2,nppnl(ippnl)
      WRITE (fmtstr(3:4),"(I2)") 15*i + 6
      WRITE (UNIT=6,FMT=fmtstr) (cppnl(ippnl,i,j),j=i,nppnl(ippnl))
    END DO
  END DO

END PROGRAM cpmd_to_qs
