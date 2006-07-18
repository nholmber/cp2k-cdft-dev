PROGRAM cpmd_to_qs

  ! Purpose: Convert the output file XX of the program pseudo.x for the
  !          generation of Goedecker-Teter-Hutter (GTH) pseudo potentials which
  !          is written in CPMD-format to the Quickstep potential database
  !          format.

  ! History: - Creation (12.12.2003,MK)

  ! ***************************************************************************

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(14,200),&
                        maxl = 3,&    ! Max. angular momentum quantum number
                        max_ppl = 4,& ! Max. number of local projectors
                        max_ppnl = 4  ! Max. number of non-local projectors

  REAL(KIND=wp), PARAMETER :: eps_zero = 1.0E-10_wp

  CHARACTER(LEN=200) :: input_file1,input_file2,line,output_file
  CHARACTER(LEN=80)  :: elec_string,string
  CHARACTER(LEN=17)  :: fmtstr1
  CHARACTER(LEN=16)  :: fmtstr3
  CHARACTER(LEN=12)  :: fmtstr2,xc_string
  REAL(KIND=wp)      :: nelec,rloc,q,z,zeff
  INTEGER            :: i,idigits,ippnl,istat,iz,j,l,lmax,n,narg,ncore,nppl,&
                        nppnl_max,nvalence,xc_code

  CHARACTER(LEN=2), DIMENSION(103) :: elesym =&
    (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg","Al","Si",&
      "P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni",&
      "Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr","Nb","Mo",&
      "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe","Cs","Ba",&
      "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",&
      "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",&
      "At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf",&
      "Es","Fm","Md","No","Lr"/)

  REAL(KIND=wp), DIMENSION(maxl+1)       :: elec_conf
  REAL(KIND=wp), DIMENSION(max_ppl)      :: cppl
  REAL(KIND=wp), DIMENSION(max_ppnl)     :: rppnl
  REAL(KIND=wp), DIMENSION(max_ppnl,3,3) :: cppnl

  INTEGER, DIMENSION(max_ppnl) :: nppnl

  INTEGER :: iargc
  LOGICAL :: frac_elec

  ! ---------------------------------------------------------------------------

  input_file1 = ""
  input_file2 = ""
  output_file = ""

  z = 0.0_wp
  zeff = 0.0_wp

  nppnl(:) = 0

  cppl(:) = 0.0_wp
  rppnl(:) = 0.0_wp
  cppnl(:,:,:) = 0.0_wp

  ! Check the number of arguments

  narg = iargc()
  IF (narg /= 3) THEN
    PRINT*,"ERROR: Expected three valid file names as arguments, found",narg
    PRINT*,"       <input_file1> <input_file2> <output_file>"
    PRINT*,"       e.g. XX atom.dat QS"
    STOP
  END IF

  ! Open I/O units

  CALL getarg(1,input_file1)
  OPEN (UNIT=1,&
        FILE=TRIM(input_file1),&
        STATUS="OLD",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="READ",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the first input file "//TRIM(input_file1)
    STOP
  END IF

  CALL getarg(2,input_file2)
  OPEN (UNIT=2,&
        FILE=TRIM(input_file2),&
        STATUS="OLD",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="READ",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the second input file "//TRIM(input_file2)
    STOP
  END IF

  CALL getarg(3,output_file)
  OPEN (UNIT=3,&
        FILE=TRIM(output_file),&
        STATUS="REPLACE",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="WRITE",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the output file "//TRIM(output_file)
    STOP
  END IF

  OPEN (UNIT=4,&
        FILE="INFO",&
        STATUS="REPLACE",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="WRITE",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the output file INFO"
    STOP
  END IF

  OPEN (UNIT=10,&
        FILE="TEXTAB",&
        STATUS="REPLACE",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="WRITE",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the output file TEXTAB"
    STOP
  END IF

  ! Read first input file (CPMD format)

  DO
    READ (UNIT=1,FMT="(A)",IOSTAT=istat) line
    IF (istat /= 0) THEN
      PRINT*,"ERROR reading the first input file "//TRIM(input_file1)
      STOP
    END IF
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

  IF (z == 0.0_wp) THEN
    PRINT*,"ERROR: Z is zero in the first input file "//TRIM(input_file1)
    STOP
  END IF
  IF (zeff == 0.0_wp) THEN
    PRINT*,"ERROR: Z(eff) is zero in the first input file "//TRIM(input_file1)
    STOP
  END IF

  SELECT CASE (xc_code)
  CASE (1312)
    xc_string = "BLYP"
  CASE (1111)
    xc_string = "BP"
  CASE (0900)
    xc_string = "PADE"
  CASE (0055)
    xc_string = "HCTH" ! CPMD version of HCTH120
  CASE (0066)
    xc_string = "HCTH93"
  CASE (0077)
    xc_string = "HCTH120"
  CASE (0088)
    xc_string = "HCTH147"
  CASE (0099)
    xc_string = "HCTH407"
  CASE (1134)
    xc_string = "PBE"
  CASE (0302)
    xc_string = "OLYP"
  CASE DEFAULT
    PRINT*,"ERROR: Invalid XC code found in the first input file "//TRIM(input_file1)
    STOP
  END SELECT

  READ (UNIT=1,FMT="(A)") line
  IF (INDEX(line,"GOEDECKER") == 0) THEN
    PRINT*,"ERROR: No keyword GOEDECKER found in the first input file "//TRIM(input_file1)
    STOP
  END IF
  READ (UNIT=1,FMT=*) nppnl_max
  READ (UNIT=1,FMT=*) rloc
  READ (UNIT=1,FMT=*) nppl,(cppl(i),i=1,nppl)
  DO ippnl=1,nppnl_max
    READ (UNIT=1,FMT=*) rppnl(ippnl),nppnl(ippnl),&
                        ((cppnl(ippnl,i,j),j=i,nppnl(ippnl)),i=1,nppnl(ippnl))
  END DO

  iz = NINT(z)
  WRITE (UNIT=4,FMT="(2(A,/,/),A,I3)")&
    " "//REPEAT("*",64),&
    " Atomic symbol                       : "//TRIM(elesym(iz)),&
    " Atomic number                       : ",iz
  IF (ABS(zeff - INT(zeff)) > eps_zero) THEN
    WRITE (UNIT=4,FMT="(A,F8.3)")&
      " Effective core charge               : ",zeff
  ELSE
    WRITE (UNIT=4,FMT="(A,I3)")&
      " Effective core charge               : ",INT(zeff)
  END IF

  ! Read second input file (electronic configuration)

  READ (UNIT=2,FMT="(A)") string
  IF (TRIM(ADJUSTL(string)) /= TRIM(elesym(iz))) THEN
    PRINT*,"ERROR: Mismatching element symbols found in the input files"
    STOP
  END IF
  READ (UNIT=2,FMT="(A)") string
  IF (TRIM(ADJUSTL(string)) /= xc_string) THEN
    PRINT*,"ERROR: Mismatching XC functionals found in the input files"
    STOP
  END IF
  READ (UNIT=2,FMT=*) line ! dummy read
  READ (UNIT=2,FMT=*) line ! dummy read
  READ (UNIT=2,FMT=*) line ! dummy read
  READ (UNIT=2,FMT=*) ncore,nvalence
  WRITE (UNIT=4,FMT="(A,I3,/,A,I3)")&
    " Number of core states               : ",ncore,&
    " Number of valence states            : ",nvalence

  ! Build electronic configuration

  frac_elec = .FALSE.
  elec_conf(:) = 0.0_wp

  DO i=1,nvalence
    READ (UNIT=2,FMT=*) n,l,q
    elec_conf(l+1) = elec_conf(l+1) + q
    IF (ABS(q - INT(q)) > eps_zero) frac_elec = .TRUE. ! we have a non-integer number of electrons
  END DO

  ! This is a hack for the elements from Hf to Pt, since
  ! the 4f electrons are not treated as valence electrons

  IF ((72 <= iz).AND.(iz <= 78)) elec_conf(4) = 0.0_wp

  nelec = SUM(elec_conf) 

  ! Find the occupied orbitals with the largest l value

  DO i=maxl+1,1,-1
    IF (elec_conf(i) > eps_zero) THEN
      lmax = i - 1
      EXIT
    END IF
  END DO

  ! Build a printable string with the electronic configuration

  elec_string = ""

  DO i=0,lmax
    IF (frac_elec) THEN
      WRITE (UNIT=elec_string(8*i+1:),FMT="(F8.3)") elec_conf(i+1)
    ELSE
      WRITE (UNIT=elec_string(5*i+1:),FMT="(I5)") INT(elec_conf(i+1))
    END IF
  END DO

  WRITE (UNIT=4,FMT="(A,A)")&
    " Electronic configuration (s,p,d,...): ",TRIM(elec_string)

  IF (ABS(zeff - nelec) > eps_zero) THEN
    WRITE (UNIT=4,FMT="(A,F8.3)")&
      " Ionic charge                        : ",zeff - nelec
  END IF

  ! TeX tabular format

  WRITE (UNIT=10,FMT="(T2,A5,I3,A3,F10.6,4(A3,F13.6),A)")&
    elesym(iz)//" & ",NINT(zeff)," & ",rloc,(" & ",cppl(i),i=1,4),"\\\\"

  DO ippnl=1,nppnl_max
    IF (nppnl(ippnl) > 0) THEN
      WRITE (UNIT=10,FMT="(T2,A11,F10.6,4(A3,F13.6),A)")&
        "   &     & ",rppnl(ippnl),(" & ",cppnl(ippnl,1,j),j=1,4),"\\\\"
      DO i=2,nppnl(ippnl)
        WRITE (UNIT=10,FMT="(T5,A1,T11,A1,T23,3(A3,F13.6),T72,A1,T87,A)")&
          "&","&",(" & ",cppnl(ippnl,i,j),j=1,3),"&","\\\\"
      END DO
    END IF
  END DO

  ! Quickstep database format

  string = ""

  WRITE (UNIT=string,FMT="(I5)") NINT(zeff)
  IF (xc_string == "PADE") THEN
  WRITE (UNIT=3,FMT="(A,1X,A)")&
    TRIM(elesym(iz)),"GTH-"//TRIM(xc_string)//"-q"//TRIM(ADJUSTL(string))//&
    " GTH-LDA-q"//TRIM(ADJUSTL(string))
  ELSE
    WRITE (UNIT=3,FMT="(A,1X,A)")&
      TRIM(elesym(iz)),"GTH-"//TRIM(xc_string)//"-q"//TRIM(ADJUSTL(string))
  END IF

  WRITE (UNIT=3,FMT="(A)") TRIM(elec_string)

  ! Set output precision for the QS database files

  idigits = 8

  fmtstr1 = "(F15. ,I5,4F15. )"
  WRITE (UNIT=fmtstr1(6:6),FMT="(I1)") idigits
  WRITE (UNIT=fmtstr1(16:16),FMT="(I1)") idigits

  WRITE (UNIT=3,FMT=fmtstr1) rloc,nppl,(cppl(i),i=1,nppl)

  WRITE (UNIT=3,FMT="(I5)") nppnl_max
  DO ippnl=1,nppnl_max
    WRITE (UNIT=3,FMT=fmtstr1)&
      rppnl(ippnl),nppnl(ippnl),(cppnl(ippnl,1,j),j=1,nppnl(ippnl))
    fmtstr2 = "(T  ,4F15. )"
    WRITE (UNIT=fmtstr2(11:11),FMT="(I1)") idigits
    DO i=2,nppnl(ippnl)
      WRITE (UNIT=fmtstr2(3:4),FMT="(I2)") 15*i + 6
      WRITE (UNIT=3,FMT=fmtstr2) (cppnl(ippnl,i,j),j=i,nppnl(ippnl))
    END DO
  END DO

  ! Write the data set also to INFO file

  WRITE (UNIT=4,FMT="(/,A)")&
    " Exchange-correlation functional     : "//TRIM(xc_string)

  ! Set output precision for the CPMD INFO section

  idigits = 6

  fmtstr3 = "(/,F13. ,4F13. )"
  WRITE (UNIT=fmtstr3(8:8),FMT="(I1)") idigits
  WRITE (UNIT=fmtstr3(15:15),FMT="(I1)") idigits

  WRITE (UNIT=4,FMT="(/,T8,A)",ADVANCE="NO") "r(loc)"
  DO i=1,nppl
    WRITE (UNIT=4,FMT="(9X,A2,I1,A1)",ADVANCE="NO") "C(",i,")"
  END DO
  WRITE (UNIT=4,FMT=fmtstr3) rloc,(cppl(i),i=1,nppl)

  DO ippnl=1,nppnl_max
    IF (nppnl(ippnl) > 0) THEN
      IF (ippnl == 1) THEN
        WRITE (UNIT=4,FMT="(/,9X,A2,I1,A1,5X,A,I1)",ADVANCE="NO")&
          "r(",ippnl-1,")","h(i,j)^",ippnl-1
      ELSE
        WRITE (UNIT=4,FMT="(9X,A2,I1,A1,5X,A,I1)",ADVANCE="NO")&
          "r(",ippnl-1,")","h(i,j)^",ippnl-1
      END IF
      WRITE (UNIT=4,FMT=fmtstr3)&
        rppnl(ippnl),(cppnl(ippnl,1,j),j=1,nppnl(ippnl))
      fmtstr2 = "(T  ,4F13. )"
      WRITE (UNIT=fmtstr2(11:11),FMT="(I1)") idigits
      DO i=2,nppnl(ippnl)
        WRITE (UNIT=fmtstr2(3:4),FMT="(I2)") 13*i + 1
        WRITE (UNIT=4,FMT=fmtstr2) (cppnl(ippnl,i,j),j=i,nppnl(ippnl))
      END DO
    END IF
  END DO

  WRITE (UNIT=4,FMT="(/,A,/,6(/,A),/,/,A)")&
    " Please cite:",&
    " - S. Goedecker, M. Teter, and J. Hutter,",&
    "   Phys. Rev. B 54, 1703 (1996)",&
    " - C. Hartwigsen, S. Goedecker, and J. Hutter,",&
    "   Phys. Rev. B 58, 3641 (1998)",&
    " - M. Krack,",&
    "   Theor. Chem. Acc. 114, 145 (2005)",&
    " "//REPEAT("*",64)

END PROGRAM cpmd_to_qs
