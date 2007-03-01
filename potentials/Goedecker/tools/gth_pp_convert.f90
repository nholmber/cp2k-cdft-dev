PROGRAM gth_pp_convert

  ! Purpose: Convert the output file psp.par of the program pseudo.x for the
  !          generation of Goedecker-Teter-Hutter (GTH) pseudo potentials to
  !          the Quickstep and Abinit potential database format.

  ! History: - Creation (12.12.2003,MK)
  !          - ABINIT format output added (27.02.07,MK)

  ! ***************************************************************************

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(14,200),&
                        maxl = 3,&      ! Max. angular momentum quantum number
                        max_ppl = 4,&   ! Max. number of local projectors
                        max_ppnl = 4,&  ! Max. number of non-local projectors
                        unit_abinit   = 10,& ! Units ...
                        unit_atom_dat = 11,&
                        unit_info     = 12,&
                        unit_psp_par  = 13,&
                        unit_cp2k     = 14,&
                        unit_textab   = 15,&
                        unit_xx       = 16

  REAL(KIND=wp), PARAMETER :: eps_zero = 1.0E-8_wp

  CHARACTER(LEN=1), DIMENSION(7), PARAMETER :: llabel = (/"s","p","d","f",&
                                                          "g","h","i"/)

  CHARACTER(LEN=200) :: input_file1,input_file2,input_file3,line
  CHARACTER(LEN=80)  :: elec_string,string
  CHARACTER(LEN=17)  :: fmtstr1
  CHARACTER(LEN=16)  :: fmtstr3
  CHARACTER(LEN=12)  :: fmtstr2,xc_string,xcf_psp_par
  REAL(KIND=wp)      :: nelec,rloc,q,xx_rloc,xx_z,xx_zeff,z,zeff
  INTEGER            :: i,idigits,ippnl,istat,iz,izeff,j,l,lmax,lppnl_max,n,&
                        narg,ncore,nppl,nppnl_max,nvalence,xc_code,xc_code_abinit

  CHARACTER(LEN=2), DIMENSION(103) :: elesym =&
    (/"H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg","Al","Si",&
      "P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni",&
      "Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y ","Zr","Nb","Mo",&
      "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe","Cs","Ba",&
      "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",&
      "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",&
      "At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf",&
      "Es","Fm","Md","No","Lr"/)

  REAL(KIND=wp), DIMENSION(maxl+1)         :: elec_conf
  REAL(KIND=wp), DIMENSION(max_ppl)        :: cppl,xx_cppl
  REAL(KIND=wp), DIMENSION(max_ppnl)       :: rppnl,xx_rppnl
  REAL(KIND=wp), DIMENSION(3,3,max_ppnl)   :: hppnl,kppnl,xx_hppnl
  REAL(KIND=wp), DIMENSION(3,3,max_ppnl,2) :: cppnl

  INTEGER, DIMENSION(max_ppnl) :: nppnl

  INTEGER :: iargc
  LOGICAL :: frac_elec

  ! ---------------------------------------------------------------------------

  ! Initialization

  input_file1 = ""
  input_file2 = ""
  input_file3 = ""

  iz = 0
  izeff = 0
  xc_code_abinit = 0

  z = 0.0_wp
  xx_z = 0.0_wp
  zeff = 0.0_wp
  xx_zeff = 0.0_wp
  cppl(:) = 0.0_wp
  xx_cppl(:) = 0.0_wp
  rppnl(:) = 0.0_wp
  xx_rppnl(:) = 0.0_wp
  hppnl(:,:,:) = 0.0_wp
  xx_hppnl(:,:,:) = 0.0_wp
  kppnl(:,:,:) = 0.0_wp
  cppnl(:,:,:,:) = 0.0_wp

  nppnl(:) = 0

  ! Check the number of arguments

  narg = iargc()
  IF (narg /= 3) THEN
    PRINT*,"ERROR: Expected four valid file names as arguments, found",narg
    PRINT*,"       <input_file1> <input_file2> <input_file3>"
    PRINT*,"       e.g. XX atom.dat psp.par"
    STOP
  END IF

  ! Open I/O units

  CALL getarg(1,input_file1)
  OPEN (UNIT=unit_xx,&
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
  OPEN (UNIT=unit_atom_dat,&
        FILE=TRIM(input_file2),&
        STATUS="OLD",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="READ",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the first input file "//TRIM(input_file2)
    STOP
  END IF

  CALL getarg(3,input_file3)
  OPEN (UNIT=unit_psp_par,&
        FILE=TRIM(input_file3),&
        STATUS="OLD",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="READ",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the second input file "//TRIM(input_file3)
    STOP
  END IF

  OPEN (UNIT=unit_cp2k,&
        FILE="CP2K",&
        STATUS="REPLACE",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="WRITE",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the output file CP2K"
    STOP
  END IF

  OPEN (UNIT=unit_abinit,&
        FILE="ABINIT",&
        STATUS="REPLACE",&
        ACCESS="SEQUENTIAL",&
        FORM="FORMATTED",&
        POSITION="REWIND",&
        ACTION="WRITE",&
        IOSTAT=istat)
  IF (istat /= 0) THEN
    PRINT*,"ERROR: Could not open the output file ABINIT"
    STOP
  END IF

  OPEN (UNIT=unit_info,&
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

  OPEN (UNIT=unit_textab,&
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

  ! Read first input file XX (CPMD format)

  DO
    READ (UNIT=unit_xx,FMT="(A)",IOSTAT=istat) line
    IF (istat /= 0) THEN
      PRINT*,"ERROR reading the first input file "//TRIM(input_file1)
      STOP
    END IF
    IF (INDEX(line,"Z  =") > 0) THEN
      READ (UNIT=line(7:),FMT=*) xx_z
    ELSE IF (INDEX(line,"ZV =") > 0) THEN
      READ (UNIT=line(7:),FMT=*) xx_zeff
    ELSE IF (INDEX(line,"XC =") > 0) THEN
      READ (UNIT=line(7:),FMT=*) xc_code
    ELSE IF (INDEX(line,"&POTENTIAL") > 0) THEN
      EXIT
    END IF
  END DO

  IF (xx_z == 0.0_wp) THEN
    PRINT*,"ERROR: Z is zero in the first input file "//TRIM(input_file1)
    STOP
  END IF
  IF (xx_zeff == 0.0_wp) THEN
    PRINT*,"ERROR: Z(eff) is zero in the first input file "//TRIM(input_file1)
    STOP
  END IF

  SELECT CASE (xc_code)
  CASE (1312)
    xc_string = "BLYP"
    xc_code_abinit = 0
  CASE (1111)
    xc_string = "BP"
    xc_code_abinit = 0
  CASE (0900)
    xc_string = "PADE"
    xc_code_abinit = 1
  CASE (0055)
    xc_string = "HCTH" ! CPMD version of HCTH120
    xc_code_abinit = 16
  CASE (0066)
    xc_string = "HCTH93"
    xc_code_abinit = 0
  CASE (0077)
    xc_string = "HCTH120"
    xc_code_abinit = 16
  CASE (0088)
    xc_string = "HCTH147"
    xc_code_abinit = 0
  CASE (0099)
    xc_string = "HCTH407"
    xc_code_abinit = 0
  CASE (1134)
    xc_string = "PBE"
    xc_code_abinit = 11
  CASE (0302)
    xc_string = "OLYP"
    xc_code_abinit = 0
  CASE DEFAULT
    PRINT*,"ERROR: Invalid XC code found in the first input file "//TRIM(input_file1)
    STOP
  END SELECT

  READ (UNIT=unit_xx,FMT="(A)") line
  IF (INDEX(line,"GOEDECKER") == 0) THEN
    PRINT*,"ERROR: No keyword GOEDECKER found in the first input file "//TRIM(input_file1)
    STOP
  END IF
  READ (UNIT=unit_xx,FMT=*) nppnl_max
  READ (UNIT=unit_xx,FMT=*) xx_rloc
  READ (UNIT=unit_xx,FMT=*) nppl,(xx_cppl(i),i=1,nppl)
  DO ippnl=1,nppnl_max
    READ (UNIT=unit_xx,FMT=*) xx_rppnl(ippnl),nppnl(ippnl),&
      ((xx_hppnl(i,j,ippnl),j=i,nppnl(ippnl)),i=1,nppnl(ippnl))
  END DO

  ! Read input file psp.par and compare to the XX data

  READ (UNIT=unit_psp_par,FMT="(A)",IOSTAT=istat) line ! dummy read
  READ (UNIT=unit_psp_par,FMT="(A)",IOSTAT=istat) line ! dummy read
  READ (UNIT=unit_psp_par,FMT="(A)",IOSTAT=istat) line ! dummy read
  READ (UNIT=unit_psp_par,FMT=*,IOSTAT=istat) xcf_psp_par
  READ (UNIT=unit_psp_par,FMT=*,IOSTAT=istat) z,zeff,rloc,(cppl(i),i=1,max_ppl)
  IF (z /= xx_z) THEN
    PRINT*,"ERROR: nuclear charges do not match:",z,xx_z
    STOP
  END IF
  IF (zeff /= xx_zeff) THEN
    PRINT*,"ERROR: effective nuclear charges do not match:",zeff,xx_zeff
    STOP
  END IF
  IF (ABS(rloc - xx_rloc) > eps_zero) THEN
    PRINT*,"ERROR: rloc values do not match:",rloc,xx_rloc
    STOP
  END IF
  DO i=1,max_ppl
    IF (ABS(cppl(i) - xx_cppl(i)) > eps_zero) THEN
      PRINT*,"ERROR: cppl(",i,") values do not match:",cppl(i),xx_cppl(i)
      STOP
    END IF
  END DO
  READ (UNIT=unit_psp_par,FMT=*,IOSTAT=istat) lppnl_max
  IF (lppnl_max + 1 /= nppnl_max) THEN
    PRINT*,"ERROR: number of non-local projectors does not match:",lppnl_max+1,nppnl_max
    STOP
  END IF
  DO ippnl=1,nppnl_max
    READ (UNIT=unit_psp_par,FMT=*) rppnl(ippnl),&
      ((cppnl(j,i,ippnl,1),j=1,i),i=1,3)
    IF (ABS(rppnl(ippnl) - xx_rppnl(ippnl)) > eps_zero) THEN
      PRINT*,"ERROR: rppnl(",ippnl,") values do not match:",rppnl(ippnl),xx_rppnl(ippnl)
      STOP
    END IF
    IF (ippnl > 1) THEN
      READ (UNIT=unit_psp_par,FMT=*)&
        ((cppnl(j,i,ippnl,2),j=1,i),i=1,3)
    END IF
    DO i=1,3
      DO j=i,3
        hppnl(i,j,ippnl) = (REAL(ippnl,KIND=wp)*cppnl(i,j,ippnl,1) +&
                            REAL(ippnl-1,KIND=wp)*cppnl(i,j,ippnl,2))/&
                           REAL(2*ippnl-1,KIND=wp)
        IF (ABS(hppnl(i,j,ippnl) - xx_hppnl(i,j,ippnl)) > eps_zero) THEN
          PRINT*,"ERROR: hppnl values do not match:",hppnl(i,j,ippnl),xx_hppnl(i,j,ippnl)
          STOP
        END IF
        kppnl(i,j,ippnl) = 2.0_wp*(cppnl(i,j,ippnl,1) - cppnl(i,j,ippnl,2))/&
                           REAL(2*ippnl-1,KIND=wp)
      END DO
    END DO
    ! Symmetrize matrices
    DO i=1,3
      DO j=i+1,3
        cppnl(j,i,ippnl,1) = cppnl(i,j,ippnl,1)
        cppnl(j,i,ippnl,2) = cppnl(i,j,ippnl,2)
        hppnl(j,i,ippnl) = hppnl(i,j,ippnl)
        kppnl(j,i,ippnl) = kppnl(i,j,ippnl)
      END DO
    END DO
  END DO

  ! Write INFO file (INFO section for CPMD format)

  iz = NINT(z)
  izeff = NINT(zeff)
  WRITE (UNIT=unit_info,FMT="(2(A,/,/),A,I3)")&
    " "//REPEAT("*",64),&
    " Atomic symbol                       : "//TRIM(elesym(iz)),&
    " Atomic number                       : ",iz
  IF (ABS(zeff - izeff) > eps_zero) THEN
    WRITE (UNIT=unit_info,FMT="(A,F8.3)")&
      " Effective core charge               : ",zeff
  ELSE
    WRITE (UNIT=unit_info,FMT="(A,I3)")&
      " Effective core charge               : ",izeff
  END IF

  ! Read second input file atom.dat (electronic configuration)

  READ (UNIT=unit_atom_dat,FMT="(A)") string
  IF (TRIM(ADJUSTL(string)) /= TRIM(elesym(iz))) THEN
    PRINT*,"ERROR: Mismatching element symbols found in the input files"
    STOP
  END IF
  READ (UNIT=unit_atom_dat,FMT="(A)") string
  IF (TRIM(ADJUSTL(string)) /= xc_string) THEN
    PRINT*,"ERROR: Mismatching XC functionals found in the input files"
    STOP
  END IF
  READ (UNIT=unit_atom_dat,FMT=*) line ! dummy read
  READ (UNIT=unit_atom_dat,FMT=*) line ! dummy read
  READ (UNIT=unit_atom_dat,FMT=*) line ! dummy read
  READ (UNIT=unit_atom_dat,FMT=*) ncore,nvalence
  WRITE (UNIT=unit_info,FMT="(A,I3,/,A,I3)")&
    " Number of core states               : ",ncore,&
    " Number of valence states            : ",nvalence

  ! Build electronic configuration

  frac_elec = .FALSE.
  elec_conf(:) = 0.0_wp

  DO i=1,nvalence
    READ (UNIT=unit_atom_dat,FMT=*) n,l,q
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

  WRITE (UNIT=unit_info,FMT="(A,A)")&
    " Electronic configuration (s,p,d,...): ",TRIM(elec_string)

  IF (ABS(zeff - nelec) > eps_zero) THEN
    WRITE (UNIT=unit_info,FMT="(A,F8.3)")&
      " Ionic charge                        : ",zeff - nelec
  END IF

  ! TeX tabular format (TEXTAB)

  WRITE (UNIT=unit_textab,FMT="(T2,A5,I3,A3,F10.6,4(A3,F13.6),A)")&
    elesym(iz)//" & ",izeff," & ",rloc,(" & ",cppl(i),i=1,4),"\\\\"

  DO ippnl=1,nppnl_max
    IF (nppnl(ippnl) > 0) THEN
      WRITE (UNIT=unit_textab,FMT="(T2,A11,F10.6,4(A3,F13.6),A)")&
        "   &     & ",rppnl(ippnl),(" & ",hppnl(1,j,ippnl),j=1,4),"\\\\"
      DO i=2,nppnl(ippnl)
        WRITE (UNIT=unit_textab,FMT="(T5,A1,T11,A1,T23,3(A3,F13.6),T72,A1,T87,A)")&
          "&","&",(" & ",hppnl(i,j,ippnl),j=1,3),"&","\\\\"
      END DO
    END IF
  END DO

  ! CP2K/Quickstep database format

  string = ""

  WRITE (UNIT=string,FMT="(I5)") izeff
  IF (xc_string == "PADE") THEN
  WRITE (UNIT=unit_cp2k,FMT="(A,1X,A)")&
    TRIM(elesym(iz)),"GTH-"//TRIM(xc_string)//"-q"//TRIM(ADJUSTL(string))//&
    " GTH-LDA-q"//TRIM(ADJUSTL(string))
  ELSE
    WRITE (UNIT=unit_cp2k,FMT="(A,1X,A)")&
      TRIM(elesym(iz)),"GTH-"//TRIM(xc_string)//"-q"//TRIM(ADJUSTL(string))
  END IF

  WRITE (UNIT=unit_cp2k,FMT="(A)") TRIM(elec_string)

  ! Set output precision for the database files

  idigits = 8

  fmtstr1 = "(F15. ,I5,4F15. )"
  WRITE (UNIT=fmtstr1(6:6),FMT="(I1)") idigits
  WRITE (UNIT=fmtstr1(16:16),FMT="(I1)") idigits
  fmtstr2 = "(T  ,4F15. )"
  WRITE (UNIT=fmtstr2(11:11),FMT="(I1)") idigits

  WRITE (UNIT=unit_cp2k,FMT=fmtstr1) rloc,nppl,(cppl(i),i=1,nppl)
  WRITE (UNIT=unit_cp2k,FMT="(I5)") nppnl_max
  DO ippnl=1,nppnl_max
    WRITE (UNIT=unit_cp2k,FMT=fmtstr1)&
      rppnl(ippnl),nppnl(ippnl),(hppnl(1,j,ippnl),j=1,nppnl(ippnl))
    DO i=2,nppnl(ippnl)
      WRITE (UNIT=fmtstr2(3:4),FMT="(I2)") 15*i + 6
      WRITE (UNIT=unit_cp2k,FMT=fmtstr2) (hppnl(i,j,ippnl),j=i,nppnl(ippnl))
    END DO
  END DO

  ! Write ABINIT database file

  WRITE (UNIT=unit_abinit,FMT="(A)")&
    "Goedecker pseudopotential for "//TRIM(elesym(iz))
  WRITE (UNIT=unit_abinit,FMT="(I5,I4,A)")&
    iz,izeff,"  070301 zatom,zion,pspdat"
  WRITE (UNIT=unit_abinit,FMT="(I2,2I3,I2,I5,I2,A)")&
    3,xc_code_abinit,lppnl_max,MAX(0,2*(nppl-1)),2001,0,&
    "  pspcod,pspxc,lmax,lloc,mmax,r2well"
  line = " "
  WRITE (UNIT=line,FMT=fmtstr1)&
    rloc,nppl,(cppl(i),i=1,nppl)
  WRITE (UNIT=line(70:),FMT="(A,4(A2,I1))")&
    "rloc nloc",(" c",i,i=1,nppl)
  WRITE (UNIT=unit_abinit,FMT="(A)") TRIM(line)
  WRITE (UNIT=unit_abinit,FMT="(I5,T70,A)")&
    nppnl_max,"nnonloc"
  DO ippnl=1,nppnl_max
    line = " "
    WRITE (UNIT=line,FMT=fmtstr1)&
      rppnl(ippnl),nppnl(ippnl),(hppnl(1,j,ippnl),j=1,nppnl(ippnl))
    WRITE (UNIT=line(70:),FMT="(A,3(A3,I2))")&
      "r"//llabel(ippnl)//" n"//llabel(ippnl),&
      (" h"//llabel(ippnl),10+j,j=1,nppnl(ippnl))
    WRITE (UNIT=unit_abinit,FMT="(A)") TRIM(line)
    DO i=2,nppnl(ippnl)
      WRITE (UNIT=fmtstr2(3:4),FMT="(I2)") 15*i + 6
      line = " "
      WRITE (UNIT=line,FMT=fmtstr2)&
        (hppnl(i,j,ippnl),j=i,nppnl(ippnl))
      WRITE (UNIT=line(70+5*i:),FMT="(3(A3,I2))")&
        (" h"//llabel(ippnl),10*i+j,j=i,nppnl(ippnl))
      WRITE (UNIT=unit_abinit,FMT="(A)") TRIM(line)
    END DO
    IF (ippnl > 1) THEN
      DO i=1,nppnl(ippnl)
        WRITE (UNIT=fmtstr2(3:4),FMT="(I2)") 15*i + 6
        line = " "
        WRITE (UNIT=line,FMT=fmtstr2)&
          (kppnl(i,j,ippnl),j=i,nppnl(ippnl))
        WRITE (UNIT=line(70+5*i:),FMT="(3(A3,I2))")&
          (" k"//llabel(ippnl),10*i+j,j=i,nppnl(ippnl))
        WRITE (UNIT=unit_abinit,FMT="(A)") TRIM(line)
      END DO
    END IF
  END DO

  ! Write the data set also to INFO file

  WRITE (UNIT=unit_info,FMT="(/,A)")&
    " Exchange-correlation functional     : "//TRIM(xc_string)

  ! Set output precision for the CPMD INFO section

  idigits = 6

  fmtstr3 = "(/,F13. ,4F13. )"
  WRITE (UNIT=fmtstr3(8:8),FMT="(I1)") idigits
  WRITE (UNIT=fmtstr3(15:15),FMT="(I1)") idigits

  WRITE (UNIT=unit_info,FMT="(/,T8,A)",ADVANCE="NO") "r(loc)"
  DO i=1,nppl
    WRITE (UNIT=unit_info,FMT="(9X,A2,I1,A1)",ADVANCE="NO") "C(",i,")"
  END DO
  WRITE (UNIT=unit_info,FMT=fmtstr3) rloc,(cppl(i),i=1,nppl)

  DO ippnl=1,nppnl_max
    IF (nppnl(ippnl) > 0) THEN
      IF (ippnl == 1) THEN
        WRITE (UNIT=unit_info,FMT="(/,9X,A2,I1,A1,5X,A,I1)",ADVANCE="NO")&
          "r(",ippnl-1,")","h(i,j)^",ippnl-1
      ELSE
        WRITE (UNIT=unit_info,FMT="(9X,A2,I1,A1,5X,A,I1)",ADVANCE="NO")&
          "r(",ippnl-1,")","h(i,j)^",ippnl-1
      END IF
      WRITE (UNIT=unit_info,FMT=fmtstr3)&
        rppnl(ippnl),(hppnl(1,j,ippnl),j=1,nppnl(ippnl))
      fmtstr2 = "(T  ,4F13. )"
      WRITE (UNIT=fmtstr2(11:11),FMT="(I1)") idigits
      DO i=2,nppnl(ippnl)
        WRITE (UNIT=fmtstr2(3:4),FMT="(I2)") 13*i + 1
        WRITE (UNIT=unit_info,FMT=fmtstr2) (hppnl(i,j,ippnl),j=i,nppnl(ippnl))
      END DO
    END IF
  END DO

  WRITE (UNIT=unit_info,FMT="(/,A,/,6(/,A),/,/,A)")&
    " Please cite:",&
    " - S. Goedecker, M. Teter, and J. Hutter,",&
    "   Phys. Rev. B 54, 1703 (1996)",&
    " - C. Hartwigsen, S. Goedecker, and J. Hutter,",&
    "   Phys. Rev. B 58, 3641 (1998)",&
    " - M. Krack,",&
    "   Theor. Chem. Acc. 114, 145 (2005)",&
    " "//REPEAT("*",64)

END PROGRAM gth_pp_convert
