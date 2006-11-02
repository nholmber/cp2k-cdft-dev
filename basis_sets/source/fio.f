C=======================================================================
C                                                                      =
C    F R E E - F O R M A T  F O R T R A N  I N P U T  R O U T I N E S  =
C                                                                      =
C                                                                      =
C MELDF - QUANTUM CHEMISTRY GROUP - INDIANA UNIVERSITY - AUGUST 21,1987=
C                                                                      =
C                                                                      =
C Subroutines:   -CFIELD                                               =
C                -INTGFL                                               =
C                -RDCARD                                               =
C                -REALFL                                               =
C                -STRTFL                                               =
C                                                                      =
C=======================================================================
C
C************************************
      SUBROUTINE CFIELD(STROUT,LIMIT)
C************************************
C
C       EXTRACT A CHARACTER FIELD AND PLACE IT IN STROUT
C
C    NOTE: FIELDS ARE TERMINATED BY END-OF-CARD, BLANK, COMMA,
C          =, OR * UNLESS STARTED BY *. FIELDS BEGINNNING WITH
C          "*" ARE TERMINATED ONLY BY ANOTHER "*" OR AN END-OF-CARD.
C
C  ON INPUT:
C  --------
C  BLANK     - .TRUE. IF A BLANK FIELD HAS BEEN DETECTED
C  STRING    - CHARACTER STRING IMAGE OF THE ENTIRE CARD
C  NCOL      - COLUMN POINTER TO THE START OF THE CURRENT FIELD
C  NCMAX     - COLUMN POINTER TO THE RIGHTMOST NONBLANK CHARACTER
C
C  INTERNAL:
C  --------
C  STARFL    - .TRUE. IF A "*" CHARACTER HAS BEEN DETECTED
C
C  ON EXIT:
C  -------
C  STROUT(*) - CHARACTER STRING TO BE RETURNED
C  LIMIT     - NUMBER CHARACTERS DESIRED IN STROUT
C  IV        - CHARACTER STRING CONTAINING ALL CHARACTERS BETWEEN THE
C              START OF THE FIELD AND THE NEXT * OR DELIMITER
C  LENGTH    - THE NUMBER OF NONBLANK CHARACTERS IN IV
C
      CHARACTER*200 STRING,IV
      CHARACTER STROUT*(*)
      LOGICAL STARFL,BLANK
      COMMON /CARDC/ IV
      LOGICAL LEOF
      COMMON /CARDI/ LENGTH,NCOL,NCMAX,NCOLP,NCMAXP,LEOF
      COMMON /IOFILE/IR,IW,IP
      SAVE 
      STARFL=.FALSE.
      LENGTH=0
C
C      --- LOCATE THE START OF THE NEXT NON-BLANK FIELD OR ---
C      --- READ A COMPLETE BLANK FIELD                     ---
C
      CALL STRTFL(BLANK,STRING)
C
C      --- TEST IF BLANK FIELD ---
C
      IF(BLANK) GO TO 115
C
C      --- TEST IF THE FIELD BEGINS WITH A STAR "*" ---
C
         IF(STRING(NCOL:NCOL).EQ.'*') THEN
            NCOL=NCOL+1
            STARFL=.TRUE.
         END IF
100      CONTINUE
C<<<<<<<<<<<<<<<<<<<<<<< BEGIN PROGRAM LOOP TO SCAN FOR THE END OF THE F
C                        EXIT ON ENCOUTERING A TERMINATOR
C
C      --- TEST FOR END OF CARD ---
C
            IF(NCOL.GT.NCMAX) GO TO 115
            IF(STARFL) THEN
C
C      --- FIELD STARTED BY "*" SO FIND TERMINAL QUOTE * ---
C
               IF(STRING(NCOL:NCOL).NE.'*') GO TO 110
               NCOL=NCOL+1
               GO TO 115
            ELSE
C
C      --- FIELD NOT STARTED BY "*", TEST FOR DELIMITER ---
C
               IF(STRING(NCOL:NCOL).EQ.' ') GO TO 115
               IF(STRING(NCOL:NCOL).EQ.',') GO TO 115
               IF(STRING(NCOL:NCOL).EQ.'=') GO TO 115
               IF(STRING(NCOL:NCOL).EQ.'$') GO TO 115
               IF(STRING(NCOL:NCOL).EQ.'*') GO TO 115
            END IF
110         CONTINUE
C
C      --- STORE THE CHARACTER IN THE SCRATCH STRING IV            ---
C      --- IT WILL CONTAIN ALL CHARACTERS UP TO THE NEXT DELIMITER ---
C
            LENGTH=LENGTH+1
            IV(LENGTH:LENGTH)=STRING(NCOL:NCOL)
            NCOL=NCOL+1
            GO TO 100
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END PROGRAM LOOP 100
115   CONTINUE
C
C      --- IF MORE CHARACTERS WERE ASKED FOR THAN WERE AVAILABLE ---
C      --- FILL WITH BLANKS AND THEN PACK INTO STRING IV         ---
C
      LL=LIMIT-LENGTH
      IF(LL.GT.0) THEN
         DO 120 I=1,LL
            LIM=LIMIT-I+1
            IV(LIM:LIM)=' '
120      CONTINUE
      END IF
C
C      --- PACK CHARACTERS INTO STROUT ---
C
      DO 130 I=1,LIMIT
         STROUT(I:I) = IV(I:I)
130   CONTINUE
      RETURN
C...END OF CFIELD
      END
C***************************************
      SUBROUTINE INTGFL(INTGER,ERROR,IX)
C***************************************
      USE basic_data_types, ONLY: dp
      IMPLICIT REAL(dp) (A-H,O-Z)
C
C     EXTRACT AN INTEGER FIELD
C
C  ON ENTRY:
C  --------
C
C  INTERNAL:
C  --------
C  FPTNUM  - A FLOATING POINT NUMBER WHICH WILL TRANSLATED INTO AN INTEG
C
C  ON EXIT:
C  --------
C  IV      - THE INTEGER TO BE RETURNED
C  ERROR   - .TRUE. IF A NONVALID INTEGER WAS DETECTED
C  IX      - A CHARACTER VARIABLE CONTAINING THE FIRST 4 CHARACTERS OF
C            THE STRING WHICH WAS SUPPOSED TO HAVE BEEN AN INTEGER
C
      CHARACTER*4 IX
      LOGICAL ERROR
      CALL REALFL(FPTNUM,ERROR,IX)
      INTGER = FPTNUM
      X      = INTGER
C
C      --- TEST IF FIELD CONTAINS AN INTEGER ---
C
      IF(X.NE.FPTNUM) ERROR=.TRUE.
      RETURN
C...END OF INTGFL
      END
C******************************
      SUBROUTINE RDCARD(STRING)
C******************************
C
C            READ THE NEXT PHYSICAL CARD.
C
C  ON ENTRY:
C  --------
C  NOTHING
C
C  ON EXIT:
C  -------
C  STRING  - STRING CONTAINING THE CARD IMAGE JUST READ
C  NCOL    - CURRENT COLUMN POINTER TO THE START OF A FIELD. SET TO 1 ON
C  NCOLP   - PREVIOUS COLUMN POINTER. SET TO 1.
C  NCMAX   - THE COLUMN NUMBER OF THE LAST NONBLANK CHARACTER FROM THE L
C            (IE. THIS IS THE END-OF-CARD MARKER)
C
      CHARACTER*200 STRING
      LOGICAL LEOF
      COMMON /CARDI/ LENGTH,NCOL,NCMAX,NCOLP,NCMAXP,LEOF
      COMMON /IOFILE/IR,IW,IP
      SAVE
C
C      --- READ A NEW PHYSICAL CARD OF (UP TO) 200 CHARACTERS ---
C
      J   = 0
      LEOF=.TRUE.
      READ(IR,900,END=999,ERR=999) STRING
      LEOF=.FALSE.
C
C      --- BEGIN SCANNING FROM RIGHT TO LEFT FOR THE FIRST NONBLANK CHAR
C      --- SET NCMAX = THE COLUMN CONTAINING THE LAST NONBLANK CHARACTER
C
      J=1
      DO 100 I=1,200
         N = 200-I+1
         IF(STRING(N:N).NE.' ') GO TO 105
100   CONTINUE
      I = 200
      J = 0
105   CONTINUE
      NCMAX = 200+1-I
C
C      --- CONVERT ANY LOWER CASE CHARACTERS TO UPPER CASE ---
C
      DO 120 I=1,NCMAX
         IF (STRING(I:I) .EQ. 'a') THEN
            STRING(I:I) = 'A'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'b') THEN
            STRING(I:I) = 'B'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'c') THEN
            STRING(I:I) = 'C'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'd') THEN
            STRING(I:I) = 'D'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'e') THEN
            STRING(I:I) = 'E'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'f') THEN
            STRING(I:I) = 'F'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'g') THEN
            STRING(I:I) = 'G'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'h') THEN
            STRING(I:I) = 'H'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'i') THEN
            STRING(I:I) = 'I'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'j') THEN
            STRING(I:I) = 'J'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'k') THEN
            STRING(I:I) = 'K'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'l') THEN
            STRING(I:I) = 'L'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'm') THEN
            STRING(I:I) = 'M'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'n') THEN
            STRING(I:I) = 'N'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'o') THEN
            STRING(I:I) = 'O'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'p') THEN
            STRING(I:I) = 'P'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'q') THEN
            STRING(I:I) = 'Q'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'r') THEN
            STRING(I:I) = 'R'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 's') THEN
            STRING(I:I) = 'S'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 't') THEN
            STRING(I:I) = 'T'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'u') THEN
            STRING(I:I) = 'U'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'v') THEN
            STRING(I:I) = 'V'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'w') THEN
            STRING(I:I) = 'W'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'x') THEN
            STRING(I:I) = 'X'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'y') THEN
            STRING(I:I) = 'Y'
            GO TO 120
         END IF
         IF (STRING(I:I) .EQ. 'z') THEN
            STRING(I:I) = 'Z'
            GO TO 120
         END IF
120   CONTINUE
C
C      --- PRINT AN ECHO OF THE CARD IMAGE ON THE OUTPUT FILE ---
C
C     WRITE(0,950) (STRING(I:I),I=1,NCMAX)
C
C      --- RESET COLUMN MARKER AND PREVIOUS COLUMN MARKER TO 1 ---
C
 999  NCOL  = 1
      NCOLP = 1
C
C      --- TREAT A BLANK CARD AS TWO DELIMITERS IN A ROW BY ---
C      --- PLACING TWO ASTERISKS IN COLUMNS 1 AND 2         ---
C
      IF(J.EQ.0) THEN
         NCMAX=2
         STRING(1:1) = '*'
         STRING(2:2) = '*'
      END IF
900   FORMAT(A200)
C950  FORMAT(1X,')',200A1)
C...END OF RDCARD
      RETURN
      END
C***************************************
      SUBROUTINE REALFL(FPTNUM,ERROR,IX)
C***************************************
      USE basic_data_types, ONLY: dp
      IMPLICIT REAL(dp) (A-H,O-Z)
C
C     READ A FULL PRECISION (64 BIT) FLOATING POINT NUMBER.
C
C  NOTE: ACCEPTED FORM IS  
C        (SIGN)___.___E(SIGN)__/(SIGN)___.___E(SIGN)__X__R
C        WHERE R IS A REPITITION FACTOR. THE EXPONENT PORTION IS OPTIONAL
C        AS IS THE DECIMAL POINT. A FIELD STARTED BY E IS ASSUMED TO HAVE
C        MANTISSA OF 1.
C
C  ON ENTRY:
C  --------
C
C  INTERNAL:
C  --------
C    POINT   - LOGICAL VARIABLE
C    EXP     - LOGICAL VARIABLE
C    REP     - .TRUE. IF THE REPITITION CHARACTER "X" IS DETECTED
C    IRC     - THE REPITITION COUNT
C
C  ON EXIT:
C  -------
C   FPTNUM - WORD TO CONTAINING THE FLOATING POINT NUMBER
C   ERROR  - LOGICAL VARIABLE, = .TRUE. IF ERROR ENCOUNTERED
C   IX     - CHARACTER VARIABLE, IF ERROR ENCOUNTERED THIS WILL CONTAIN
C            THE FIRST FOUR CHARACTERS OF THE FIELD WHICH WAS SUPPOSED
C            TO BE INTERPRETED AS A FLOATING POINT NUMBER
C
      CHARACTER*200 IV
      CHARACTER*4  IX
      LOGICAL POINT,EXP,ERROR,REP
      COMMON /CARDC/ IV
      LOGICAL LEOF
      COMMON /CARDI/ LENGTH,NCOL,NCMAX,NCOLP,NCMAXP,LEOF
      COMMON /IOFILE/IR,IW,IP
      DATA  IRC/0/,ZERO,ONE,TEN/0.D0,1.D0,10.D0/,IFOUR/4/
      SAVE
C
      ERROR=.FALSE.
      IF(IRC.GT.0) GO TO 200
      INV    = 0
      FPTNUM = ZERO
C
C      --- OBTAIN A COMPLETE FIELD IN STRING IV ---
C
      CALL CFIELD(IX,IFOUR)
      IF(LENGTH.EQ.0) RETURN
      V1=ZERO
      SIGN=ONE
      SCALE=ONE
      EXP   = .FALSE.
      REP   = .FALSE.
      POINT = .FALSE.
      NCH   = 15
      ISEXP = 1
      NEXP  = 0
      DO 150 J=1,LENGTH
C
C      --- IDENTIFY THE J-TH CHARACTER ---
C
         IINDEX=0
         IF(IV(J:J).EQ.'0') GO TO 105
         IINDEX=1
         IF(IV(J:J).EQ.'1') GO TO 105
         IINDEX=2
         IF(IV(J:J).EQ.'2') GO TO 105
         IINDEX=3
         IF(IV(J:J).EQ.'3') GO TO 105
         IINDEX=4
         IF(IV(J:J).EQ.'4') GO TO 105
         IINDEX=5
         IF(IV(J:J).EQ.'5') GO TO 105
         IINDEX=6
         IF(IV(J:J).EQ.'6') GO TO 105
         IINDEX=7
         IF(IV(J:J).EQ.'7') GO TO 105
         IINDEX=8
         IF(IV(J:J).EQ.'8') GO TO 105
         IINDEX=9
         IF(IV(J:J).EQ.'9') GO TO 105
         IF(REP) GO TO 100
         IINDEX=10
         IF(IV(J:J).EQ.'/') GO TO 105
         IINDEX=11
         IF(IV(J:J).EQ.'.') GO TO 105
         IINDEX=12
         IF(IV(J:J).EQ.'+') GO TO 105
         IINDEX=13
         IF(IV(J:J).EQ.'-') GO TO 105
         IINDEX=14
         IF(IV(J:J).EQ.'E') GO TO 105
         IINDEX=15
         IF(IV(J:J).EQ.'D') GO TO 105
100      CONTINUE
         IF(IV(J:J).EQ.'X') GO TO 175
C
C      --- ILLEGAL CHARACTER FOUND ---
C
            ERROR=.TRUE.
            RETURN
105      CONTINUE
C
C      --- TEST IF THIS CHARACTER IS A NON-NUMERAL ---
C
         IF(IINDEX.GT.9) GO TO 115
C
C      --- TEST IF PART OF REPETITION ---
C
         IF(REP) GO TO 180
C
C      --- TEST IF PART OF EXPONENT ---
C
            IF(EXP) GO TO 110
C
C      --- ADD DIGIT TO MANTISSA ---
C
               INV=INV+1
               IF(POINT) THEN
               SCALE=SCALE/TEN
               FPTNUM=FPTNUM+FLOAT(IINDEX)*SCALE
               ELSE
               FPTNUM=TEN*FPTNUM+FLOAT(IINDEX)
               END IF
            GO TO 150
110         CONTINUE
C
C      --- ADD DIGIT TO EXPONENT ---
C
               NEXP=10*NEXP+IINDEX
            GO TO 150
115      CONTINUE
C
C      --- PROCESS NON-NUMERALS CHARACTERS ---
         JJ=IINDEX-9
         GO TO (120,125,130,135,145,145),JJ
120      CONTINUE
C
C      --- SLASH DETECTED (FIELD EXPRESSED AS A FRACTION) NUMERATOR COMP
C
      V1=FPTNUM*SIGN*(TEN**(ISEXP*NEXP))
      FPTNUM=ZERO
      INV=0
      IF(V1.EQ.ZERO) GO TO 155
            SIGN=ONE
            ISEXP=1
            EXP=.FALSE.
            SCALE=ONE
            POINT=.FALSE.
            NEXP=0
         GO TO 150
125      CONTINUE
C
C      --- DECIMAL POINT DETECTED ---
C
            POINT=.TRUE.
         GO TO 150
130      CONTINUE
C
C      --- PLUS SIGN DETECTED, TEST IF START OF MANTISSA OR EXPONENT ---
C
            IF(FPTNUM) 145,150,145
135      CONTINUE
C
C      --- MINUS SIGN, TEST IF START OF MANTISSA OF EXPONENT ---
C
            IF(FPTNUM.EQ.ZERO) GO TO 140
               ISEXP=-1
            GO TO 145
140         CONTINUE
               SIGN=-ONE
            GO TO 150
145      CONTINUE
C
C      --- E, D OR EMBEDDED + OR - STARTS EXPONENT FIELD ---
C
            IF(INV.EQ.0) FPTNUM=ONE
            EXP=.TRUE.
            GO TO 150
175         CONTINUE
            REP=.TRUE.
            GO TO 150
180         CONTINUE
            IRC=10*IRC+IINDEX
            GO TO 150
150   CONTINUE
      FPTNUM=FPTNUM*SIGN*(TEN**(ISEXP*NEXP))
C
C      --- THE NUMERATOR IS FINISHED, TEST FOR NO DENOMINATOR ---
C
      IF(V1.EQ.ZERO) GO TO 155
      IF(FPTNUM.EQ.ZERO) FPTNUM=ONE
      FPTNUM=V1/FPTNUM
155   CONTINUE
      V1=FPTNUM
      IF((REP).AND.(IRC.EQ.0)) ERROR=.TRUE.
      IF(IRC.NE.0) IRC=IRC-1
      RETURN
200   CONTINUE
      ERROR=.FALSE.
      IRC=IRC-1
      FPTNUM=V1
      RETURN
C...END OF REALFL
      END
C************************************
      SUBROUTINE STRTFL(BLANK,STRING)
C************************************
C
C              SCAN FOR THE START OF A FIELD
C
C    NOTE: THIS ROUTINE CONSIDERS A STRING OF BLANKS (WITH OR WITHOUT AN
C          EMBEDDED COMMA OR EQUALS SIGN AS A DELIMITER. A DELIMITER FOL
C          BY A COMMENT FIELD $...$ IS TREATED AS A BLANK STRING. COMMEN
C          FIELDS MAY ALSO BE TERMINATED BY AN END-OF-CARD. THE "*" IN A
C          *...* CHARACTER FIELD IS TREATED AS A BLANK DELIMITER. EMBEDD
C          BLANKS IN A *..* FIELD ARE TREATED AS TEXT.
C
C  ON ENTRY:
C  --------
C
C  ON EXIT:
C  -------
C   BLANK  - .TRUE. IF A BLANK FIELD HAS BEEN DETECTED.
C   STRING - A CHARACTER STRING CONTAINING THE CONTENTS OF CARD
C   NCOL   - COLUMN POINTER TO THE START OF THE CURRENT FIELD
C   NCOLP  - COLUMN POINTER TO THE START OF PREVIOUS FIELD
C   NCMAX  - POINTER TO THE COLUMN CONTAINING THE RIGHTMOST NONBLANK CHA
C   NCMAXP - POINTER TO THE END OF THE PREVIOUS CARD
C
      LOGICAL BLANK,DELIM,FIRST
      CHARACTER*200 STRING
      LOGICAL LEOF
      COMMON /CARDI/ LENGTH,NCOL,NCMAX,NCOLP,NCMAXP,LEOF
      COMMON /IOFILE/IR,IW,IP
C
      DATA FIRST/.TRUE./
      SAVE
C
C      --- INITIALIZE VARIABLES ON FIRST CALL TO THIS ROUTINE ---
C
      IF(FIRST) THEN
         FIRST=.FALSE.
         NCOL  = 1
         NCMAX = 0
      END IF
C
C      --- RETAIN START OF PREVIOUS FIELD IN NCOLP ---
C
      NCOLP  = NCOL
      NCMAXP = NCMAX
      DELIM=.FALSE.
      BLANK=.TRUE.
100   CONTINUE
C<<<<< BEGIN PROGRAM LOOP SCANNING FROM LEFT TO RIGHT FOR A NONBLANK CHA
C   EXIT ON ENCOUTERING EITHER A COMPLETE BLANK FIELD OR THE START OF NO
C
C      --- READ ANOTHER PHYSICAL CARD INTO "STRING" IF NECESSARY ---
C
         IF(NCOL.GT.NCMAX) CALL RDCARD(STRING)
         IF(STRING(NCOL:NCOL).EQ.' ') GO TO 150
C
C      --- TEST FOR COMMENT FIELD ---
C
            IF(STRING(NCOL:NCOL).EQ.'$') GO TO 115
               IF(DELIM) GO TO 105
C
C      --- TEST FOR FIRST EQUALS SIGN OR COMMA ---
C
                  DELIM=.TRUE.
                  IF(STRING(NCOL:NCOL).EQ.'=') GO TO 150
                  IF(STRING(NCOL:NCOL).EQ.',') GO TO 150
                  GO TO 110
105            CONTINUE
C
C      --- TEST FOR SECOND = OR , INDICATING BLANK FIELD (EXIT IF YES)--
C
               IF(STRING(NCOL:NCOL).EQ.'=') GO TO 160
               IF(STRING(NCOL:NCOL).EQ.',') GO TO 160
110            CONTINUE
C
C      --- A NON-BLANK FIELD HAS BEEN FOUND (EXIT)---
C
               BLANK=.FALSE.
               GO TO 160
115         CONTINUE
C
C      --- COMMENT FIELD DETECTED (STARTED BY $) ---
C
               DELIM=.FALSE.
               NCOL=NCOL+1
140            CONTINUE
C<<<<<<<<<<<<<< BEGIN PROGRAM LOOP TO SKIP TO THE END OF THE COMMENT FIE
C               EXIT ON ENCOUNTERING THE SECOND $ OR AN END-OF-CARD
               IF(NCOL.GT.NCMAX)            GO TO 150
               IF(STRING(NCOL:NCOL).EQ.'$') GO TO 150
               NCOL=NCOL+1
               GO TO 140
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END PROGRAM LOOP 1
150   CONTINUE
C
C      --- INCREMENT COLUMN POINTER AND GO ON TO THE NEXT CHARACTER ---
C
      NCOL=NCOL+1
      GO TO 100
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END PROGRAM LOOP 1
160   CONTINUE
      RETURN
C...END OF STRTFL
      END
C     ==================================================================
      FUNCTION INSCAN(IUNIT,LABEL)
      USE basic_data_types, ONLY: dp
      IMPLICIT REAL(dp) (A-H,O-Z)
      COMMON /IOFILE/IR,IW,IP
      CHARACTER LABEL*(*),LINE*200
      LOGICAL EXISTS
C     ==--------------------------------------------------------------==
C     ==  Scans unit 'iunit' for the section header 'label'           ==
C     ==--------------------------------------------------------------==
      INQUIRE(UNIT=IUNIT,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
        WRITE(0,'(A,I3,A)') ' INSCAN: Unit',iunit,' does not exist'
        STOP
      ENDIF
      INSCAN=1
      REWIND(IUNIT)
      IR=IUNIT
      IW=6
      IP=6
C     ==--------------------------------------------------------------==
   10 CONTINUE
      READ(IUNIT,END=20,ERR=20,FMT='(A200)') LINE
      IFIND=INDEX(LINE,LABEL)
      IF(IFIND.NE.0) THEN
        INSCAN=0
        GOTO 20
      ENDIF
      GOTO 10
   20 CONTINUE
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

C************************************
      SUBROUTINE NEWLINE
C************************************
C
C    Jump to the beginning of a new line and forget
C    the rest of the old line (trial version by G.Lippert)
C  -------
C   NCOL   - COLUMN POINTER TO THE START OF THE CURRENT FIELD
C   NCMAX  - POINTER TO THE COLUMN CONTAINING THE RIGHTMOST NONBLANK CHA
C
      COMMON /CARDI/ LENGTH,NCOL,NCMAX,NCOLP,NCMAXP,LEOF
      NCOL=NCMAX+1
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================

