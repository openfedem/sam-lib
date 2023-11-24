C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE FIIOPN (NCHLIN,NCMTCH,LUNRD,LUNER,LUNECO,IERR)
C
C **********************************************************************
C
C   S A M  library routine :  FIIOPN                  GROUP 1 / PUBLIC
C
C     TASK :  To initialize key parameters and character strings used
C             by the FII-routines, i.e., to "open" the FII-package
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-22 / 1.0
C                       90-01-24 / 1.1   K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           IERR,LUNECO,LUNER,LUNRD,NCHLIN,NCMTCH
C
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           I
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IERR = 0
C                                               ! input parameters
C
      IF (NCHLIN.LT.1 .OR. NCHLIN.GT.132)  GO TO 90
      MAXCHL = NCHLIN
      IF (NCMTCH.LT.1)  GO TO 90
      NCHREC = NCMTCH
      IF (LUNRD.LT.0)   GO TO 90
      LRU = LUNRD
      LPU = LUNER
      LEU = LUNECO
C                                               ! initiate alphabets
      CHALPL = 'abcdefghijklmnopqrstuvwxyz'
      CHALPU = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
C                                               ! defaults
      KLCV   = 1
      KEOL   = 0
      CHCOM1 = '*'
      CHCOM2 = '!'
      CHSEP  = ','
C                                               ! zero and blank
      NEOL   = 0
      IER    = 0
      LCNT   = 0
      LDGT   = 0
      LEXP   = 0
      LUIC   = 0
      LNUM   = 0
      NCHNAM = 0
      DPVAL  = 0.0D0
      CHBLK  = ' '
      CHLET  = CHBLK
      DO 10 I=1,16
         CHNAM(I:I) = CHBLK
   10 CONTINUE
C
      GO TO 100
C ----------------------------------------------- error exit
   90 IERR =-1
      IF (LPU.EQ.0) THEN
         WRITE ( * ,600) NCHLIN,NCMTCH,LUNRD
      ELSEIF (LPU.GT.0) THEN
         WRITE (LPU,600) NCHLIN,NCMTCH,LUNRD
      ENDIF
C -----------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  600 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE FIIOPN'/
     +            5X,'ILLEGAL PARAMETER(S) :  NCHLIN =',I5 /
     +           29X,                        'NCMTCH =',I5 /
     +           29X,                        'LUNRD  =',I5 )
C
      END
      SUBROUTINE FIIRLS
C
C **********************************************************************
C
C   S A M  library routine :  FIIRLS                  GROUP 1 / PUBLIC
C
C     TASK :  The same as FIIRLN, i.e., to read a new data line.
C             The only difference is that this subroutine always goes
C             to STOP on error and on attempted read through an end of
C             file.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIIRLN      (SAM - 1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-12-20 / 1.0
C                       90-01-24 / 1.1    K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C                                               ! Local variable
      INTEGER           KOPY
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIIRLN
C ----------------------------------------------------------------------
C
      KOPY = KEOL
      KEOL =-1
      CALL FIIRLN
      IF (IER .LT. 0)   STOP
      KEOL = KOPY
C
      RETURN
      END
      LOGICAL FUNCTION LNMFII()
C
C **********************************************************************
C
C   S A M  library routine :  LNMFII                  GROUP 1 / PUBLIC
C
C     TASK : To look for a valid name in the next character positions
C            of the current input line, and, if found, place the name
C            in the leftmost positions of CHNAM
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIISKP, LLTFII,
C                                   LDGFII and LDLFII        (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           L,LD,LP,LLPP
      CHARACTER         CHLL*1,CHLN*16
      LOGICAL           LDGFII,LLTFII,LDLFII
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIISKP,LDGFII,LLTFII,LDLFII
C ----------------------------------------------------------------------
      LNMFII =.FALSE.
      IF (LPP.GT.NCHL) GO TO 100
      CALL FIISKP
      LLPP = LPP
C
      LP   = 1
      LD   = LDGT
      CHLL = CHLET
      IF (LLTFII()) THEN
         CHLN(LP:LP) = CHLET
   10    LP = LP+1
         IF (LP.GT.16) THEN
            IF (LDLFII()) THEN
C                                               ! exit - true
               LNMFII =.TRUE.
               CHNAM  = CHLN
               NCHNAM = 16
               GO TO 100
            ELSE
C                                               ! exit - false
               LPP   = LLPP
               CHLET = CHLL
               LDGT  = LD
               GO TO 100
            ENDIF
         ENDIF
         IF (LLTFII()) THEN
            CHLN(LP:LP) = CHLET
         ELSEIF (LDGFII())THEN
            L = LPP-1
            CHLN(LP:LP) = CHLCOP(L:L)
         ELSEIF (LDLFII()) THEN
C                                               ! exit - true
            LNMFII =.TRUE.
            NCHNAM = LP-1
            CHNAM(1:NCHNAM) = CHLN(1:NCHNAM)
            GO TO 100
         ELSE
C                                               ! exit - false
            LPP   = LLPP
            CHLET = CHLL
            LDGT  = LD
            GO TO 100
         ENDIF
C                                               ! continue loop
         GO TO 10
      ENDIF
C
  100 RETURN
      END
      LOGICAL FUNCTION LFLFII()
C
C **********************************************************************
C
C   S A M  library routine :  LFLFII                  GROUP 1 / PUBLIC
C
C     TASK : To look for a valid floating point number in the next
C            positions in the current input line, and, if found, to
C            decode/compute its value and store it in DPVAL.
C            This function has been modified to also accept simple ex-
C            pressions, involving floating point constants, trigono-
C            metric functions sine and cosine and the operators
C            +, - , * , / and **.  Arguments of sine and cosine must be
C            enclosed in parentheses, which is the only legal use of
C            parentheses.
C            The expression, which must not contain blank characters,
C            is evaluated from left to right.
C
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  LNXFII, LXXFII, LDLFII   (SAM-1)
C                                   LSCFII, LUFFII           (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C                       87-12-07 / 1.1    K.Bell      (LXXFII)
C                       90-01-28 / 2.0    K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
C                                               ! Local variables
      INTEGER           ISIGN,LD,LE,LLPP,LU
      LOGICAL           LNXFII,LSCFII,LDLFII,LUFFII,LXXFII
      DOUBLE PRECISION  DPC,DPV
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LNXFII,LSCFII,LDLFII,LUFFII,LXXFII
C ----------------------------------------------------------------------
      LFLFII =.FALSE.
      IF (LPP.GT.NCHL)  GO TO 200
      LLPP = LPP
      LU   = LUIC
      LD   = LDGT
      LE   = LEXP
      DPC  = DPVAL
C
      ISIGN = 1
      IF (LNXFII('+')) THEN
         ISIGN = 1
      ELSEIF (LNXFII('-')) THEN
         ISIGN =-1
      ENDIF
C
      IF (LUFFII() .OR. LSCFII()) THEN
         DPV = DPVAL
      ELSE
         GO TO 90
      ENDIF
C
   10 IF (LDLFII()) THEN
         IF (ISIGN .LT. 0)  DPV =-DPV
         DPVAL  = DPV
         LFLFII =.TRUE.
         GO TO 100
      ELSEIF (LXXFII('+')) THEN
         IF (LUFFII() .OR. LSCFII()) THEN
            DPV = DPV + DPVAL
         ELSE
            GO TO 90
         ENDIF
      ELSEIF (LXXFII('-')) THEN
         IF (LUFFII() .OR. LSCFII()) THEN
            DPV = DPV - DPVAL
         ELSE
            GO TO 90
         ENDIF
      ELSEIF (LXXFII('*')) THEN
         IF (LXXFII('*')) THEN
            IF (LUFFII() .OR. LSCFII()) THEN
               DPV = DPV**DPVAL
            ELSE
               GO TO 90
            ENDIF
         ELSEIF (LUFFII() .OR. LSCFII()) THEN
            DPV = DPV*DPVAL
         ELSE
            GO TO 90
         ENDIF
      ELSEIF (LXXFII('/')) THEN
         IF (LUFFII() .OR. LSCFII()) THEN
            IF (DPVAL .NE. 0.0D0) THEN
               DPV = DPV/DPVAL
            ELSE
               DPV = 0.0D0
            ENDIF
         ELSE
            GO TO 90
         ENDIF
      ELSE
         GO TO 90
      ENDIF
C
      GO TO 10
C
   90 LPP   = LLPP
      DPVAL = DPC
C
  100 LUIC  = LU
      LDGT  = LD
      LEXP  = LE
C
  200 RETURN
      END
      LOGICAL FUNCTION LSCFII()
C
C **********************************************************************
C
C   S A M  library routine :  LSCFII                 GROUP 1 / PUBLIC
C
C     TASK : To look for an unsigned trigonometric function expression
C            in the next positions in the current input line, and, if
C            found, to decode/compute its value and store it in DPVAL.
C            Only sine (SIN) and cosine (COS) are considered, and the
C            argument, which may be a floating point expression of con-
C            stants and operators (with no blanks), must appear within
C            parentheses.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  LUFFII, LXXFII   (SAM-1)
C                                   SIN,    COS      (Fortran library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-28 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
C
      INTEGER           ISIGN,K,LD,LE,LLPP,LU,NTMP
      DOUBLE PRECISION  DPC,DPV
      LOGICAL           LUFFII,LXXFII
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LUFFII,LXXFII
C ----------------------------------------------------------------------
      LSCFII =.FALSE.
      IF (LPP.GT.NCHL)  GO TO 200
      LLPP   = LPP
      NTMP   = NCHREC
      NCHREC = 4
C
      IF (LXXFII('SIN(')) THEN
         K = 1
      ELSEIF (LXXFII('COS(')) THEN
         K = 2
      ELSE
         GO TO 150
      ENDIF
C
      LU  = LUIC
      LD  = LDGT
      LE  = LEXP
      DPC = DPVAL
C
      ISIGN = 1
      IF (LXXFII('+')) THEN
         ISIGN = 1
      ELSEIF (LXXFII('-')) THEN
         ISIGN =-1
      ENDIF
C
      IF (LUFFII()) THEN
         DPV = DPVAL
   10    IF (LXXFII(')')) THEN
            IF (ISIGN .LT. 0)  DPV   =-DPV
            IF (K .EQ. 1)      DPVAL = SIN(DPV)
            IF (K .EQ. 2)      DPVAL = COS(DPV)
            LSCFII =.TRUE.
            GO TO 100
         ELSEIF (LXXFII('+')) THEN
            IF (LUFFII()) THEN
               DPV = DPV + DPVAL
            ELSE
               GO TO 90
            ENDIF
         ELSEIF (LXXFII('-')) THEN
            IF (LUFFII()) THEN
               DPV = DPV - DPVAL
            ELSE
               GO TO 90
            ENDIF
         ELSEIF (LXXFII('*')) THEN
            IF (LXXFII('*')) THEN
               IF (LUFFII()) THEN
                  DPV = DPV**DPVAL
               ELSE
                  GO TO 90
               ENDIF
            ELSEIF (LUFFII()) THEN
               DPV = DPV*DPVAL
            ELSE
               GO TO 90
            ENDIF
         ELSEIF (LXXFII('/')) THEN
            IF (LUFFII()) THEN
               IF (DPVAL .NE. 0.0D0) THEN
                  DPV = DPV/DPVAL
               ELSE
                  DPV = 0.0D0
               ENDIF
            ELSE
               GO TO 90
            ENDIF
         ELSE
            GO TO 90
         ENDIF
C
         GO TO 10
      ELSE
         GO TO 90
      ENDIF
C
   90 LPP    = LLPP
      DPVAL  = DPC
  100 LUIC   = LU
      LDGT   = LD
      LEXP   = LE
  150 NCHREC = NTMP
C
  200 RETURN
      END
      LOGICAL FUNCTION LUFFII()
C
C **********************************************************************
C
C   S A M  library routine :  LUFFII                  GROUP 1 / PUBLIC
C
C     TASK : To look for an unsigned floating point constant in the next
C            positions in the current input line, and, if found, to de-
C            code its value and store it in DPVAL.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  LICFII, LDGFII, LXPFII   (SAM-1)
C                                   LXXFII                   (SAM-1)
C
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-27 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
C                                               ! Local variables
      INTEGER           LD,LE,LLPP,LU,NTMP
      DOUBLE PRECISION  DPF,DPI,DPT,DPV,DPX
      LOGICAL           LDGFII,LICFII,LXPFII,LXXFII
C
      PARAMETER         ( DPI = 3.141592653589793D0 )
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LDGFII,LICFII,LXPFII,LXXFII
C ----------------------------------------------------------------------
      LUFFII =.FALSE.
      LLPP   = LPP
      LU     = LUIC
      LD     = LDGT
      LE     = LEXP
      NTMP   = NCHREC
      NCHREC = 2
C
      DPV = 0.0D0
      DPF = 0.0D0
      DPX = 1.0D0
      DPT = 1.0D1
C
      IF (LICFII()) THEN
         DPV = LUIC
         IF (LXXFII('.')) THEN
            IF (LDGFII()) THEN
   10          DPX = DPX/DPT
               DPF = DPF + LDGT*DPX
               IF (LDGFII())  GO TO 10
               DPV = DPV + DPF
               IF (LXPFII())  DPV = DPV*(DPT**LEXP)
            ELSEIF (LXPFII()) THEN
               DPV = DPV*(DPT**LEXP)
            ENDIF
         ELSEIF (LXPFII()) THEN
            DPV = DPV*(DPT**LEXP)
         ENDIF
         LUFFII =.TRUE.
         DPVAL  = DPV
      ELSEIF (LXXFII('.')) THEN
         IF (LDGFII()) THEN
   20       DPX = DPX/DPT
            DPF = DPF + LDGT*DPX
            IF (LDGFII())  GO TO 20
            DPV = DPV + DPF
            IF (LXPFII())  DPV = DPV*(DPT**LEXP)
            LUFFII =.TRUE.
            DPVAL  = DPV
         ELSE
            LPP = LLPP
         ENDIF
      ELSEIF (LXXFII('PI')) THEN
         LUFFII =.TRUE.
         DPVAL  = DPI
      ENDIF
C
      LUIC   = LU
      LDGT   = LD
      LEXP   = LE
      NCHREC = NTMP
C
      RETURN
      END
      LOGICAL FUNCTION LINFII()
C
C **********************************************************************
C
C   S A M  library routine :  LINFII                  GROUP 1 / PUBLIC
C
C     TASK :  To look for a valid integer number in the next positions
C             in the current input line, and, if found, to decode its
C             value and store it in LNUM
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  LNXFII, LDLFII, LICFII    (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           ISIGN,LLPP,LP
      LOGICAL           LNXFII,LDLFII,LICFII
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LNXFII,LDLFII,LICFII
C ----------------------------------------------------------------------
      LINFII =.FALSE.
      IF (LPP.GT.NCHL)  GO TO 100
      LLPP  = LPP
C
      ISIGN = 1
      IF (LNXFII('+')) THEN
         ISIGN = 1
      ELSEIF (LNXFII('-')) THEN
         ISIGN =-1
      ENDIF
      LP = LUIC
      IF (LICFII()) THEN
         IF (LDLFII()) THEN
            LINFII =.TRUE.
            LNUM   = ISIGN*LUIC
         ELSE
            LUIC  = LP
            LPP   = LLPP
         ENDIF
      ELSE
         LPP = LLPP
      ENDIF
C
  100 RETURN
      END
      LOGICAL FUNCTION LICFII()
C
C **********************************************************************
C
C   S A M  library routine :  LICFII                  GROUP 1 / PUBLIC
C
C     TASK :  To look for an unsigned integer constant (uic) in the next
C             character position(s) in the current input line, and if
C             found, to decode its value and store it in  LUIC
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :   LDGFII         (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      LOGICAL           LDGFII
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LDGFII
C ----------------------------------------------------------------------
      LICFII =.FALSE.
C
      IF (LDGFII()) THEN
         LICFII =.TRUE.
         LUIC   = LDGT
         IF (LDGFII()) THEN
   10       LUIC = 10*LUIC + LDGT
            IF (LDGFII()) GO TO 10
         ENDIF
      ENDIF
C
      RETURN
      END
      LOGICAL FUNCTION LXPFII()
C
C **********************************************************************
C
C   S A M  library routine :  LXPFII                  GROUP 1 / PUBLIC
C
C     TASK :  To look for a valid exponent in the next positions in the
C             current input line, and, if found, to decode it and store
C             its value in LEXP
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  LDGFII and LXXFII   (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C                       87-12-07 / 1.1    K.Bell  (LXXFII)
C                       88-10-15 / 1.2    K.Bell  (4 digits)
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           ISIGN,LLPP
      LOGICAL           LDGFII,LXXFII
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LDGFII,LXXFII
C ----------------------------------------------------------------------
      LXPFII =.FALSE.
      LLPP  = LPP
C
      IF (LXXFII('E')) THEN
         ISIGN = 1
         IF (LXXFII('+')) THEN
            ISIGN = 1
         ELSEIF (LXXFII('-')) THEN
            ISIGN =-1
         ENDIF
         IF (LDGFII()) THEN
            LEXP  = LDGT
            IF (LDGFII())  LEXP = 10*LEXP + LDGT
            IF (LDGFII())  LEXP = 10*LEXP + LDGT
            IF (LDGFII())  LEXP = 10*LEXP + LDGT
            LEXP   = ISIGN*LEXP
            LXPFII =.TRUE.
         ELSE
            LPP = LLPP
         ENDIF
      ENDIF
C
      RETURN
      END
      LOGICAL FUNCTION LSPFII()
C
C **********************************************************************
C
C   S A M  library routine :  LSPFII                 GROUP 1 / PUBLIC
C
C     TASK :  To look for a valid separator in the next character
C             positions of the current input line
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIISKP      (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-03-23 / 1.0
C                       96-05-19 / 1.1    K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIISKP
C ----------------------------------------------------------------------
      LSPFII =.FALSE.
C
      CALL FIISKP
      IF (CHLCOP(LPP:LPP) .EQ. CHSEP) THEN
         LPP = LPP+1
         LSPFII =.TRUE.
         CALL FIISKP
      ENDIF
C
      RETURN
      END
      LOGICAL FUNCTION LDLFII()
C
C **********************************************************************
C
C   S A M  library routine :  LDLFII                 GROUP 1 / PUBLIC
C
C     TASK :  To look for a valid delimiter in the next character
C             positions of the current input line.
C
C             This funtion replaces the "old" version of LSPFII.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIISKP and LELFII   (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C                       90-01-24 / 2.0     K.Bell
C                       96-05-19 / 2.1     K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      LOGICAL           LELFII
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIISKP,LELFII
C ----------------------------------------------------------------------
      LDLFII =.TRUE.
C
      IF (CHLCOP(LPP:LPP) .EQ. CHBLK) THEN
         CALL FIISKP
         IF (CHLCOP(LPP:LPP) .EQ. CHSEP) THEN
            LPP = LPP+1
            CALL FIISKP
         ENDIF
         GO TO 100
      ENDIF
C
      IF (CHLCOP(LPP:LPP) .EQ. CHSEP) THEN
         LPP = LPP+1
         CALL FIISKP
         GO TO 100
      ENDIF
C
      IF (LELFII())  GO TO 100
C
      LDLFII =.FALSE.
C
  100 RETURN
      END
      LOGICAL FUNCTION LELFII()
C
C **********************************************************************
C
C   S A M  library routine :  LELFII                  GROUP 1 / PUBLIC
C
C     TASK :  To look for end of line, i.e., to check that there are
C             no more non-blank characters in the current input line
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :   FIISKP    (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIISKP
C ----------------------------------------------------------------------
      LELFII =.TRUE.
      CALL FIISKP
      IF (LPP.GT.NCHL)  GO TO 100
      LELFII =.FALSE.
C
  100 RETURN
      END
      LOGICAL FUNCTION LEFFII()
C
C **********************************************************************
C
C   S A M  library routine :  LEFFII                 GROUP 1 / PUBLIC
C
C     TASK :  To look for end of file.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-28 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      LEFFII =.FALSE.
      IF (NEOL .GT. 0)  LEFFII =.TRUE.
C
      RETURN
      END
      LOGICAL FUNCTION LLTFII()
C
C **********************************************************************
C
C   S A M  library routine :  LLTFII                  GROUP 0 / PUBLIC
C
C     TASK :  To look for a letter in the current position of the
C             current input line
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  INDEX    (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-10-02 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           I
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      LLTFII =.FALSE.
      IF (LPP.GT.NCHL)  GO TO 100
C
      I = INDEX(CHALPU,CHLCOP(LPP:LPP))
      IF (I.NE.0) THEN
         LLTFII =.TRUE.
         CHLET = CHLCOP(LPP:LPP)
         LPP   = LPP+1
      ENDIF
C
  100 RETURN
      END
      LOGICAL FUNCTION LDGFII()
C
C **********************************************************************
C
C   S A M  library routine :  LDGFII                  GROUP 1 / PUBLIC
C
C     TASK :  To look for a digit in the current position of the current
C             input line
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :   INDEX    (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           I
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      LDGFII =.FALSE.
      IF (LPP.GT.NCHL)  GO TO 100
C
      I = INDEX('0123456789',CHLCOP(LPP:LPP))
      IF (I.EQ.0)       GO TO 100
      LDGFII =.TRUE.
      LDGT  = I-1
      LPP   = LPP+1
C
  100 RETURN
      END
      LOGICAL FUNCTION LNXFII(CHSTR)
C
C **********************************************************************
C
C   S A M  library routine :  LNXFII                  GROUP 1 / PUBLIC
C
C     TASK :  To compare the character string CHSTR with the characters
C             in the current data line, starting with the first non-
C             blank character at or after position LPP.
C             The string is recognized if all characters in CHSTR
C             match or if at least the NCHREC first characters match
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIISKP and FIIPOS  (SAM-1)
C                                   LEN                (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-19 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER*(*)     CHSTR
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           LLPP,LP,NC,NCR
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIISKP,FIIPOS
C ----------------------------------------------------------------------
      LNXFII =.FALSE.
C
      NC   = LEN(CHSTR)
      NCR  = 0
      CALL FIISKP
      IF (LPP.GT.NCHL)  GO TO 100
C
      LLPP = LPP
      DO 10 LP=1,NC
         IF (CHSTR(LP:LP).EQ.CHLCOP(LLPP:LLPP)) THEN
            NCR = NCR+1
         ELSE
            GO TO 20
         ENDIF
         LLPP = LLPP+1
   10 CONTINUE
C
   20 IF (NCR.LT.NC.AND.NCR.LT.NCHREC)  GO TO 100
C
      LNXFII =.TRUE.
      LPP    = LLPP
      IF (NCR.LT.NC) CALL FIIPOS(0)
C
  100 RETURN
      END
      LOGICAL FUNCTION LXXFII(CHSTR)
C
C **********************************************************************
C
C   S A M  library routine :  LXXFII                  GROUP 1 / PUBLIC
C
C     TASK :  To compare the character string CHSTR with the characters
C             in the current data line, starting with the character in
C             the current line position LPP.
C             The string is recognized if all characters in CHSTR
C             match or if at least the NCHREC first characters match
C
C     Note :  This routine is similar to LNXFII except that it does not
C             skip leading blank characters.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIIPOS             (SAM-1)
C                                   LEN                (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-12-07 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER*(*)     CHSTR
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           LLPP,LP,NC,NCR
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIIPOS
C ----------------------------------------------------------------------
      LXXFII =.FALSE.
C
      NC   = LEN(CHSTR)
      NCR  = 0
      IF (LPP.GT.NCHL)  GO TO 100
C
      LLPP = LPP
      DO 10 LP=1,NC
         IF (CHSTR(LP:LP).EQ.CHLCOP(LLPP:LLPP)) THEN
            NCR = NCR+1
         ELSE
            GO TO 20
         ENDIF
         LLPP = LLPP+1
   10 CONTINUE
C
   20 IF (NCR.LT.NC.AND.NCR.LT.NCHREC)  GO TO 100
C
      LXXFII =.TRUE.
      LPP    = LLPP
      IF (NCR.LT.NC) CALL FIIPOS(0)
C
  100 RETURN
      END
      LOGICAL FUNCTION LSTFII(CH)
C
C **********************************************************************
C
C   S A M  library routine :  LSTFII                  GROUP 1 / PUBLIC
C
C     TASK :  To compare the character CH with the last non-blank
C             character in the current data line.
C
C     Note :  The line position pointer (LPP) is not changed, regardless
C             of the outcome of the comparison.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-12-20 / 1.0
C                       90-01-24 / 1.1   K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER         CH
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      LSTFII =.FALSE.
C
      IF (CHLCOP(NCHL:NCHL) .EQ. CH)  LSTFII =.TRUE.
C
      RETURN
      END
      LOGICAL FUNCTION LNLFII(CHSTR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LNLFII                  GROUP 1 / PUBLIC
C
C     TASK :  To find the next data line, starting with the current
C             line, that has the character string CHSTR as its first
C             significant (non-blank) characters.
C             If found, LNLFII returns with the value .TRUE., and the
C             line position pointer is moved past the last character
C             of CHSTR.
C             If not found before end of file, LNLFII returns with the
C             value .FALSE.
C             Note that this function does not terminate execution due
C             to end-of-file, regardless of KEOL. However, if a Fortran
C             READ error is encountered, execution is terminated.
C
C     Member of the FII-package.
C
C     ROUTINES CALLED/REFERENCED :    LNXFII      (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-29 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER*(*)     CHSTR
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           KOPY
      LOGICAL           LNXFII
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          LNXFII
C ----------------------------------------------------------------------
      KOPY = KEOL
      KEOL = 2
C
   10 IF (LNXFII(CHSTR)) THEN
         LNLFII =.TRUE.
         GO TO 100
      ELSE
         CALL FIIRLN
         IF (NEOL .GT. 0) THEN
            LNLFII =.FALSE.
            GO TO 100
         ELSEIF (IER .LT. 0) THEN
            STOP
         ENDIF
      ENDIF
      GO TO 10
C
  100 KEOL = KOPY
      RETURN
C
      END
      CHARACTER FUNCTION CLTFII()
C
C **********************************************************************
C
C   S A M  library routine :  CLTFII                  GROUP 1 / PUBLIC
C
C     TASK :  To retrieve the last letter identified by a reference to
C             LLTFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-22 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
C
      SAVE              /FIICHR/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C ----------------------------------------------------------------------
      CLTFII = CHLET
      RETURN
      END
      INTEGER FUNCTION IDGFII()
C
C **********************************************************************
C
C   S A M  library routine :  IDGFII                  GROUP 1 / PUBLIC
C
C     TASK :  To retrieve the value of the last digit identified by a
C             reference to LDGFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-22 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IDGFII = LDGT
      RETURN
      END
      INTEGER FUNCTION IICFII()
C
C **********************************************************************
C
C   S A M  library routine :  IICFII                  GROUP 1 / PUBLIC
C
C     TASK :  To retrieve the value of last unsigned integer constant
C             (uic) identified by a reference to LICFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IICFII = LUIC
      RETURN
      END
      INTEGER FUNCTION IXPFII()
C
C **********************************************************************
C
C   S A M  library routine :  IXPFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the value of the last exponent identified by
C            a reference to LXPFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IXPFII = LEXP
      RETURN
      END
      INTEGER FUNCTION INNFII()
C
C **********************************************************************
C
C   S A M  library routine :  INNFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the value of the last integer number identified
C            by a reference to LINFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      INNFII = LNUM
      RETURN
      END
      REAL FUNCTION SFLFII()
C
C **********************************************************************
C
C   S A M  library routine :  SFLFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the (single precision) value of the last
C            floating point number identified by a reference to LFLFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      SFLFII = REAL(DPVAL)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DFLFII()
C
C **********************************************************************
C
C   S A M  library routine :  DFLFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the double precision value of the last floating
C            point number identified by a reference to LFLFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      DFLFII = DPVAL
      RETURN
      END
      CHARACTER*(*) FUNCTION CNMFII()
C
C **********************************************************************
C
C   S A M  library routine :  CNMFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the last name identified by a reference to
C            LNMFII
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  LEN
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C                       95-03-08 / 1.1    K.Bell
C
C **********************************************************************
C
      IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                                ! local variables
      INTEGER           I,L
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      L = LEN(CNMFII)
      IF (L .LT. NCHNAM) THEN
         CNMFII(1:L) = CHNAM(1:L)
      ELSEIF (L .GT. NCHNAM) THEN
         CNMFII(1:NCHNAM) = CHNAM(1:NCHNAM)
         DO 10 I=NCHNAM+1,L
            CNMFII(I:I) = ' '
   10    CONTINUE
      ELSE
         CNMFII(1:NCHNAM) = CHNAM(1:NCHNAM)
      ENDIF
C
      RETURN
      END
      CHARACTER*(*) FUNCTION CSTFII(LPS,LPE)
C
C **********************************************************************
C
C   S A M  library routine :  CSTFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the portion of the current input line starting
C            with the character in position LPS and ending with the
C            character in position LPE
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           LPE,LPS
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C                                               ! Local variables
      INTEGER           LP,N
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      LP = LPE
      IF (LP .GT. MAXCHL)       LP = MAXCHL
      N = LP-LPS+1
      IF (N.GT.0) THEN
         CSTFII(1:N) = CHLINE(LPS:LP)
      ENDIF
      RETURN
      END
      SUBROUTINE FIICHC (CH,N)
C
C **********************************************************************
C
C   S A M  library routine :  FIICHC                 GROUP 1 / PUBLIC
C
C     TASK : To reset, depending on N,  one of the three one-character
C            variables, CHCOM1, CHCOM2 or CHSEP.
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-28 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           N
      CHARACTER         CH*1
C                                               ! COMMON variables
C
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
C
      SAVE              /FIICHR/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C ----------------------------------------------------------------------
      IF (N .EQ. 1)  CHCOM1 = CH
      IF (N .EQ. 2)  CHCOM2 = CH
      IF (N .EQ. 3)  CHSEP  = CH
C
      RETURN
      END
      SUBROUTINE FIIPAR (IPAR,N)
C
C **********************************************************************
C
C   S A M  library routine :  FIIPAR                 GROUP 1 / PUBLIC
C
C     TASK : To reset one of the 7 integer parameters initially set by
C            FIIOPN
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-28 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           IPAR,N
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IF (N .EQ. 1) THEN
         LRU = IPAR
      ELSEIF (N .EQ. 2) THEN
         LPU = IPAR
      ELSEIF (N .EQ. 3) THEN
         LEU = IPAR
      ELSEIF (N .EQ. 11) THEN
         KLCV = IPAR
      ELSEIF (N .EQ. 12) THEN
         KEOL = IPAR
      ELSEIF (N .EQ. 21) THEN
         MAXCHL = IPAR
         IF (MAXCHL .GT. 132)  MAXCHL = 132
      ELSEIF (N .EQ. 22) THEN
         NCHREC = IPAR
         IF (NCHREC .LT. 1)    NCHREC = 1
      ENDIF
C
      RETURN
      END
      SUBROUTINE FIISST (CHC)
C
C **********************************************************************
C
C   S A M  library routine :  FIISST                  GROUP 1 / PUBLIC
C
C     TASK : To reset the one-character variable, CHCOM1, identifying
C            comment lines (when appearing in the first line position)
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER         CHC
C                                               ! COMMON variables
C
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
C
      SAVE              /FIICHR/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C ----------------------------------------------------------------------
      CHCOM1 = CHC
      RETURN
      END
      SUBROUTINE FIISLN (LUN,N)
C
C **********************************************************************
C
C   S A M  library routine :  FIISLN                  GROUP 1 / PUBLIC
C
C     TASK : To reset the value of one of the three logical I/O units
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           LUN,N
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IF (N.EQ.1)  LRU = LUN
      IF (N.EQ.2)  LPU = LUN
      IF (N.EQ.3)  LEU = LUN
C
      RETURN
      END
      SUBROUTINE FIISPR (IVAL,N)
C
C **********************************************************************
C
C   S A M  library routine :  FIISPR                  GROUP 1 / PUBLIC
C
C     TASK : To reset the letter conversion flag (KLCV) or the para-
C            meter defining string recognition (NCHREC)
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           IVAL,N
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IF (N.EQ.1) THEN
         IF (IVAL.EQ.0) KLCV = 0
         IF (IVAL.EQ.1) KLCV = 1
      ELSEIF (N.EQ.2) THEN
         NCHREC = IVAL
      ENDIF
C
      RETURN
      END
      INTEGER FUNCTION ILPFII(N)
C
C **********************************************************************
C
C   S A M  library routine :  ILPFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the current value of the
C            - line position pointer LPP     (N=1)
C            - parameter MAXCHL              (N=2)
C            - parameter NCHREC              (N=3)
C            - parameter KLCV                (N=4)
C            - logical read unit LRU         (N=5)
C            - logical error print unit LPU  (N=6)
C            - logical echo print unit  LEU  (N=7)
C            - parameter KEOL                (N=8)
C
C            or return the
C            - number of the current input
C              line (relative to start of
C              file)                         (N=9)
C            - pointer to the last signi-
C              ficant non-blank character
C              in the current input line     (N=10)
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   88-01-17 / 1.0
C                       90-01-24 / 1.1   K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           N
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
C
      IF (N .EQ. 1) THEN
         ILPFII = LPP
      ELSEIF (N .EQ. 2) THEN
         ILPFII = MAXCHL
      ELSEIF (N .EQ. 3) THEN
         ILPFII = NCHREC
      ELSEIF (N .EQ. 4) THEN
         ILPFII = KLCV
      ELSEIF (N .EQ. 5) THEN
         ILPFII = LRU
      ELSEIF (N .EQ. 6) THEN
         ILPFII = LPU
      ELSEIF (N .EQ. 7) THEN
         ILPFII = LEU
      ELSEIF (N .EQ. 8) THEN
         ILPFII = KEOL
      ELSEIF (N .EQ. 9) THEN
         ILPFII = LCNT
      ELSEIF (N .EQ. 10) THEN
         ILPFII = NCHL
      ELSE
         ILPFII = 0
      ENDIF
C
      RETURN
      END
      INTEGER FUNCTION IERFII()
C
C **********************************************************************
C
C   S A M  library routine :  IERFII                  GROUP 1 / PUBLIC
C
C     TASK : To retrieve the value of the error flag
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IERFII = IER
      RETURN
      END
      SUBROUTINE FIICTU (CH)
C
C **********************************************************************
C
C   S A M  library routine :  FIICTU                  GROUP 0 / PUBLIC
C
C     TASK :  To check the variable  CH*1  to see whether it contains a
C             lower case letter and if so convert it to upper case
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  INDEX    (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-10-02 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      CHARACTER         CH*1
C                                               ! COMMON variables
C
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
C                                               ! Local variables
      INTEGER           I
C
      SAVE              /FIICHR/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C ----------------------------------------------------------------------
      I = INDEX(CHALPL,CH)
      IF (I.NE.0)  CH = CHALPU(I:I)
C
      RETURN
      END
      SUBROUTINE FIISKP
C
C **********************************************************************
C
C   S A M  library routine :  FIISKP                  GROUP 1 / PUBLIC
C
C     TASK :  To skip blank characters in the input line (CHLINE), i.e.,
C             to move the current line position pointer (LPP) forward to
C             the next non-blank character
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C                       90-01-24 / 1.1   K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IF (CHLCOP(LPP:LPP) .EQ. CHBLK) THEN
   10    LPP = LPP+1
         IF (CHLCOP(LPP:LPP) .EQ. CHBLK)  GO TO 10
      ENDIF
C
      RETURN
      END
      SUBROUTINE FIIPOS (NC)
C
C **********************************************************************
C
C   S A M  library routine :  FIIPOS                  GROUP 1 / PUBLIC
C
C     TASK :  To move the line position pointer (LPP)
C                 - forward if  NC > 0
C                 - back    if  NC < 0
C                 - forward to the next blank character if  NC = 0
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-20 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      INTEGER           NC
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      IF (NC .EQ. 0) THEN
   10    IF (LPP .GT. NCHL)               GO TO 100
         IF (CHLCOP(LPP:LPP) .EQ. CHBLK)  GO TO 100
         LPP = LPP+1
         GO TO 10
      ELSE
         LPP = LPP+NC
         IF (LPP.LT.1)     LPP = 1
         IF (LPP.GT.NCHL)  LPP = NCHL+1
      ENDIF
C
  100 RETURN
      END
      SUBROUTINE FIIWLN (CHSTR,LPU)
C
C **********************************************************************
C
C   S A M  library routine :  FIIWLN                  GROUP 1 / PUBLIC
C
C     TASK :  To print the character string CHSTR on unit LPU
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-19 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER*(*)     CHSTR
      INTEGER           LPU
C ----------------------------------------------------------------------
      IF (LPU.EQ.0) THEN
         WRITE ( * ,600) CHSTR
      ELSEIF (LPU.GT.0) THEN
         WRITE (LPU,600) CHSTR
      ENDIF
C
      RETURN
C
  600 FORMAT (1X,A)
C
      END
      SUBROUTINE FIIWST (CHSTR,LPS,LPE,LPOS,LPU)
C
C **********************************************************************
C
C   S A M  library routine :  FIIWST                  GROUP 1 / PUBLIC
C
C     TASK :  To print the character substring  CHSTR(LPS:LPE) on unit
C             LPU, such that the first character is printed in line
C             position LPOS
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIIWLN      (SAM-1)
C                                   LEN         (FORTRAN Library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-19 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER*(*)     CHSTR
      INTEGER           LPE,LPOS,LPS,LPU
C                                               ! local variables
      CHARACTER         CHLN*132
      INTEGER           L,LAST,LL,NC
C
      EXTERNAL          FIIWLN
C ----------------------------------------------------------------------
      NC = LEN(CHSTR)
      LL = LPE
      IF (LL.GT.NC)      LL = NC
      IF (LPS.GT.LL)     GO TO 100
C
      LAST = LPOS+LL-LPS
      IF (LAST.GT.132)   LAST =132
      IF (LPOS.GT.LAST)  GO TO 100
C
      IF (LPOS.GT.1) THEN
         DO 10 L=1,LPOS-1
            CHLN(L:L) = ' '
   10    CONTINUE
      ENDIF
      CHLN(LPOS:LAST) = CHSTR(LPS:LL)
C
      CALL FIIWLN (CHLN(1:LAST),LPU)
C
  100 RETURN
      END
      SUBROUTINE FIIRWD
C
C **********************************************************************
C
C   S A M  library routine :  FIIRWD                 GROUP 1 / PUBLIC
C
C     TASK : To rewind the input file (associated with LRU)
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-01-29 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIINMB/
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C ----------------------------------------------------------------------
      REWIND LRU
C
      RETURN
      END
      SUBROUTINE FIIERR (CHMSG)
C
C **********************************************************************
C
C   S A M  library routine :  FIIERR                  GROUP 1 / PUBLIC
C
C     TASK :  To write the (error) message contained in the string
C             CHMSG, followed by the current input line with a question
C             mark underneath the current (offending) character
C
C     Member of the FII-package
C
C     ROUTINES CALLED/REFERENCED :  FIIWLN and FIIWST  (SAM-1)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-19 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
C
      CHARACTER*(*)     CHMSG
C                                               ! COMMON variables
C
      INTEGER           IER,KEOL,KLCV,LCNT,LDGT,LEU,LEXP,LNUM,LPP,
     +                  LPU,LRU,LUIC,MAXCHL,NCHNAM,NCHL,NCHREC,NEOL
      CHARACTER         CHLINE*133,CHLCOP*133,CHALPL*26,CHALPU*26,
     +                  CHNAM*16,CHCOM1*1,CHCOM2*1,CHBLK*1,CHLET*1,
     +                  CHSEP*1
      DOUBLE PRECISION  DPVAL
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIIWLN,FIIWST
C ----------------------------------------------------------------------
      IF (LPU.EQ.0)     WRITE ( * ,600)
      IF (LPU.GT.0)     WRITE (LPU,600)
C
      CALL FIIWLN (CHMSG,LPU)
C
      IF (LPU.EQ.0)     WRITE ( * ,600)
      IF (LPU.GT.0)     WRITE (LPU,600)
C
      CALL FIIWLN (CHLCOP(1:NCHL),LPU)
      CALL FIIWST ('?',1,1,LPP,LPU)
C
      RETURN
C
  600 FORMAT (' ')
C
      END
