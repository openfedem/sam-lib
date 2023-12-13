C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE IMINT (IM,M,N,KONST)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  IMINT                   GROUP 0 / PUBLIC
C
C     IMINT  SETS ALL ELEMENTS OF AN INTEGER MATRIX EQUAL TO A CONSTANT
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C
      INTEGER           KONST,M,N,     IM(*)
      INTEGER           I,K,MN
C ----------------------------------------------------------------------
      MN = M*N
      IF (MN.LT.1)      GO TO 100
      K  = KONST
C
      DO 50 I=1,MN
         IM(I) = K
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE RMINT (RM,M,N,SCALAR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RMINT                   GROUP 0 / PUBLIC
C
C     RMINT  SETS ALL ELEMENTS OF A REAL MATRIX EQUAL TO A CONSTANT
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           M,N
      DOUBLE PRECISION  SCALAR,    RM(*)
C
      INTEGER           I,MN
      DOUBLE PRECISION  S
C ----------------------------------------------------------------------
      MN = M*N
      IF (MN.LT.1)      GO TO 100
      S  = SCALAR
C
      DO 50 I=1,MN
         RM(I) = S
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE JCOPY (IAR,JAR,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  JCOPY                   GROUP 0 / PUBLIC
C
C     THIS SUBROUTINE COPIES N WORDS (1 WORD = 1 INTEGER STORAGE
C     LOCATION) FROM ARRAY IAR, STATING WITH IAR(1), INTO ARRAY JAR,
C     STARTING WITH JAR(1)
C     HERE COPY IMPLIES AN EXACT COPY OF THE WORD, BIT BY BIT.  THIS
C     DFINITION MAY, ON SOME SYSTEMS, REQUIRE ASSEMBLER (AND THUS
C     MACHINE DEPENDENT) CODE
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   81-03-31 / 1.0
C
C **********************************************************************
C
      INTEGER           N
      INTEGER           IAR(N),JAR(N)
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (N .LE. 0)     GO TO 100
C
      DO 10 I=1,N
         JAR(I) = IAR(I)
   10 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE RCOPY (A,B,M,N,ISGN)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RCOPY                   GROUP 0 / PUBLIC
C
C     RCOPY  MAKES A COPY A REAL MATRIX, WITH OR WITHOUT CHANGE OF SIGN
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           ISGN,M,N
      DOUBLE PRECISION  A(*),B(*)
C
      INTEGER           I,MN
C ----------------------------------------------------------------------
      MN = M*N
      IF (MN.LT.1)      GO TO 100
      IF (ISGN.LT.0)    GO TO  50
C
      DO 25 I=1,MN
         B(I) = A(I)
   25 CONTINUE
      GO TO 100
C
   50 DO 75 I=1,MN
         B(I) =-A(I)
   75 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE SCALE (A,B,M,N,SCALAR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SCALE                   GROUP 0 / PUBLIC
C
C     SCALE  MULTIPLIES A RECTANGULAR MATRIX BY A SCALAR
C
C     ROUTINES CALLED/REFERENCED :  MOD    (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           M,N
      DOUBLE PRECISION  SCALAR,   A(*),B(*)
C
      INTEGER           I,L,LP1,MN
      DOUBLE PRECISION  S
C ----------------------------------------------------------------------
      MN = M*N
      IF(MN.LT.1)       GO TO 100
      S  = SCALAR
      L  = MOD(MN,5)
      IF (L.EQ.0)       GO TO  50
C                                               ** CLEAN-UP CODE
      DO 25 I=1,L
         B(I) = S*A(I)
   25 CONTINUE
      IF (MN.LT.5)      GO TO 100
C                                               ** MAIN LOOP
   50 LP1 = L+1
      DO 75 I=LP1,MN,5
         B(I  ) = S*A(I  )
         B(I+1) = S*A(I+1)
         B(I+2) = S*A(I+2)
         B(I+3) = S*A(I+3)
         B(I+4) = S*A(I+4)
   75 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE SCADD (A,B,C,M,N,SCALAR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SCADD                   GROUP 0 / PUBLIC
C
C     SCADD  MULTIPLIES A MATRIX BY A SCALAR AND ADDS THE RESULTING
C     MATRIX TO ANOTHER MATRIX  (SCALE AND ADD)
C
C     ROUTINES CALLED/REFERENCED :  MOD    (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           M,N
      DOUBLE PRECISION  SCALAR,    A(*),B(*),C(*)
C
      INTEGER           I,L,LP1,MN
      DOUBLE PRECISION  S
C ----------------------------------------------------------------------
      MN = M*N
      IF (MN.LT.1)      GO TO 100
      S  = SCALAR
      L  = MOD(MN,5)
      IF (L.EQ.0)       GO TO 50
C                                               ** CLEAN-UP CODE
      DO 25 I=1,L
         C(I) = A(I) + S*B(I)
   25 CONTINUE
      IF (MN.LT.5)      GO TO 100
C                                               ** MAIN LOOP
   50 LP1 = L+1
      DO 75 I=LP1,MN,5
         C(I  ) = A(I  ) + S*B(I  )
         C(I+1) = A(I+1) + S*B(I+1)
         C(I+2) = A(I+2) + S*B(I+2)
         C(I+3) = A(I+3) + S*B(I+3)
         C(I+4) = A(I+4) + S*B(I+4)
   75 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE VASCA (X,Y,INCX,INCY,SCALAR,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  VASCA                   GROUP 0 / PUBLIC
C
C     VASCA  SCALES A VECTOR AND ADDS THE RESULT TO ANOTHER VECTOR
C     (CONSTANT TIMES A VECTOR PLUS A VECTOR)
C     UNROLLED LOOPS ARE USED FOR INCREMENTS EQUAL TO ONE
C
C     ROUTINES CALLED/REFERENCED :  MOD     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 /1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           INCX,INCY,N
      DOUBLE PRECISION  SCALAR,       X(*),Y(*)
C
      INTEGER           I,IX,IY,M,MP1
      DOUBLE PRECISION  S,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IF (N.LT.1)                   GO TO 100
      IF (SCALAR.EQ.ZERO)           GO TO 100
      S = SCALAR
      IF (INCX.EQ.1.AND.INCY.EQ.1)  GO TO 50
C ----------------------------------------------------------------------
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1 (ONE)
C ----------------------------------------------------------------------
      IX = 1
      IY = 1
      DO 20 I=1,N
         X(IX) = X(IX) + S*Y(IY)
         IX    = IX + INCX
         IY    = IY + INCY
   20 CONTINUE
      GO TO 100
C ----------------------------------------------------------------------
C  CODE FOR BOTH INCREMENTS EQUAL TO 1 (ONE)
C ----------------------------------------------------------------------
   50 M = MOD(N,4)
      IF (M.EQ.0)  GO TO 70
C                                               ** CLEAN-UP CODE
      DO 60 I=1,M
         X(I) = X(I) + S*Y(I)
   60 CONTINUE
      IF (N.LT.4)  GO TO 100
C
   70 MP1 = M+1
      DO 80 I=MP1,N,4
         X(I)   = X(I)   + S*Y(I)
         X(I+1) = X(I+1) + S*Y(I+1)
         X(I+2) = X(I+2) + S*Y(I+2)
         X(I+3) = X(I+3) + S*Y(I+3)
   80 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE MATADD (A,B,C,M,N,ISGN)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MATADD                  GROUP 0 / PUBLIC
C
C     MATADD  ADDS/SUBTRACTS TWO RECTANGULAR MATRICES
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           ISGN,M,N
      DOUBLE PRECISION  A(*),B(*),C(*)
C
      INTEGER           I,MN
C ----------------------------------------------------------------------
      MN = M*N
      IF (MN.LT.1)      GO TO 100
      IF (ISGN.LT.0)    GO TO  50
C
      DO 25 I=1,MN
         C(I) = A(I)+B(I)
   25 CONTINUE
      GO TO 100
C
   50 DO 75 I=1,MN
         C(I) = A(I)-B(I)
   75 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE COLMOD (B,V,M,N,ISGN)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  COLMOD                  GROUP 0 / PUBLIC
C
C     COLMOD  ADDS (IF ISGN.GE.0) OR SUBTRACTS (IF ISGN.LT.0) THE
C     ELEMENTS OF VECTOR  V  TO/FROM EACH COLUMN OF MATRIX  B
C
C     ROUTINES CALLED/REFERENCED :  MATADD        (SAM-3)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 /1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
      INTEGER           ISGN,M,N
      DOUBLE PRECISION  B(M,N),V(M)
C
      INTEGER           J
C
      EXTERNAL          MATADD
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)  GO TO 100
C
      DO 10 J=1,N
         CALL MATADD (B(1,J),V,B(1,J),M,1,ISGN)
   10 CONTINUE
C
  100 RETURN
      END
      DOUBLE PRECISION FUNCTION PRACC(A,B,INCA,INCB,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRACC                   GROUP 0 / PUBLIC
C
C     FORMS THE INNER (DOT) PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL  (AFTER  J.DONGARRA)
C     DATE/VERSION  :   83-12-28 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           INCA,INCB,N
      DOUBLE PRECISION  A(*),B(*)
C
      INTEGER           I,IA,IB,M,MP1
      DOUBLE PRECISION  ACC,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      ACC = ZERO
      IF (N.LE.0)                   GO TO 100
      IF (INCA.EQ.1.AND.INCB.EQ.1)  GO TO 20
C ----------------------------------------------------------------------
C     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
C ----------------------------------------------------------------------
      IA  = 1
      IB  = 1
      DO 10 I=1,N
         ACC = ACC + A(IA)*B(IB)
         IA  = IA  + INCA
         IB  = IB  + INCB
   10 CONTINUE
      GO TO 100
C ----------------------------------------------------------------------
C     CODE FOR BOTH INCREMENTS EQUAL TO 1
C ----------------------------------------------------------------------
   20 M = MOD(N,5)
      IF (M .EQ. 0)     GO TO 40
C                                               ** CLEAN-UP LOOP
      DO 30 I=1,M
         ACC = ACC + A(I)*B(I)
   30 CONTINUE
      IF (N .LT. 5)     GO TO 100
C                                               ** MAIN LOOP
   40 MP1 = M+1
      DO 50 I=MP1,N,5
         ACC = ACC + A(I  )*B(I  ) + A(I+1)*B(I+1) + A(I+2)*B(I+2)
     +             + A(I+3)*B(I+3) + A(I+4)*B(I+4)
   50 CONTINUE
C
  100 PRACC = ACC
      RETURN
      END
      REAL FUNCTION CPUSEC (T)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  CPUSEC                  GROUP 0 / PUBLIC
C
C     THIS MACHINE DEPENDENT FUNCTION RETURNS AS ITS VALUE ACCUMULATED
C     CPU-TIME, IN DECIMAL SECONDS, RELATIVE TO JOB START, MINUS THE
C     VALUE OF INPUT ARGUMENT T
C
C     ROUTINES CALLED/REFERENCED :  CPU_TIME / ETIME
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ............
C     DATE/VERSION  :   THIS IS A  D U M M Y VERSION RETURNING
C                       CPUSEC = 0.
C
C **********************************************************************
C
      REAL              T
C
#if defined(win32) || defined(win64) || defined(aix) || defined(aix64)
      REAL TIM
      CALL CPU_TIME (TIM)
      CPUSEC = TIM - T
#else
      REAL TOTALT,TARRAY(2)
      CALL ETIME(TARRAY,TOTALT)
      CPUSEC = TARRAY(1) - T
#endif
C
      RETURN
      END
      SUBROUTINE DATIME (IYEAR,IMONTH,IDAY,IHOUR,IMIN,ISEC)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  DTIME                   GROUP 0 / PUBLIC
C
C     RETURNS THE DATE AND WALL CLOCK TIME
C
C     ROUTINES CALLED/REFERENCED :  DATE_AND_TIME / IDATE, ITIME
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ........
C     DATE/VERSION  :   THIS IS A  D U M M Y VERSION RETURNING
C                       ZERO IN ALL ARGUMENTS
C
C **********************************************************************
C
      INTEGER           IDAY,IHOUR,IMIN,IMONTH,ISEC,IYEAR
C
#if defined(win32) || defined(win64) || defined(sunos)
      INTEGER VALUES(8)
      CALL DATE_AND_TIME (VALUES=VALUES)
      IYEAR  = VALUES(1)
      IMONTH = VALUES(2)
      IDAY   = VALUES(3)
      IHOUR  = VALUES(5)
      IMIN   = VALUES(6)
      ISEC   = VALUES(7)
#elif defined(aix) || defined(aix64)
      CALL IDATE (IYEAR,IMONTH,IDAY,IHOUR,IMIN,ISEC)
#else
      integer, parameter :: i4 = selected_int_kind(9) ! 4-byte integer
      INTEGER(kind=i4) VALUES(3)
c      CALL IDATE (VALUES(1),VALUES(2),VALUES(3))
c     IDATE call updated for linux/gnu gfortran
      CALL IDATE (VALUES)
      IMONTH = VALUES(1)
      IDAY   = VALUES(2)
      IYEAR  = 2000 + MOD(VALUES(3),100_i4)
      CALL ITIME (VALUES)
      IHOUR  = VALUES(1)
      IMIN   = VALUES(2)
      ISEC   = VALUES(3)
#endif
C
      RETURN
      END
      SUBROUTINE PUTREC (LFU,KEY,INCAD,NWR,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PUTREC                  GROUP 0 / PUBLIC
C
C     THIS MACHINE DEPENDENT SUBROUTINE WRITES A RECORD OF NWR WORDS,
C     STORED IN ARRAY INCAD, INTO FILE LFU AT ADDRESS KEY.
C     THE TRANSFER IS ACCOMPLISHED, SECTOR BY SECTOR, BY MEANS OF THE
C     DIRECT ACCESS I/O-FEATURE OF FORTRAN-77 (IMPLEMENTED IN SUB-
C     ROUTINE DIRT).  THE FILE ADDRESS, KEY, IS UPDATED
C
C     ROUTINES CALLED/REFERENCED :  DIRT  (SAM - 0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ........   NO.OF WORDS PER SECTOR (NWS) = 128
C                                  NWS SHOULD BE SET FOR EACH PARTICULAR
C                                  COMPUTER
C     DATE/VERSION  :   84-01-07 / 1.0
C
C **********************************************************************
C
      INTEGER           IERR,KEY,LFU,NWR
      INTEGER           INCAD(*)
C
      INTEGER           IP,NW,NWS
C
      PARAMETER         ( NWS=128 )
C
      EXTERNAL          DIRT
C ----------------------------------------------------------------------
      IERR = 0
      IF (NWR.LT.1)               GO TO 100
      IF (LFU.LT.1.OR.LFU.GT.99)  GO TO  91
      IF (KEY.LT.1)               GO TO  92
C
      IP   = 1
      NW   = NWS
      IF (NW.GT.NWR)              NW = NWR
C
C     DO-WHILE MORE SECTORS LEFT
   10    CALL DIRT (LFU,KEY,INCAD(IP),NW,1,IERR)
         IF (IERR.NE.0)           GO TO 93
         KEY = KEY+1
         IP  = IP+NW
         IF (IP.GT.NWR)           GO TO 100
         IF ((NWR-IP).LT.NW)      NW = NWR-IP+1
         GO TO 10
C     END-DO
C ------------------------------------------------ ERROR EXIT
   91 IERR =-991
      GO TO 100
   92 IERR =-992
      GO TO 100
   93 IERR =-IERR
C ------------------------------------------------
  100 RETURN
      END
      SUBROUTINE GETREC (LFU,KEY,INCAD,NWR,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  GETREC                  GROUP 0 / PUBLIC
C
C     THIS MACHINE DEPENDENT SUBROUTINE READS A RECORD OF NWR WORDS FROM
C     FILE LFU, STARTING AT ADDRESS KEY, INTO THE CORE ARRAY INCAD.
C     THE TRANSFER IS ACCOMPLISHED, SECTOR BY SECTOR, BY MEANS OF THE
C     DIRECT ACCESS I/O-FEATURE OF FORTRAN-77 (IMPLEMENTED IN SUB-
C     ROUTINE DIRT).  THE FILE ADDRESS, KEY, IS UPDATED.
C
C     ROUTINES CALLED/REFERENCED :  DIRT    (SAM - 0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ........   NO.OF WORDS PER SECTOR (NWS) = 128
C                                  NWS SHOULD BE SET FOR EACH PARTICULAR
C                                  COMPUTER
C     DATE/VERSION  :   84-01-07 / 1.0
C
C **********************************************************************
C
      INTEGER           IERR,KEY,LFU,NWR
      INTEGER           INCAD(*)
C
      INTEGER           IP,NW,NWS
C
      PARAMETER         ( NWS=128 )
C
      EXTERNAL          DIRT
C ----------------------------------------------------------------------
      IERR = 0
      IF (NWR.LT.1)               GO TO 100
      IF (LFU.LT.1.OR.LFU.GT.99)  GO TO  91
      IF (KEY.LT.1)               GO TO  92
C
      IP   = 1
      NW   = NWS
      IF (NW.GT.NWR)              NW = NWR
C
C     DO-WHILE MORE SECTORS LEFT
C
   10    CALL DIRT (LFU,KEY,INCAD(IP),NW,2,IERR)
         IF (IERR.NE.0)           GO TO  93
         KEY = KEY+1
         IP  = IP+NW
         IF (IP.GT.NWR)           GO TO 100
         IF ((NWR-IP).LT.NW)      NW = NWR-IP+1
         GO TO 10
C     END-DO
C ------------------------------------------------ ERROR EXIT
   91 IERR =-991
      GO TO 100
   92 IERR =-992
      GO TO 100
   93 IERR =-IERR
C ------------------------------------------------
  100 RETURN
      END
      SUBROUTINE DIRT (LFU,KEY,INCAD,NWT,IOP,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  DIRT                    GROUP 0 / PRIVATE
C
C     THIS IS AN AUXILIARY ROUTINE FOR THE FORTRAN-77 VERSION OF
C     PUTREC/GETREC TRANSFER OF DATA BETWEEN PRIMARY AND SECONDARY
C     STORAGE.
C     DIRT WRITES (IOP=1) OR READS (IOP=2) ONE SECTOR (OR PART OF A
C     SECTOR), TO/FROM FILE UNIT LFU AT FILE ADDRESS KEY, FROM/TO
C     PRIMARY STORAGE INCAD(NWT), BY MEANS OF THE DIRECT ACCESS
C     READ/WRITE FEATURE OF FORTRAN-77  (NOTE:  SECTOR DESIGNATES THE
C     SAME STORAGE UNIT AS THE TERM RECORD USED BY THE FORTRAN-77
C     STANDARD).  NWT MUST BE GREATER THAN ZERO, BUT IT MUST NOT EXCEED
C     THE 'RECORD' LENGTH DEFINED FOR THE FILE WHEN OPENED.
C     DIRT ENABLES DIRECT TRANSFER WITHOUT BUFFERING OR IMPLIED DO.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-01-09 / 1.0
C                       88-03-07 / 1.1
C
C **********************************************************************
C
      INTEGER           IERR,IOP,LFU,KEY,NWT
      INTEGER           INCAD(NWT)
C ----------------------------------------------------------------------
      IF (IOP.EQ.1)     WRITE (LFU,REC=KEY,IOSTAT=IERR,ERR=100) INCAD
      IF (IOP.EQ.2)     READ  (LFU,REC=KEY,IOSTAT=IERR,ERR=100) INCAD
C
  100 RETURN
      END
      INTEGER FUNCTION NEWKEY (KEY,NWR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NEWKEY                  GROUP 0 / PUBLIC
C
C     THIS MACHINE DEPENDENT FUNCTION DETERMINES, AND RETURNS AS ITS
C     VALUE,THE NEXT AVAILABLE RECORD KEY (FILE ADDRESS IN TERMS OF
C     'SECTORS'), ASSUMING A RECORD WITH NWR WORDS IS TO BE STORED/READ
C     FROM THE ADDRESS HELD BY THE CURRENT VALUE OF ARGUMENT KEY (KEY
C     IS NOT ALTERED)
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ........   ONE WORD = 32 BITS
C                                  NO. OF WORDS PER SECTOR (NWS) = 128
C                                  NWS SHOULD BE SET FOR EACH PARTICULAR
C                                  COMPUTER
C     DATE/VERSION  :   80-11-16 / 1.0
C
C **********************************************************************
C
      INTEGER           KEY,NWR
C
      INTEGER           NS,NWS
C
      PARAMETER         ( NWS=128 )
C ----------------------------------------------------------------------
      NS     = NWR/NWS
      IF (NS*NWS.LT.NWR)  NS = NS+1
      NEWKEY = KEY+NS
      RETURN
      END
      SUBROUTINE ERREC (LFU,KEY,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ERREC                   GROUP 0 / PUBLIC
C
C     THIS SUBROUTINE PRINTS AN ERROR MESSAGE FOR AN ERROR SITUATION
C     OCCURING DURING DATA TRANSFER USING PUTREC OR GETREC.
C     ERREC MUST BE CALLED WITH ARGUMENTS LFU, KEY AND IERR AS RE-
C     TURNED FROM  PUTREC/GETREC.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-01-09 / 1.0
C
C **********************************************************************
C
      INTEGER           IERR,KEY,LFU,LPU
      INTEGER           LERR
C ----------------------------------------------------------------------
      IF (IERR.EQ.0)       GO TO 100
      IF (LPU.LT.1)        GO TO 100
      LERR = -IERR
      WRITE (LPU,600)      LFU,KEY
      IF (IERR.EQ.(-991).OR.IERR.EQ.(-992))  LERR = 0
      IF (IERR.EQ.(-991))  WRITE (LPU,610)
      IF (IERR.EQ.(-992))  WRITE (LPU,620)
      IF (LERR.NE.0)       WRITE (LPU,630)   LERR
C
  100 RETURN
C ----------------------------------------------------------------------
  600 FORMAT (///' *** ERROR DURING DATA TRANSFER TO/FROM :'
     +        /5X,    'FILE UNIT',I5,'  (AT ADDRESS',I6,'  )')
  610 FORMAT ( 5X,    'ILLEGAL FILE UNIT')
  620 FORMAT ( 5X,    'ILLEGAL FILE ADDRESS')
  630 FORMAT ( 5X,    'ERROR CODE =',I5,3X,'(SYSTEM PARAMETER)')
C
      END
      SUBROUTINE OPNFIL (CHNAME,LFU,ISTAT,IACC,IFORM,IBLANK,
     +                   LREC,IPSW,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  OPNFIL                  GROUP 0 / PUBLIC
C
C     THIS MACHINE DEPENDENT ROUTINE OPENS A FILE WITH NAME  CHNAME
C     ON LOGICAL UNIT  LFU  WITH THE FOLLOWING CHARACTERISTICS :
C
C     STATUS CODE   ISTAT =  1 :  UNKNOWN  (DEFAULT)
C                         =  2 :  OLD
C                         =  3 :  NEW
C                         =  4 :  SCRATCH
C     ACCESS CODE   IACC  = 10 :  SEQUENTIAL / READ AND WRITE (DEFAULT)
C                         = 20 :  DIRECT / READ AND WRITE
C                         = 30 :  DIRECT - PUTREC/GETREC TRANSFER
C     FORMAT CODE   IFORM =  1 :  FORMATTED READ/WRITE
C                                 (DEFAULT IF SEQUENTIAL)
C                         =  2 :  UNFORMATTED READ/WRITE
C                                 (DEFAULT IF DIRECT)
C     RECORD LENGTH       = LREC, IN TERMS OF WORDS (NO DEFAULT)
C                           DUMMY IF IACC=10 OR 30
C     BLANKS CODE  IBLANK =  0 : ALL BLANKS IN NUMERIC INPUT FIELDS
C                                ARE IGNORED  (DEFAULT)
C                         =  1 : ALL NON-LEADING BLANKS IN NUMERIC
C                                INPUT FIELDS ARE TREATED AS ZEROS
C                           DUMMY IF IFORM=2
C
C     ROUTINES CALLED/REFERENCED :  ABS     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ........
C     DATE/VERSION  :   83-12-30 / 1.0
C
C **********************************************************************
C
      INTEGER           IACC,IBLANK,IERR,IFORM,IPSW,ISTAT,LFU,LPU,LREC
      CHARACTER         CHNAME*(*)
C
      INTEGER           LENREC
      CHARACTER         CHACC*10,CHBLNK*4,CHFORM*11,CHSTAT*7
C ----------------------------------------------------------------------
      IERR = 0
C ----------------------------------------------------------------------
C     CHECK INPUT AND PREPARE PARAMETERS FOR OPEN STATEMENT
C ----------------------------------------------------------------------
      IF (LFU.LT.1.OR.LFU.GT.99)   GO TO 910
C                                               ** DEFAULTS
      CHSTAT = 'UNKNOWN'
      CHACC  = 'SEQUENTIAL'
      CHFORM = 'FORMATTED'
      CHBLNK = 'NULL'
C                                               ** FILE STATUS
      IF (ISTAT.EQ.2)    CHSTAT = 'OLD'
      IF (ISTAT.EQ.3)    CHSTAT = 'NEW'
      IF (ISTAT.EQ.4)    CHSTAT = 'SCRATCH'
C                                               ** FORMAT
      IF (IACC.EQ.20)    CHFORM = 'UNFORMATTED'
      IF (IACC.EQ.30)    CHFORM = 'UNFORMATTED'
      IF (IFORM.EQ.1)    CHFORM = 'FORMATTED'
      IF (IFORM.EQ.2)    CHFORM = 'UNFORMATTED'
C                                               ** RECORD LENGTH
      LENREC = 0
      IF (IACC.EQ.20)    GO TO 10
      IF (IACC.EQ.30)    GO TO 20
      GO TO 50
   10 IF (LREC.LT.1)     GO TO 920
      LENREC = LREC
      IF (IFORM.EQ.1)    LENREC = 4*LENREC
      GO TO 50
   20 LENREC = 512
C                                               ** ACCESS MODE
   50 IF (IACC.EQ.20)    CHACC  = 'DIRECT'
      IF (IACC.EQ.30)    CHACC  = 'DIRECT'
C                                               ** BLANKS HANDLING
      IF (CHFORM.EQ.'UNFORMATTED')  GO TO 100
      IF (IBLANK.EQ.1)              CHBLNK = 'ZERO'
C ----------------------------------------------------------------------
C     FILE OPENING
C ----------------------------------------------------------------------
  100 IF (CHSTAT.EQ.'SCRATCH')      GO TO 400
      IF (CHACC.EQ.'SEQUENTIAL')    GO TO 300
C
      IF (CHFORM.EQ.'FORMATTED')    GO TO 210
C                                               ** UNFORMATTED DIRECT
C
      OPEN (UNIT=LFU,FILE=CHNAME,ACCESS=CHACC,STATUS=CHSTAT,
     +      RECL=LENREC,FORM=CHFORM,ERR=930,IOSTAT=IERR)
      GO TO 500
C                                               ** FORMATTED DIRECT
C
  210 OPEN (UNIT=LFU,FILE=CHNAME,ACCESS=CHACC,STATUS=CHSTAT,
     +      RECL=LENREC,FORM=CHFORM,BLANK=CHBLNK,ERR=930,IOSTAT=IERR)
      GO TO 500
  300 IF (CHFORM.EQ.'FORMATTED')    GO TO 310
C                                               ** UNFORMATTED
C                                                  SEQUENTIAL
C
      OPEN (UNIT=LFU,FILE=CHNAME,ACCESS=CHACC,STATUS=CHSTAT,
     +      FORM=CHFORM,ERR=930,IOSTAT=IERR)
      GO TO 500
C                                               ** FORMATTED SEQUENTIAL
C
  310 OPEN (UNIT=LFU,FILE=CHNAME,ACCESS=CHACC,STATUS=CHSTAT,
     +      FORM=CHFORM,BLANK=CHBLNK,ERR=930,IOSTAT=IERR)
      GO TO 500
C
  400 IF (CHACC.EQ.'DIRECT')        GO TO 420
      IF (CHFORM.EQ.'FORMATTED')    GO TO 410
C                                               ** UNFORMATTED SCRATCH
C                                                  SEQUENTIAL
      OPEN (UNIT=LFU,ACCESS=CHACC,STATUS=CHSTAT,
     +      FORM=CHFORM,ERR=930,IOSTAT=IERR)
      GO TO 500
C                                               ** FORMATTED SCRATCH
C                                                  SEQUENTIAL
  410 OPEN (UNIT=LFU,ACCESS=CHACC,STATUS=CHSTAT,
     +      FORM=CHFORM,BLANK=CHBLNK,ERR=930,IOSTAT=IERR)
      GO TO 500
C
  420 IF (CHFORM.EQ.'FORMATTED')    GO TO 430
C
C                                               ** UNFORMATTED SCRATCH
C                                                  DIRECT
      OPEN (UNIT=LFU,ACCESS=CHACC,STATUS=CHSTAT,
     +      RECL=LENREC,FORM=CHFORM,ERR=930,IOSTAT=IERR)
      GO TO 500
C                                               ** FORMATTED SCRATCH
C                                                  DIRECT
  430 OPEN (UNIT=LFU,ACCESS=CHACC,STATUS=CHSTAT,
     +      RECL=LENREC,FORM=CHFORM,BLANK=CHBLNK,ERR=930,IOSTAT=IERR)
C
  500 IF (IPSW.LT.1.OR.LPU.LT.1)    GO TO 1000
C ----------------------------------------------------------------------
C     PRINT FILE OPENING INFORMATION
C ----------------------------------------------------------------------
      IF (CHSTAT.EQ.'SCRATCH')      GO TO 510
      WRITE (LPU,6010)  CHNAME,LFU
      GO TO 520
  510 WRITE (LPU,6020)  LFU
  520 WRITE (LPU,6030)  CHSTAT,CHACC,CHFORM,LENREC
      IF (CHFORM.EQ.'FORMATTED')    WRITE (LPU,6040) CHBLNK
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IF (LPU .GT. 0)   WRITE (LPU,6900) LFU
      IF (LPU .GT. 0)   WRITE (LPU,6910)
      IERR = -901
      GO TO 1000
  920 IF (LPU .GT. 0)   WRITE (LPU,6900) LFU
      IF (LPU .GT. 0)   WRITE (LPU,6920) LREC
      IERR = -902
      GO TO 1000
  930 IF (LPU .GT. 0)   WRITE (LPU,6900) LFU
      IF (LPU .GT. 0)   WRITE (LPU,6930) IERR
      IERR = -ABS(IERR)
C ------------------------------------------------
 1000 RETURN
C ------------------------------------------------ FORMATS
C
 6010 FORMAT (///5X,'FILE  ',A
     +          /5X,'IS OPENED ON UNIT NO.',I4)
 6020 FORMAT (///5X,'A SCRATCH FILE IS OPENED ON UNIT NO.',I4)
 6030 FORMAT (   5X,'FILE CHARACTERISTICS :'
     +         /10X,'STATUS =',A7
     +         /10X,'ACCESS =',A10
     +         /10X,'FORM   =',A11
     +         /10X,'LENGTH =',I5,' WORDS/BYTES PER RECORD (SECTOR)')
 6040 FORMAT (  10X,'BLANK  =',A4)
C
 6900 FORMAT (///' *** ERROR DURING OPENING OF FILE UNIT',I5)
 6910 FORMAT (5X,'ILLEGAL FILE UNIT')
 6920 FORMAT (5X,'ILLEGAL RECORD LENGTH =',I5)
 6930 FORMAT (5X,'ERROR CODE =',I5,'  (SYSTEM PARAMETER)')
C
      END
      SUBROUTINE CLSFIL (LFU,KSTAT,IPSW,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  CLSFIL                  GROUP 0 / PUBLIC
C
C     THIS (FORTRAN-77) ROUTINE CLOSES THE FILE CONNECTED TO UNIT NO.
C     LFU  -  THE DISPOSITION OF THE FILE IS :
C             KSTAT = 1 : KEEP
C             KSTAT = 2 : DELETE
C             KSTAT = ANY OTHER VALUE : SYSTEM DEFAULT WHICH IS
C                                       - DELETE FOR SCRATCH FILES
C                                       - KEEP FOR ALL OTHER FILES
C
C     ROUTINES CALLED/REFERENCED :  ABS   (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-01-03 / 1.0
C
C **********************************************************************
C
      INTEGER           IERR,IPSW,KSTAT,LFU,LPU
C ----------------------------------------------------------------------
      IERR = 0
      IF (KSTAT .LT. 1)  CLOSE (LFU,ERR=90,IOSTAT=IERR)
      IF (KSTAT .GT. 2)  CLOSE (LFU,ERR=90,IOSTAT=IERR)
      IF (KSTAT .EQ. 1)  CLOSE (LFU,STATUS='KEEP',ERR=90,IOSTAT=IERR)
      IF (KSTAT .EQ. 2)  CLOSE (LFU,STATUS='DELETE',ERR=90,IOSTAT=IERR)
      IF (IPSW  .LT. 1)  GO TO 100
      IF (LPU   .LT. 0)  GO TO 100
      WRITE (LPU,600)    LFU
      GO TO 100
C ------------------------------------------------ ERROR EXIT
   90 IF (LPU .GT. 0)    WRITE (LPU,690) LFU,IERR
      IERR = -ABS(IERR)
C ------------------------------------------------
  100 RETURN
C ------------------------------------------------ FORMATS
  600 FORMAT (///5X,'FILE CONNECTED TO UNIT NO.',I5,'  IS CLOSED')
C
  690 FORMAT (///' *** ERROR DURING CLOSING OF FILE UNIT',I5
     +       /5X,'ERROR CODE =',I5,'  (SYSTEM PARAMETER)')
      END
      INTEGER FUNCTION IMP (NUM)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  IMP                     GROUP 0 / PUBLIC
C
C     RETURNS AS ITS VALUE THE FOLLOWING MACHINE-DEPENDENT PARAMETERS :
C
C     NUM = 1 :  STANDARD (PERMANENTLY OPEN) READ UNIT NO.
C         = 2 :  STANDARD (PERMANENTLY OPEN) WRITE UNIT NO.
C         = 3 :  NUMBER OF SIGNIFICANT (DECIMAL) DIGITS IN FLOATING-
C                POINT ARITHMETIC
C         = 4 :  LARGEST (POSITIVE) INTEGER NUMBER IN STANDARD
C                (DEFAULT) INTEGER MODE
C         = 5 :  MAXIMUM (DECIMAL) EXPONENT
C         = 6 :  MINIMUM (DECIMAL) EXPONENT
C         = 7 :  NO. OF BITS PER (STANDARD INTEGER) WORD
C         = 8 :  NO. OF CHARACTERS PER (STANDARD INTEGER) WORD
C         = 9 :  NO. OF WORDS PER FLOATING POINT DATUM
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     FOR USE ON    :   ............
C     DATE/VERSION  :   THIS IS A  D U M M Y  VERSION RETURNING ZERO
C                       FOR ALL VALUES OF NUM, EXCEPT FOR NUM=3
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           NUM
      INTEGER           I,NVAL
      DOUBLE PRECISION  AUX,EPS,TEN,ONE,   DUM
C
      PARAMETER         ( ONE = 1.0D0 , TEN = 10.0D0 )
C
      COMMON            /IMPDUM/ DUM
C ----------------------------------------------------------------------
C
      NVAL = 0
      IF (NUM.EQ.3)     GO TO 10
      GO TO 50
C
   10 NVAL = 0
      EPS  = ONE
      DO 20 I=1,50
         EPS = EPS/TEN
         AUX = ONE + EPS
         DUM = AUX
         IF (AUX.GT.ONE)  GO TO 20
         NVAL = I-1
         GO TO 50
   20 CONTINUE
C
   50 IMP = NVAL
C
      RETURN
      END
      SUBROUTINE FIIRLN
C
C **********************************************************************
C
C   S A M  library routine :  FIIRLN                  GROUP 0 / PUBLIC
C
C     TASK :  To read a new data (input) line from unit LRU into the
C             character string CHLINE and echo the line on unit LEU.
C             If the line is a comment line or if it is a blank line,
C             a new line is read.
C             A copy is made in CHLCOP, from character position 1 up
C             to and including the last significant (non-blank/non-tab)
C             character preceeding an in-line comment marker (if any)
C             or end of line.  If the "letter conversion flag" is on
C             (KLCV=1) all lower case letters in CHLCOP are converted
C             to upper case.  Also all occurences of the tab character
C             (ASCII value 9) are replaced by the blank character
C             (ASCII value 32).
C             The line position pointer (LPP) is set to point at the
C             first non-blank character, and the line length "pointer"
C             NCHL is set to include the last significant character on
C             the line (i.e., the last non-blank, non-comment character.
C             End-of-file condition is handled according to the EOL-flag
C             KEOL.
C
C     Member of the FII-package (but maintained in the MDR package, see
C     below).
C
C     NOTE :  This subroutine is "machine dependent" to the extent that
C             it assumes that the ASCII character set is used.
C
C
C     ROUTINES CALLED/REFERENCED :  FIIWLN           (SAM-1)
C                                   CHAR,   ICHAR    (Fortran library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   87-09-19 / 1.0
C                       90-01-20 / 2.0   K.Bell        ### ASCII ###
C                       96-06-27 / 2.1   K.Bell
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
C                                               ! Local variables
      INTEGER           IOS,L,LP,NA
C
      SAVE              /FIICHR/ , /FIINMB/
C
      COMMON /FIICHR/   CHLINE,CHLCOP,CHALPL,CHALPU,CHNAM,
     +                  CHCOM1,CHCOM2,CHBLK,CHLET,CHSEP
C
      COMMON /FIINMB/   DPVAL,MAXCHL,NCHL,NCHREC,KEOL,KLCV,LRU,LPU,LEU,
     +                  IER,LCNT,LDGT,LEXP,LUIC,LNUM,NCHNAM,NEOL,LPP
C
      EXTERNAL          FIIWLN
C ----------------------------------------------------------------------
      IER  = 0
      LPP  = 1
C
   10 LCNT = LCNT+1
      IF (LRU .EQ. 0) THEN
         READ( * ,500,IOSTAT=IOS,ERR=90,END=91)  CHLINE(1:MAXCHL)
      ELSE
         READ(LRU,500,IOSTAT=IOS,ERR=90,END=91)  CHLINE(1:MAXCHL)
      ENDIF
C
      CALL FIIWLN (CHLINE(1:MAXCHL),LEU)
C                                               ! comment line ?
      IF (CHLINE(1:1) .EQ. CHCOM1)  GO TO 10
C                                               ! set marker at end+1
      CHLINE(MAXCHL+1 : MAXCHL+1) = 'S'
C                                               ! move past blanks and
C                                                 tabs
      LP = 0
   20 LP = LP+1
      IF (CHLINE(LP:LP) .EQ. CHBLK)     GO TO 20
      IF (ICHAR(CHLINE(LP:LP)) .EQ. 9)  GO TO 20
C
C                                               ! comment line ?
      IF (CHLINE(LP:LP) .EQ. CHCOM2)  GO TO 10
C                                               ! blank line ?
      IF (LP .GT. MAXCHL)             GO TO 10
C ----------------------------------------------------------------------
C     Copy input line (CHLINE) into CHLCOP and convert to upper case
C     (if KLCV=1) and convert tabs to blanks - also look for in-line
C     comment marker
C ----------------------------------------------------------------------
      LPP = LP
      IF (LPP .GT. 1) THEN
         DO 25 LP=1,LPP-1
            CHLCOP(LP:LP) = CHBLK
   25    CONTINUE
         LP = LPP
      ENDIF
C
   30 NA  = ICHAR(CHLINE(LP:LP))
C
      IF (NA.GT.96 .AND. NA.LT.123) THEN
         IF (KLCV .EQ. 1) THEN
C                                               ! convert to upper case
            CHLCOP(LP:LP) = CHAR(NA-32)
         ELSE
            CHLCOP(LP:LP) = CHLINE(LP:LP)
         ENDIF
         NCHL = LP
         LP   = LP+1
C
      ELSEIF (NA.EQ.32 .OR. NA.EQ.9) THEN
C                                               ! blanks and tabs
   40    LP = LP+1
         IF (CHLINE(LP:LP) .EQ. CHBLK)     GO TO 40
         IF (ICHAR(CHLINE(LP:LP)) .EQ. 9)  GO TO 40
         IF (LP .GT. MAXCHL)               GO TO 80
         DO 45 L=NCHL+1,LP-1
            CHLCOP(L:L) = CHBLK
   45    CONTINUE
C
      ELSEIF (CHLINE(LP:LP) .EQ. CHCOM2) THEN
C                                               ! in-line comment
         IF (LP .EQ. LPP)  GO TO 10
         GO TO 80
C
      ELSE
C                                               ! other character
         CHLCOP(LP:LP) = CHLINE(LP:LP)
         NCHL = LP
         LP   = LP+1
C
      ENDIF
      IF (LP .LE. MAXCHL)  GO TO 30
C                                               ! set marker at end+1
   80 CHLCOP(NCHL+1:NCHL+1) = 'S'
      GO TO 100
C ----------------------------------------------- error condition
   90 IER =-IOS
      IF (LPU.EQ.0) THEN
         WRITE( * ,600) LCNT,LRU
         WRITE( * ,620) IOS
      ELSEIF (LPU.GT.0) THEN
         WRITE(LPU,600) LCNT,LRU
         WRITE(LPU,620) IOS
      ENDIF
      GO TO 100
C                                               ! end of file
   91 NEOL = NEOL+1
      IF (KEOL .GT. 0) THEN
         IF (KEOL .GE. NEOL)  GO TO 100
      ENDIF
      IER =-999
      IF (LPU.EQ.0) THEN
         WRITE( * ,600) LCNT,LRU
         WRITE( * ,610)
         IF (KEOL .GT. 0) WRITE ( * ,630) KEOL
      ELSEIF (LPU.GT.0) THEN
         WRITE(LPU,600) LCNT,LRU
         WRITE(LPU,610)
         IF (KEOL .GT. 0) WRITE ( * ,630) KEOL
      ENDIF
      IF (KEOL .EQ. 0)  GO TO 100
      STOP
C -----------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  500 FORMAT (A)
  600 FORMAT (/' *** ERROR DURING READING OF LINE',I6,'  ON UNIT',I4)
  610 FORMAT (   5X,'ATTEMPT TO READ THROUGH AN END OF FILE')
  620 FORMAT (   5X,'IOSTAT =',I4)
  630 FORMAT (   5X,'(',I5,' TIMES )')
C
      END
