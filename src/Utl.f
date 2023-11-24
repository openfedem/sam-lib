C     SPDX-License-Identifier: Apache-2.0
C
      INTEGER FUNCTION MINDX(V,IS,IE,LORD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MINDX                   GROUP 8 / PUBLIC
C
C     T A S K :  TO RETURN (THROUGH THE FUNCTION NAME) THE INDEX IN
C                ARRAY  V , BETWEEN IS AND IE, CONTAINING
C                - THE ALGEBRAICALLY SMALLEST   (LORD=1),
C                - THE NUMERICALLY SMALLEST     (LORD=2),
C                - THE ALGEBRAICALLY LARGEST    (LORD=3), OR
C                - THE NUMERICALLY LARGEST      (LORD=4)
C                ELEMENT OF THE ARRAY.
C
C     ROUTINES CALLED/REFERENCED :  ABS     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-06-08 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IE,IS,LORD
      DOUBLE PRECISION  V(IE)
C
      INTEGER           I,ISP1,K
      DOUBLE PRECISION  X
C ----------------------------------------------------------------------
      K    = 0
      IF (IS.GT.IE)     GO TO 500
      K    = IS
      IF (IS.EQ.IE)     GO TO 500
      ISP1 = IS+1
C
      IF (LORD.EQ.1) THEN
C                                               ** ALGEBRAIC. SMALLEST
         X = V(IS)
         DO 100 I=ISP1,IE
            IF (V(I).LT.X) THEN
               X = V(I)
               K = I
            ENDIF
  100    CONTINUE
C
      ELSEIF (LORD.EQ.2) THEN
C                                               ** NUMERICALLY SMALLEST
         X = ABS(V(IS))
         DO 200 I=ISP1,IE
            IF (ABS(V(I)).LT.X) THEN
               X = ABS(V(I))
               K = I
            ENDIF
  200    CONTINUE
C
      ELSEIF (LORD.EQ.3) THEN
C                                               ** ALGERAIC. LARGEST
         X = V(IS)
         DO 300 I=ISP1,IE
            IF (V(I).GT.X) THEN
               X = V(I)
               K = I
            ENDIF
  300    CONTINUE
C
      ELSEIF (LORD.EQ.4) THEN
C                                               ** NUMERICALLY LARGEST
         X = ABS(V(IS))
         DO 400 I=ISP1,IE
            IF (ABS(V(I)).GT.X) THEN
               X = ABS(V(I))
               K = I
            ENDIF
  400    CONTINUE
C
      ELSE
         K = 0
      ENDIF
C
  500 MINDX = K
C
      RETURN
      END
      SUBROUTINE RANVEC (V,RAN,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RANVEC                  GROUP 8 / PUBLIC
C
C     T A S K :  TO INTRODUCE PSEUDO-RANDOM NUMBERS, BETWEEN -0.5 AND
C                0.5, INTO ALL N LOCATIONS OF VECTOR  V
C                THE NUMBERS ARE BASED ON THE INPUT VALUE OF PARAMETER
C                RAN  -  RAN RETURNS WITH THE VALUE OF V(N)
C
C     ROUTINES CALLED/REFERENCED :  INT, REAL  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-06-14 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           N
      DOUBLE PRECISION  RAN,V(N)
C
      INTEGER           I,K
      DOUBLE PRECISION  HALF,ONE,PI,TEN,X,ZERO
C
      PARAMETER         (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TEN=10.0D0)
      PARAMETER         ( PI = 3.141592653589793D0 )
C ----------------------------------------------------------------------
C
      IF (N.GT.0) THEN
         X = RAN
         DO 100 I=1,N
            IF (X.EQ.ZERO)   X = PI/TEN
            X = X*X/REAL(N)
   10       X = TEN*X
            IF (X.LT.ONE)    GO TO 10
            K = INT(X)
            X = X - REAL(K)
            IF (X.GT.HALF)   X = X-ONE
            V(I) = X
  100    CONTINUE
         RAN = X
      ENDIF
C
      RETURN
      END
      SUBROUTINE REARRV (V,N,LORD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE : REARRV                  GROUP 8 / PUBLIC
C
C     T A S K : TO REARRANGE THE ELEMENTS OF VECTOR  V  IN AN ORDER OF
C               - INCREASING ALGEBRAIC VALUE  (LORD = 1)
C                 (ASCENDING ORDER)
C               - INCREASING NUMERICAL VALUE  (LORD = 2)
C               - DECREASING ALGEBRAIC VALUE  (LORD = 3)
C                 (DESCENDING ORDER)
C               - DECREASING NUMERICAL VALUE  (LORD = 4)
C
C     ROUTINES CALLED/REFERENCED :  MINDX   (SAM-8)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-06-14 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           LORD,N
      DOUBLE PRECISION  V(*)
C
      INTEGER           J,K,NM1,       MINDX
      DOUBLE PRECISION  E
C
      EXTERNAL          MINDX
C ----------------------------------------------------------------------
      IF (LORD.LT.1)    GO TO 100
      IF (LORD.GT.4)    GO TO 100
      IF (N.LT.2)       GO TO 100
C
      NM1 = N-1
      DO 50 J=1,NM1
         K = MINDX(V,J,N,LORD)
         IF (K.GT.J) THEN
            E    = V(J)
            V(J) = V(K)
            V(K) = E
         ENDIF
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE CONVRT (DPA,SPA,N,KON)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  CONVRT                  GROUP 8 / PUBLIC
C
C     T A S K :  TO CONVERT A DOUBLE PRECISION ARRAY  DPA  INTO A
C                SINGLE PRECISION ARRAY  SPA  (KON=21) OR VICE VERSA
C                (KON=12)
C     THE TWO ARRAYS MAY START AT THE SAME ADDRESS, IN WHICH CASE THE
C     `ORIGINAL` ARRAY IS, OF COURSE, OVERWRITTEN
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :  KOLBEIN BELL
C     DATE/VERSION  :  86-06-14 / 1.0
C
C **********************************************************************
      IMPLICIT          NONE
      INTEGER           KON,N,
     +                  I
      REAL              SPA(N)
      DOUBLE PRECISION  DPA(N)
C ----------------------------------------------------------------------
      IF (N.LT.1)    GO TO 50
C
      IF (KON.EQ.21) THEN
C                                               ** FROM DP TO SP
         DO 20 I=1,N
            SPA(I) = REAL(DPA(I))
   20    CONTINUE
      ELSEIF (KON.EQ.12) THEN
C                                               ** FROM SP TO DP
         DO 40 I=N,1,-1
            DPA(I) = DBLE(SPA(I))
   40    CONTINUE
C
      ENDIF
C
   50 RETURN
      END
      SUBROUTINE PRITAB (IARR,NEL,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRITAB                  GROUP 8 / PUBLIC
C
C     THIS SUBROUTINE PRINTS A ONE-DIMENSIONAL, INTEGER ARRAY IN TABULAR
C     (TWO-DIMENSIONAL) FORM (A4-FORMAT)
C     A FIELD WIDTH OF ONLY 6 IS USED - EACH NUMBER SHOULD THEREFORE
C     CONTAIN NO MORE THAN 5 DIGITS
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   79-12-06 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           LPU,NEL
      INTEGER           IARR(NEL)
C
      INTEGER           K,KE,KK,KS
C ----------------------------------------------------------------------
      IF (LPU.LT.1)     GO TO 100
C
      WRITE (LPU,6000)
      KK = 0
   10 KS = KK+1
      KE = KS+9
      IF (KE.GT.NEL)    KE = NEL
      WRITE (LPU,6100)  KK,(IARR(K),K=KS,KE)
      IF (KE.EQ.NEL)    GO TO 20
      KK = KK+10
      GO TO 10
   20 WRITE (LPU,6200)
C
  100 RETURN
C ------------------------------------------------ FORMATS
C
 6000 FORMAT(/6X,72('-') / 14X,'I' / 7X,'INDEX',2X,'I',
     +        5X,'1',5X,'2',5X,'3',5X,'4',5X,'5',5X,'6',5X,'7',
     +        5X,'8',5X,'9',4X,'10' /
     +       14X,'I'     / 6X,72('-') / 14X,'I' )
 6100 FORMAT(I12,2X,'I',I7,9I6  / 14X,'I' )
 6200 FORMAT(6X,72('-'))
C
      END
      SUBROUTINE PREMAT (A,NROW,NCOL,NRP,MAXNCR,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PREMAT                  GROUP 8 / PUBLIC
C
C     THIS SUBROUTINE PRINTS THE FIRST NRP ROWS AND ALL NCOL COLUMNS
C     OF THE REAL MATRIX A, USING E-FORMAT
C     MAXIMUM NO. OF COLUMNS PER ROW (MAXNCR) MAY BE
C              4    (4E15.7),
C              5    (5E13.4),
C              8    (8E15.7)  OR
C             10    (10E12.3)       (10 IS DEFAULT)
C
C     ROUTINES CALLED/REFERENCED :  MOD     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   80-11-23 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           LPU,MAXNCR,NCOL,NROW,NRP
      REAL              A(NROW,NCOL)
C
      INTEGER           I,IB,J,JC,JE,JS,NC,NR,    ICN(10)
C ----------------------------------------------------------------------
      IF (LPU.LT.1)     GO TO 100
C
      NR = NRP
      IF (NRP.GT.NROW)  NR = NROW
      NC = 10
      IF (MAXNCR.EQ.4)  NC = 4
      IF (MAXNCR.EQ.5)  NC = 5
      IF (MAXNCR.EQ.8)  NC = 8
      JS = 1
      JC = 1
C
   20 JE = JS+NC-1
      IF (JE.GT.NCOL)   JE = NCOL
C                                               ** PRINT COL. NUMBERS
      DO 40 J=1,NC
         ICN(J) = JS+J-1
         IF (ICN(J).LT.JE)  GO TO 40
         JC = J
         GO TO 50
   40 CONTINUE
   50 IF (NC.EQ.4)      WRITE (LPU,6104) (ICN(J),J=1,JC)
      IF (NC.EQ.5)      WRITE (LPU,6105) (ICN(J),J=1,JC)
      IF (NC.EQ.8)      WRITE (LPU,6108) (ICN(J),J=1,JC)
      IF (NC.EQ.10)     WRITE (LPU,6110) (ICN(J),J=1,JC)
      WRITE (LPU,6000)
C                                               ** PRINT ARRAY ELEMENTS
      DO 60 I=1,NR
         IF (NC.EQ.4)   WRITE (LPU,6204) I,(A(I,J),J=JS,JE)
         IF (NC.EQ.5)   WRITE (LPU,6205) I,(A(I,J),J=JS,JE)
         IF (NC.EQ.8)   WRITE (LPU,6208) I,(A(I,J),J=JS,JE)
         IF (NC.EQ.10)  WRITE (LPU,6210) I,(A(I,J),J=JS,JE)
         IB = MOD(I,10)
         IF (IB.EQ.0)   WRITE (LPU,6000)
   60 CONTINUE
      JS = JE+1
C                                               ** MORE COLUMNS ?
      IF (JS.LE.NCOL)  GO TO 20
C
  100 RETURN
C ------------------------------------------------ FORMATS
 6000 FORMAT (' ')
 6104 FORMAT (///7X,4I15)
 6105 FORMAT (///9X,5I13)
 6108 FORMAT (///2X,8I15)
 6110 FORMAT (///5X,10I12)
 6204 FORMAT (I10,3X,1P4E15.7)
 6205 FORMAT (I10,3X,1P5E13.4)
 6208 FORMAT (I5 ,3X,1P8E15.7)
 6210 FORMAT (I5 ,3X,1P10E12.3)
C
      END
      SUBROUTINE PRDMAT (DPA,NROW,NCOL,NRP,MAXNCR,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRDMAT                  GROUP 8 / PUBLIC
C
C     THIS SUBROUTINE PRINTS THE FIRST NRP ROWS AND ALL NCOL COLUMNS
C     OF THE DOUBLE PRECISION MATRIX  DPA
C     MAXIMUM NO. OF COLUMNS PER ROW (MAXNCR) MAY BE
C              2   (2D25.15),
C              5   (5D13.4),
C              6   (6D20.11)  OR
C             10   (10D12.3)        (10 IS DEFAULT)
C
C     ROUTINES CALLED/REFERENCED :  MOD     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   81-03-16 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           LPU,MAXNCR,NCOL,NROW,NRP
      DOUBLE PRECISION  DPA(NROW,NCOL)
C
      INTEGER           I,IB,J,JC,JE,JS,NC,NR,     ICN(10)
C ----------------------------------------------------------------------
      IF (LPU.LT.1)     GO TO 100
C
      NR = NRP
      IF (NRP.GT.NROW)  NR = NROW
      NC = 10
      IF (MAXNCR.EQ.2)  NC = 2
      IF (MAXNCR.EQ.5)  NC = 5
      IF (MAXNCR.EQ.6)  NC = 6
      JS = 1
      JC = 1
C
   20 JE = JS+NC-1
      IF (JE.GT.NCOL)   JE = NCOL
C                                               ** PRINT COL. NUMBERS
      DO 40 J=1,NC
         ICN(J) = JS+J-1
         IF (ICN(J).LT.JE)  GO TO 40
         JC = J
         GO TO 50
   40 CONTINUE
   50 IF (NC.EQ.2)      WRITE (LPU,6102) (ICN(J),J=1,JC)
      IF (NC.EQ.5)      WRITE (LPU,6105) (ICN(J),J=1,JC)
      IF (NC.EQ.6)      WRITE (LPU,6106) (ICN(J),J=1,JC)
      IF (NC.EQ.10)     WRITE (LPU,6110) (ICN(J),J=1,JC)
      WRITE (LPU,6000)
C                                               ** PRINT ARRAY ELEMENTS
      DO 60 I=1,NR
         IF (NC.EQ.2)   WRITE (LPU,6202) I,(DPA(I,J),J=JS,JE)
         IF (NC.EQ.5)   WRITE (LPU,6205) I,(DPA(I,J),J=JS,JE)
         IF (NC.EQ.6)   WRITE (LPU,6206) I,(DPA(I,J),J=JS,JE)
         IF (NC.EQ.10)  WRITE (LPU,6210) I,(DPA(I,J),J=JS,JE)
         IB = MOD(I,10)
         IF (IB.EQ.0)   WRITE (LPU,6000)
   60 CONTINUE
      JS = JE+1
C                                               ** MORE COLUMNS ?
      IF (JS.LE.NCOL)  GO TO 20
C
  100 RETURN
C ------------------------------------------------ FORMATS
 6000 FORMAT (' ')
 6102 FORMAT (///5X,2I25)
 6105 FORMAT (///9X,5I13)
 6106 FORMAT (///6I20)
 6110 FORMAT (///5X,10I12)
 6202 FORMAT (I10,5X,1P2D25.15)
 6205 FORMAT (I10,3X,1P5D13.4)
 6206 FORMAT ( I5,3X,1P6D20.11)
 6210 FORMAT ( I5,3X,1P10D12.3)
C
      END
      SUBROUTINE PRIMAT (MAT,NROW,NCOL,NRP,MAXNCR,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRIMAT                  GROUP 8 / PUBLIC
C
C     THIS SUBROUTINE PRINTS THE FIRST NRP ROWS AND ALL NCOL COLUMNS
C     OF THE INTEGER MATRIX MAT
C     MAXIMUM NO. OF COLUMNS PER ROW MAY BE
C             MAXNCR =  5  (FIELD WIDTH = 10  -  A4)
C             MAXNCR = 10  (FIELD WIDTH =  6  -  A4)
C             MAXNCR = 20  (FIELD WIDTH =  6)
C     FOR ALL OTHER VALUES OF MAXNCR A MAXIMUM OF
C             10 COLUMNS AND A FIELD WIDTH OF 10
C     ARE PRINTED PER ROW
C
C     ROUTINES CALLED/REFERENCED :  MOD     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-01-03 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           LPU,MAXNCR,NCOL,NROW,NRP
      INTEGER           MAT(NROW,NCOL)
C
      INTEGER           I,IB,J,JC,JE,JS,NC,NR,NW,     ICN(20)
C ----------------------------------------------------------------------
      IF (LPU.LT.1)      GO TO 100
C
      NR = NRP
      IF (NRP.GT.NROW)   NR = NROW
      NC = 10
      NW = 10
      IF (MAXNCR.EQ. 5)  NC =  5
      IF (MAXNCR.EQ.20)  NC = 20
      IF (MAXNCR.EQ.10)  NW =  6
      JS = 1
      JC = 1
C
   10 JE = JS+NC-1
      IF (JE.GT.NCOL)    JE = NCOL
C                                               ** PRINT COL. NUMBERS
      DO 20 J=1,NC
         ICN(J) = JS+J-1
         IF (ICN(J).LT.JE)  GO TO 20
         JC = J
         GO TO 30
   20 CONTINUE
   30 IF (NC.EQ.10)      GO TO 40
      IF (NC.EQ. 5)      WRITE (LPU,6100) (ICN(J),J=1,JC)
      IF (NC.EQ.20)      WRITE (LPU,6120) (ICN(J),J=1,JC)
      GO TO 50
   40 IF (NW.EQ. 6)      WRITE (LPU,6140) (ICN(J),J=1,JC)
      IF (NW.EQ.10)      WRITE (LPU,6160) (ICN(J),J=1,JC)
   50 WRITE (LPU,6000)
C                                               ** PRINT ARRAY ELEMENTS
      DO 80 I=1,NR
         IF (NC.EQ.10)   GO TO 60
         IF (NC.EQ. 5)   WRITE (LPU,6200) I,(MAT(I,J),J=JS,JE)
         IF (NC.EQ.20)   WRITE (LPU,6220) I,(MAT(I,J),J=JS,JE)
         GO TO 70
   60    IF (NW.EQ. 6)   WRITE (LPU,6240) I,(MAT(I,J),J=JS,JE)
         IF (NW.EQ.10)   WRITE (LPU,6260) I,(MAT(I,J),J=JS,JE)
   70    IB = MOD(I,10)
         IF (IB.EQ.0)    WRITE (LPU,6000)
   80 CONTINUE
      JS = JE+1
C                                               ** MORE COLUMNS ?
      IF (JS.LE.NCOL)    GO TO 10
C
  100 RETURN
C ------------------------------------------------ FORMATS
 6000 FORMAT (' ')
 6100 FORMAT (///3X,'COLUMN :',3X,5I10)
 6120 FORMAT (///1X,'COLUMN :',20I6)
 6140 FORMAT (///1X,'COLUMN :',10I6)
 6160 FORMAT (///3X,'COLUMN :',3X,10I10)
 6200 FORMAT (3X,'ROW',I6,3X,5I10)
 6220 FORMAT (1X,'ROW',I4,2X,20I6)
 6240 FORMAT (1X,'ROW',I4,2X,10I6)
 6260 FORMAT (3X,'ROW',I6,3X,10I10)
C
      END
      SUBROUTINE SKYCNV (ASKY,MSKY,ASQR,N,KSA)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SKYCNV                  GROUP 8 / PUBLIC
C
C     T A S K :  To convert a symmetric "skyline-stored" (ABS(KSA)=1)
C                or diagonal (ABS(KSA)=2) matrix, ASKY, to a full,
C                square N by N matrix, ASQR.
C
C
C     ROUTINES CALLED/REFERENCED :  RMINT     (SAM-0)
C                                   ABS, MOD  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-01-25 / 1.0
C                       03-07-30 / 1.1    K.M.Okstad  (KSA<0 and KSA>10)
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           KSA,N,    MSKY(N)
      DOUBLE PRECISION  ASKY(*),ASQR(N,N)
C                                                ! Local variables
      INTEGER           I,IS,J,JP
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IF (N .LT. 1)  GO TO 100
C
      IF (ABS(KSA) .LT. 10) CALL RMINT (ASQR,N,N,ZERO)
      IF (KSA .GE. 0) THEN
         ASQR(1,1) = ASKY(1)
      ELSE
         ASQR(1,1) = -ASKY(1)
      ENDIF
C
      IF (MOD(ABS(KSA),10) .EQ. 2) THEN
C                                                ! Diagonal matrix
         DO 10 J=2,N
            IF (KSA .GE. 0) THEN
               ASQR(J,J) = ASKY(J)
            ELSE
               ASQR(J,J) = -ASKY(J)
            ENDIF
   10    CONTINUE
      ELSE
C                                                ! Skyline matrix
         JP = 2
         DO 30 J=2,N
            IS = J - MSKY(J) + MSKY(J-1) + 1
            DO 20 I=IS,J
               IF (KSA .GE. 0) THEN
                  ASQR(I,J) = ASKY(JP)
                  ASQR(J,I) = ASKY(JP)
               ELSE
                  ASQR(I,J) = -ASKY(JP)
                  ASQR(J,I) = -ASKY(JP)
               ENDIF
               JP = JP+1
   20       CONTINUE
   30    CONTINUE
      ENDIF
C
  100 RETURN
      END
      SUBROUTINE MMXMAT (A,EMM,MR,MC,N,KODE,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MMXMAT                  GROUP 8 / PUBLIC
C
C     T A S K :  To determine largest and smallest diagonal and off-
C                diagonal element of a square matrix A.
C                The elements are returned in array EMM as:
C                   EMM(1) = max. diagonal element
C                   EMM(2) = min. diagonal element
C                   EMM(3) = max. off-diagonal element
C                   EMM(4) = min. off-diagonal element
C                The corresponding row and column indices are stored in
C                MR and MC, respectively.
C                If  KODE = 1 :  test is for absolute (numerical) value
C                If  KODE = 2 :  test is for algebraic value
C                If  LPU > 0  the values are printed on unit no. LPU.
C
C
C     ROUTINES CALLED/REFERENCED :  ABS     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-03-16 / 1.0
C                       93-10-22 / 1.1    K.Bell
C                       97-02-16 / 1.2    K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           KODE,LPU,N,    MC(4),MR(4)
      DOUBLE PRECISION  A(N,N),EMM(4)
C                                                ! local variables
      INTEGER           I,J
      DOUBLE PRECISION  BIG,ZERO
C
      PARAMETER         ( BIG = 1.0D10 , ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      DO 10 I=1,4
         MR(I) = 1
         MC(I) = 1
   10 CONTINUE
      EMM(2) = BIG
      EMM(4) = BIG
C
      IF (KODE .EQ. 2) THEN
         EMM(1) = -BIG
         EMM(3) = -BIG
         DO 100 J=1,N
            IF (A(J,J) .GT. EMM(1)) THEN
               MR(1)  = J
               EMM(1) = A(J,J)
            ENDIF
            IF (A(J,J) .LT. EMM(2)) THEN
               MR(2)  = J
               EMM(2) = A(J,J)
            ENDIF
            DO 50 I=1,N
               IF (I .NE. J) THEN
                  IF (A(I,J) .GT. EMM(3)) THEN
                     MR(3)  = I
                     MC(3)  = J
                     EMM(3) = A(I,J)
                  ENDIF
                  IF (A(I,J) .LT. EMM(4)) THEN
                     MR(4)  = I
                     MC(4)  = J
                     EMM(4) = A(I,J)
                  ENDIF
               ENDIF
   50       CONTINUE
  100    CONTINUE
C
      ELSE
         EMM(1) = ZERO
         EMM(3) = ZERO
         DO 200 J=1,N
            IF (ABS(A(J,J)) .GT. EMM(1)) THEN
               MR(1)  = J
               EMM(1) = ABS(A(J,J))
            ENDIF
            IF (ABS(A(J,J)) .LT. EMM(2)) THEN
               MR(2)  = J
               EMM(2) = ABS(A(J,J))
            ENDIF
            DO 150 I=1,N
               IF (I .NE. J) THEN
                  IF (ABS(A(I,J)) .GT. EMM(3)) THEN
                     MR(3)  = I
                     MC(3)  = J
                     EMM(3) = ABS(A(I,J))
                  ENDIF
                  IF (ABS(A(I,J)) .LT. EMM(4)) THEN
                     MR(4)  = I
                     MC(4)  = J
                     EMM(4) = ABS(A(I,J))
                  ENDIF
               ENDIF
  150       CONTINUE
  200    CONTINUE
      ENDIF
      MC(1) = MR(1)
      MC(2) = MR(2)
C
      IF (LPU .GT. 0) THEN
C                                                ! print
         IF (KODE .EQ. 1) THEN
            WRITE (LPU,610) N
         ELSE
            WRITE (LPU,620) N
         ENDIF
         WRITE (LPU,625)
         WRITE (LPU,630) EMM(1),MR(1),MC(1)
         WRITE (LPU,640) EMM(2),MR(2),MC(2)
         WRITE (LPU,650) EMM(3),MR(3),MC(3)
         WRITE (LPU,660) EMM(4),MR(4),MC(4)
      ENDIF
C
      RETURN
C
  610 FORMAT (///6X,'EXTREME (NUMERICALLY) ELEMENTS OF MATRIX (DIM = ',
     &           I5,')')
  620 FORMAT (///6X,'EXTREME (ALGEBRAICALLY) ELEMENTS OF MATRIX (DIM = '
     &          ,I5,')')
  625 FORMAT (/44X,'VALUE',14X,'ROW',3X,'COLUMN'/)
  630 FORMAT (6X,'Largest diagonal element   :',1PE24.15,I8,I7)
  640 FORMAT (6X,'Smallest diagonal element  :',1PE24.15,I8,I7)
  650 FORMAT (6X,'Largest off-diag. element  :',1PE24.15,I8,I7)
  660 FORMAT (6X,'Smallest off-diag. element :',1PE24.15,I8,I7 //)
C
      END
      SUBROUTINE SASCO (MTARG,MPASC,N,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SASCO                  GROUP 8 / PUBLIC
C
C     T A S K :  To sort, in ascending order, the elements of the target
C                array MTARG in such a way that MPASC(1) contains the
C                pointer or index in the target array (MTARG) containing
C                the lowest valued entry, MPASC(2) points to the second
C                lowest entry, and so on.
C
C
C     ROUTINES CALLED/REFERENCED :  INDXTV  (SAM-8)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   97-02-15 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,N
      INTEGER           MTARG(*),MPASC(*)
C
C                                                ! Local variables
      INTEGER           I,J,K,MV,     INDXTV
      EXTERNAL          INDXTV
C ----------------------------------------------------------------------
C
      IERR = 0
C                                                ! Copy, check and
C                                                  sort  MTARG
      IF (N .EQ. 1) THEN
         MPASC(1) = 1
         GOTO 100
      ELSEIF (N .LT. 1) THEN
         IERR = IERR-1
         GOTO 100
      ENDIF
C
      DO 10 I=1,N
         MPASC(I) = MTARG(I)
   10 CONTINUE
C
      DO 40 I=1,N-1
         K  = I
         MV = MPASC(I)
         DO 20 J=I+1,N
            IF (MPASC(J) .LT. MV) THEN
               MV = MPASC(J)
               K  = J
            ELSEIF (MPASC(J) .EQ. MV) THEN
               IERR = IERR-1
               GOTO 100
            ENDIF
   20    CONTINUE
         IF (K .GT. I) THEN
            J        = MPASC(I)
            MPASC(I) = MPASC(K)
            MPASC(K) = J
         ENDIF
   40 CONTINUE
C                                                ! Establish MPASC
      DO 50 I=1,N
         K = MPASC(I)
         MPASC(I) = INDXTV(MTARG,K,N)
   50 CONTINUE
C
  100 RETURN
      END
      INTEGER FUNCTION INDXTV(IARR,IVAL,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  INDXTV                 GROUP 8 / PUBLIC
C
C     T A S K :  To return the index in array IARR where IVAL is stored.
C                If IVAL is not found, INDXTV returns the value 0.
C                If IVAL is present more than once, the first occurence
C                is reported.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   97-02-15 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IVAL,N,    IARR(N)
C                                                ! Local variables
      INTEGER           I
C ----------------------------------------------------------------------
      INDXTV = 0
      IF (N .LE. 0)  GOTO 100
C
      DO 10 I=1,N
         IF (IARR(I) .EQ. IVAL) THEN
            INDXTV = I
            GOTO 100
         ENDIF
   10 CONTINUE
C
  100 RETURN
      END
