C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SKYSOL (A,B,MSKY,EPS,NEQ,NEQ1,NRHS,
     +                   LPU,IFLAG,NND,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  SKYSOL                 GROUP 4 / PUBLIC
C
C     THIS SUBROUTINE SOLVES (NEQ1=NEQ) OR REDUCES (NEQ1.LT.NEQ) A
C     SYSTEM OF EQUATIONS    A*X = B    ( A  IS SYMMETRIC)  THROUGH
C     THREE STEPS:
C        1) FACTORIZATION OF A         -  A = L*D*(LT)
C        2) FORWARD REDUCTION OF B     -  Y = (D-1)*(L-1)*B
C        3) BACKSUBSTITUTION TO GIVE X -  X = (L-T)*Y
C     STORAGE-WISE  LT AND D REPLACES A, Y REPLACES B AND X REPLACES Y.
C     THE ACTIVE COLUMNS OF THE UPPER TRIANGULAR PART OF  A  ARE STORED
C     CONSECUTIVELY IN ARRAY A  (LOCATION OF DIAGONAL ELEMENT NO. I IS
C     STORED IN MSKY(I))
C     DEPENDING ON THE VALUE OF THE OPERATION FLAG (IFLAG), ONE OF THE
C     THREE STEPS, THE FIRST TWO, THE LAST TWO, OR ALL THREE STEPS
C     MAY BE CARRIED OUT IN ONE REFERENCE TO SKYSOL
C
C     ROUTINES CALLED/REFERENCED :  PRACC        (SAM-0)
C                                   ABS AND MIN  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL   (AFTER E.L.WILSON AND H.H.DOVEY)
C     DATE/VERSION  :   80-04-29 / 1.0
C                       91-11-05 / 1.1   K.Bell
C                       96-12-10 / 1.2   K.Bell
C                       04-05-06 / 1.3   K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,LPU,NEQ,NEQ1,NND,NRHS,  MSKY(*)
      DOUBLE PRECISION  EPS, A(*),B(NEQ,*)
C
      INTEGER           I,ID,IE,II,IJ,IS,J,JS,K,L,NJ,NN
      DOUBLE PRECISION  R,S,T,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
C  PRELIMINARIES  -  CHECK PARAMETERS AND DIRECT EXECUTION
C ----------------------------------------------------------------------
      IERR = 0
      IF (NEQ .LT.1)         GO TO 910
      IF (NEQ1.LT.1)         GO TO 910
      IF (NEQ1.GT.NEQ)       GO TO 910
      IF (IFLAG.EQ.1)        GO TO 100
      IF (IFLAG.EQ.2)        GO TO 100
      IF (IFLAG.EQ.3)        GO TO 100
      IF (IFLAG.EQ.4)        GO TO 400
      IF (IFLAG.EQ.5)        GO TO 400
      IF (IFLAG.EQ.6)        GO TO 600
      GO TO 910
C ----------------------------------------------------------------------
C  STEP 1  -  FACTORIZATION
C ----------------------------------------------------------------------
  100 J = 1
      IF (A(1).EQ.ZERO)      GO TO 930
      IF (NEQ.EQ.1)          GO TO 350
C
      DO 300 J=2,NEQ
         NJ = MSKY(J) - MSKY(J-1)
         IF (NJ. EQ. 1) THEN
            ID = MSKY(J)
            IF (A(ID) .EQ. ZERO) THEN
               GO TO 930
            ELSE
               GO TO 300
            ENDIF
         ENDIF
         K  = J-NJ+1
         IF (K.LT.1)         GO TO 920
C                                               ** FORM COLUMN J OF  U
         I  = K-1
         IF (I.EQ.0)         I=1
C
  200    I  = I+1
         NN = MIN((NJ-J+I),(MSKY(I)-MSKY(I-1))) - 1
         IS = MSKY(I)-NN
         IJ = MSKY(J)-J+I
C
         IF (I.EQ.J)         GO TO 250
C
         JS = IJ-NN
         IF (I.GT.NEQ1)      NN = NN+NEQ1-I+1
         IF (NN.GT.0)        A(IJ) = A(IJ) - PRACC(A(IS),A(JS),1,1,NN)
         GO TO 200
C                                               ** FORM COLUMN J OF  LT
C                                                  AND  D(J,J)
  250    IE = MSKY(I)-1
         IF (I.GT.NEQ1)      IE = IE+NEQ1-I+1
         IF (IE.LT.IS)       GO TO 300
         R  = ZERO
         DO 275 II=IS,IE
            ID    = MSKY(K)
            T     = A(II)
            A(II) = A(II)/A(ID)
            R     = R + A(II)*T
            K     = K + 1
  275    CONTINUE
         T     = ABS(A(IJ))
         A(IJ) = A(IJ)-R
         IF (J.GT.NEQ1)      GO TO 300
C                                               ** CHECK FOR ILL COND.
         IF (A(IJ).EQ.ZERO)  GO TO 930
         S  = ABS(A(IJ))
         IF (S.LT.(EPS*T))   GO TO 940
  300 CONTINUE
C
  350 IF (IFLAG.EQ.1)        GO TO 800
C ----------------------------------------------------------------------
C  STEP 2  -  FORWARD REDUCTION
C ----------------------------------------------------------------------
  400 IF (NRHS.LT.1)         GO TO 910
      DO 500 L=1,NRHS
         IF (NEQ.EQ.1)       GO TO 450
         DO 425 I=2,NEQ
            NN = MSKY(I) - MSKY(I-1) - 1
            IS = MSKY(I) - NN
            K  = I-NN
            IF (I.GT.NEQ1)  NN = NN+NEQ1-I+1
            IF (NN.GT.0)    B(I,L) = B(I,L) - PRACC(A(IS),B(K,L),1,1,NN)
  425    CONTINUE
  450    DO 475 I=1,NEQ1
            ID = MSKY(I)
            B(I,L) = B(I,L)/A(ID)
  475    CONTINUE
  500 CONTINUE
C
      IF (IFLAG.EQ.2)        GO TO 800
      IF (IFLAG.EQ.5)        GO TO 1000
C ----------------------------------------------------------------------
C  STEP 3  -  BACKSUBSTITUTION
C ----------------------------------------------------------------------
  600 IF (NRHS.LT.1)         GO TO 910
      DO 700 L=1,NRHS
         I  = NEQ
  625    IF (I.EQ.1)         GO TO 700
         K  = I - MSKY(I) + MSKY(I-1) + 1
         IS = MSKY(I-1) + 1
         IE = MSKY(I) - 1
         IF (I.GT.NEQ1)      IE = IE+NEQ1-I+1
         IF (IE.LT.IS)       GO TO 675
         DO 650 II=IS,IE
            B(K,L) = B(K,L) - A(II)*B(I,L)
            K = K+1
  650    CONTINUE
  675    I = I-1
         GO TO 625
  700 CONTINUE
C
      IF (IFLAG.GT.3)        GO TO 1000
C ----------------------------------------------------------------------
C  CHECK FOR NEGATIVE DIAGONAL ELEMENTS
C ----------------------------------------------------------------------
  800 NND = 0
      DO 850 I=1,NEQ
         J = MSKY(I)
         IF (A(J).LT.ZERO)   NND = NND+1
  850 CONTINUE
      IF (NND.GT.0)          IERR = 3
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      NND = J
      GO TO 950
  930 IERR =-3
      NND = J
      GO TO 950
  940 IERR =-4
      NND = J
  950 IF (LPU.LT.0)          GO TO 1000
C
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))      WRITE (LPU,6910) NEQ,NEQ1,NRHS,IFLAG
      IF (IERR.EQ.(-2))      WRITE (LPU,6920) J,MSKY(J-1),MSKY(J)
      IF (IERR.EQ.(-3))      WRITE (LPU,6930) J
      IF (IERR.EQ.(-4))      WRITE (LPU,6940) J,EPS,T,S
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE SKYSOL')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  NEQ    NEQ1    NRHS   IFLAG'/
     +                                    I32,     I8,     I8,     I8)
 6920 FORMAT(/5X,'INCORRECT ENTRIES IN MSKY :  J   MSKY(J-1)  MSKY(J)'/
     +                                       I35,       I11,      I9 )
 6930 FORMAT(/5X,'ZERO PIVOT ELEMENT IN EQUATION',I6)
 6940 FORMAT(/5X,'EQUATION',I6,'  FAILED THE STABILITY TEST'
     +       /5X,'EPS =',1PE12.4,5X,'A(I,I) =',1PE12.4,5X,'D(I,I) =',
     +        1PE12.4 )
C
      END
      SUBROUTINE MINV (A,C,WA,IWA,EPS,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  MINV                   GROUP 4 / PUBLIC
C
C     THIS SUBROUTINE INVERTS A NON-SYMMETRIC MATRIX  A  AND STORES THE
C     INVERSE IN  C
C     CROUT FACTORIZATION WITH PARTIAL PIVOTING IS USED
C
C     ROUTINES CALLED/REFERENCED :  PRACC         (SAM-0)
C                                   ABS AND SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   80-02-20 / 1.0
C                       93-03-18 / 1.1    K.Bell
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LPU,N
      INTEGER           IWA(N)
      DOUBLE PRECISION  EPS,  A(N,N),C(N,N),WA(N)
C
      INTEGER           I,IP1,J,K,KM1,KP1,L,NM1
      DOUBLE PRECISION  S,T,TOL,ONE,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
C  PRELIMINARIES
C ----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1)          GO TO 910
      IF (N.GT.1)          GO TO  20
      IF (A(1,1).EQ.ZERO)  GO TO 920
      C(1,1) = ONE/A(1,1)
      GO TO 1000
C                                                ! initiate IWA and find
C                                                  threshold value
   20 TOL = ZERO
      DO 60 J=1,N
         IWA(J) = J
         DO 40 I=1,N
            IF (ABS(A(I,J)) .GT. TOL)  TOL = ABS(A(I,J))
   40    CONTINUE
   60 CONTINUE
      TOL = TOL*EPS
C                                               ** SCALING FACTORS FOR
C                                                  THE PIVOT SEARCH
      DO 80 I=1,N
         S = PRACC(A(I,1),A(I,1),N,N,N)
         IF (S.LE.ZERO)    GO TO 920
         WA(I) = SQRT(S)
   80 CONTINUE
C ----------------------------------------------------------------------
C  CROUT FACTORIZATION
C ----------------------------------------------------------------------
      DO 300 K=1,N
         KM1 = K-1
         KP1 = K+1
         L   = K
         IF (K.EQ.1)       GO TO 120
C                                               ** K'TH COLUMN OF L
         DO 100 I=K,N
            A(I,K) = A(I,K) - PRACC(A(I,1),A(1,K),N,1,KM1)
  100    CONTINUE
C
  120    S = ABS(A(K,K))*WA(K)
         IF (K.EQ.N)       GO TO 160
C                                               ** SEARCH FOR K'TH PIVOT
         DO 140 I=KP1,N
            T = ABS(A(I,K))*WA(I)
            IF (T.LE.S)    GO TO 140
            S = T
            L = I
  140    CONTINUE
C                                               ** CHECK FOR ILL COND.
  160    IF (S.EQ.ZERO)    GO TO 920
         I = IWA(L)
         IF (ABS(A(L,K)) .LT. TOL)  GO TO 920
         IF (K.EQ.N)       GO TO 300
         IF (L.EQ.K)       GO TO 200
C                                               ** INTERCHANGE ROWS
         S      = WA(K)
         WA(K)  = WA(L)
         WA(L)  = S
         IWA(L) = IWA(K)
         IWA(K) = I
         DO 180 J=1,N
            S      = A(L,J)
            A(L,J) = A(K,J)
            A(K,J) = S
  180    CONTINUE
C                                               ** THE K'TH ROW OF U
  200    S = ONE/A(K,K)
         IF (K.GT.1)       GO TO 240
         DO 220 J=KP1,N
            A(K,J) = S*A(K,J)
  220    CONTINUE
         GO TO 300
  240    DO 260 J=KP1,N
            A(K,J) = S*(A(K,J) - PRACC(A(K,1),A(1,J),N,1,KM1))
  260    CONTINUE
  300 CONTINUE
C ----------------------------------------------------------------------
C  SOLVE FOR THE INVERSE MATRIX  -  COLUMN BY COLUMN
C ----------------------------------------------------------------------
      NM1 = N-1
      DO 600 K=1,N
         KP1 = K+1
         KM1 = K-1
         J   = IWA(K)
         IF (K.EQ.1)       GO TO 340
         DO 320 I=1,KM1
            C(I,J) = ZERO
  320    CONTINUE
C                                               ** FORWARD REDUCTION IN
C                                                  COLUMN J
  340    C(K,J) = ONE/A(K,K)
         IF (K.EQ.N)       GO TO 420
         DO 400 I=KP1,N
            L      = I-K
            C(I,J) = -PRACC(A(I,K),C(K,J),N,1,L)/A(I,I)
  400    CONTINUE
C                                               ** BACKSUBSTITUTION IN
C                                                  COLUMN J
  420    DO 500 L=1,NM1
            I      = N-L
            IP1    = I+1
            C(I,J) = C(I,J) - PRACC(A(I,IP1),C(IP1,J),N,1,L)
  500    CONTINUE
  600 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
  950 IF (LPU.LT.0)        GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))    WRITE (LPU,6910) N
      IF (IERR.EQ.(-2))    WRITE (LPU,6920)
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE MINV')
 6910 FORMAT(/5X,'ILLEGAL MATRIX DIMENSION  N =',I6 )
 6920 FORMAT(/5X,'SINGULAR OR NEAR-SINGULAR MATRIX')
C
      END
      SUBROUTINE MINVIP (A,EPS,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  MINVIP                 GROUP 4 / PUBLIC
C
C     THIS SUBROUTINE INVERTS A NON-SYMMETRIC MATRIX 'IN PLACE', THAT IS
C     WITHOUT USING EXTRA STORAGE, BY THE GAUSS-JORDAN METHOD
C
C     ROUTINES CALLED/REFERENCED :  ABS AND SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   JAMES M.NELSON AND KOLBEIN BELL
C     DATE/VERSION  :   80-02-21 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,LPU,N
      DOUBLE PRECISION  EPS, A(N,N)
C
      INTEGER           I,J,K
      DOUBLE PRECISION  S,ONE,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C ----------------------------------------------------------------------
C
      IERR = 0
      IF (N.LT.1)               GO TO 910
      S = ZERO
      DO 10 I=1,N
         S = S + A(I,I)*A(I,I)
   10 CONTINUE
      S = SQRT(S)*EPS
C
      DO 60 I=1,N
         IF (ABS(A(I,I)).LT.S)  GO TO 920
         A(I,I) = ONE/A(I,I)
         DO 20 J=1,N
            IF (J.EQ.I)         GO TO 20
            A(I,J) = A(I,J)*A(I,I)
   20    CONTINUE
         DO 40 K=1,N
            IF (K.EQ.I)         GO TO 40
            DO 30 J=1,N
               IF (J.EQ.I)      GO TO 30
               A(K,J) = A(K,J) - A(I,J)*A(K,I)
   30       CONTINUE
   40    CONTINUE
         DO 50 K=1,N
            IF (K.EQ.I)         GO TO 50
            A(K,I) = -A(K,I)*A(I,I)
   50    CONTINUE
   60 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
  950 IF (LPU.LT.0)             GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))         WRITE (LPU,6910) N
      IF (IERR.EQ.(-2))         WRITE (LPU,6920)
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE MINVIP')
 6910 FORMAT(/5X,'HILLEGAL MATRIX DIMENSION  N =',I6)
 6920 FORMAT(/5X,'SINGULAR OR NEAR-SINGULAR MATRIX')
C
      END
      SUBROUTINE GJMINV (A,MC,MP,MR,EPS,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  GJMINV                 GROUP 4 / PUBLIC
C
C     TASK :  To invert a non-symmetric matrix "in-place", that is with-
C             out using extra storage, by the Gauss-Jordan method.
C             Full pivoting is used, and the code is based on subroutine
C             GAUSSJ in "NUMERICAL RECIPES" by W.H.Press et. al.
C             (Cambridge Univ. Press, 1986).
C
C     ROUTINES CALLED/REFERENCED :  ABS       (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-03-16 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LPU,N,    MC(N),MP(N),MR(N)
      DOUBLE PRECISION  EPS, A(N,N)
C                                                ! local variables
      INTEGER           I,ICOL,IROW,J,K
      DOUBLE PRECISION  BIG,DUM,PINV,ONE,TOL,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
      IF (N .LT. 1)  GO TO 910
      BIG  = ZERO
      DO 20 J=1,N
         MP(J) = 0
         DO 10 I=1,N
            IF (ABS(A(I,J)) .GT. BIG)  BIG = ABS(A(I,J))
   10    CONTINUE
   20 CONTINUE
      IF (BIG .LE. ZERO)  GO TO 920
C
      IF (N .EQ. 1) THEN
C                                                !  N = 1
         A(1,1) = ONE / A(1,1)
         GO TO 1000
      ENDIF
C                                                ! threshold value
      TOL  = BIG*EPS
      ICOL = 1
      IROW = 1
C
      DO 100 J=1,N
         BIG = ZERO
C                                                ! find pivot
         DO 40 I=1,N
            IF (MP(I) .NE. 1) THEN
               DO 30 K=1,N
                  IF (MP(K) .EQ. 0) THEN
                     IF (ABS(A(I,K)) .GE. BIG) THEN
                        BIG  = ABS(A(I,K))
                        IROW = I
                        ICOL = K
                     ENDIF
                  ELSEIF (MP(K) .GT. 1) THEN
                     GO TO 930
                  ENDIF
   30          CONTINUE
            ENDIF
   40    CONTINUE
C
         MP(ICOL) = MP(ICOL) + 1
         IF (IROW .NE. ICOL) THEN
C                                                ! row interchange
            DO 50 K=1,N
               DUM       = A(IROW,K)
               A(IROW,K) = A(ICOL,K)
               A(ICOL,K) = DUM
   50       CONTINUE
         ENDIF
         MR(J) = IROW
         MC(J) = ICOL
         IF (ABS(A(ICOL,ICOL)) .LE. TOL)  GO TO 930
         PINV  = ONE / A(ICOL,ICOL)
         A(ICOL,ICOL) = ONE
C                                                ! divide pivot row
         DO 60 K=1,N
            A(ICOL,K) = A(ICOL,K)*PINV
   60    CONTINUE
C                                                ! reduce rows
         DO 80 I=1,N
C                                                ! except pivot row
            IF (I .NE. ICOL) THEN
               DUM       = A(I,ICOL)
               A(I,ICOL) = ZERO
               DO 70 K=1,N
                  A(I,K) = A(I,K) - A(ICOL,K)*DUM
   70          CONTINUE
            ENDIF
   80    CONTINUE
  100 CONTINUE
C
C --- unscramble by column interchange in reverse order
C
      DO 200 J=N,1,-1
         IF (MR(J) .NE. MC(J)) THEN
            DO 150 I=1,N
               DUM        = A(I,MR(J))
               A(I,MR(J)) = A(I,MC(J))
               A(I,MC(J)) = DUM
  150       CONTINUE
         ENDIF
  200 CONTINUE
      GO TO 1000
C ------------------------------------------------ error exit
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      GO TO 950
  930 IERR =-3
  950 IF (LPU .GT. 0) THEN
         WRITE (LPU,6900)
         IF (IERR .EQ. (-1))  WRITE (LPU,6910) N
         IF (IERR .EQ. (-2))  WRITE (LPU,6920)
         IF (IERR .EQ. (-3))  WRITE (LPU,6930)
      ENDIF
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN from  S A M  LIBRARY routine GJMINV')
 6910 FORMAT(/5X,'Illegal matrix dimension:   N =',I6)
 6920 FORMAT(/5X,'Cannot invert the null matrix !')
 6930 FORMAT(/5X,'Singular or near-singular matrix')
C
      END
      SUBROUTINE SYMINV (A,EPS,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  SYMINV                 GROUP 4 / PUBLIC
C
C     THIS SUBROUTINE INVERTS A SYMMETRIC MATRIX  A  AND STORES THE
C     INVERSE IN  A
C
C     ROUTINES CALLED/REFERENCED :  PRACC   (SAM-0)
C                                   ABS     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   80-02-22 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,LPU,N
      DOUBLE PRECISION  EPS, A(N,N)
C
      INTEGER           I,IM1,IP1,J,JM1,JP1,K,NM1
      DOUBLE PRECISION  S,T,ONE,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1)                    GO TO 910
      J = 1
      IF (A(1,1).EQ.ZERO)            GO TO 920
      IF (N.GT.1)                    GO TO 100
      A(1,1) = ONE/A(1,1)
      GO TO 1000
C ----------------------------------------------------------------------
C  FACTORIZE A INTO   L*D*LT  (LT AND D ARE STORED)
C ----------------------------------------------------------------------
  100 DO 200 J=2,N
         JM1 = J-1
         IF (JM1.EQ.1)               GO TO 140
C                                               ** FORM COLUMN J OF U
         DO 120 I=2,JM1
            IM1    = I-1
            A(I,J) = A(I,J) - PRACC(A(1,I),A(1,J),1,1,IM1)
  120    CONTINUE
C                                               ** FORM COLUMN J OF LT
C                                                  AND D(J,J)
  140    S = ZERO
         DO 160 I=1,JM1
            T      = A(I,J)
            A(I,J) = A(I,J)/A(I,I)
            S      = S + T*A(I,J)
  160    CONTINUE
         T      = ABS(A(J,J))
         A(J,J) = A(J,J) - S
C                                               ** CHECK FOR ILL COND.
         IF (A(J,J).EQ.ZERO)         GO TO 920
         IF (ABS(A(J,J)).LT.(T*EPS)) GO TO 920
  200 CONTINUE
C ----------------------------------------------------------------------
C  SOLVE FOR THE LOWER TRIANGULAR PART OF  A-INVERSE
C ----------------------------------------------------------------------
      NM1 = N-1
      DO 400 J=1,NM1
         JP1    = J+1
         S      = A(J,J)
         A(J,J) = ONE
C                                               ** FORWARD REDUCTION IN
C                                                  COLUMN J
         DO 320 I=JP1,N
            K      = I-J
            A(I,J) =-PRACC(A(J,I),A(J,J),1,1,K)
  320    CONTINUE
         A(J,J) = ONE/S
         DO 340 I=JP1,N
            A(I,J) = A(I,J)/A(I,I)
  340    CONTINUE
C                                               ** BACKSUBSTITUTION IN
C                                                  COLUMN J
         I      = N
  350    I      = I-1
         IP1    = I+1
         K      = N-I
         A(I,J) = A(I,J) - PRACC(A(I,IP1),A(IP1,J),N,1,K)
         IF (I.GT.J)  GO TO 350
  400 CONTINUE
      I      = N
      A(I,I) = ONE/A(I,I)
C ------------------------------------------------ BY SYMMRTY
      DO 600 I=1,N
         DO 550 J=I,N
            A(I,J) = A(J,I)
  550    CONTINUE
  600 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
  950 IF (LPU.LT.0)                  GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))              WRITE (LPU,6910) N
      IF (IERR.EQ.(-2))              WRITE (LPU,6920) J
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE SYMINV')
 6910 FORMAT(/5X,'ILLEGAL MATRIX DIMENSION  N =',I6)
 6920 FORMAT(/5X,'SINGULAR OR NEAR-SINGULAR MATRIX  (ROW',I5,' )')
C
      END
      SUBROUTINE BANSOL (A,B,EPS,NEQ,NBW,NRHS,LPU,KODE,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  BANSOL                 GROUP 4 / PUBLIC
C
C     BANSOL SOLVES A SYSTEM OF LINEAR EQUATIONS  A*X = B  THROUGH
C     THREE STEPS :
C        1) FACTORIZATION OF  A  ( A = L*D*LT )
C        2) FORWARD REDUCTION    ( Z = (L-1)*B )
C        3) BACKSUBSTITUTION     ( X = (L-T)*(D-1)*Z )
C     A  IS A SYMMETRIC BAND MATRIX (THE UPPER HALF BAND IS STORED IN
C     ARRAY A(NBW,NEQ)).   L  IS A LOWER TRIANGULAR BAND MATRIX WITH
C     UNIT ELEMENTS ON THE DIAGONAL, AND  D  IS A DIAGONAL MATRIX
C     ( U = D*LT  IS AN UPPER TRIANGULAR MATRIX ).
C     DEPENDING ON THE VALUE OF ARGUMENT  KODE,
C        ONE OF THE THREE STEPS,
C        THE LAST TWO STEPS OR
C        ALL THREE STEPS
C     MAY BE CARRIED OUT IN ONE REFERENCE TO BANSOL
C
C     ROUTINES CALLED/REFERENCED :  PRACC             (SAM-0)
C                                   ABS, MAX AND MIN  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY  :  KOLBEIN BELL
C     DATE/VERSION   :  84-12-05 / 1.0
C                       86-11-09 / 2.0   K.BELL
C                       87 01-11 / 2.1   K.BELL
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C                       D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KODE,LPU,NBW,NEQ,NRHS
      DOUBLE PRECISION  A(NBW,NEQ),B(NEQ,NRHS),EPS
C
      INTEGER           I,IS,J,JMNBW,K,L,LMIN,N,NBWM1,NBWP1
      DOUBLE PRECISION  S,T,ZERO,    PRACC
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
C                                               ** PRELIMINARIES
      IERR = 0
      IF (NEQ.LT.1)     GO TO 910
      IF (NBW.LT.1)     GO TO 910
      IF (NBW.GT.NEQ)   GO TO 910
      NBWM1 = NBW-1
      NBWP1 = NBW+1
C                                               ** BRANCHING
      IF (KODE.EQ.1)    GO TO 100
      IF (KODE.EQ.2)    GO TO 100
      IF (KODE.EQ.3)    GO TO 500
      IF (KODE.EQ.4)    GO TO 500
      IF (KODE.EQ.5)    GO TO 600
      GO TO 910
C
C ----------------------------------------------------------------------
C  STEP 1  -  FACTORIZATION   ( KODE = 1 AND 2)
C ----------------------------------------------------------------------
C
  100 IF (NBW.EQ.1) THEN
C                                               ** DIAGONAL MATRIX
         DO 125 J=1,NEQ
            IF (A(1,J).EQ.ZERO)  GO TO 920
  125    CONTINUE
C
      ELSE
C                                               ** BAND MATRIX
         J = 1
         IF (A(NBW,1).EQ.ZERO)   GO TO 920
C
         DO 200 J=2,NEQ
            LMIN = MAX(1,NBWP1-J)
            IF (J.GT.2) THEN
C                                               ** FORM COLUMN J OF  U
               IS    = LMIN+1
               JMNBW = J-NBW
               DO 150 I=IS,NBWM1
                  K = JMNBW+I
                  L = MAX(NBWP1-K,NBWP1-I)
                  N = NBW-L
                  A(I,J) = A(I,J) - PRACC(A(L,K),A(LMIN,J),1,1,N)
  150          CONTINUE
            ENDIF
C                                               ** FORM COLUMN J OF  LT
C                                                  AND  D(J,J)
            K = MAX(1,J-NBWM1)
            S = ZERO
            DO 175 I=LMIN,NBWM1
               T      = A(I,J)
               A(I,J) = A(I,J)/A(NBW,K)
               S      = S + A(I,J)*T
               K      = K+1
  175       CONTINUE
            T        = A(NBW,J)
            A(NBW,J) = T-S
C                                               ** CHECK FOR ILL COND.
C
            IF (A(NBW,J).EQ.ZERO)               GO TO 920
            IF (ABS(A(NBW,J)).LT.(EPS*ABS(T)))  GO TO 930
  200    CONTINUE
      ENDIF
C                                               ** BRANCHING
      IF (KODE.EQ.1)     GO TO 1000
C
C ----------------------------------------------------------------------
C  STEP 2  -  FORWARD REDUCTION    (KODE = 2, 3 AND 4)
C ----------------------------------------------------------------------
C
  500 IF (NRHS.LT.1)     GO TO 910
      IF (NBW.GT.1) THEN
C
         DO 550 J=1,NRHS
            DO 525 I=2,NEQ
               L = MAX(1,NBWP1-I)
               K = MAX(1,I-NBWM1)
               N = NBW-L
               B(I,J) = B(I,J) - PRACC(A(L,I),B(K,J),1,1,N)
  525       CONTINUE
  550    CONTINUE
      ENDIF
C                                               ** BRANCHING
      IF (KODE.EQ.4)     GO TO 1000
C
C ----------------------------------------------------------------------
C  STEP 3  -  BACKSUBSTITUTION    (KODE = 2, 3 AND 5)
C ----------------------------------------------------------------------
C
  600 IF (NRHS.LT.1)     GO TO 910
C
      DO 650 J=1,NRHS
         DO 625 I=1,NEQ
            B(I,J) = B(I,J)/A(NBW,I)
  625    CONTINUE
  650 CONTINUE
C
      IF (NBW.GT.1) THEN
C
         DO 700 J=1,NRHS
            DO 675 L=1,NEQ-1
               I = NEQ-L
               N = MIN(NBWM1,L)
               B(I,J) = B(I,J)
     +                  - PRACC(A(NBWM1,I+1),B(I+1,J),NBWM1,1,N)
  675       CONTINUE
  700    CONTINUE
      ENDIF
      GO TO 1000
C
C ------------------------------------------------ ERROR EXIT
C
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      GO TO 950
  930 IERR =-3
C
  950 IF (LPU.LE.0)      GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))  WRITE (LPU,6910) NEQ,NBW,NRHS,KODE
      IF (IERR.EQ.(-2))  WRITE (LPU,6920) J
      IF (IERR.EQ.(-3))  WRITE (LPU,6930) J,T,A(NBW,J),EPS
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE BANSOL')
 6910 FORMAT(5X,'ILLEGAL PARAMETER(S) NEQ, NBW, NRHS AND/OR KODE'
     +      /4I10)
 6920 FORMAT(5X,'ZERO PIVOT ELEMENT IN EQUATION ',I6)
 6930 FORMAT(5X,'NEAR-SINGULAR SYSTEM DETECTED IN EQ. I =',I6 /
     +       5X,'A(I,I), D(I,I) AND EPS =',1P3E14.4)
C
      END
      SUBROUTINE COMSOL (A,B,MSKY,EPS,NEQ2,NRHS,LPU,IFLAG,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  COMSOL                 GROUP 4 / PUBLIC
C
C     THIS SUBROUTINE SOLVES A SYSTEM OF COMPLEX EQUATIONS,  A*X = B ,
C     WHERE  A  IS SYMMETRIC, THROUGH THREE STEPS :
C        1) FACTORIZATION OF A         -  A = L*D*(LT)
C        2) FORWARD REDUCTION OF B     -  Y = (D-1)*(L-1)*B
C        3) BACKSUBSTITUTION TO GIVE X -  X = (L-T)*Y
C     STORAGE-WISE  LT AND D REPLACES A, Y REPLACES B AND X REPLACES Y.
C     THE ACTIVE COLUMNS OF THE UPPER TRIANGULAR PART OF  A  ARE STORED
C     CONSECUTIVELY IN ARRAY A  (LOCATION OF DIAGONAL ELEMENT NO. I IS
C     STORED IN MSKY(I))
C     DEPENDING ON THE VALUE OF THE OPERATION FLAG (IFLAG), ONE OF THE
C     THREE STEPS, THE FIRST TWO, THE LAST TWO, OR ALL THREE STEPS
C     MAY BE CARRIED OUT IN ONE REFERENCE TO COMSOL
C
C     NOTE 1 :  FORTRAN COMPLEX IS  N O T  USED, BUT STORAGE-WISE
C               ARRAYS  A  AND  B  ARE ARRANGED EXACTLY AS THEY WOULD
C               BE IF OF TYPE COMPLEX.
C
C     NOTE 2 :  ARGUMENT  NEQ2  IS TWICE THE NUMBER OF EQUATIONS
C
C     NOTE 3 :  COMSOL IS A MODIFIED VERSION OF  SKYSOL
C
C     ROUTINES CALLED/REFERENCED :  COMPAC        (SAM-4)
C                                   MIN AND SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   85-10-20 / 1.0
C                       86-06-04 / 1.1    K.BELL
C                       96-12-10 / 1.2    K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,LPU,NEQ2,NRHS,       MSKY(*)
      DOUBLE PRECISION  EPS, A(*),B(NEQ2,*)
C
      INTEGER           I,ID,IE,II,IIC,IJ,IJC,IS,J,JS,K,L,NEQ,NJ,NN
      DOUBLE PRECISION  CI,CR,RI,RR,S,T,ZERO
      DOUBLE PRECISION  COMPAC
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          COMPAC
C ----------------------------------------------------------------------
C  PRELIMINARIES  -  CHECK PARAMETERS AND DIRECT EXECUTION
C ----------------------------------------------------------------------
      IERR = 0
      IF (NEQ2.LT.2)         GO TO 910
      NEQ  = NEQ2/2
      IF (IFLAG.EQ.1)        GO TO 100
      IF (IFLAG.EQ.2)        GO TO 100
      IF (IFLAG.EQ.3)        GO TO 100
      IF (IFLAG.EQ.4)        GO TO 400
      IF (IFLAG.EQ.5)        GO TO 400
      IF (IFLAG.EQ.6)        GO TO 600
      GO TO 910
C ----------------------------------------------------------------------
C  STEP 1  -  FACTORIZATION
C ----------------------------------------------------------------------
  100 J = 1
      T = SQRT(A(1)**2 + A(2)**2)
      IF (T.EQ.ZERO)         GO TO 930
      IF (NEQ.EQ.1)          GO TO 350
C
      DO 300 J=2,NEQ
         NJ = MSKY(J) - MSKY(J-1)
         IF (NJ .EQ. 1) THEN
            ID = MSKY(J)
            T  = SQRT(A(2*ID-1)**2 + A(2*ID)**2)
            IF (T .EQ. ZERO) THEN
               GO TO 930
            ELSE
               GO TO 300
            ENDIF
         ENDIF
         K  = J-NJ+1
         IF (K.LT.1)         GO TO 920
C                                               ** FORM COLUMN J OF  U
         I  = K-1
         IF (I.EQ.0)         I=1
C
  200    I  = I+1
         NN = MIN((NJ-J+I),(MSKY(I)-MSKY(I-1))) - 1
         IF (NN.LT.1)        GO TO 200
         IS = MSKY(I)-NN
         IJ = MSKY(J)-J+I
         IJC = 2*IJ-1
C
         IF (I.EQ.J)         GO TO 250
C
         JS = IJ-NN
         A(IJC)   = A(IJC)   - COMPAC(A(2*IS-1),A(2*JS-1),1,1,NN,1)
         A(IJC+1) = A(IJC+1) - COMPAC(A(2*IS-1),A(2*JS-1),1,1,NN,2)
         GO TO 200
C                                               ** FORM COLUMN J OF  LT
C                                                  AND  D(J,J)
  250    IE = MSKY(I)-1
         IF (IE.LT.IS)       GO TO 300
         RR = ZERO
         RI = ZERO
         DO 275 II=IS,IE
            ID       = MSKY(K)
            IIC      = 2*II-1
            CR       = A(IIC)
            CI       = A(IIC+1)
            T        = A(2*ID-1)**2 + A(2*ID)**2
            A(IIC)   = (CR*A(2*ID-1) + CI*A(2*ID))/T
            A(IIC+1) = (CI*A(2*ID-1) - CR*A(2*ID))/T
            RR       = RR + CR*A(IIC) - CI*A(IIC+1)
            RI       = RI + CR*A(IIC+1) + CI*A(IIC)
            K        = K + 1
  275    CONTINUE
         T        = SQRT(A(IJC)**2 + A(IJC+1)**2)
         A(IJC)   = A(IJC)   - RR
         A(IJC+1) = A(IJC+1) - RI
C                                               ** CHECK FOR ILL COND.
         S = SQRT(A(IJC)**2 + A(IJC+1)**2)
         IF (S.EQ.ZERO)      GO TO 930
         IF (S.LT.(EPS*T))   GO TO 940
  300 CONTINUE
C
  350 IF (IFLAG.EQ.1)        GO TO 1000
C ----------------------------------------------------------------------
C  STEP 2  -  FORWARD REDUCTION
C ----------------------------------------------------------------------
  400 IF (NRHS.LT.1)         GO TO 910
      DO 500 L=1,NRHS
         IF (NEQ.EQ.1)       GO TO 450
         DO 425 I=2,NEQ
            NN = MSKY(I) - MSKY(I-1) - 1
            IF (NN.LT.1)     GO TO 425
            IS = MSKY(I) - NN
            K  = I-NN
            B(2*I-1,L) = B(2*I-1,L) - COMPAC(A(2*IS-1),B(2*K-1,L),
     +                                                     1,1,NN,1)
            B(2*I  ,L) = B(2*I  ,L) - COMPAC(A(2*IS-1),B(2*K-1,L),
     +                                                     1,1,NN,2)
  425    CONTINUE
  450    DO 475 I=1,NEQ
            ID = MSKY(I)
            RR = B(2*I-1,L)
            RI = B(2*I  ,L)
            T  = A(2*ID-1)**2 + A(2*ID)**2
            B(2*I-1,L) = (RR*A(2*ID-1) + RI*A(2*ID))/T
            B(2*I  ,L) = (RI*A(2*ID-1) - RR*A(2*ID))/T
  475    CONTINUE
  500 CONTINUE
C
      IF (IFLAG.EQ.2)        GO TO 1000
      IF (IFLAG.EQ.5)        GO TO 1000
C ----------------------------------------------------------------------
C  STEP 3  -  BACKSUBSTITUTION
C ----------------------------------------------------------------------
  600 IF (NRHS.LT.1)         GO TO 910
      DO 700 L=1,NRHS
         I  = NEQ
  625    IF (I.EQ.1)         GO TO 700
         K  = I - MSKY(I) + MSKY(I-1) + 1
         IS = MSKY(I-1) + 1
         IE = MSKY(I) - 1
         IF (IE.LT.IS)       GO TO 675
         DO 650 II=IS,IE
            B(2*K-1,L) = B(2*K-1,L) - A(2*II-1)*B(2*I-1,L)
     +                              + A(2*II)*B(2*I,L)
            B(2*K  ,L) = B(2*K  ,L) - A(2*II-1)*B(2*I,L)
     +                              - A(2*II)*B(2*I-1,L)
            K = K+1
  650    CONTINUE
  675    I = I-1
         GO TO 625
  700 CONTINUE
C
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      GO TO 950
  930 IERR =-3
      GO TO 950
  940 IERR =-4
  950 IF (LPU.LT.0)          GO TO 1000
C
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))      WRITE (LPU,6910) NEQ2,NRHS,IFLAG
      IF (IERR.EQ.(-2))      WRITE (LPU,6920) J,MSKY(J-1),MSKY(J)
      IF (IERR.EQ.(-3))      WRITE (LPU,6930) J
      IF (IERR.EQ.(-4))      WRITE (LPU,6940) J,EPS,T,S
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE COMSOL')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  NEQ2    NRHS   IFLAG'/
     +                                     I33,     I8,     I8)
 6920 FORMAT(/5X,'INCORRECT ENTRIES IN MSKY :  J   MSKY(J-1)  MSKY(J)'/
     +                                       I35,       I11,      I9 )
 6930 FORMAT(/5X,'ZERO PIVOT ELEMENT IN EQUATION',I6)
 6940 FORMAT(/5X,'EQUATION',I6,'  FAILED THE STABILITY TEST'
     +       /5X,'EPS =',1PE12.4,5X,'ABS(A(I,I)) =',1PE12.4,5X,
     +        'ABS(D(I,I)) =',1PE12.4 )
C
      END
      DOUBLE PRECISION FUNCTION COMPAC(A,B,INCA,INCB,N,IRI)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  COMPAC                  GROUP 4 / PRIVATE
C
C     FORMS THE REAL (IRI=1) OR THE IMAGINARY (IRI=2) PART OF THE INNER
C     (DOT) PRODUCT OF TWO COMPLEX ARRAYS  A  AND  B.
C     FORTRAN COMPLEX IS  N O T  USED.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   85-10-10 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           INCA,INCB,IRI,N
      DOUBLE PRECISION  A(*),B(*)
C
      INTEGER           I,IA,IB,LINCA,LINCB
      DOUBLE PRECISION  ACC,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      ACC   = ZERO
      IF (N.LT.1)       GO TO 100
C
      LINCA = 2*INCA
      LINCB = 2*INCB
      IA    = 1
      IB    = 1
      IF (IRI.EQ.2)     GO TO 20
C                                               ** REAL PART
      DO 10 I=1,N
         ACC = ACC + A(IA)*B(IB) - A(IA+1)*B(IB+1)
         IA  = IA + LINCA
         IB  = IB + LINCB
   10 CONTINUE
      GO TO 100
C                                               ** IMAGINARY PART
   20 DO 30 I=1,N
         ACC = ACC + A(IA)*B(IB+1) + B(IB)*A(IA+1)
         IA  = IA + LINCA
         IB  = IB + LINCB
   30 CONTINUE
C
  100 COMPAC = ACC
C
      RETURN
      END
      SUBROUTINE SOLDLT (A,B,EPS,NEQ,NRHS,LPU,KODE,IERR)
C
C **********************************************************************
C
C   S A M  library routine :  SOLDLT                 GROUP 4 / PUBLIC
C
C     TASK :  To solve a system of linear equations  A*X = B through
C             three steps:
C             1) Factorization of A  ( A = L*D*Lt = L*U )
C             2) Forward reduction   ( Z = (L-1)*B )
C             3) Backsubstitution    ( X = (L-t)*(D-1)*Z = (U-1)*Z )
C     A is a symmetric matrix (stored in full), L is a lower triangular
C     matrix with unit elements on the diagonal and D is a diagonal
C     matrix  ( U = D*Lt is an upper triangular matrix).
C     Depending on the value of argument KODE,
C             - one of the three steps,
C             - the last two steps, or
C             - all three steps
C     may be carried out in one reference to SOLDLT.
C
C     ROUTINES CALLED/REFERENCED :  ABS   (Fortran library)
C
C
C     PROGRAMMED BY  :  Kolbein Bell
C     DATE/VERSION   :  90-03-27 / 1.0
C
C **********************************************************************
C     Conditionals   :  S - Single precision     D - Double precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,KODE,LPU,NEQ,NRHS
      DOUBLE PRECISION  A(NEQ,NEQ),B(NEQ,*),EPS
C
      INTEGER           I,J,K,L
      DOUBLE PRECISION  X,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
C
C --- check input parameters
C
      IF (KODE.LT.1 .OR.  KODE.GT.5)  GO TO 91
      IF (KODE.GT.1 .AND. NRHS.LT.1)  GO TO 91
      IF (NEQ .LT. 1)                 GO TO 91
C
      IF (KODE .LT. 3) THEN
C
C ------ step 1  -  factorization (KODE = 1 or 2)
C
         J = 1
         IF (A(1,1) .EQ. ZERO)  GO TO 92
         IF (NEQ .GT. 1) THEN
            J = 2
            X = A(2,2)
            A(2,1) = A(1,2) / A(1,1)
            A(2,2) = A(2,2) - A(2,1)*A(1,2)
C                                                ! ill conditioned ?
C
            IF (A(2,2) .EQ. ZERO)               GO TO 92
            IF (ABS(A(2,2)) .LT. (EPS*ABS(X)))  GO TO 93
C
            IF (NEQ .GT. 2) THEN
               DO 40 J=3,NEQ
C                                                ! form column J of U
                  DO 10 I=2,J-1
                     DO 5 K=1,I-1
                        A(I,J) = A(I,J) - A(I,K)*A(K,J)
    5                CONTINUE
   10             CONTINUE
C                                                ! form row J of L
                  DO 20 I=1,J-1
                     A(J,I) = A(I,J) / A(I,I)
   20             CONTINUE
C                                                ! form and check Djj
                  X = A(J,J)
                  DO 30 K=1,J-1
                     A(J,J) = A(J,J) - A(J,K)*A(K,J)
   30             CONTINUE
                  IF (A(J,J) .EQ. ZERO)               GO TO 92
                  IF (ABS(A(J,J)) .LT. (EPS*ABS(X)))  GO TO 93
   40          CONTINUE
            ENDIF
         ENDIF
      ENDIF
C
      IF (KODE.GT.1 .AND. KODE.LT.5) THEN
C
C ------ step 2  -  forward reduction  (KODE = 2, 3 or 4)
C
         IF (NEQ .GT. 1) THEN
            DO 60 J=1,NRHS
               DO 55 I=2,NEQ
                  DO 50 K=1,I-1
                     B(I,J) = B(I,J) - A(I,K)*B(K,J)
   50             CONTINUE
   55          CONTINUE
   60       CONTINUE
         ENDIF
      ENDIF
C
      IF (KODE.EQ.2 .OR. KODE.EQ.3 .OR. KODE.EQ.5) THEN
C
C ------ step 3  -  backsubstitution  (KODE = 2, 3 or 5)
C
         DO 80 J=1,NRHS
            B(NEQ,J) = B(NEQ,J) / A(NEQ,NEQ)
            IF (NEQ .GT. 1) THEN
               DO 75 L=1,NEQ-1
                  I = NEQ-L
                  DO 70 K=I+1,NEQ
                     B(I,J) = B(I,J) - A(I,K)*B(K,J)
   70             CONTINUE
                  B(I,J) = B(I,J) / A(I,I)
   75          CONTINUE
            ENDIF
   80    CONTINUE
      ENDIF
C
      GO TO 100
C
C ------------------------------------------------ error exit
C
   91 IERR =-1
      GO TO 95
   92 IERR =-2
      GO TO 95
   93 IERR =-3
C
   95 IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
         IF (IERR.EQ.(-1) )  WRITE (LPU,691) NEQ,NRHS,KODE
         IF (IERR.EQ.(-2))   WRITE (LPU,692) J
         IF (IERR.EQ.(-3))   WRITE (LPU,693) J,X,A(J,J),EPS
      ENDIF
C ----------------------------------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  690 FORMAT(///' *** ERROR return from  S A M  library routine SOLDLT')
  691 FORMAT(5X,'Illegal parameter(s) NEQ, NRHS and/or KODE:',3I8)
  692 FORMAT(5X,'Zero pivot element in equation',I6)
  693 FORMAT(5X,'Near-singular system detected in eq. I =',I6 /
     +       5X,'A(I,I), D(I,I) and EPS =',1P3E14.4)
C
      END
      SUBROUTINE FACHOL (A,EPS,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  FACHOL                 GROUP 4 / PUBLIC
C
C     TASK :  To factorize a symmetric, positive definite N by N matrix
C             A into its Cholesky factors L*LT.
C             On successful exit A stores LT on and above the diagonal.
C
C
C     ROUTINES CALLED/REFERENCED :   PRACC           (SAM-0)
C                                    SQRT            (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-06 / 1.0
C                       96-06-27 / 1.1  K.Bell
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LPU,N
      DOUBLE PRECISION  EPS, A(N,N)
C                                                ! local variables
      INTEGER           I,I1,J
      DOUBLE PRECISION  X,ZERO,          PRACC
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IERR = 0
      IF (EPS .LT. ZERO) THEN
         EPS = 1.0D-10
      ENDIF
      IF (N .GT. 0) THEN
         IF (A(1,1) .GT. ZERO) THEN
            A(1,1) = SQRT(A(1,1))
            IF (N .GT. 1) THEN
               DO 50 J=2,N
                  IF (A(J,J) .LE. ZERO) THEN
                     IERR = -2
                     GO TO 90
                  ENDIF
                  A(1,J) = A(1,J) / A(1,1)
                  DO 25 I=2,J-1
                     I1 = I-1
                     X = A(I,J) - PRACC(A(1,I),A(1,J),1,1,I1)
                     A(I,J) = X / A(I,I)
   25             CONTINUE
                  I1 = J-1
                  X = A(J,J) - PRACC(A(1,J),A(1,J),1,1,I1)
                  IF (X .LE. ZERO) THEN
                     IERR = -2
                     GO TO 90
                  ELSEIF (X .LT. EPS*A(J,J)) THEN
                     IERR = -3
                     GO TO 90
                  ELSE
                     A(J,J) = SQRT(X)
                  ENDIF
   50          CONTINUE
            ENDIF
         ELSE
            J    = 1
            IERR =-2
            GO TO 90
         ENDIF
      ELSE
         IERR =-1
         GO TO 90
      ENDIF
      GO TO 100
C ----------------------------------------------- error exit
   90 IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
         IF (IERR .EQ. (-1))  WRITE (LPU,691) N
         IF (IERR .EQ. (-2))  WRITE (LPU,692) J
         IF (IERR .EQ. (-3))  WRITE (LPU,693) J
      ENDIF
C ----------------------------------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  690 FORMAT(///' *** ERROR RETURN from  S A M  LIBRARY routine FACHOL')
  691 FORMAT(/5X,'Illegal matrix dimension:   N =',I6)
  692 FORMAT(/5X,'Matrix is not pos.def.  (eqn.',I5,' )')
  693 FORMAT(/5X,'Singular or near-singular matrix  (eqn.',I5,' )')
C
      END
      SUBROUTINE FSCHOL (A,B,N,NRHS,KODE)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  FSCHOL                 GROUP 4 / PUBLIC
C
C     TASK :  For KODE .ge. 0 :
C
C                Solve a lower triangular system (forward substitution)
C                     L X = B   >  X = (L-1) B
C
C             For KODE .lt. 0 and NRHS=N (B is square):
C
C                Solve a lower triangular system (forward substitution)
C                such that
C                     X = XT = (L-1) BT
C                In other words, the right-hand side is B-transpose, and
C                the solution matrix (stored in B) is symmetric.
C             A contains the transpose of L (which is a Cholesky factor)
C             on and above the diagonal (as returned from FACHOL).
C
C
C     ROUTINES CALLED/REFERENCED :   PRACC      (SAM-0)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-07 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           KODE,N,NRHS
      DOUBLE PRECISION  A(N,N),B(N,NRHS)
C                                                ! local variables
      INTEGER           I,I1,J
      DOUBLE PRECISION  X,     PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
C
      IF (N .GT. 0) THEN
         IF (KODE .GE. 0) THEN
C
C --------  Ordinary solution (forward substitution)
C
            DO 100 J=1,NRHS
               B(1,J) = B(1,J) / A(1,1)
               DO 50 I=2,N
                  I1 = I-1
                  X      = B(I,J) - PRACC(A(1,I),B(1,J),1,1,I1)
                  B(I,J) = X / A(I,I)
   50          CONTINUE
  100       CONTINUE
         ELSEIF (KODE.LT.0 .AND. NRHS.EQ.N) THEN
C
C --------  Right-hand side is square (=BT) and solution is symmetric
C
            B(1,1) = B(1,1) / A(1,1)
            DO 200 J=2,N
               B(J,1) = B(J,1) / A(1,1)
               DO 150 I=2,J
                  I1 = I-1
                  X      = B(J,I) - PRACC(A(1,I),B(J,1),1,N,I1)
                  B(J,I) = X / A(I,I)
  150          CONTINUE
  200       CONTINUE
            DO 300 I=1,N
               DO 250 J=I,N
                  B(I,J) = B(J,I)
  250          CONTINUE
  300       CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE BSCHOL (A,B,N,NRHS,KODE)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  BSCHOL                 GROUP 4 / PUBLIC
C
C     TASK :  For KODE .ge. 0 :
C
C               Solve an upper triangular system (backward substitution)
C                     LT X = B   >  X = (L-T) B
C
C             For KODE .lt. 0 and NRHS=N (B is square):
C
C               Solve an upper triangular system (backward substitution)
C               such that
C                     X = XT = (L-T) B
C                In other words, the right-hand side is a square matrix
C                B and the solution matrix (stored in B) is symmetric.
C             A contains the transpose of L (which is a Cholesky factor)
C             on and above the diagonal (as returned from FACHOL).
C
C
C     ROUTINES CALLED/REFERENCED :   PRACC      (SAM-0)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-07 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           KODE,N,NRHS
      DOUBLE PRECISION  A(N,N),B(N,NRHS)
C                                                ! local variables
      INTEGER           I,I1,J,NI
      DOUBLE PRECISION  X,     PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
C
      IF (N .GT. 0) THEN
         IF (KODE .GE. 0) THEN
C
C --------  Ordinary solution (backward substitution)
C
            DO 100 J=1,NRHS
               B(N,J) = B(N,J) / A(N,N)
               DO 50 I=N-1,1,-1
                  I1 = I+1
                  NI = N-I
                  X      = B(I,J) - PRACC(A(I,I1),B(I1,J),N,1,NI)
                  B(I,J) = X / A(I,I)
   50          CONTINUE
  100       CONTINUE
         ELSEIF (KODE.LT.0 .AND. NRHS.EQ.N) THEN
C
C --------  Right-hand side is square and solution is symmetric
C
            DO 200 J=1,N
               B(N,J) = B(N,J) / A(N,N)
               DO 150 I=N-1,J,-1
                  I1 = I+1
                  NI = N-I
                  X      = B(I,J) - PRACC(A(I,I1),B(I1,J),N,1,NI)
                  B(I,J) = X / A(I,I)
  150          CONTINUE
  200       CONTINUE
            DO 300 I=1,N
               DO 250 J=I,N
                  B(I,J) = B(J,I)
  250          CONTINUE
  300       CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE FACSKY (A,MSKY,EPS,N,LPU,IEQ,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  FACSKY                 GROUP 4 / PUBLIC
C
C     TASK :  To factorize a symmetric, positive definite N by N matrix
C             A into its Cholesky factors UT*U.
C             A is stored in "skyline" format, i.e., the active columns
C             of the upper triangular part of A are stored consecutively
C             in array A - location of diagonal element no. I is stored
C             in MSKY(I).
C             On successful exit A stores U.
C
C
C     ROUTINES CALLED/REFERENCED :   PRACC           (SAM-0)
C                                    MIN, SQRT       (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-14 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IEQ,LPU,N,     MSKY(*)
      DOUBLE PRECISION  EPS, A(*)
C                                                ! local variables
C
      INTEGER           I,IJ,IS,J,JJ,JS,K,NJ,NN
      DOUBLE PRECISION  X,ZERO,                 PRACC
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IERR = 0
      IEQ  = 0
C
      IF (N .GT. 0) THEN
         IF (A(1) .GT. ZERO) THEN
            A(1) = SQRT(A(1))
            IF (N .GT. 1) THEN
               DO 50 J=2,N
                  JJ = MSKY(J)
                  IF (A(JJ) .LE. ZERO) THEN
                     IERR = -2
                     IEQ  = J
                     GO TO 90
                  ENDIF
                  NJ = JJ - MSKY(J-1)
                  IF (NJ .GT. 1) THEN
                     K = J-NJ+1
                     IF (K .GT. 0) THEN
C
C --------------------- form column J of U
C
                        I = K-1
C                                                ! start row loop
   10                   I = I+1
                        IF (I .GT. 1) THEN
                           NN = MIN((NJ-J+I),(MSKY(I)-MSKY(I-1))) - 1
                           IS = MSKY(I)-NN
                           IJ = MSKY(J)-J+I
                           JS = IJ-NN
C
                           IF (NN .GT. 0) THEN
                              X = A(IJ) - PRACC(A(IS),A(JS),1,1,NN)
                           ELSE
                              X = A(IJ)
                           ENDIF
C
                           IF (I .LT. J) THEN
                              A(IJ) = X / A(MSKY(I))
                           ELSE
                              IF (X .GT. EPS*A(IJ)) THEN
                                 A(IJ) = SQRT(X)
                              ELSEIF (X .LE. ZERO) THEN
                                 IERR = -2
                                 IEQ  = J
                                 GO TO 90
                              ELSE
                                 IERR = -3
                                 IEQ  = J
                                 GO TO 90
                              ENDIF
C                                                ! exit row loop
                              GO TO 50
                           ENDIF
                        ELSE
                           IJ    = MSKY(J-1) + 1
                           A(IJ) = A(IJ) / A(1)
                        ENDIF
                        GO TO 10
C                                                ! end row loop
                     ELSE
                        IERR =-4
                        IEQ  = J
                        GO TO 90
                     ENDIF
                  ELSE
                     A(JJ) = SQRT(A(JJ))
                  ENDIF
   50          CONTINUE
            ENDIF
         ELSE
            IERR = -2
            IEQ  = 1
            GO TO 90
         ENDIF
      ELSE
         IERR = -1
         GO TO 90
      ENDIF
      GO TO 100
C ----------------------------------------------- error exit
   90 IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
         IF (IERR .EQ. (-1))  WRITE (LPU,691) N
         IF (IERR .EQ. (-2))  WRITE (LPU,692) IEQ
         IF (IERR .EQ. (-3))  WRITE (LPU,693) IEQ
         IF (IERR .EQ. (-4))  WRITE (LPU,694)
      ENDIF
C ----------------------------------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  690 FORMAT(///' *** ERROR RETURN from  S A M  LIBRARY routine FACSKY')
  691 FORMAT(/5X,'Illegal matrix dimension:   N =',I6)
  692 FORMAT(/5X,'Matrix is not pos.def.  (eqn.',I5,' )')
  693 FORMAT(/5X,'Singular or near-singular matrix  (eqn.',I5,' )')
  694 FORMAT(/5X,'Incorrect entries in MSKY :  J   MSKY(J-1)  MSKY(J)'/
     &       22X,                            I10,       I11,      I9)
C
      END
      SUBROUTINE FSBSKY (A,B,MSKY,N,NRHS)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  FSBSKY                 GROUP 4 / PUBLIC
C
C     TASK :  To solve a lower triangular system (forward substitution)
C
C                     UT X = B   >  X = (U-T) B
C
C             for NRHS right-hand sides.
C             A contains the Cholesky factor U in "skyline" storage
C             format (as returned from FACSKY).
C
C
C     ROUTINES CALLED/REFERENCED :   PRACC      (SAM-0)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-16 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           N,NRHS,        MSKY(*)
      DOUBLE PRECISION  A(*),B(N,NRHS)
C                                                ! local variables
      INTEGER           I,ID,IS,K,L,NN
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
C
      IF (N .GT. 0) THEN
         DO 100 L=1,NRHS
            B(1,L) = B(1,L) / A(1)
            DO 50 I=2,N
               ID = MSKY(I)
               NN = ID - MSKY(I-1) - 1
               IF (NN .GT. 0) THEN
                  IS = ID - NN
                  K  = I-NN
                  B(I,L) = B(I,L) - PRACC(A(IS),B(K,L),1,1,NN)
               ENDIF
               B(I,L) = B(I,L) / A(ID)
   50       CONTINUE
  100    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE BCKSKY (A,B,MSKY,N,NRHS)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  BCKSKY                 GROUP 4 / PUBLIC
C
C     TASK :  To solve an upper triangular system (backsubstitution)
C
C                     U X = B   >  X = (U-1) B
C
C             for NRHS right-hand sides.
C             A contains the Cholesky factor U stored in "skyline"
C             format (as returned from FACSKY).
C
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-16 / 1.0
C
C **********************************************************************
C     CONDITIONALS   :  S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           N,NRHS,         MSKY(*)
      DOUBLE PRECISION  A(*),B(N,NRHS)
C                                                ! local variables
      INTEGER           I,IE,IS,J,K,L
C ----------------------------------------------------------------------
C
      IF (N .GT. 0) THEN
         DO 100 L=1,NRHS
            DO 50 J=N,2,-1
               B(J,L) = B(J,L) / A(MSKY(J))
               K      = J - MSKY(J) + MSKY(J-1) + 1
               IS     = MSKY(J-1) + 1
               IE     = MSKY(J) - 1
               DO 25 I=IS,IE
                  B(K,L) = B(K,L) - A(I)*B(J,L)
                  K      = K+1
   25          CONTINUE
   50       CONTINUE
            B(1,L) = B(1,L) / A(1)
  100    CONTINUE
      ENDIF
C
      RETURN
      END
