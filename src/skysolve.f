C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SKYSOLVE (A,B,MSKY,EPS,SCALE,NEQ,NEQ1,NRHS,LPU,IFLAG,
     +                     DCOR,NND,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE  :  SKYSOLVE                 GROUP 4 / PUBLIC
C
C     THIS SUBROUTINE SOLVES (NEQ1=NEQ) OR REDUCES (NEQ1.LT.NEQ) A
C     SYSTEM OF EQUATIONS A*X = B (A IS SYMMETRIC) THROUGH THREE STEPS:
C        1) FACTORIZATION OF A         -  A = L*D*(LT)
C        2) FORWARD REDUCTION OF B     -  Y = (D-1)*(L-1)*B
C        3) BACKSUBSTITUTION TO GIVE X -  X = (L-T)*Y
C     STORAGE-WISE LT AND D REPLACES A, Y REPLACES B AND X REPLACES Y.
C     THE ACTIVE COLUMNS OF THE UPPER TRIANGULAR PART OF A ARE STORED
C     CONSECUTIVELY IN ARRAY A (LOCATION OF DIAGONAL ELEMENT NO. I IS
C     STORED IN MSKY(I))
C     DEPENDING ON THE VALUE OF THE OPERATION FLAG (IFLAG), ONE OF THE
C     THREE STEPS, THE FIRST TWO, THE LAST TWO, OR ALL THREE STEPS MAY
C     BE CARRIED OUT IN ONE REFERENCE TO SKYSOLVE
C
C     ROUTINES CALLED/REFERENCED :  DDOT, DAXPY  (BLAS)
C                                   SPRF3P       (SAM)
C
C     PROGRAMMED BY :   KOLBEIN BELL   (AFTER E.L.WILSON AND H.H.DOVEY)
C     DATE/VERSION  :   80-04-29 / 1.0
C                       03-07-29 / 2.0   K.M.Okstad (added DCOR, etc).
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IERR,IFLAG,LPU,NEQ,NEQ1,NND,NRHS,MSKY(NEQ)
      DOUBLE PRECISION  EPS,SCALE,A(*),B(NEQ,NRHS),DCOR(2,NEQ)
C
      INTEGER           I,ID,IE,II,IJ,IS,J,JS,K,L,NJ,NN
      DOUBLE PRECISION  R,T,ZERO,ONE, DDOT
C
      PARAMETER       ( ZERO = 0.0D0, ONE = 1.0D0 )
C
      EXTERNAL          DDOT, DAXPY, SPRF3P
C
C ----------------------------------------------------------------------
C  PRELIMINARIES  -  CHECK PARAMETERS AND DIRECT EXECUTION
C ----------------------------------------------------------------------
      IERR = 0
      IF (NEQ  .LT. 1)                    GO TO 910
      IF (NEQ1 .LT. 1 .OR. NEQ1 .GT. NEQ) GO TO 910
C
      GO TO (100,100,100,400,400,600) ABS(IFLAG)
      GO TO 910
C
C ----------------------------------------------------------------------
C  STEP 1  -  FACTORIZATION
C ----------------------------------------------------------------------
  100 DO 150 J = 1, NEQ
         DCOR(1,1) = ZERO
         DCOR(2,1) = ZERO
         ID = MSKY(J)
         IF (A(ID) .EQ. ZERO) THEN
C     TODO: Check that the entire column (and row) actually is zero.
            DCOR(2,J) = ONE
            A(ID) = DCOR(2,J)
            IERR = IERR - 1
         END IF
  150 CONTINUE
C
      DO 300 J = 2, NEQ
         NJ = MSKY(J) - MSKY(J-1)
         IF (NJ .EQ. 1) THEN
            GO TO 300
         ELSE IF (NJ .GT. J .OR. NJ .LT. 1) THEN
            GO TO 920
         END IF
C                                               ** FORM COLUMN J OF  U
         K = J - NJ + 1
         IF (K .LE. 2) THEN
            I = 1
         ELSE
            I = K - 1
         END IF
C
  200    I  = I + 1
         NN = MIN(NJ-J+I,MSKY(I)-MSKY(I-1)) - 1
         IS = MSKY(I) - NN
         IJ = MSKY(J) - J+I
         IF (I.EQ.J)         GO TO 250
C
         JS = IJ - NN
         IF (I .GT. NEQ1)    NN = NN + NEQ1-I+1
         IF (NN .GT. 0)      A(IJ) = A(IJ) - DDOT(NN,A(IS),1,A(JS),1)
         GO TO 200
C                                               ** FORM COLUMN J OF  LT
C                                                  AND  D(J,J)
  250    IE = MSKY(I)-1
         IF (I .GT. NEQ1)    IE = IE + NEQ1-I+1
         IF (IE .LT. IS)     GO TO 300
         R = ZERO
         DO 275 II = IS, IE
            ID    = MSKY(K)
            T     = A(II)
            A(II) = A(II)/A(ID)
            R     = R + A(II)*T
            K     = K + 1
  275    CONTINUE
         T     = ABS(A(IJ))
         A(IJ) = A(IJ)-R
         IF (J .GT. NEQ1)    GO TO 300
C                                               ** CHECK FOR ILL COND.
         IF (ABS(A(IJ)) .LT. EPS*T) THEN
            DCOR(1,J) = A(IJ)
            DCOR(2,J) = SIGN(EPS*MAX(ONE,SCALE*T),A(IJ))
C     TODO: Insert DIAG = ONE instead and zero-out the whole column/row.
            A(IJ) = DCOR(2,J)
            IERR = IERR - 1
         END IF
  300 CONTINUE
C
      IF (ABS(IFLAG) .EQ. 1) GO TO 800
C ----------------------------------------------------------------------
C  STEP 2  -  FORWARD REDUCTION
C ----------------------------------------------------------------------
  400 IF (NRHS .LT. 1)       GO TO 910
      DO 500 L = 1, NRHS
         IF (IERR .LT. 0 .AND. IFLAG .LT. 0) THEN
C                                               ** SUPPRESS SINGULAR EQS
            DO 410 I = 1, NEQ1
               IF (DCOR(2,I) .NE. ZERO) B(I,L) = ZERO
  410       CONTINUE
         END IF
         DO 425 I = 2, NEQ
            IS = MSKY(I-1) + 1
            NN = MSKY(I) - IS
            K  = I - NN
            IF (I .GT. NEQ1) NN = NN+NEQ1-I+1
            IF (NN .GT. 0)   B(I,L) = B(I,L) - DDOT(NN,A(IS),1,B(K,L),1)
  425    CONTINUE
         DO 475 I = 1, NEQ1
            B(I,L) = B(I,L)/A(MSKY(I))
  475    CONTINUE
  500 CONTINUE
C
      IF (ABS(IFLAG) .EQ. 2) GO TO 800
      IF (ABS(IFLAG) .EQ. 5) GO TO 900
C ----------------------------------------------------------------------
C  STEP 3  -  BACKSUBSTITUTION
C ----------------------------------------------------------------------
  600 IF (NRHS .LT. 1)       GO TO 910
      DO 700 L = 1, NRHS
         DO 650 I = NEQ, 2, -1
            IS = MSKY(I-1) + 1
            NN = MSKY(I) - IS
            K  = I - NN
            IF (I .GT. NEQ1) NN = NN + NEQ1-I+1
            IF (NN .GT. 0)   CALL DAXPY (NN,-B(I,L),A(IS),1,B(K,L),1)
  650    CONTINUE
  700 CONTINUE
C
      IF (ABS(IFLAG) .GT. 3) GO TO 900
C ----------------------------------------------------------------------
C  CHECK FOR NEGATIVE DIAGONAL ELEMENTS
C ----------------------------------------------------------------------
  800 NND = 0
      DO 850 I = 1, NEQ
         IF (A(MSKY(I)) .LT. ZERO) NND = NND+1
  850 CONTINUE
      IF (NND.GT.0 .AND. IERR.EQ.0) IERR = 3
  900 IF (IERR .LT. 0) THEN
         NN = -IERR
         IERR = -3
         GO TO 950
      END IF
      RETURN
C ------------------------------------------------ ERROR EXIT
  910 IERR = -1
      GO TO 950
  920 IERR = -2
      NND = J
  950 IF (LPU .LT. 0)        RETURN
C
      WRITE (LPU,6900)
      IF (IERR .EQ. -1)      WRITE (LPU,6910) NEQ,NEQ1,NRHS,IFLAG
      IF (IERR .EQ. -2)      WRITE (LPU,6920) J,MSKY(J-1),MSKY(J)
      IF (IERR .EQ. -3)      THEN
         WRITE (LPU,6930)
         CALL SPRF3P (NEQ,NN,EPS,DCOR,LPU)
         IF (IFLAG .LT. 0)   WRITE (LPU,6931)
      END IF
      RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE SKYSOL')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  NEQ    NEQ1    NRHS   IFLAG'/
     +                                    I32,     I8,     I8,     I8)
 6920 FORMAT(/5X,'INCORRECT ENTRIES IN MSKY :  J   MSKY(J-1)  MSKY(J)'/
     +                                       I35,       I11,      I9 )
 6930 FORMAT(/5X,'SINGULAR OR POORLY CONDITIONED EQUATION SYSTEM')
 6931 FORMAT( 5X,'THE ABOVE EQUATIONS ARE AUTOMATICALLY SUPPRESSED'/)
C
      END
