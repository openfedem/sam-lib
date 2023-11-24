C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE CMVY18 (M,N,Y,APNT,A)
C
C***********************************************************************
C
C     Created: 11 Jan 1994 (ACD)
C
C     Purpose: Perform a matrix-vector multiplication Y = Y + A*X,
C              assuming data structures used in sparse CC' codes.
C              Uses level 8 loop unrolling.
C
C     Input:
C        M      - Number of rows in A
C        N      - Number of columns in A
C        Y      - Vector to which A*X will be added
C        APNT   - Index vector for A
C                 APNT(I) points to the first nonzero in column I of A
C        A      - Nonzero terms of A
C
C     Output:
C        Y      - Contains Y = Y + A*X
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER            M, N, APNT(*)
      DOUBLE PRECISION   Y(*), A(*)
C
      INTEGER            I, I1, I2, I3, I4, I5, I6, I7, I8, J, M1, REM1
      DOUBLE PRECISION   A1, A2, A3, A4, A5, A6, A7, A8
C
C***********************************************************************
C
      M1   = M + 1
      REM1 = MOD(N,8) + 1
C
      GOTO (800,100,200,300,400,500,600,700) REM1
C
  100 CONTINUE
      I1 = APNT(2) - M1
      A1 = -A(I1+1)
      DO 150 I = 1, M
         Y(I) = Y(I) + A1*A(I1+I)
  150 CONTINUE
      GOTO 800
C
  200 CONTINUE
      I1 = APNT(2) - M1
      I2 = APNT(3) - M1
      A1 = -A(I1+1)
      A2 = -A(I2+1)
      DO 250 I = 1, M
         Y(I) = (Y(I) + A1*A(I1+I)) + A2*A(I2+I)
  250 CONTINUE
      GOTO 800
C
  300 CONTINUE
      I1 = APNT(2) - M1
      I2 = APNT(3) - M1
      I3 = APNT(4) - M1
      A1 = -A(I1+1)
      A2 = -A(I2+1)
      A3 = -A(I3+1)
      DO 350 I = 1, M
         Y(I) = ((Y(I) + A1*A(I1+I)) + A2*A(I2+I)) + A3*A(I3+I)
  350 CONTINUE
      GOTO 800
C
  400 CONTINUE
      I1 = APNT(2) - M1
      I2 = APNT(3) - M1
      I3 = APNT(4) - M1
      I4 = APNT(5) - M1
      A1 = -A(I1+1)
      A2 = -A(I2+1)
      A3 = -A(I3+1)
      A4 = -A(I4+1)
      DO 450 I = 1, M
         Y(I) = (((Y(I) + A1*A(I1+I)) + A2*A(I2+I)) + A3*A(I3+I))
     +                  + A4*A(I4+I)
  450 CONTINUE
      GOTO 800
C
  500 CONTINUE
      I1 = APNT(2) - M1
      I2 = APNT(3) - M1
      I3 = APNT(4) - M1
      I4 = APNT(5) - M1
      I5 = APNT(6) - M1
      A1 = -A(I1+1)
      A2 = -A(I2+1)
      A3 = -A(I3+1)
      A4 = -A(I4+1)
      A5 = -A(I5+1)
      DO 550 I = 1, M
         Y(I) = ((((Y(I) + A1*A(I1+I)) + A2*A(I2+I)) + A3*A(I3+I))
     +                   + A4*A(I4+I)) + A5*A(I5+I)
  550 CONTINUE
      GOTO 800
C
  600 CONTINUE
      I1 = APNT(2) - M1
      I2 = APNT(3) - M1
      I3 = APNT(4) - M1
      I4 = APNT(5) - M1
      I5 = APNT(6) - M1
      I6 = APNT(7) - M1
      A1 = -A(I1+1)
      A2 = -A(I2+1)
      A3 = -A(I3+1)
      A4 = -A(I4+1)
      A5 = -A(I5+1)
      A6 = -A(I6+1)
      DO 650 I = 1, M
         Y(I) = (((((Y(I) + A1*A(I1+I)) + A2*A(I2+I)) + A3*A(I3+I))
     +                    + A4*A(I4+I)) + A5*A(I5+I)) + A6*A(I6+I)
  650 CONTINUE
      GOTO 800
C
  700 CONTINUE
      I1 = APNT(2) - M1
      I2 = APNT(3) - M1
      I3 = APNT(4) - M1
      I4 = APNT(5) - M1
      I5 = APNT(6) - M1
      I6 = APNT(7) - M1
      I7 = APNT(8) - M1
      A1 = -A(I1+1)
      A2 = -A(I2+1)
      A3 = -A(I3+1)
      A4 = -A(I4+1)
      A5 = -A(I5+1)
      A6 = -A(I6+1)
      A7 = -A(I7+1)
      DO 750 I = 1, M
         Y(I) = ((((((Y(I) + A1*A(I1+I)) + A2*A(I2+I)) + A3*A(I3+I))
     +                     + A4*A(I4+I)) + A5*A(I5+I)) + A6*A(I6+I))
     +                     + A7*A(I7+I)
  750 CONTINUE
      GOTO 800
C
  800 CONTINUE
      DO 1000 J = REM1, N, 8
         I1 = APNT(J+1) - M1
         I2 = APNT(J+2) - M1
         I3 = APNT(J+3) - M1
         I4 = APNT(J+4) - M1
         I5 = APNT(J+5) - M1
         I6 = APNT(J+6) - M1
         I7 = APNT(J+7) - M1
         I8 = APNT(J+8) - M1
         A1 = -A(I1+1)
         A2 = -A(I2+1)
         A3 = -A(I3+1)
         A4 = -A(I4+1)
         A5 = -A(I5+1)
         A6 = -A(I6+1)
         A7 = -A(I7+1)
         A8 = -A(I8+1)
         DO 900 I = 1, M
            Y(I) = (((((((Y(I) + A1*A(I1+I)) + A2*A(I2+I)) + A3*A(I3+I))
     +                         + A4*A(I4+I)) + A5*A(I5+I)) + A6*A(I6+I))
     +                         + A7*A(I7+I)) + A8*A(I8+I)
  900    CONTINUE
 1000 CONTINUE
C
      RETURN
      END
