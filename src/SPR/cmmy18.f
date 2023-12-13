C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE CMMY18 (M,N,Q,XPNT,X,Y,LDY)
C
C***********************************************************************
C
C     Created: 11 Jan 1994 (ACD)
C
C     Purpose: Perform a matrix-matrix multiplication, Y = Y + X*A,
C              assuming data structures used in the sparse CC' codes.
C              Uses level 8 loop unrolling.
C
C     Input:
C        M      - Number of rows in X and Y
C        N      - Number of columns in X and rows in A
C        Q      - Number of columns in A and Y
C        XPNT   - XPNT(J+1) points to one location beyond the end of the
C                 J-th column of X, also used to access the rows of A
C        X      - Contains the columns of X and the rows of A
C        Y      - Matrix to which X*A will be added
C        LDY    - Length of first column of Y
C
C     Output:
C        Y      - Contains Y = Y + X*A
C
C***********************************************************************
C
      IMPLICIT NONE
C
      INTEGER            M, N, Q, XPNT(*), LDY
      DOUBLE PRECISION   X(*), Y(*)
C
      INTEGER            I, I1, I2, I3, I4, I5, I6, I7, I8,
     +                   IYL, IY1, IY2, J, LY, M1, M2, REM1, YCOL
      DOUBLE PRECISION   A1, A2, A3, A4, A5, A6, A7, A8
C
C***********************************************************************
C
      M1   = M
      LY   = LDY
      IYL  = 0
      REM1 = MOD(N,8) + 1
C
C     Loop over each column of Y
      DO 2000 YCOL = 1, Q
C
C        Update index range variables
         IY1 = IYL + 1
         IY2 = IY1 + M1 - 1
         IYL = IYL + LY
         M2  = IY1 + M1
C
C        Perform approprate matrix vector multiplication X*A(*,YCOL)
         GOTO (800,100,200,300,400,500,600,700) REM1
C
  100    CONTINUE
         I1 = XPNT(2) - M2
         A1 = -X(I1+IY1)
         DO 150 I = IY1, IY2
            Y(I) = Y(I) + A1*X(I1+I)
  150    CONTINUE
         GOTO 800
C
  200    CONTINUE
         I1 = XPNT(2) - M2
         I2 = XPNT(3) - M2
         A1 = -X(I1+IY1)
         A2 = -X(I2+IY1)
         DO 250 I = IY1, IY2
            Y(I) = (Y(I) + A1*X(I1+I)) + A2*X(I2+I)
  250    CONTINUE
         GOTO 800
C
  300    CONTINUE
         I1 = XPNT(2) - M2
         I2 = XPNT(3) - M2
         I3 = XPNT(4) - M2
         A1 = -X(I1+IY1)
         A2 = -X(I2+IY1)
         A3 = -X(I3+IY1)
         DO 350 I = IY1, IY2
            Y(I) = ((Y(I) + A1*X(I1+I)) + A2*X(I2+I)) + A3*X(I3+I)
  350    CONTINUE
         GOTO 800
C
  400    CONTINUE
         I1 = XPNT(2) - M2
         I2 = XPNT(3) - M2
         I3 = XPNT(4) - M2
         I4 = XPNT(5) - M2
         A1 = -X(I1+IY1)
         A2 = -X(I2+IY1)
         A3 = -X(I3+IY1)
         A4 = -X(I4+IY1)
         DO 450 I = IY1, IY2
            Y(I) = (((Y(I) + A1*X(I1+I)) + A2*X(I2+I)) + A3*X(I3+I))
     +                     + A4*X(I4+I)
  450    CONTINUE
         GOTO 800
C
  500    CONTINUE
         I1 = XPNT(2) - M2
         I2 = XPNT(3) - M2
         I3 = XPNT(4) - M2
         I4 = XPNT(5) - M2
         I5 = XPNT(6) - M2
         A1 = -X(I1+IY1)
         A2 = -X(I2+IY1)
         A3 = -X(I3+IY1)
         A4 = -X(I4+IY1)
         A5 = -X(I5+IY1)
         DO 550 I = IY1, IY2
            Y(I) = ((((Y(I) + A1*X(I1+I)) + A2*X(I2+I)) + A3*X(I3+I))
     +                      + A4*X(I4+I)) + A5*X(I5+I)
  550    CONTINUE
         GOTO 800
C
  600    CONTINUE
         I1 = XPNT(2) - M2
         I2 = XPNT(3) - M2
         I3 = XPNT(4) - M2
         I4 = XPNT(5) - M2
         I5 = XPNT(6) - M2
         I6 = XPNT(7) - M2
         A1 = -X(I1+IY1)
         A2 = -X(I2+IY1)
         A3 = -X(I3+IY1)
         A4 = -X(I4+IY1)
         A5 = -X(I5+IY1)
         A6 = -X(I6+IY1)
         DO 650 I = IY1, IY2
            Y(I) = (((((Y(I) + A1*X(I1+I)) + A2*X(I2+I)) + A3*X(I3+I))
     +                       + A4*X(I4+I)) + A5*X(I5+I)) + A6*X(I6+I)
  650    CONTINUE
         GOTO 800
C
  700    CONTINUE
         I1 = XPNT(2) - M2
         I2 = XPNT(3) - M2
         I3 = XPNT(4) - M2
         I4 = XPNT(5) - M2
         I5 = XPNT(6) - M2
         I6 = XPNT(7) - M2
         I7 = XPNT(8) - M2
         A1 = -X(I1+IY1)
         A2 = -X(I2+IY1)
         A3 = -X(I3+IY1)
         A4 = -X(I4+IY1)
         A5 = -X(I5+IY1)
         A6 = -X(I6+IY1)
         A7 = -X(I7+IY1)
         DO 750 I = IY1, IY2
            Y(I) = ((((((Y(I) + A1*X(I1+I)) + A2*X(I2+I)) + A3*X(I3+I))
     +                        + A4*X(I4+I)) + A5*X(I5+I)) + A6*X(I6+I))
     +                        + A7*X(I7+I)
  750    CONTINUE
         GOTO 800
C
  800    CONTINUE
         DO 1000 J = REM1, N, 8
            I1 = XPNT(J+1) - M2
            I2 = XPNT(J+2) - M2
            I3 = XPNT(J+3) - M2
            I4 = XPNT(J+4) - M2
            I5 = XPNT(J+5) - M2
            I6 = XPNT(J+6) - M2
            I7 = XPNT(J+7) - M2
            I8 = XPNT(J+8) - M2
            A1 = -X(I1+IY1)
            A2 = -X(I2+IY1)
            A3 = -X(I3+IY1)
            A4 = -X(I4+IY1)
            A5 = -X(I5+IY1)
            A6 = -X(I6+IY1)
            A7 = -X(I7+IY1)
            A8 = -X(I8+IY1)
            DO 900 I = IY1, IY2
               Y(I) = (((((((Y(I) + A1*X(I1+I)) + A2*X(I2+I))
     +                            + A3*X(I3+I)) + A4*X(I4+I))
     +                            + A5*X(I5+I)) + A6*X(I6+I))
     +                            + A7*X(I7+I)) + A8*X(I8+I)
  900       CONTINUE
 1000    CONTINUE
C
         M1 = M1 - 1
         LY = LY - 1
 2000 CONTINUE
C
      RETURN
      END
