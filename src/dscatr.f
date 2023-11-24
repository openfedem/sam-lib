C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE DSCATR (N,IX,A,B,INCA)
      INTEGER            N,INCA,IX(N)
      DOUBLE PRECISION   A(*),B(*)
C
C     DSCATR : Scatter operator for large arrays.
C
      INTEGER :: I, IA, IB, INCB, M
C
      IF (N .LT. 1) THEN
         RETURN
      ELSE IF (N .EQ. 1) THEN
         B(IX(1)) = A(1)
      ELSE IF (INCA .NE. 1) THEN
C
C     Code for increments larger than 1 (no loop unrolling)
C
         IF (INCA .LT. 0) THEN
            INCB = -1
            IA = 1 + (1-N)*INCA
            IB = N
         ELSE
            INCB = 1
            IA = 1
            IB = 1
         END IF
         DO I = 1, N
            B(IX(IB)) = A(IA)
            IA = IA + INCA
            IB = IB + INCB
         END DO
C
      ELSE
C
C     Code for increment equal to 1 (use loop unrolling)
C
         M = MOD(N,7)
         DO I = 1, M
            B(IX(I)) = A(I)
         END DO
         DO I = M+1, N, 7
            B(IX(I  )) = A(I  )
            B(IX(I+1)) = A(I+1)
            B(IX(I+2)) = A(I+2)
            B(IX(I+3)) = A(I+3)
            B(IX(I+4)) = A(I+4)
            B(IX(I+5)) = A(I+5)
            B(IX(I+6)) = A(I+6)
         END DO
C
      END IF
C
      RETURN
      END
