C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE DGATHR (N,IX,A,B,INCB)
      INTEGER            N,INCB,IX(N)
      DOUBLE PRECISION   A(*),B(*)
C
C     DGATHR : Gather operator for large arrays.
C
      INTEGER :: I, IA, IB, INCA, M
C
      IF (N .LT. 1) THEN
         RETURN
      ELSE IF (N .EQ. 1) THEN
         B(1) = A(IX(1))
      ELSE IF (INCB .NE. 1) THEN
C
C     Code for increments larger than 1 (no loop unrolling)
C
         IF (INCB .LT. 0) THEN
            INCA = -1
            IA = N
            IB = 1 + (1-N)*INCB
         ELSE
            INCA = 1
            IA = 1
            IB = 1
         END IF
         DO I = 1, N
            B(IB) = A(IX(IA))
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
            B(I) = A(IX(I))
         END DO
         DO I = M+1, N, 7
            B(I  ) = A(IX(I  ))
            B(I+1) = A(IX(I+1))
            B(I+2) = A(IX(I+2))
            B(I+3) = A(IX(I+3))
            B(I+4) = A(IX(I+4))
            B(I+5) = A(IX(I+5))
            B(I+6) = A(IX(I+6))
         END DO
C
      END IF
C
      RETURN
      END
