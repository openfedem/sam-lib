C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE DMVV (N,A,INCA,B,INCB,C)
      INTEGER          N,INCA,INCB
      DOUBLE PRECISION A(*),B(*),C(*)
C
C     DMVV : Element-wise vector multiplication.
C
      INTEGER :: I, IA, IB, M
C
      IF (N .LT. 1) THEN
         RETURN
      ELSE IF (N .EQ. 1) THEN
         C(1) = A(1)*B(1)
      ELSE IF (INCA .EQ. 1 .AND. INCB .EQ. 1) THEN
C
C     Code for increment equal to 1 (use loop unrolling)
C
         M = MOD(N,7)
         DO I = 1, M
            C(I) = A(I)*B(I)
         END DO
         DO I = M+1, N, 7
            C(I  ) = A(I  )*B(I  )
            C(I+1) = A(I+1)*B(I+1)
            C(I+2) = A(I+2)*B(I+2)
            C(I+3) = A(I+3)*B(I+3)
            C(I+4) = A(I+4)*B(I+4)
            C(I+5) = A(I+5)*B(I+5)
            C(I+6) = A(I+6)*B(I+6)
         END DO
C
      ELSE IF (INCA .EQ. 1 .AND. INCB .EQ. 0) THEN
C
C     Code for increment equal to 1 (use loop unrolling)
C
         M = MOD(N,7)
         DO I = 1, M
            C(I) = A(I)*B(1)
         END DO
         DO I = M+1, N, 7
            C(I  ) = A(I  )*B(1)
            C(I+1) = A(I+1)*B(1)
            C(I+2) = A(I+2)*B(1)
            C(I+3) = A(I+3)*B(1)
            C(I+4) = A(I+4)*B(1)
            C(I+5) = A(I+5)*B(1)
            C(I+6) = A(I+6)*B(1)
         END DO
C
      ELSE
C
C     Code for increments larger than 1 (no loop unrolling)
C
         IF (INCA .LT. 0) THEN
            IA = 1 + (1-N)*INCA
         ELSE
            IA = 1
         END IF
         IF (INCB .LT. 0) THEN
            IB = 1 + (1-N)*INCB
         ELSE
            IB = 1
         END IF
         DO I = 1, N
            C(I) = A(IA)*B(IB)
            IA = IA + INCA
            IB = IB + INCB
         END DO
C
      END IF
C
      RETURN
      END
