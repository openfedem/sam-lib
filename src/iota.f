C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE IOTA (N, IX, INCX, I0)
C
C    Description of Parameters
C
C     --Input--
C        N  number of elements to generate
C     INCX  storage spacing between elements of IX
C
C     --Output--
C       IX  integer vector with N elements
C
C     Make a sequence of increasing integers in array IX.
C     For I = 0 to N-1, copy  I0+I to IX(LX+I*INCX),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX.
C     Basically, with INCX=1 it works like std::iota() of the C++ library.
C
      INTEGER :: N, IX(*), INCX, I0, IIX, I, M
C
      IF (N .LE. 0) THEN
         RETURN
      ELSE IF (N .EQ. 1) THEN
         IX(1) = I0
      ELSE IF (INCX .NE. 1) THEN
C
C     Code for unequal or nonpositive increments.
C
         IF (INCX .LT. 0) THEN
            IIX = 1 + (1-N)*INCX
         ELSE
            IIX = 1
         END IF
         DO I = 1, N
            IX(IIX) = I0+I-1
            IIX = IIX + INCX
         END DO
C
      ELSE
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
         M = MOD(N,7)
         DO I = 1, M
            IX(I) = I0+I-1
         END DO
C
         DO I = M, N-1, 7
            IX(I+1) = I0+I
            IX(I+2) = I0+I+1
            IX(I+3) = I0+I+2
            IX(I+4) = I0+I+3
            IX(I+5) = I0+I+4
            IX(I+6) = I0+I+5
            IX(I+7) = I0+I+6
         END DO
C
      END IF
C
      RETURN
      END
