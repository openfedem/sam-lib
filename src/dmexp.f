C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE DMEXP (M,N,MM,NN,A)
      INTEGER           M,N,MM,NN
      DOUBLE PRECISION  A(*)
C
C     DMEXP  : Expand the M by N matrix A into a MM by NN matrix.
C
      INTEGER            I,J,JM,JN
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D0 )
C
      IF (M  .LE. 0 .OR.  N  .LE. 0) RETURN
      IF (MM .LT. M .OR.  NN .LT. N) RETURN
      IF (MM .EQ. M .AND. NN .EQ. N) RETURN
C
      DO 10 I = M*N+1,MM*NN
         A(I) = ZERO
   10 CONTINUE
C
      IF (MM .EQ. M) RETURN
C
      DO 200 J = N,2,-1
         JN =  M*J
         JM = MM*J - MM + M
         DO 100 I = 1,M
            A(JM) = A(JN)
            A(JN) = ZERO
            JN = JN - 1
            JM = JM - 1
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
