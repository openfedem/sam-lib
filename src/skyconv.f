C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SKYCONV (ASKY,MSKY,ASQR,LDA,N,KSA)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SKYCONV                  GROUP 8 / PUBLIC
C
C     T A S K :  To convert a symmetric "skyline-stored" (ABS(KSA)=1)
C                or diagonal (ABS(KSA)=2) matrix, ASKY, to a full,
C                square N by N matrix, ASQR.
C
C
C     ROUTINES CALLED/REFERENCED :  ABS, MOD, SIGN  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-01-25 / 1.0
C                       03-07-30 / 2.0    K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           LDA,N,KSA,MSKY(N)
      DOUBLE PRECISION  ASKY(*),ASQR(LDA,N)
C                                                ! Local variables
      INTEGER           I,IS,J,JP
      DOUBLE PRECISION  ZERO
C
      PARAMETER       ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IF (N .LT. 1 .OR. N .GT. LDA) RETURN
C
      IF (ABS(KSA) .LT. 10) THEN
         DO 15 J = 1, N
            DO 10 I = 1, N
               ASQR(I,J) = ZERO
  10        CONTINUE
  15     CONTINUE
      ENDIF
      ASQR(1,1) = ASKY(1)*SIGN(1,KSA)
C
      IF (MOD(ABS(KSA),10) .EQ. 2) THEN
C                                                ! Diagonal matrix
         DO 100 J = 2, N
            ASQR(J,J) = ASKY(J)*SIGN(1,KSA)
  100    CONTINUE
C
      ELSE
C                                                ! Skyline matrix
         JP = 2
         DO 300 J = 2, N
            IS = J - MSKY(J) + MSKY(J-1) + 1
            DO 200 I = IS, J
               ASQR(I,J) = ASKY(JP)*SIGN(1,KSA)
               ASQR(J,I) = ASKY(JP)*SIGN(1,KSA)
               JP = JP + 1
  200       CONTINUE
  300    CONTINUE
C
      ENDIF
C
      RETURN
      END
