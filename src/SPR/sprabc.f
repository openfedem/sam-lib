C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRABC (MSPAR,MSIFA,MEQNA,MEQNB,NDOF,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRABC                    GROUP 9 / PUBLIC
C
C     T A S K :  To determine permutation arrays between two different
C                equation orderings, as given by MEQNA and MEQNB. It is
C                used when, e.g., we have an eigenvalue problem where
C                both the stiffness and mass matrix are on sparse form,
C                but the mass matrix has a different sparsity pattern
C                than the stiffness matrix (when using lumped mass).
C                The permutation arrays are placed first in MSIFA.
C
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Knut Morten Okstad
C     DATE/VERSION  :   03-02-27 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           NDOF,LPU,IERR
      INTEGER           MSPAR(*),MSIFA(*),MEQNA(NDOF),MEQNB(NDOF)
      INTEGER           I,IEQNA,IEQNB,NEQ,PERE2I,PERI2E
C ----------------------------------------------------------------------
C
      IERR   = 0
      NEQ    = MSPAR(8)
      PERI2E = 0
      PERE2I = PERI2E + NEQ
C
      DO 100 I = 1, NEQ
         MSIFA(PERI2E+I) = 0
         MSIFA(PERE2I+I) = 0
  100 CONTINUE
C
C     MEQNA  - gives the internal equation ordering for each DOF
C     MEQNB  - gives the external equation ordering for each DOF
C     PERE2I - permutation from external to internal ordering
C     PERI2E - permutation from internal to external ordering
C
      DO 200 I = 1, NDOF
         IEQNA = MEQNA(I)
         IEQNB = MEQNB(I)
         IF (IEQNA .GT. 0 .AND. IEQNB .GT. 0) THEN
            IF (IEQNA .GT. NEQ .OR. IEQNB .GT. NEQ) GOTO 900
            MSIFA(PERI2E+IEQNA) = IEQNB
            MSIFA(PERE2I+IEQNB) = IEQNA
         ELSE IF (IEQNA .GT. 0 .OR. IEQNB .GT. 0) THEN
            GOTO 900
         END IF
  200 CONTINUE
C
      DO 300 I = 1, NEQ
         IF (MSIFA(PERI2E+I) .EQ. 0) GOTO 910
         IF (MSIFA(PERE2I+I) .EQ. 0) GOTO 910
  300 CONTINUE
C
      MSPAR(58) = 2
C
      RETURN
C
  900 CONTINUE
      IERR = -1
      WRITE(LPU,6000)
      WRITE(LPU,6010) I
      RETURN
C
  910 CONTINUE
      IERR = -2
      WRITE(LPU,6000)
      WRITE(LPU,6020) I
      RETURN
C
 6000 FORMAT(/' *** ERROR RETURN FROM SPRABC'
     $       /'     INCONSISTENT EQUATION ORDERING DATA.' )
 6010 FORMAT( '     ENCOUNTERED WHILE PROCESSING DOF NO.',I8)
 6020 FORMAT( '     ENCOUNTERED WHILE CHECKING EQUATION NO.',I8)
C
      END
