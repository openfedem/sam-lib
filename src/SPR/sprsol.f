C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRSOL (IOP, MSPAR, MTREES, MSIFA, SM, B,
     +                   LDB, NRHS, TOL, IWORK, RWORK, LPU, IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRSOL                  GROUP 9 / PUBLIC
C
C     THIS SUBROUTINE SOLVES A SYSTEM OF EQUATIONS  A*X = B
C     (A  IS SYMMETRIC) THROUGH THREE STEPS:
C        1) FACTORIZATION OF A         -  A = L*D*(LT)
C        2) FORWARD REDUCTION OF B     -  Y = (D-1)*(L-1)*B
C        3) BACKSUBSTITUTION TO GIVE X -  X = (L-T)*Y
C     DEPENDING ON THE VALUE OF THE OPERATION FLAG (IOP), ONE OF THE
C     THREE STEPS, THE FIRST TWO, THE LAST TWO, OR ALL THREE STEPS
C     MAY BE CARRIED OUT IN ONE REFERENCE TO SPRSOL
C
C **********************************************************************
C
      IMPLICIT NONE
C
      INTEGER           IOP, LDB, NRHS, LPU, IERR
      INTEGER           MSPAR(*), MTREES(*), MSIFA(*), IWORK(*)
      DOUBLE PRECISION  SM(*), B(*), TOL(3), RWORK(*)
C
      EXTERNAL          LMVY18, LMMY18
C
      IERR = 0
      IF (IOP .LE. 3) THEN
C
         CALL SPRFC1 (MSPAR, MTREES, MSIFA, SM, TOL,
     +                IWORK, RWORK, LPU, IERR, LMVY18, LMMY18)
C
      END IF
      IF (IOP .GE. 2 .AND. IOP .LE. 5 .AND. IERR .GE. 0) THEN
C
         CALL SPRFS1 (MSPAR, MTREES, MSIFA, SM, B, LDB, NRHS,
     +                IWORK, RWORK, LPU, IERR)
C
      END IF
      IF (IOP .GE. 3 .AND. IOP .NE. 5 .AND. IERR .GE. 0) THEN
C
         CALL SPRBS1 (MSPAR, MTREES, MSIFA, SM, B, LDB, NRHS,
     +                RWORK, LPU, IERR)
C
      END IF
C
      RETURN
      END
