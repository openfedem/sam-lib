C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRSEQ (MADOF,MINEX,MMNPC,MPMCEQ,MMCEQ,
     &                   NANOD,NDOF,NMMNPC,LPU,MSC,MPAR,
     &                   MASTER,NUMNOD,MEQN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRSEQ                 GROUP 9 / PRIVATE
C
C     TASK :  To check control matrices MSC and MMCEQ, to determine
C             parameters NDOF1, NDOF2, NSPDOF, NSDOF, NPDOF, NDDOF, NEQ
C             and NMMCEQ and store them in MPAR and, finally, to
C             determine control matrix MEQN.
C
C
C     ROUTINES CALLED/REFERENCED :  SPRDOF and SPRER1    (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-01-11 / 1.0
C                       02-09-04 / 1.1   K.M.Okstad
C                       05-10-10 / 1.2   K.M.Okstad
C                       18-05-25 / 2.0   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           NANOD,NDOF,NMMNPC,LPU,IERR
      INTEGER           MADOF(NANOD+1),MINEX(NANOD),MMNPC(NMMNPC),
     &                  MPMCEQ(*),MMCEQ(*),MSC(NDOF),MPAR(50),
     &                  MASTER(NDOF),NUMNOD(NANOD),MEQN(NDOF)
C
C                                                ! local variables
      INTEGER           I,IDOF,INOD,IP,ISC,J,JDOF,LDOF,
     &                  NCEQ,NDDOF,NDF(2),NM,NPDOF,NSPDOF
C ----------------------------------------------------------------------
      IERR = 0
      NCEQ = MPAR(7)
C
C --- Mark all active nodes (i.e., those that are connected to elements)
C
      DO INOD = 1, NANOD
         NUMNOD(INOD) = 0
      END DO
      DO IP = 1, NMMNPC
         INOD = MMNPC(IP)
         IF (INOD .GT. 0 .AND. INOD .LE. NANOD) THEN
            NUMNOD(INOD) = NUMNOD(INOD) + 1
         ENDIF
      END DO
C
C --- Mark all DOFs that are masters in constraint equations
C
      DO IDOF = 1, NDOF
         MASTER(IDOF) = 0
      END DO
      DO I = 1, NCEQ
         DO IP = MPMCEQ(I)+1, MPMCEQ(I+1)-1
            IDOF = MMCEQ(IP)
            IF (IDOF .GT. 0 .AND. IDOF .LE. NDOF) THEN
               MASTER(IDOF) = MASTER(IDOF) + 1
            ENDIF
         END DO
      END DO
C
      NDF(1) = 0
      NDF(2) = 0
      NSPDOF = 0
      NPDOF  = 0
      NDDOF  = 0
C
C --- Check MSC and determine NDOF1, NDOF2 and NSPDOF
C
      DO INOD = 1, NANOD
         DO IDOF = MADOF(INOD), MADOF(INOD+1)-1
            ISC  = MSC(IDOF)
            IF (ISC .EQ. 1 .OR. ISC .EQ. 2) THEN
               IF (NUMNOD(INOD) .GT. 0 .OR. MASTER(IDOF) .GT. 0) THEN
                  NDF(ISC) = NDF(ISC) + 1
               ELSE ! Unused free DOF, deactivate it (specified to zero)
                  MSC(IDOF) = 0
                  MEQN(IDOF) = 0
                  NSPDOF = NSPDOF + 1
               ENDIF
            ELSEIF (ISC .EQ. 0)  THEN
               MEQN(IDOF) = 0
               NSPDOF = NSPDOF + 1
            ELSE
               IERR = IERR - 1
               LDOF = IDOF - MADOF(INOD) + 1
               IP   = MINEX(INOD)
               CALL SPRER1 (1,'SPRSEQ',ISC,LDOF,IP,LPU,IERR)
            ENDIF
         END DO
      END DO
C
      MPAR( 4) = NDF(1)
      MPAR( 5) = NDF(2)
      MPAR( 6) = NSPDOF
      MPAR( 8) = NSPDOF - NCEQ
      MPAR(11) = NDF(1) + NDF(2)
C
C --- Check MMCEQ, determine NPDOF and NDDOF, and
C     set negative values in MEQN for prescribed and dependent DOFs
C
      DO I = 1, NCEQ
         IP = MPMCEQ(I)
         IDOF = MMCEQ(IP)
         IF (IDOF .LT. 1 .OR. IDOF .GT. NDOF) THEN
            IERR = IERR - 1
            CALL SPRER1 (3,'SPRSEQ',IDOF,I,I,LPU,IERR)
         ELSEIF (MSC(IDOF) .NE. 0) THEN
            IERR = IERR - 1
            CALL SPRDOF (MADOF,MINEX,NANOD,IDOF,LDOF,INOD)
            CALL SPRER1 (4,'SPRSEQ',LDOF,INOD,I,LPU,IERR)
         ELSE
            MEQN(IDOF) = -I
            NM = MPMCEQ(I+1) - IP - 1
            IF (NM .EQ. 0) THEN                  ! prescribed DOF
               NPDOF = NPDOF + 1
            ELSEIF (NM .GT. 0) THEN              ! slave DOF
               NDDOF = NDDOF + 1
               DO J = 1, NM
                  IP = IP + 1
                  JDOF = MMCEQ(IP)
                  IF (JDOF .GT. 0 .AND. JDOF .LE. NDOF) THEN
                     IF (MSC(JDOF) .LT. 1) THEN
                        IERR = IERR - 1
                        CALL SPRDOF (MADOF,MINEX,NANOD,JDOF,LDOF,INOD)
                        CALL SPRER1 (6,'SPRSEQ',LDOF,INOD,I,LPU,IERR)
                     ENDIF
                  ENDIF
               END DO
            ELSE
               IERR = IERR - 1
               CALL SPRER1 (5,'SPRSEQ',I,I,I,LPU,IERR)
            ENDIF
         ENDIF
      END DO
C
      MPAR( 9) = NPDOF
      MPAR(10) = NDDOF
      MPAR(16) = MPMCEQ(NCEQ+1) - 1
C
C --- Complete MEQN
C
      I = 1
      J = NDF(1) + 1
      DO IDOF = 1, NDOF
         IF (MSC(IDOF) .EQ. 1) THEN
            MEQN(IDOF) = I
            I = I + 1
         ELSEIF (MSC(IDOF) .EQ. 2) THEN
            MEQN(IDOF) = J
            J = J + 1
         ENDIF
      END DO
C
      RETURN
      END
