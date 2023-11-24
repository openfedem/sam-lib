C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRPRP (MADOF,MINEX,
     $                   MPMNPC,MMNPC,
     $                   MPMCEQ,MMCEQ,
     $                   MSC,NSPAR,LPU,
     $                   MPAR,MSPAR,
     $                   MEQN,MWORK,
     $                   IERR        )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRPRP                  GROUP 9 / PUBLIC
C
C     T A S K :  To prepare (check and determine) control information
C                for the "sparse" assembly process ("pre-assembly"),
C                and to determine a safe estimate for the size of the
C                control matrix MSICA, that is, estimate NMSICA.
C
C
C     ROUTINES CALLED/REFERENCED: SPRSEQ, SPRMXD, SPRSEN, SPRER1 (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-01-11 / 1.0
*                       98-08-24 / 1.1   A.C.Damhaug
*                                        MSPAR(50) -> MSPAR(*), and
*                                        added superelement option.
C                       03-01-15 / 1.2   K.M.Okstad
C                                        Moved initialization of
C                                        MSPAR(35:36) to SPRSEN.
C                       03-03-11 / 1.3   K.M.Okstad
C                                        Added call to SPRMXD.
C                       04-09-21 / 1.4   K.M.Okstad
C                                        MSPAR(*) -> MSPAR(NSPAR).
C                                        NSPAR is new input argument.
C                       18-05-25 / 1.5   K.M.Okstad
C                                        Added MMNPC,NANOD,NDOF,NMMNPC,
C                                        MWORK(1) and MWORK(IP1) as
C                                        arguments in call to SPRSEQ.
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,LPU,NSPAR
      INTEGER           MADOF(*),MEQN(*),MINEX(*),MMCEQ(*),MMNPC(*),
     $                  MPAR(50),MPMCEQ(*),MPMNPC(*),MSC(*),
     $                  MSPAR(NSPAR),MWORK(*)
C                                                ! Local variables
C
      INTEGER           IP1,IP2,NANOD,NDOF,NEL,NMMNPC
C ----------------------------------------------------------------------
      IERR  = 0
      NANOD = MPAR(1)
      NEL   = MPAR(2)
      NDOF  = MADOF(NANOD+1) - 1
      IF (MPAR(3) .NE. NDOF) THEN
         IERR = -1
         CALL SPRER1 (7,'SPRPRP',NDOF,NDOF,NDOF,LPU,IERR)
         RETURN
      ENDIF
      NMMNPC   = MPMNPC(NEL+1) - 1
      MPAR(15) = NMMNPC
C
C --- Calculate MEQN
      IP1 = 1   + NDOF
      IP2 = IP1 + NANOD
      CALL SPRSEQ (MADOF,MINEX,MMNPC,MPMCEQ,MMCEQ,
     $             NANOD,NDOF,NMMNPC,LPU,MSC,MPAR,
     $             MWORK(1),MWORK(IP1),MEQN,IERR)
      IF (IERR .LT. 0) RETURN
C
C --- Calculate maximum element dofs and store in MPAR(20)
      CALL SPRMXD (MADOF,MPMNPC,MMNPC,MPAR)
C
C --- Estimate NMSICA
      CALL SPRSEN (MPAR,MPMNPC,MMNPC,MADOF,MSC,
     $             MPMCEQ,MMCEQ,MEQN,NEL,NANOD,NDOF,NMMNPC,NSPAR,
     $             MSPAR,MWORK(1),MWORK(IP1),MWORK(IP2))

      RETURN
      END
