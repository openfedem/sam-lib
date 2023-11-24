C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE RGD2D (MINEX,MADOF,MEKN,
     &                  MPMNPC,MMNPC,MPMCEX,MMCEX,
     &                  TXC,TYZC,CH,TTCCX,
     &                  MPAR,MSC,MPMCEQ,MMCEQ,
     &                  TTCC,NCEX,KRIGID,MAXSIZ,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGD2D                 GROUP 9 / PUBLIC
C
C     TASK :  To generate new versions of the constraint information
C             arrays, MPMCEQ,MMCEQ and TTCC, that account for both
C             implicitly defined constraints via rigid elements and the
C             NCEX explicitly defined constraints recorded in arrays
C             MPMCEX, MMCEX and TTCCX.  Rigid elements have an entry in
C             MEKN equal to KRIGID.
C             Table TYZC contains the y-coordinates of the nodal points
C             if CH = 'Y' or 'y', or the z-coordinates if CH = 'Z' or
C             'z'.
C             Status codes in MSC are updated and the final number of
C             constraint equations (NCEQ) is stored in MPAR(7).
C             Arrays MMCEQ and TTCC contain MAXSIZ elements on input;
C             the actual number of elements containing relevant infor-
C             formation (NMMCEQ) is determined by the routine and
C             returned in MPAR(16).
C             On input to the routine the number of (active) nodal
C             points (NANOD) and the total number of elements (NEL)
C             are assumed to be stored in MPAR(1) and MPAR(2), re-
C             spectively, and the total number of dofs (NDOF) is assumed
C             to be stored in MPAR(3).
C             For MPMCEQ a storage size of NDOF (total number of dofs),
C             which is the theoretical maximum, is recommended.  It
C             should be noted that no check is made on "overflow" of
C             MPMCEQ.
C
C
C     ROUTINES CALLED/REFERENCED :   RGDC0,  RGDCE2, RGDXPL   (SAM-9)
C                                    RGDSWP, RGDERR           (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-06-29 / 1.0
C                       92-04-28 / 2.0    K.Bell
C                       99-10-14 / 2.1    K.Bell
C                       02-07-09 / 2.2    K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,KRIGID,LPU,MAXSIZ,NCEX
      INTEGER           MADOF(*),MEKN(*),MINEX(*),MMCEQ(*),MMCEX(*),
     &                  MMNPC(*),MPAR(*),MPMCEQ(*),MPMCEX(*),
     &                  MPMNPC(*),MSC(*)
      DOUBLE PRECISION  TTCC(*),TTCCX(*),TXC(*),TYZC(*)
      CHARACTER*1       CH
C
C                                                ! local variables
C
      INTEGER           I,ICEQ,IDOFM1,IDOFM2,IDOFS,IEL,
     &                  INODM,INODS,IPCUR,ISC,JP,NM,NNDOF
      DOUBLE PRECISION  C,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
      IF (CH.EQ.'Y' .OR. CH.EQ.'y') THEN
         CONTINUE
      ELSEIF (CH.EQ.'Z' .OR. CH.EQ.'z') THEN
         CONTINUE
      ELSE
         CALL RGDERR (8,I,I,I,LPU,IERR)
         GO TO 1000
      ENDIF
C
C --- set status codes (in MSC) for explicitly constrained dofs to -1
C     for prescribed dofs, and to -2 for slave dofs
C
      IF (NCEX .GT. 0) THEN
         DO 10 I=1,NCEX
            JP    = MPMCEX(I)
            IDOFS = MMCEX(JP)
            IF (IDOFS.LT.1 .OR. IDOFS.GT.MPAR(3)) THEN
               CALL RGDERR (1,I,IDOFS,IDOFS,LPU,IERR)
               GO TO 1000
            ENDIF
            NM    = MPMCEX(I+1) - JP - 1
            IF (NM .EQ. 0) THEN
               MSC(IDOFS) = -1
            ELSEIF (NM .GT. 0) THEN
               MSC(IDOFS) = -2
            ELSE
               CALL RGDERR (7,I,NM,NM,LPU,IERR)
               GO TO 1000
            ENDIF
   10    CONTINUE
      ENDIF
C
C --- For each rigid element check:
C     - status of each master dof
C     - status of each slave dof
C     and increase legal positive status code by 10
C
      DO 100 IEL=1,MPAR(2)
         IF (MEKN(IEL) .EQ. KRIGID) THEN
            JP    = MPMNPC(IEL)
            INODM = MMNPC(JP)
            IDOFM1 = MADOF(INODM)
            DO 30 I=1,3
               IF (MSC(IDOFM1).GT.0 .AND. MSC(IDOFM1).LT.11) THEN
                  MSC(IDOFM1) = MSC(IDOFM1) + 10
               ELSEIF (MSC(IDOFM1) .EQ. (-2)) THEN
                  CALL RGDERR (5,I,MINEX(INODM),I,LPU,IERR)
                  GO TO 1000
               ENDIF
               IDOFM1 = IDOFM1+1
   30       CONTINUE
            DO 50 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
               INODS = MMNPC(JP)
               DO 40 IDOFS=MADOF(INODS),MADOF(INODS+1)-1
                  ISC = MSC(IDOFS)
                  IF (ISC .GT. 0) THEN
                     MSC(IDOFS) = MSC(IDOFS) + 10
                  ELSEIF (ISC.EQ.(-1)) THEN
                     CALL RGDERR (10,1,MINEX(INODS),I,LPU,IERR)
                     GO TO 1000
                  ENDIF
   40          CONTINUE
   50       CONTINUE
         ENDIF
  100 CONTINUE
C
C --- initialize current slave (ICEQ) and pointer in MMCEQ/TTCC (IPCUR)
C
      ICEQ  = 1
      IPCUR = 1
C
C ======================================================================
C     For each dof of each slave node of each rigid element:
C     establish its constraint equation in MPMCEQ/MMCEQ/TTCC
C ======================================================================
C
      DO 200 IEL=1,MPAR(2)
         IF (MEKN(IEL) .EQ. KRIGID) THEN
            JP    = MPMNPC(IEL)
            INODM = MMNPC(JP)
            DO 150 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
               INODS  = MMNPC(JP)
               NNDOF  = MADOF(INODS+1) - MADOF(INODS)
               IDOFS  = MADOF(INODS)
               ISC    = MSC(IDOFS)
               IDOFM1 = MADOF(INODM)
               IDOFM2 = IDOFM1+2
               IF (CH.EQ.'Y' .OR. CH.EQ.'y') THEN
                  C = TYZC(INODM) - TYZC(INODS)
               ELSE
                  C = TYZC(INODS) - TYZC(INODM)
               ENDIF
               CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                      MPMCEQ,MMCEQ,TTCC,
     &                      C,IDOFS,IDOFM1,IDOFM2,
     &                      NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
               IF (IERR .LT. 0) THEN
                  CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                  GO TO 1000
               ENDIF
               IF (ISC .EQ. (-2)) THEN
                  CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                         NCEX,INODS,IDOFS,LPU,IERR)
                  IF (IERR .LT. 0)  GO TO 1000
               ENDIF
C
               IDOFS  = IDOFS+1
               ISC    = MSC(IDOFS)
               IDOFM1 = IDOFM1+1
               IF (CH.EQ.'Y' .OR. CH.EQ.'y') THEN
                  C = TXC(INODS) - TXC(INODM)
               ELSE
                  C = TXC(INODM) - TXC(INODS)
               ENDIF
               CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                      MPMCEQ,MMCEQ,TTCC,
     &                      C,IDOFS,IDOFM1,IDOFM2,
     &                      NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
               IF (IERR .LT. 0) THEN
                  CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                  GO TO 1000
               ENDIF
               IF (ISC .EQ. (-2)) THEN
                  CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                         NCEX,INODS,IDOFS,LPU,IERR)
                  IF (IERR .LT. 0)  GO TO 1000
               ENDIF
C
               IF (NNDOF .EQ. 2) THEN
                  CONTINUE
               ELSEIF (NNDOF .EQ. 3) THEN
                  IDOFS  = IDOFS+1
                  ISC    = MSC(IDOFS)
                  IDOFM1 = IDOFM2
                  IDOFM2 = 0
                  C      = ZERO
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
               ELSE
                  CALL RGDERR (3,MINEX(INODS),NNDOF,NNDOF,LPU,IERR)
                  GO TO 1000
               ENDIF
  150       CONTINUE
         ENDIF
  200 CONTINUE
C
C ======================================================================
C     Append explicitly defined constraint equations to MPMCEQ/MMCEQ/
C     TTCC, and expand/contract on any specified master dofs such that
C     all master dofs included in the final version of the constraint
C     array MMCEQ are free dofs
C ======================================================================
C
      IF (NCEX .GT. 0) THEN
         CALL RGDXPL (MPAR,MINEX,MADOF,
     &                MPMCEX,MMCEX,TTCCX,
     &                MSC,MPMCEQ,MMCEQ,TTCC,
     &                NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
         IF (IERR .LT. 0)  GO TO 1000
      ENDIF
C
C --- restore MSC, set NCEQ and NMMCEQ and complete MPMCEQ
C
      DO 300 I=1,MPAR(3)
         IF (MSC(I) .LT. 0) THEN
            MSC(I) = 0
         ELSEIF (MSC(I) .GT. 10) THEN
            MSC(I) = MSC(I) - 10
         ENDIF
  300 CONTINUE
      MPAR(7)  = ICEQ-1
      MPAR(16) = IPCUR-1
      IF (ICEQ .GT. 1) THEN
         MPMCEQ(ICEQ) = IPCUR
      ENDIF
C
 1000 RETURN
      END
      SUBROUTINE RGD3D (MINEX,MADOF,MEKN,
     &                  MPMNPC,MMNPC,MPMCEX,MMCEX,
     &                  TXC,TYC,TZC,TTCCX,
     &                  MPAR,MSC,MPMCEQ,MMCEQ,
     &                  TTCC,NCEX,KRIGID,KRXY,
     &                  KRXZ,KRYZ,MAXSIZ,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGD3D                 GROUP 9 / PUBLIC
C
C     TASK :  To generate new versions of the constraint information
C             arrays, MPMCEQ,MMCEQ and TTCC, that account for both
C             implicitly defined constraints via rigid elements and the
C             NCEX explicitly defined constraints recorded in arrays
C             MPMCEX, MMCEX and TTCCX.
C             Four types of rigid elements may be accounted for, each
C             of which is recognized by a special "code value" in the
C             appropriate entry of MEKN:
C             KRIGID - a completely rigid 3-dim element with 6 dofs
C             KRXY   - an element that is completely rigid in the xy-
C                      plane, but has no stiffness out of this plane;
C                      the element has 3 significant dofs
C             KRXZ   - an element that is completely rigid in the xz-
C                      plane, but has no stiffness out of this plane;
C                      the element has 3 significant dofs
C             KRYZ   - an element that is completely rigid in the yz-
C                      plane, but has no stiffness out of this plane;
C                      the element has 3 significant dofs
C             A code value of zero (=0) indicates that the particular
C             type of element is not present in the model.
C             Status codes in MSC are updated and the final number of
C             constraint equations (NCEQ) is stored in MPAR(7).
C             Arrays MMCEQ and TTCC contain MAXSIZ (insignificant)
C             elements on input;  the actual number of elements con-
C             taining relevant information (NMMCEQ) is determined by
C             the routine and returned in MPAR(16).
C             On input to the routine the total number of elements (NEL)
C             is assumed to be stored in MPAR(2), and the total number
C             of dofs (NDOF) is assumed to be stored in MPAR(3).
C             For MPMCEQ a storage size of NDOF (total number of dofs),
C             which is the theoretical maximum, is recommended.  It
C             should be noted that no check is made on "overflow" of
C             MPMCEQ.
C
C
C     ROUTINES CALLED/REFERENCED :   RGDC0,  RGDCE2, RGDCE3   (SAM-9)
C                                    RGDERR, RGDXPL, RGDSWP   (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-07-06 / 1.0
C                       92-05-09 / 2.0   K.Bell
C                       99-10-14 / 2.1   K.Bell
C                       00-11-19 / 2.2   B.Haugen
C                       02-07-09 / 2.3   K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,KRIGID,KRXY,KRXZ,KRYZ,LPU,MAXSIZ,NCEX
      INTEGER           MADOF(*),MEKN(*),MINEX(*),MMCEQ(*),MMCEX(*),
     &                  MMNPC(*),MPAR(*),MPMCEQ(*),MPMCEX(*),
     &                  MPMNPC(*),MSC(*)
      DOUBLE PRECISION  TTCC(*),TTCCX(*),TXC(*),TYC(*),TZC(*)
C
C                                                ! local variables
C
      INTEGER           I,ICEQ,IDOFM1,IDOFM2,IDOFM3,IDOFS,IEL,
     &                  INODM,INODS,IPCUR,ISC,JP,KIND,NM,NNDOF
      DOUBLE PRECISION  C1,C2,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
C
C --- set status codes (in MSC) for explicitly constrained dofs to -1
C     for prescribed dofs, and to -2 for slave dofs
C
      IF (NCEX .GT. 0) THEN
         DO 10 I=1,NCEX
            JP    = MPMCEX(I)
            IDOFS = MMCEX(JP)
            IF (IDOFS.LT.1 .OR. IDOFS.GT.MPAR(3)) THEN
               CALL RGDERR (1,I,IDOFS,IDOFS,LPU,IERR)
               GO TO 1000
            ENDIF
            NM    = MPMCEX(I+1) - JP - 1
            IF (NM .EQ. 0) THEN
               MSC(IDOFS) = -1
            ELSEIF (NM .GT. 0) THEN
               MSC(IDOFS) = -2
            ELSE
               CALL RGDERR (7,I,NM,NM,LPU,IERR)
               GO TO 1000
            ENDIF
   10    CONTINUE
      ENDIF
C
C --- For each rigid element check:
C     - status of each master dof
C     - status of each slave dof
C     and increase legal positive status code by 10
C
      DO 100 IEL=1,MPAR(2)
         KIND = MEKN(IEL)
         IF (KIND.EQ.KRIGID .OR. KIND.EQ.KRXY .OR.
     &       KIND.EQ.KRYZ   .OR. KIND.EQ.KRXZ )     THEN
            JP    = MPMNPC(IEL)
            INODM = MMNPC(JP)
            NNDOF = MADOF(INODM+1) - MADOF(INODM)
            IF (NNDOF .NE. 6) THEN
               CALL RGDERR (13,NNDOF,MINEX(INODM),IEL,LPU,IERR)
               GO TO 1000
            ENDIF
            IDOFM1 = MADOF(INODM)
            DO 30 I=1,6
               IF (MSC(IDOFM1).GT.0 .AND. MSC(IDOFM1).LT.11) THEN
                  MSC(IDOFM1) = MSC(IDOFM1) + 10
               ELSEIF (MSC(IDOFM1) .EQ. (-2)) THEN
                  CALL RGDERR (5,I,MINEX(INODM),I,LPU,IERR)
                  GO TO 1000
               ENDIF
               IDOFM1 = IDOFM1+1
   30       CONTINUE
            DO 50 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
               INODS = MMNPC(JP)
               NNDOF = MADOF(INODS+1) - MADOF(INODS)
               IF (NNDOF.EQ.3 .OR. NNDOF.EQ.6) THEN
                  CONTINUE
               ELSE
                  CALL RGDERR (14,NNDOF,MINEX(INODS),IEL,LPU,IERR)
                  GO TO 1000
               ENDIF
               DO 40 IDOFS=MADOF(INODS),MADOF(INODS+1)-1
                  ISC = MSC(IDOFS)
                  IF (ISC .GT. 0) THEN
                     MSC(IDOFS) = MSC(IDOFS) + 10
                  ELSEIF (ISC.EQ.(-1)) THEN
                     CALL RGDERR (10,1,MINEX(INODS),I,LPU,IERR)
                     GO TO 1000
                  ENDIF
   40          CONTINUE
   50       CONTINUE
         ENDIF
  100 CONTINUE
C
C --- initialize current slave (ICEQ) and pointer in MMCEQ/TTCC (IPCUR)
C
      ICEQ  = 1
      IPCUR = 1
C
C ======================================================================
C     For each dof of each slave node of each rigid element:
C     establish its constraint equation in MPMCEQ/MMCEQ/TTCC
C ======================================================================
C
      IF (KRIGID .NE. 0) THEN
C ------------------------------------------------ completely rigid 3-D
C                                                  element
         DO 200 IEL=1,MPAR(2)
            IF (MEKN(IEL) .EQ. KRIGID) THEN
               JP    = MPMNPC(IEL)
               INODM = MMNPC(JP)
               DO 150 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
                  INODS  = MMNPC(JP)
                  NNDOF  = MADOF(INODS+1) - MADOF(INODS)
                  IDOFS  = MADOF(INODS)
                  if (nndof .eq. 6) then
                     isc = msc(idofs+3) + msc(idofs+4) + msc(idofs+5)
                     if (isc .eq. 0) nndof = 3
                  end if
                  ISC    = MSC(IDOFS)
                  IDOFM1 = MADOF(INODM)
                  IDOFM2 = IDOFM1+4
                  IDOFM3 = IDOFM1+5
                  C1     = TZC(INODS) - TZC(INODM)
                  C2     = TYC(INODM) - TYC(INODS)
                  CALL RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IDOFS  = IDOFS+1
                  ISC    = MSC(IDOFS)
                  IDOFM1 = IDOFM1+1
                  IDOFM2 = IDOFM1+2
                  IDOFM3 = IDOFM1+4
                  C1     = TZC(INODM) - TZC(INODS)
                  C2     = TXC(INODS) - TXC(INODM)
                  CALL RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
               IF (ISC .EQ. (-2)) THEN
                  CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                         NCEX,INODS,IDOFS,LPU,IERR)
                  IF (IERR .LT. 0)  GO TO 1000
               ENDIF
C
                  IDOFS  = IDOFS+1
                  ISC    = MSC(IDOFS)
                  IDOFM1 = IDOFM1+1
                  IDOFM2 = IDOFM1+1
                  IDOFM3 = IDOFM1+2
                  C1     = TYC(INODS) - TYC(INODM)
                  C2     = TXC(INODM) - TXC(INODS)
                  CALL RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IF (NNDOF .EQ. 3) THEN
                     CONTINUE
C
                  ELSEIF (NNDOF .EQ. 6) THEN
                     IDOFS  = IDOFS+1
                     ISC    = MSC(IDOFS)
                     IDOFM1 = IDOFM1+1
                     IDOFM2 = 0
                     IDOFM3 = 0
                     C1     = ZERO
                     C2     = ZERO
                     CALL RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                            MPMCEQ,MMCEQ,TTCC,
     &                            C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                            NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),
     &                               I,LPU,I)
                        GO TO 1000
                     ENDIF
                     IF (ISC .EQ. (-2)) THEN
                        CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                               NCEX,INODS,IDOFS,LPU,IERR)
                        IF (IERR .LT. 0)  GO TO 1000
                     ENDIF
C
                     IDOFS  = IDOFS+1
                     ISC    = MSC(IDOFS)
                     IDOFM1 = IDOFM1+1
                     IDOFM2 = 0
                     IDOFM3 = 0
                     C1     = ZERO
                     C2     = ZERO
                     CALL RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                            MPMCEQ,MMCEQ,TTCC,
     &                            C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                            NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),
     &                               I,LPU,I)
                        GO TO 1000
                     ENDIF
                     IF (ISC .EQ. (-2)) THEN
                        CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                               NCEX,INODS,IDOFS,LPU,IERR)
                        IF (IERR .LT. 0)  GO TO 1000
                     ENDIF
C
                     IDOFS  = IDOFS+1
                     ISC    = MSC(IDOFS)
                     IDOFM1 = IDOFM1+1
                     IDOFM2 = 0
                     IDOFM3 = 0
                     C1     = ZERO
                     C2     = ZERO
                     CALL RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                            MPMCEQ,MMCEQ,TTCC,
     &                            C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                            NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),
     &                               I,LPU,I)
                        GO TO 1000
                     ENDIF
                     IF (ISC .EQ. (-2)) THEN
                        CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                               NCEX,INODS,IDOFS,LPU,IERR)
                        IF (IERR .LT. 0)  GO TO 1000
                     ENDIF
                  ELSE
                     CALL RGDERR (3,MINEX(INODS),NNDOF,NNDOF,LPU,IERR)
                     GO TO 1000
                  ENDIF
  150          CONTINUE
            ENDIF
  200    CONTINUE
      ENDIF
C
      IF (KRXY .NE. 0) THEN
C ----------------------------------------------- element is rigid only
C                                                 in xy-plane
         DO 300 IEL=1,MPAR(2)
            IF (MEKN(IEL) .EQ. KRXY) THEN
               JP    = MPMNPC(IEL)
               INODM = MMNPC(JP)
               DO 250 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
                  INODS  = MMNPC(JP)
                  NNDOF  = MADOF(INODS+1) - MADOF(INODS)
                  IDOFS  = MADOF(INODS)
                  ISC    = MSC(IDOFS)
                  IDOFM1 = MADOF(INODM)
                  IDOFM2 = IDOFM1+5
                  C1     = TYC(INODM) - TYC(INODS)
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IDOFS  = IDOFS+1
                  ISC    = MSC(IDOFS)
                  IDOFM1 = IDOFM1+1
                  C1     = TXC(INODS) - TXC(INODM)
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IF (NNDOF .EQ. 3) THEN
                     CONTINUE
                  ELSEIF (NNDOF .EQ. 6) THEN
                     IDOFS  = IDOFS+4
                     ISC    = MSC(IDOFS)
                     IDOFM1 = IDOFM2
                     IDOFM2 = 0
                     C1     = ZERO
                     CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                            MPMCEQ,MMCEQ,TTCC,
     &                            C1,IDOFS,IDOFM1,IDOFM2,
     &                            NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),
     &                               I,LPU,I)
                        GO TO 1000
                     ENDIF
                     IF (ISC .EQ. (-2)) THEN
                        CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                               NCEX,INODS,IDOFS,LPU,IERR)
                        IF (IERR .LT. 0)  GO TO 1000
                     ENDIF
                  ELSE
                     CALL RGDERR (3,MINEX(INODS),NNDOF,NNDOF,LPU,IERR)
                     GO TO 1000
                  ENDIF
  250          CONTINUE
            ENDIF
  300    CONTINUE
      ENDIF
C
      IF (KRXZ .NE. 0) THEN
C ----------------------------------------------- element is rigid only
C                                                 in xz-plane
         DO 400 IEL=1,MPAR(2)
            IF (MEKN(IEL) .EQ. KRXZ) THEN
               JP    = MPMNPC(IEL)
               INODM = MMNPC(JP)
               DO 350 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
                  INODS  = MMNPC(JP)
                  NNDOF  = MADOF(INODS+1) - MADOF(INODS)
                  IDOFS  = MADOF(INODS)
                  ISC    = MSC(IDOFS)
                  IDOFM1 = MADOF(INODM)
                  IDOFM2 = IDOFM1+4
                  C1     = TZC(INODS) - TZC(INODM)
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IDOFS  = IDOFS+2
                  ISC    = MSC(IDOFS)
                  IDOFM1 = IDOFM1+2
                  C1     = TXC(INODM) - TXC(INODS)
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IF (NNDOF .EQ. 3) THEN
                     CONTINUE
                  ELSEIF (NNDOF .EQ. 6) THEN
                     IDOFS  = IDOFS+2
                     ISC    = MSC(IDOFS)
                     IDOFM1 = IDOFM2
                     IDOFM2 = 0
                     C1     = ZERO
                     CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                            MPMCEQ,MMCEQ,TTCC,
     &                            C1,IDOFS,IDOFM1,IDOFM2,
     &                            NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),
     &                               I,LPU,I)
                        GO TO 1000
                     ENDIF
                     IF (ISC .EQ. (-2)) THEN
                        CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                               NCEX,INODS,IDOFS,LPU,IERR)
                        IF (IERR .LT. 0)  GO TO 1000
                     ENDIF
                  ELSE
                     CALL RGDERR (3,MINEX(INODS),NNDOF,NNDOF,LPU,IERR)
                     GO TO 1000
                  ENDIF
  350          CONTINUE
            ENDIF
  400    CONTINUE
      ENDIF
C
      IF (KRYZ .NE. 0) THEN
C ----------------------------------------------- element is rigid only
C                                                 in yz-plane
         DO 500 IEL=1,MPAR(2)
            IF (MEKN(IEL) .EQ. KRYZ) THEN
               JP    = MPMNPC(IEL)
               INODM = MMNPC(JP)
               DO 450 JP=MPMNPC(IEL)+1,MPMNPC(IEL+1)-1
                  INODS  = MMNPC(JP)
                  NNDOF  = MADOF(INODS+1) - MADOF(INODS)
                  IDOFS  = MADOF(INODS)+1
                  ISC    = MSC(IDOFS)
                  IDOFM1 = MADOF(INODM)+1
                  IDOFM2 = IDOFM1+2
                  C1     = TZC(INODM) - TZC(INODS)
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IDOFS  = IDOFS+1
                  ISC    = MSC(IDOFS)
                  IDOFM1 = IDOFM1+1
                  C1     = TYC(INODS) - TYC(INODM)
                  CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                         MPMCEQ,MMCEQ,TTCC,
     &                         C1,IDOFS,IDOFM1,IDOFM2,
     &                         NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),I,LPU,I)
                     GO TO 1000
                  ENDIF
                  IF (ISC .EQ. (-2)) THEN
                     CALL RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                            NCEX,INODS,IDOFS,LPU,IERR)
                     IF (IERR .LT. 0)  GO TO 1000
                  ENDIF
C
                  IF (NNDOF .EQ. 3) THEN
                     CONTINUE
                  ELSEIF (NNDOF .EQ. 6) THEN
                     IDOFS  = IDOFS+1
                     ISC    = MSC(IDOFS)
                     IDOFM1 = IDOFM2
                     IDOFM2 = 0
                     C1     = ZERO
                     CALL RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                            MPMCEQ,MMCEQ,TTCC,
     &                            C1,IDOFS,IDOFM1,IDOFM2,
     &                            NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL RGDERR (-2,MINEX(INODM),MINEX(INODS),
     &                               I,LPU,I)
                        GO TO 1000
                     ENDIF
                  ELSE
                     CALL RGDERR (3,MINEX(INODS),NNDOF,NNDOF,LPU,IERR)
                     GO TO 1000
                  ENDIF
  450          CONTINUE
            ENDIF
  500    CONTINUE
      ENDIF
C
C ======================================================================
C     Append explicitly defined constraint equations to MPMCEQ/MMCEQ/
C     TTCC, and expand/contract on any specified master dofs such that
C     all master dofs included in the final version of the constraint
C     array MMCEQ are free dofs
C ======================================================================
C
      IF (NCEX .GT. 0) THEN
         CALL RGDXPL (MPAR,MINEX,MADOF,
     &                MPMCEX,MMCEX,TTCCX,
     &                MSC,MPMCEQ,MMCEQ,TTCC,
     &                NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
         IF (IERR .LT. 0)  GO TO 1000
      ENDIF
C
C --- restore MSC, set NCEQ and NMMCEQ and complete MPMCEQ
C
      DO 600 I=1,MPAR(3)
         IF (MSC(I) .LT. 0) THEN
            MSC(I) = 0
         ELSEIF (MSC(I) .GT. 10) THEN
            MSC(I) = MSC(I) - 10
         ENDIF
  600 CONTINUE
      MPAR(7)  = ICEQ-1
      MPAR(16) = IPCUR-1
      IF (ICEQ .GT. 1) THEN
         MPMCEQ(ICEQ) = IPCUR
      ENDIF
C
 1000 RETURN
      END
      SUBROUTINE RGDCE2 (MSC,MPMCEX,MMCEX,TTCCX,
     &                   MPMCEQ,MMCEQ,TTCC,
     &                   C,IDOFS,IDOFM1,IDOFM2,
     &                   NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDCE2                GROUP 9 / PRIVATE
C
C     TASK :  To establish the ICEQ'th constraint equation for a dof
C             of a slave node of a rigid, 2-D, element.
C             The slave (dof no. IDOFS) is coupled to the master node
C             dofs, IDOFM1 and IDOFM2, by the equation
C
C                   D-IDOFS = D-IDOFM1 + C*(D-IDOFM2)
C
C             The constraint equation is stored in MMCEQ and TTCC, from
C             index IPCUR onwards.
C             IPCUR and ICEQ are updated.
C
C
C     ROUTINES CALLED/REFERENCED :  RGDC0, RGDERR   (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-06-25 / 1.0
C                       92-04-27 / 1.1   K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           ICEQ,IDOFS,IDOFM1,IDOFM2,IERR,IPCUR,
     &                  LPU,MAXSIZ,NCEX
      INTEGER           MMCEQ(*),MMCEX(*),MPMCEQ(*),MPMCEX(*),MSC(*)
      DOUBLE PRECISION  C,  TTCC(*),TTCCX(*)
C                                                ! local variables
      INTEGER           JP
      DOUBLE PRECISION  C0,ONE,ZERO
C
      PARAMETER         (ZERO = 0.0D0, ONE = 1.0D0)
C ----------------------------------------------------------------------
C
      IERR         = 0
      MPMCEQ(ICEQ) = 0
      MMCEQ(IPCUR) = 0
      TTCC(IPCUR)  = ZERO
C
      IF (IPCUR+2 .GT. MAXSIZ) THEN
         CALL RGDERR (4,MAXSIZ,MAXSIZ,MAXSIZ,LPU,IERR)
         GO TO 100
      ENDIF
C                                           ! First master dof (IDOFM1)
      IF (MSC(IDOFM1) .GT. 0) THEN
         MPMCEQ(ICEQ)   = IPCUR
         MMCEQ(IPCUR)   = IDOFS
         TTCC(IPCUR)    = ZERO
         MMCEQ(IPCUR+1) = IDOFM1
         TTCC(IPCUR+1)  = ONE
         IPCUR          = IPCUR+2
         MSC(IDOFS)     = -3
      ELSEIF (MSC(IDOFM1) .EQ. 0) THEN
         MSC(IDOFS)     = 0
      ELSEIF (MSC(IDOFM1) .EQ. (-1)) THEN
         CALL RGDC0 (MPMCEX,MMCEX,TTCCX,NCEX,IDOFM1,C0)
         MPMCEQ(ICEQ)   = IPCUR
         MMCEQ(IPCUR)   = IDOFS
         TTCC(IPCUR)    = C0
         IPCUR          = IPCUR+1
         MSC(IDOFS)     = -1
      ENDIF
C
      IF (IDOFM2 .GT. 0) THEN
C                                           ! Second master dof (IDOFM2)
         IF (MSC(IDOFM2) .GT. 0) THEN
            IF (MSC(IDOFS) .EQ. 0) THEN
               MPMCEQ(ICEQ) = IPCUR
               MMCEQ(IPCUR) = IDOFS
               TTCC(IPCUR)  = ZERO
               IPCUR        = IPCUR+1
            ENDIF
            MMCEQ(IPCUR) = IDOFM2
            TTCC(IPCUR)  = C
            IPCUR        = IPCUR+1
            MSC(IDOFS)   = -3
         ELSEIF (MSC(IDOFM2) .EQ. 0) THEN
            CONTINUE
         ELSEIF (MSC(IDOFM2) .EQ. -1) THEN
            CALL RGDC0 (MPMCEX,MMCEX,TTCCX,NCEX,IDOFM2,C0)
            IF (MPMCEQ(ICEQ) .GT. 0) THEN
               JP = MPMCEQ(ICEQ)
            ELSE
               JP           = IPCUR
               MPMCEQ(ICEQ) = IPCUR
               IPCUR        = IPCUR+1
            ENDIF
            MMCEQ(JP)  = IDOFS
            TTCC(JP)   = TTCC(JP) + C*C0
            IF (MSC(IDOFS) .EQ. 0) THEN
               MSC(IDOFS) = -1
            ENDIF
         ENDIF
      ENDIF
C
      IF (MPMCEQ(ICEQ) .GT. 0)  ICEQ = ICEQ+1
C
  100 RETURN
      END
      SUBROUTINE RGDCE3 (MSC,MPMCEX,MMCEX,TTCCX,
     &                   MPMCEQ,MMCEQ,TTCC,
     &                   C1,C2,IDOFS,IDOFM1,IDOFM2,IDOFM3,
     &                   NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDCE3                GROUP 9 / PRIVATE
C
C     TASK :  To establish the ICEQ'th constraint equation for a dof
C             of a slave node of a rigid, 3-D, element.
C             The slave (dof no. IDOFS) is coupled to the master node
C             dofs, IDOFM1, IDOFM2 and IDOFM3, by the equation
C
C               D-IDOFS = D-IDOFM1 + C1*(D-IDOFM2) + C2*(D-IDOFM3)
C
C             The constraint equation is stored in MMCEQ and TTCC, from
C             index IPCUR onwards.
C             IPCUR and ICEQ are updated.
C
C
C     ROUTINES CALLED/REFERENCED :  RGDC0, RGDERR   (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-07-02 / 1.0
C                       92-05-05 / 1.1   K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           ICEQ,IDOFS,IDOFM1,IDOFM2,IDOFM3,IERR,IPCUR,
     &                  LPU,MAXSIZ,NCEX
      INTEGER           MMCEQ(*),MMCEX(*),MPMCEQ(*),MPMCEX(*),MSC(*)
      DOUBLE PRECISION  C1,C2,   TTCC(*),TTCCX(*)
C                                                ! local variables
      INTEGER           JP
      DOUBLE PRECISION  C0,ONE,ZERO
C
      PARAMETER         (ZERO = 0.0D0, ONE = 1.0D0)
C ----------------------------------------------------------------------
C
      IERR         = 0
      MPMCEQ(ICEQ) = 0
      MMCEQ(IPCUR) = 0
      TTCC(IPCUR)  = ZERO
C
      IF (IPCUR+3 .GT. MAXSIZ) THEN
         CALL RGDERR (4,MAXSIZ,MAXSIZ,MAXSIZ,LPU,IERR)
         GO TO 100
      ENDIF
C                                           ! First master dof (IDOFM1)
      IF (MSC(IDOFM1) .GT. 0) THEN
         MPMCEQ(ICEQ)   = IPCUR
         MMCEQ(IPCUR)   = IDOFS
         TTCC(IPCUR)    = ZERO
         MMCEQ(IPCUR+1) = IDOFM1
         TTCC(IPCUR+1)  = ONE
         IPCUR          = IPCUR+2
         MSC(IDOFS)     = -3
      ELSEIF (MSC(IDOFM1) .EQ. 0) THEN
         MSC(IDOFS)     = 0
      ELSEIF (MSC(IDOFM1) .EQ. (-1)) THEN
         CALL RGDC0 (MPMCEX,MMCEX,TTCCX,NCEX,IDOFM1,C0)
         MPMCEQ(ICEQ)   = IPCUR
         MMCEQ(IPCUR)   = IDOFS
         TTCC(IPCUR)    = C0
         IPCUR          = IPCUR+1
         MSC(IDOFS)     = -1
      ENDIF
C
      IF (IDOFM2 .GT. 0) THEN
C                                           ! Second master dof (IDOFM2)
         IF (MSC(IDOFM2) .GT. 0) THEN
            IF (MSC(IDOFS) .EQ. 0) THEN
               MPMCEQ(ICEQ) = IPCUR
               MMCEQ(IPCUR) = IDOFS
               TTCC(IPCUR)  = ZERO
               IPCUR        = IPCUR+1
            ENDIF
            MMCEQ(IPCUR) = IDOFM2
            TTCC(IPCUR)  = C1
            IPCUR        = IPCUR+1
            MSC(IDOFS)   = -3
         ELSEIF (MSC(IDOFM2) .EQ. 0) THEN
            CONTINUE
         ELSEIF (MSC(IDOFM2) .EQ. -1) THEN
            CALL RGDC0 (MPMCEX,MMCEX,TTCCX,NCEX,IDOFM2,C0)
            IF (MPMCEQ(ICEQ) .GT. 0) THEN
               JP = MPMCEQ(ICEQ)
            ELSE
               JP           = IPCUR
               MPMCEQ(ICEQ) = IPCUR
               IPCUR        = IPCUR+1
            ENDIF
            MMCEQ(JP)  = IDOFS
            TTCC(JP)   = TTCC(JP) + C1*C0
            MSC(IDOFS) = -1
         ENDIF
      ENDIF
C
      IF (IDOFM3 .GT. 0) THEN
C                                           ! Third master dof (IDOFM3)
         IF (MSC(IDOFM3) .GT. 0) THEN
            IF (MSC(IDOFS) .EQ. 0) THEN
               MPMCEQ(ICEQ) = IPCUR
               MMCEQ(IPCUR) = IDOFS
               TTCC(IPCUR)  = ZERO
               IPCUR        = IPCUR+1
            ENDIF
            MMCEQ(IPCUR) = IDOFM3
            TTCC(IPCUR)  = C2
            IPCUR        = IPCUR+1
            MSC(IDOFS)   = -3
         ELSEIF (MSC(IDOFM3) .EQ. 0) THEN
            CONTINUE
         ELSEIF (MSC(IDOFM3) .EQ. -1) THEN
            CALL RGDC0 (MPMCEX,MMCEX,TTCCX,NCEX,IDOFM3,C0)
            IF (MPMCEQ(ICEQ) .GT. 0) THEN
               JP = MPMCEQ(ICEQ)
            ELSE
               JP           = IPCUR
               MPMCEQ(ICEQ) = IPCUR
               IPCUR        = IPCUR+1
            ENDIF
            MMCEQ(JP)  = IDOFS
            TTCC(JP)   = TTCC(JP) + C2*C0
            MSC(IDOFS) = -1
         ENDIF
      ENDIF
C
      IF (MPMCEQ(ICEQ) .GT. 0)  ICEQ = ICEQ+1
C
  100 RETURN
      END
      SUBROUTINE RGDC0 (MPMCEQ,MMCEQ,TTCC,NCEQ,IDOFS,C0)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RGDC0                  GROUP 9 / PRIVATE
C
C     TASK :  To return the value, C0, of prescribed dof IDOFS
C
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-06-20 / 1.0
C                       92-04-27 / 2.0   K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IDOFS,NCEQ,         MMCEQ(*),MPMCEQ(*)
      DOUBLE PRECISION  C0,      TTCC(*)
C                                                ! local variables
      INTEGER           I,IP
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      C0 = ZERO
C
      DO 10 I=1,NCEQ
         IP = MPMCEQ(I)
         IF (MMCEQ(IP) .EQ. IDOFS) THEN
            C0 = TTCC(IP)
            GO TO 100
         ENDIF
   10 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE RGDCNV (MINEX,MADOF,MSC,
     &                   MPMCEX,MMCEX,TTCCX,
     &                   MMCEQ,TTCC,C,IDOFS,NCEX,
     &                   NANOD,MAXSIZ,IPSLV,LPU,IPCUR,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDCNV                GROUP 9 / PRIVATE
C
C     TASK :  To resolve the situation in which an explicitly defined
C             "master dof" in a particular constraint equation is also
C             defined as a "slave dof" in another constraint equation.
C             The slave dof (IDOFS) appearing as a master dof with
C             coefficient C in MMCEX/TTCCX is replaced by its
C             appropriate master dofs (which must all be free dofs) with
C             new coefficients in MMCEQ/TTCC.
C
C     ROUTINES CALLED/REFERENCED :  RGDNDN and RGDERR   (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-06-24 / 1.0
C                       92-04-26 / 2.0    K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IDOFS,IERR,IPCUR,IPSLV,LPU,MAXSIZ,NANOD,NCEX
      INTEGER           MADOF(*),MINEX(*),MMCEQ(*),MMCEX(*),MPMCEX(*),
     &                  MSC(*)
      DOUBLE PRECISION  C,  TTCC(*),TTCCX(*)
C                                                ! local variables
      INTEGER           I,INODS,JDOF,JP,LDOFS
C ----------------------------------------------------------------------
C
      DO 100 I=1,NCEX
         JP = MPMCEX(I)
         IF (MMCEX(JP) .EQ. IDOFS) THEN
            TTCC(IPSLV) = TTCC(IPSLV) + C*TTCCX(JP)
            IF (MPMCEX(I+1)-MPMCEX(I) .GT. 1) THEN
               DO 80 JP=MPMCEX(I)+1,MPMCEX(I+1)-1
                  JDOF = MMCEX(JP)
                  IF (MSC(JDOF) .GT. 0) THEN
                     IF (IPCUR .GT. MAXSIZ) THEN
                        CALL RGDERR (4,MAXSIZ,MAXSIZ,MAXSIZ,LPU,IERR)
                        GO TO 200
                     ENDIF
                     MMCEQ(IPCUR) = JDOF
                     TTCC(IPCUR)  = C*TTCCX(JP)
                     IPCUR        = IPCUR+1
                  ELSEIF (MSC(JDOF) .LT. 0) THEN
                     CALL RGDNDN (MINEX,MADOF,NANOD,IDOFS,INODS,LDOFS)
                     CALL RGDERR (9,LDOFS,INODS,INODS,LPU,IERR)
                     GO TO 200
                  ENDIF
   80          CONTINUE
            ENDIF
            GO TO 200
         ENDIF
  100 CONTINUE
C
      CALL RGDNDN (MINEX,MADOF,NANOD,IDOFS,INODS,LDOFS)
      CALL RGDERR (6,LDOFS,INODS,INODS,LPU,IERR)
C
  200 RETURN
      END
      SUBROUTINE RGDNDN (MINEX,MADOF,NANOD,IDOF,INOD,LDOF)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDNDN                GROUP 9 / PRIVATE
C
C     TASK :  To determine the (external) node number (INOD) and local
C             dof-number (within node) (LDOF) corresponding to
C             (internal) dof-number IDOF.
C
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   92-04-26 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IDOF,INOD,LDOF,NANOD
      INTEGER           MADOF(*),MINEX(*)
C                                                ! local variables
C
      INTEGER           I
C ----------------------------------------------------------------------
      INOD = 0
      LDOF = 0
      DO 10 I=2,NANOD+1
         IF (IDOF .LT. MADOF(I)) THEN
            INOD = MINEX(I-1)
            LDOF = IDOF - MADOF(I-1) + 1
            GO TO 100
         ENDIF
   10 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE RGDSWP (MINEX,MSC,MPMCEX,MMCEX,TTCCX,
     &                   NCEX,INODS,IDOFS,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDSWP                GROUP 9 / PRIVATE
C
C     TASK :  To swap the slave with an appropriate member of its master
C             dofs in an explicit constraint equation for which the
C             original slave is also a slave dof of a rigid element.
C
C
C     ROUTINES CALLED/REFERENCED :   RGDERR   (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   92-05-01 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IDOFS,IERR,INODS,LPU,NCEX
      INTEGER           MINEX(*),MMCEX(*),MPMCEX(*),MSC(*)
      DOUBLE PRECISION  TTCCX(*)
C
C                                                ! local variables
C
      INTEGER           I,J,JE,JJ,JP,LDOF
      DOUBLE PRECISION  C,ONE,ZERO
C
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
C
      DO 100 I=1,NCEX
         JP    = MPMCEX(I)
         IF (MMCEX(JP) .EQ. IDOFS) THEN
            JE = MPMCEX(I+1) - 1
            JJ = 0
            DO 10 J=JP+1,JE
               LDOF = MMCEX(J)
               IF (MSC(LDOF) .GT. 10) THEN
                  CALL RGDERR (11,MINEX(INODS),J,J,LPU,IERR)
                  GO TO 500
               ENDIF
               IF (JJ .EQ. 0) THEN
                  IF (MSC(LDOF) .GT. 0) THEN
                     IF (TTCCX(J) .NE. ZERO) THEN
                        JJ = J
                        C  = TTCCX(J)
                     ENDIF
                  ENDIF
               ENDIF
   10       CONTINUE
            IF (JJ .GT. 0) THEN
               LDOF      = MMCEX(JJ)
               MMCEX(JJ) = IDOFS
               MMCEX(JP) = LDOF
               TTCCX(JP) = -TTCCX(JP) / C
               DO 20 J=JP+1,JE
                  IF (J .EQ. JJ) THEN
                     TTCCX(J) = ONE / C
                  ELSE
                     TTCCX(J) = -TTCCX(J) / C
                  ENDIF
   20          CONTINUE
               MSC(LDOF) = -2
            ELSE
               CALL RGDERR (12,MINEX(INODS),J,J,LPU,IERR)
               GO TO 500
            ENDIF
            GO TO 500
         ENDIF
  100 CONTINUE
C
  500 RETURN
      END
      SUBROUTINE RGDXPL (MPAR,MINEX,MADOF,
     &                   MPMCEX,MMCEX,TTCCX,
     &                   MSC,MPMCEQ,MMCEQ,TTCC,
     &                   NCEX,MAXSIZ,LPU,ICEQ,IPCUR,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDXPL                GROUP 9 / PRIVATE
C
C     TASK :  To append explicitly defined constraint equations to
C             MPMCEQ/MMCEQ/TTCC, and expand/contract on any specified
C             master dofs such that all master dofs included in the
C             final version of the constraint array MMCEQ are free dofs.
C
C
C     ROUTINES CALLED/REFERENCED :   RGDCNV and RGDERR   (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-07-02 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - Single Precision  (REAL)
C                       D - Double Precision
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           ICEQ,IERR,IPCUR,LPU,MAXSIZ,NCEX
      INTEGER           MADOF(*),MINEX(*),MMCEQ(*),MMCEX(*),MPAR(*),
     &                  MPMCEQ(*),MPMCEX(*),MSC(*)
      DOUBLE PRECISION  TTCC(*),TTCCX(*)
C
C                                                ! local variables
C
      INTEGER           I,IDOFM,IDOFS,ISC,JP,KP,N,NM
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
C
      DO 50 I=1,NCEX
         JP    = MPMCEX(I)
         KP    = IPCUR
         IDOFS = MMCEX(JP)
         IF (IDOFS .GT. 0) THEN
            IF (IPCUR .GT. MAXSIZ) THEN
               CALL RGDERR (4,MAXSIZ,MAXSIZ,MAXSIZ,LPU,IERR)
               GO TO 100
            ENDIF
            MPMCEQ(ICEQ) = IPCUR
            MMCEQ(IPCUR) = IDOFS
            TTCC(IPCUR)  = TTCCX(JP)
            ICEQ         = ICEQ+1
            IPCUR        = IPCUR+1
C
            NM = MPMCEX(I+1) - MPMCEX(I) - 1
C
            IF (NM .GT. 0) THEN
               DO 25 JP=MPMCEX(I)+1,MPMCEX(I+1)-1
                  IDOFM = MMCEX(JP)
                  ISC   = MSC(IDOFM)
                  IF (ISC .GT. 0) THEN
                     IF (IPCUR .GT. MAXSIZ) THEN
                        CALL RGDERR (4,MAXSIZ,MAXSIZ,MAXSIZ,LPU,IERR)
                        GO TO 100
                     ENDIF
                     MMCEQ(IPCUR) = IDOFM
                     TTCC(IPCUR)  = TTCCX(JP)
                     IPCUR        = IPCUR+1
                  ELSEIF (ISC.EQ.(-1) .OR. ISC.EQ.(-2)) THEN
                     CALL RGDCNV (MINEX,MADOF,MSC,MPMCEX,MMCEX,TTCCX,
     &                            MMCEQ,TTCC,TTCCX(JP),IDOFM,NCEX,
     &                            MPAR(1),MAXSIZ,KP,LPU,IPCUR,IERR)
                     IF (IERR .LT. 0)  GO TO 100
                  ELSEIF (MSC(IDOFM) .EQ. (-3)) THEN
                     N = ICEQ-2
                     CALL RGDCNV (MINEX,MADOF,MSC,MPMCEQ,MMCEQ,TTCC,
     &                            MMCEQ,TTCC,TTCCX(JP),IDOFM,N,
     &                            MPAR(1),MAXSIZ,KP,LPU,IPCUR,IERR)
                     IF (IERR .LT. 0)  GO TO 100
                  ENDIF
   25          CONTINUE
            ENDIF
C
C --------- if slave dof is suppressed, remove it
C
            IF ((IPCUR-1) .EQ. KP) THEN
               IF (TTCC(KP) .EQ. ZERO) THEN
                  ICEQ       = ICEQ-1
                  IPCUR      = IPCUR-1
                  MSC(IDOFS) = 0
               ENDIF
            ENDIF
         ENDIF
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE RGDERR (NUM,IVAL1,IVAL2,IVAL3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :   RGDERR                GROUP 9 / PRIVATE
C
C     TASK :  Print error messages and set the error flag for the RGD-
C             routines
C
C
C     ROUTINES CALLED/REFERENCED :  ABS    (Fortran library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   91-07-20 / 1.0
C                       92-05-04 / 1.1   K.Bell
C                       02-01-16 / 1.2   K.M.Okstad
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IVAL1,IVAL2,IVAL3,LPU,NUM
C
C                                                ! local variables
      INTEGER           N
C ----------------------------------------------------------------------
      IERR =-ABS(NUM)
C
      IF (LPU .GT. 0) THEN
         IF (NUM .GT. 0) THEN
            WRITE (LPU,600)
         ENDIF
         N = ABS(NUM)
C
         IF (N .EQ. 1) THEN
            WRITE (LPU,650)
            WRITE (LPU,601) IVAL1,IVAL2
         ELSEIF (N .EQ. 2) THEN
            WRITE (LPU,602) IVAL1,IVAL2
         ELSEIF (N .EQ. 3) THEN
            WRITE (LPU,603) IVAL1,IVAL2
         ELSEIF (N .EQ. 4) THEN
            WRITE (LPU,604) IVAL1
         ELSEIF (N .EQ. 5) THEN
            WRITE (LPU,605) IVAL1,IVAL2
         ELSEIF (N .EQ. 6) THEN
            WRITE (LPU,650)
            WRITE (LPU,606) IVAL1,IVAL2
         ELSEIF (N .EQ. 7) THEN
            WRITE (LPU,650)
            WRITE (LPU,607) IVAL1,IVAL2
         ELSEIF (N .EQ. 8) THEN
            WRITE (LPU,608)
         ELSEIF (N .EQ. 9) THEN
            WRITE (LPU,609) IVAL1,IVAL2
         ELSEIF (N .EQ. 10) THEN
            WRITE (LPU,610) IVAL1,IVAL2
         ELSEIF (N .EQ. 11) THEN
            WRITE (LPU,611) IVAL1
         ELSEIF (N .EQ. 12) THEN
            WRITE (LPU,612) IVAL1
         ELSEIF (N .EQ. 13) THEN
            WRITE (LPU,613) IVAL1,IVAL2
            WRITE (LPU,619) IVAL3
         ELSEIF (N .EQ. 14) THEN
            WRITE (LPU,614) IVAL1,IVAL2
            WRITE (LPU,619) IVAL3
         ENDIF
      ENDIF
#ifdef FT_DEBUG
      IF (NUM .GT. 0) IERR = IERR/0 ! Force a core dump for easy debugging
#endif
C
      RETURN
C -------------------------------------------------------------- formats
  600 FORMAT (///' *** ERROR DURING RIGID ELEMENT CONSTRAINT HANDLING')
  601 FORMAT (5X,'Explicit constraint eqn. no.',I5,' specifies non-'/
     &        5X,'existent dof no.',I8)
  602 FORMAT (5X,'Masternode',I8,';   slavenode',I8)
  603 FORMAT (5X,'Slavenode',I8,'  has',I5,'  dofs')
  604 FORMAT (5X,'Insufficient storage allocated ('
     &       ,I8,' ) for constraint information')
  605 FORMAT (5X,'Dof no.',I3,' of rigid element masternode',I8,' is'/
     &        5X,'itself defined as a slave in a constraint equation')
  606 FORMAT (5X,'Constraint eqn. for dof',I3,' of node',I8,
     &           ' not found')
  607 FORMAT (5X,'Explicit constraint eqn. no.',I5,' has',I6,
     &           ' master dofs')
  608 FORMAT (5X,'Second coordinate axis is neither y- nor z-axis')
  609 FORMAT (5X,'Illegal ("two-level") coupling via dof',I3,
     &           ' of node',I8)
  610 FORMAT (5X,'Dof no.',I3,' of rigid element slavenode',I8,' is'/
     &        5X,'illegally (explicitly) constrained')
  611 FORMAT (5X,'Illegal coupling of rigid elements at node',I8)
  612 FORMAT (5X,'Illegal (double?) constraint at node',I8)
  613 FORMAT (5X,'Illegal no. of dofs (',I3,' ) at master-node',I8)
  614 FORMAT (5X,'Illegal no. of dofs (',I3,' ) at slave-node',I8)
  619 FORMAT (5X,'Detected for element',I8)
C
  650 FORMAT (5X,'Inconsistent (corrupt?) data')
C
      END
