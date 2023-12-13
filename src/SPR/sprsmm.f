C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRSMM (A,B,X,Y,W,MIW,
     &                   MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                   N,KSA,KSB,IFLAG,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRSMM                 GROUP 9 / PRIVATE
C
C     T A S K :  To determine
C
C                   X := (A-1)*UbT*X         IFLAG = 1
C                   Y  = (A-1)*UbT*X         IFLAG = 2
C                   X := Ub*(A-1)*UbT*X      IFLAG = 3
C                   Y  = Ub*(A-1)*UbT*X      IFLAG = 4
C                   X := (Ua-T)*B*(Ua-1)*X   IFLAG = 5
C                   Y  = (Ua-T)*B*(Ua-1)*X   IFLAG = 6
C     where
C           X and Y are vectors,
C           A is a symmetric sparse (compressed) matrix, represented by
C           its factors La*Da*LaT, stored in A (KSA = 2),
C           Ub is the (Cholesky) factor of B, stored in B (KSB<0), and
C           Ua is the Cholesky factor of A, stored in A (KSA=3)
C     If  KSB=1 or -1, B / Ub is a symmetric / upper triangular sparse
C     matrix (same storage pattern as A), and
C     if  KSB=2 or -2, B / Ub is a diagonal matrix
C
C     This is the  s p a r s e  version of SPSMM (used by LANCZ2)
C
C
C     ROUTINES CALLED/REFERENCED :  SCOPY/DCOPY        (BLAS)
C                                   SPRTMM, SPRPRM   (SAM/SPR)
C                                   SPRFS1, SPRBS1   (SAM/SPR)
C                                   SPRFS2, SPRBS2   (SAM/SPR)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   96-12-22 / 1.0
*                       99-01-06 / 1.1   A.C.Damhaug
*                                        MSPAR(50) -> MSPAR(*).
C                       03-01-17 / 1.2   K.M.Okstad
C                                        Added a separate sparse
C                                        structure for matrix B.
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,IFLAG,LPU,KSA,KSB,N,MIW(*)
      INTEGER           MPARA(*),MTREEA(*),MSIFA(*)
      INTEGER           MPARB(*),MTREEB(*),MSIFB(*)
      DOUBLE PRECISION  A(*),B(*),Y(N),X(N),W(N)
C
      include 'integer_kind.inc'
C
      INTEGER           LPB,LSB
      INTEGER(kind=i4)  NN
C ----------------------------------------------------------------------
C
#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRSMM, IFLAG =',IFLAG
#endif
C
      IERR = 0
      IF (ABS(KSB) .EQ. 2) THEN
         LPB = 1
         LSB = -1
      ELSE
         LPB = MPARB(8) + 1
         LSB = 1
      ENDIF
C
      IF (IFLAG .EQ. 1) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GOTO 90
         CALL SPRTMM (B(LPB),X,X,W,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,2,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRFS1 (MPARA,MTREEA,MSIFA,A,X,N,1,MIW,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRBS1 (MPARA,MTREEA,MSIFA,A,X,N,1,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
      ELSEIF (IFLAG .EQ. 2) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GOTO 90
         CALL SPRTMM (B(LPB),X,W,Y,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,22,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRFS1 (MPARA,MTREEA,MSIFA,A,Y,N,1,MIW,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRBS1 (MPARA,MTREEA,MSIFA,A,Y,N,1,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
      ELSEIF (IFLAG .EQ. 3) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GOTO 90
         CALL SPRTMM (B(LPB),X,Y,W,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,22,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRFS1 (MPARA,MTREEA,MSIFA,A,W,N,1,MIW,X,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRBS1 (MPARA,MTREEA,MSIFA,A,W,N,1,X,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRTMM (B(LPB),W,Y,X,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,11,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
      ELSEIF (IFLAG .EQ. 4) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GOTO 90
         CALL SPRTMM (B(LPB),X,Y,W,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,22,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRFS1 (MPARA,MTREEA,MSIFA,A,W,N,1,MIW,Y,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRBS1 (MPARA,MTREEA,MSIFA,A,W,N,1,Y,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRTMM (B(LPB),W,X,Y,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,11,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
      ELSEIF (IFLAG .EQ. 5) THEN
         IF (KSA.NE.3 .OR. KSB.LT.0)  GOTO 90
         CALL SPRBS2 (MPARA,MTREEA,MSIFA,A,X,N,1,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRPRM (B(LPB),X,X,W,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,1,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRFS2 (MPARA,MTREEA,MSIFA,A,X,N,1,MIW,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
      ELSEIF (IFLAG .EQ. 6) THEN
         IF (KSA.NE.3 .OR. KSB.LT.0)  GOTO 90
         NN = int(N,i4)
         CALL DCOPY (NN,X,1_i4,W,1_i4)
         CALL SPRBS2 (MPARA,MTREEA,MSIFA,A,W,N,1,Y,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRPRM (B(LPB),W,Y,Y,MPARB,MTREEB,MSIFB,
     &                N,1,LSB,88,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
         CALL SPRFS2 (MPARA,MTREEA,MSIFA,A,Y,N,1,MIW,W,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 95
      ENDIF
      GOTO 100
C                                                ! error return
   90 IERR = -1
      IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
         WRITE (LPU,691) IFLAG,KSA,KSB
      ENDIF
      GOTO 100
   95 IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
      ENDIF
C
  100 CONTINUE
#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRSMM'
#endif
      RETURN
C
  690 FORMAT(///' *** ERROR return from  S A M library routine  SPRSMM')
  691 FORMAT(5X,'Illegal parameter(s) - IFLAG, KSA, KSB = ',3I5)
C
      END
