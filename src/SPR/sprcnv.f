C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRCNV   ( L     , A     , MSPAR , MTREES, MSIFA ,
     $                        M     , KSA   , LPU   , IERR            )

C     ------------------------------------------------------------------

      IMPLICIT          NONE

      INTEGER           M     , KSA   , LPU   , IERR
      INTEGER           MSPAR(*)              ,
     $                  MTREES(*)             , MSIFA(*)
      DOUBLE PRECISION  L(*)                  , A(M*M)

C ======================================================================
C  S A M  library routine :  SPRCNV                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To convert the compressed stored symmetric matrix L to full
C     storage in A.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 Feb. 25, 2003 (kmo)
C                 Get pointers from MSPAR instead of recalculating them.
C                 Jul. 30, 2003 (kmo)
C                 Extended the meaning of KSA to account for sign change
C                 and optional initialization.
C                 Jun. 14, 2004 (kmo)
C                 Added SPRCONV with leading dimension for matrix A.
C
C     L - DOUBLE PRECISION(*)
C      Entry : A symmetric matrix where the lower triangle is stored
C              in compressed format.
C              N O T E : If ABS(KSA)=1 then SM(NEQ+1) should be the
C                        first location of SM accessed by this routine.
C                        If ABS(KSA)=2 it is assumed that L is only
C                        provided as vector of length NEQ.
C      Exit  : Unchanged.
C     A - DOUBLE PRECISION(M*M)
C      Entry : Irrelevant.
C      Exit  : Holds the symmetric M times M matrix A. Note that L
C              is copied to both triangles of A, i.e. the symmetry
C              of A is not taken into account by only copy from L
C              to, for instance the lower triangle of A.
C     MSPAR - INTEGER(*)
C      Entry : Matrix of sparse parameters as output from SPRSMB.
C      Exit  : Unchanged.
C     MTREES - INTEGER(NTREES)
C      Entry : As output from SPRSMB.
C      Exit  : Unchanged.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : As output from SPRSMB.
C      Exit  : Unchanged.
C     M - INTEGER
C      Entry : The leading dimension of A. M .GE. NEQ, where NEQ is the
C              number of equations. If M .GT. NEQ, L is returned in A
C              with leading dimension M and in NEQ columns. If M.LT.NEQ
C              only the upper left (M,M) submatrix of A is computed.
C      Exit  : Unchanged.
C     KSA - INTEGER
C      Entry : Storage code.
C              = 1 : A is a symmetric matrix that is not diagonal.
C              = 2 : A is diagonal.
C              > 9 : Same as KSA-10, but don't initialize the matrix
C              < 0 : Same as -KSA, but negate the matrix elements
C      Exit  : Unchanged.
C     LPU - INTEGER
C      Entry : Output for error messages.
C      Exit  : Unchanged.
C     IERR - INTEGER
C      Entry : Need not be set.
C      Exit  : Error flag, it is set to zero in case of a normal return.
C
C     Working arrays
C     --------------
C
C     None
C
C     Procedures
C     ----------
C
C     SPRCNQ
C     SPRER1
C
C     Intrinsic
C     ---------
C
C     None
C
C     Include Blocks
C     --------------
C
C     None
C
C     Common Blocks
C     -------------
C
C     None
C
C     ------------------------------------------------------------------

      CALL SPRCONV  ( L     , A     , MSPAR , MTREES, MSIFA ,
     $                M     , M     , KSA   , LPU   , IERR    )

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRCNV
      END


      SUBROUTINE   SPRCONV  ( L     , A     , MSPAR , MTREES, MSIFA ,
     $                        LDA   , M     , KSA   , LPU   , IERR    )


      IMPLICIT          NONE

      INTEGER           LDA   , M     , KSA   , LPU   , IERR
      INTEGER           MSPAR(*)              ,
     $                  MTREES(*)             , MSIFA(*)
      DOUBLE PRECISION  L(*)                  , A(LDA*M)

C     ------------------------------------------------------------------

      INTEGER           I, J  , LINDX , NEQ   , NOFSUB,
     $                  NSUPER, XLINDX, XLNZ  , XSUPER
      DOUBLE PRECISION  ZERO
      PARAMETER       ( ZERO = 0.0D0 )
      EXTERNAL          SPRCNQ, SPRER1

C     ==================================================================

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -11
         IERR = -1
         CALL SPRER1 ( 30, 'SPRCNV', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         IERR = 0
      ENDIF

      I = MOD(ABS(KSA),10)
      IF ( I .LT. 1 .OR. I .GT. 2 ) THEN
         MSPAR(1) = -11
         IERR = -2
         CALL SPRER1 ( 26, 'SPRCNV', KSA, 0, 0, LPU, IERR )
         RETURN
      ENDIF

      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NOFSUB = MSPAR(15)

C     -----------------------
C     SET POINTERS TO MTREES.
C     -----------------------
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)

C     ----------------------
C     SET POINTERS TO MSIFA.
C     ----------------------
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

C     ==================================================================

      IF ( ABS(KSA) .LT. 10 ) THEN
C        ---------------------
C        INITIALIZE A TO ZERO.
C        ---------------------
         DO 110 J = 1, MIN(M,NEQ)
            DO 100 I = 1, MIN(M,NEQ)
               A(J*LDA-LDA+I) = ZERO
  100       CONTINUE
  110    CONTINUE
      ENDIF

      IF ( MOD(ABS(KSA),10) .EQ. 1 ) THEN

C        -----------------------------------
C        CONVERT L TO A BY A CALL TO SPRCNQ.
C        -----------------------------------
         CALL SPRCNQ ( NSUPER, LDA, M, NOFSUB, MTREES(XSUPER),
     $                 MSIFA(XLNZ), MTREES(XLINDX), MSIFA(LINDX),
     $                 L, A, SIGN(1,KSA) )

      ELSE

C        -----------------------
C        L IS A DIAGONAL MATRIX.
C        -----------------------
         DO 200 I = 1, MIN(M,NEQ)
            A(I*LDA-LDA+I) = L(I)*SIGN(1,KSA)
  200    CONTINUE

      ENDIF

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRCONV
      END


      SUBROUTINE   SPRCNQ   ( NSUPER, LDA,M , NOFSUB, XSUPER, XLNZ  ,
     $                        XLINDX, LINDX , LNZ   , A     , SCALE   )

      IMPLICIT          NONE

      INTEGER           NSUPER, LDA,M , NOFSUB, SCALE
      INTEGER           XSUPER(NSUPER+1)      , XLNZ(NSUPER+1)        ,
     $                  XLINDX(NSUPER+1)      , LINDX(NOFSUB)
      DOUBLE PRECISION  LNZ(*)                , A(LDA,M)

C     ----------------------------------------------------
C --- TO CONVERT A LOWER SYMMETRIC MATRIX TO FULL STORAGE.
C     ----------------------------------------------------
      INTEGER           FSTVAR, I     , II    , IPLNZ , ISTRT ,
     $                  ISTOP , J     , JS    , LSTVAR

      DO 400 JS = 1, NSUPER
         FSTVAR = XSUPER(JS)
         LSTVAR = MIN(M,XSUPER(JS+1)-1)
         ISTRT  = XLINDX(JS)
         ISTOP  = XLINDX(JS+1) - 1
         IPLNZ  = XLNZ(JS)
         DO 300 J = FSTVAR, LSTVAR
            A(J,J) = LNZ(IPLNZ)*SCALE
            IPLNZ = IPLNZ + 1
            DO 200 I = ISTRT+1, ISTOP
               II = LINDX(I)
               IF (II .LE. M) THEN
                  A(J,II) = LNZ(IPLNZ)*SCALE
                  A(II,J) = LNZ(IPLNZ)*SCALE
               ENDIF
               IPLNZ = IPLNZ + 1
  200       CONTINUE
            ISTRT = ISTRT + 1
  300    CONTINUE
  400 CONTINUE

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRCNQ
      END
