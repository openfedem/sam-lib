C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRTMM   ( U     , B     , C     , WA    , MSPAR ,
     $                        MTREES, MSIFA , M     , N     , KSU   ,
     $                        IFLAG , LPU   , IERR                    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           M     , N     , IFLAG , KSU   , LPU   , IERR
      INTEGER           MSPAR(*)              ,
     $                  MTREES(*)             , MSIFA(*)
      DOUBLE PRECISION  U(*)                  ,
     $                  B(M,N)                , C(M,N)                ,
     $                  WA(M)

C ======================================================================
C  S A M  library routine :  SPRTMM                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C                               t
C --- To compute C = UB or C = U B, where A is a compressed matrix, B
C     and C are M,N matrices.
C
C     Due to the use of the work array WA, B and C may be the same
C     array. I. e. the routine may be used to compute the products
C     C = UC or C = UtC.
C
C     Also note that on exit WA contains the last column of C, thus if
C     N = 1, C = WA.
C
C     If N=1, then we may set IFLAG = 11 or IFLAG = 22, then C is not
C     referenced and the result is returned in WA.
C     IFLAG=11 computes WA = UB, and IFLAG=22 computes WA = UtB.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
*                 Jun. 01, 1999 (acd)
*                 Modified SPRTMQ in order to account for evaluation of
*                 internal stiffness only.
C                 Feb. 25, 2003 (kmo)
C                 Get pointers from MSPAR instead of recalculating them.
C                 Added internal permutation array PERI2E.
C                 Added SPRTMR and use that instead of SPRTMQ when
C                 an internal permutation is given.
C                 Replaced the out-zeroing of external components in
C                 SPRTMQ (and now SPRTMR) by equivalent calls to
C                 SPRZR0 and SPRZR1 in SPRTMM instead.
C                 Aug. 27, 2005 (kmo)
C                 Removed calls to SPRZR0 and SPRZR1 again.
C                 Instead it is assumed that the external equations
C                 are ordered last, and that M = # internal equations
C                 when MSPAR(57) = 1.
C
C     U - DOUBLE PRECISION(*)
C      Entry : A symmetric matrix where the lower triangle is stored
C              in compressed format.
C      Exit  : Unchanged.
C     B - DOUBLE PRECISION(M,N)
C      Entry : The rectangular matrix to be premultipied by U or Ut.
C      Exit  : Unchanged, except for if B uses the same storage
C              locations as C which is permitted due to the work
C              array WA.
C     C - DOUBLE PRECISION(M,N)
C      Entry : Irrelevant.
C      Exit  : Assigned to the matrix product UB or UtB.
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
C      Entry : The leading dimension of B and C.
C      Exit  : Unchanged.
C     N - INTEGER
C      Entry : The number of columns in B and C, it is required that
C              N .GT. 0.
C      Exit  : Unchanged.
C     KSU - INTEGER
C      Entry : .GT. 0 : U holds a compressed matrix on sparse format.
C              .LE. 0 : U is diagonal and held as a vector.
C      Exit  : Unchanged.
C     IFLAG - INTEGER
C      Entry : 1 : C <- UB.
C              2 : C <- UtB.
C              The next two options require N=1.
C              11 : WA <- UB.
C              22 : WA <- UtB.
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
C     WA - DOUBLE PRECISION(M)
C      Entry : Not defined.
C      Exit  : Need not be saved, but as noted above its contents on
C              exit is the same as the last column of C.
C              Also note its special use when IFLAG = 11,22.
C
C     Procedures
C     ----------
C
C     SPRTMQ
C     SPRTMR
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

      LOGICAL           AVOIDC
      INTEGER           IN    , JN    , LINDX , NEQ   , NOFSUB, NSUPER,
     $                  NEQFUL, PERI2E, XLINDX, XLNZ  , XSUPER
      EXTERNAL          SPRTMQ, SPRTMR, SPRER1

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRTMM, IFLAG =',IFLAG
#endif

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1).LT.0 ) THEN
         MSPAR(1) = -10
         IERR = -1
         CALL SPRER1 ( 30, 'SPRTMM', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         IERR = 0
      ENDIF

* ------ Used to avoid C, note only if N=1.
      AVOIDC = (IFLAG.EQ.11 .OR. IFLAG.EQ.22) .AND. N.EQ.1

      IF ( .NOT.AVOIDC .AND. IFLAG.NE.1 .AND. IFLAG.NE.2 ) THEN
         MSPAR(1) = -10
         IERR = -2
         CALL SPRER1 ( 25, 'SPRTMM', IFLAG, N, 0, LPU, IERR )
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

C     -------------------------------------------
C     SET POINTERS TO MSIFA.
C     -------------------------------------------
      PERI2E = 1
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

C     ==================================================================

C     Use only the first M (internal) equations when M < NEQ
      NEQFUL = NEQ
      IF ( M.LT.NEQ ) THEN
         IF ( MSPAR(57).EQ.1 ) THEN
            NEQ = M
         ELSE
            MSPAR(1) = -10
            IERR = -3
            CALL SPRER1 ( 24, 'SPRTMM', 0, NEQ, M, LPU, IERR )
            RETURN
         ENDIF
      ENDIF

C     -------------------------------
C     COMPUTE C BY N CALLS TO SPRTMQ.
C     -------------------------------
      DO 200 IN = 1, N

         IF ( MSPAR(58).GE.1 ) THEN
            CALL SPRTMR ( NSUPER, NEQ, NEQFUL, NOFSUB, IFLAG,
     $                    MSIFA(PERI2E), MTREES(XSUPER), MTREES(XLINDX),
     $                    MSIFA(XLNZ), MSIFA(LINDX), U, B(1,IN), WA )
         ELSE
            CALL SPRTMQ ( NSUPER, NEQ, NOFSUB, KSU, IFLAG,
     $                    MTREES(XSUPER), MTREES(XLINDX), MSIFA(XLNZ),
     $                    MSIFA(LINDX), U, B(1,IN), WA )
         ENDIF

         IF ( .NOT.AVOIDC ) THEN
            DO 100 JN = 1, M
               C(JN,IN) = WA(JN)
  100       CONTINUE
         ENDIF

  200 CONTINUE


C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRTMM'
#endif
      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTMM
      END


      SUBROUTINE   SPRTMQ   ( NSUPER, NEQ   , NOFSUB, KSU   , IFLAG ,
     $                        XSUPER, XLINDX, XLNZ  , LINDX , LNZ   ,
     $                        B     , C       )

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , NOFSUB, KSU   , IFLAG
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  XLNZ(NSUPER+1)        , LINDX(NOFSUB)
      DOUBLE PRECISION  LNZ(*)                ,
     $                  B(NEQ)                , C(NEQ)

C     --------------------------------------------------
C --- TO PRE-MULTIPLY A COMPRESSED MATRIX WITH A VECTOR.
C     --------------------------------------------------

      INTEGER           FSTVAR, I     , II    , IPLNZ , ISTRT ,
     $                  ISTOP , J     , JS    , LSTVAR

      DOUBLE PRECISION  ZERO
      PARAMETER       ( ZERO = 0.0D0 )

      IF ( KSU.GT.0 ) THEN
C        -------------
C        INITIALIZE C.
C        -------------
         DO 100 I = 1, NEQ
            C(I) = ZERO
  100    CONTINUE

C        ------------------------------
C        SYMMETRIC NON-DIAGONAL MATRIX.
C        ------------------------------
         IF ( IFLAG.EQ.1 .OR. IFLAG.EQ.11 ) THEN
C           ---------------------------------------
C                                    T
C           COMPUTE THE PRODUCT C = L B, I.E. C=UB.
C           ---------------------------------------
            DO 400 JS = 1, NSUPER
               FSTVAR = XSUPER(JS)
               LSTVAR = XSUPER(JS+1) - 1
               ISTRT  = XLINDX(JS)
               ISTOP  = XLINDX(JS+1) - 1
               IPLNZ  = XLNZ(JS) - 1
               DO 300 J = FSTVAR, LSTVAR
                  IPLNZ = IPLNZ + 1
                  IF ( J.LE.NEQ ) THEN
                     C(J) = C(J) + LNZ(IPLNZ)*B(J)
                     DO 200 I = ISTRT+1, ISTOP
                        IPLNZ = IPLNZ + 1
                        II    = LINDX(I)
                        IF ( II.LE.NEQ ) THEN
                           C(J) = C(J) + LNZ(IPLNZ)*B(II)
                        ENDIF
  200                CONTINUE
                  ENDIF
                  ISTRT = ISTRT + 1
  300          CONTINUE
  400       CONTINUE

         ELSEIF ( IFLAG.EQ.2 .OR. IFLAG.EQ.22 ) THEN
C           ---------------------------------------
C                                               T
C           COMPUTE THE PRODUCT C = LB, I.E. C=U B.
C           ---------------------------------------
            DO 700 JS = 1, NSUPER
               FSTVAR = XSUPER(JS)
               LSTVAR = XSUPER(JS+1) - 1
               ISTRT  = XLINDX(JS)
               ISTOP  = XLINDX(JS+1) - 1
               IPLNZ  = XLNZ(JS) - 1
               DO 600 J = FSTVAR, LSTVAR
                  IPLNZ = IPLNZ + 1
                  IF ( J.LE.NEQ ) THEN
                     C(J) = C(J) + LNZ(IPLNZ)*B(J)
                     DO 500 I = ISTRT+1, ISTOP
                        IPLNZ = IPLNZ + 1
                        II    = LINDX(I)
                        IF ( II.LE.NEQ ) THEN
                           C(II) = C(II) + LNZ(IPLNZ)*B(J)
                        ENDIF
  500                CONTINUE
                  ENDIF
                  ISTRT = ISTRT + 1
  600          CONTINUE
  700       CONTINUE
         ENDIF
      ELSE
C        ----------------
C        DIAGONAL MATRIX.
C        ----------------
         DO 800 J = 1, NEQ
            C(J) = LNZ(J)*B(J)
  800    CONTINUE
      ENDIF

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRTMQ
      END


      SUBROUTINE   SPRTMR   ( NSUPER, NEQ   , NEQFUL, NOFSUB, IFLAG ,
     $                        PERI2E, XSUPER, XLINDX, XLNZ  , LINDX ,
     $                        LNZ   , B     , C       )

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , NEQFUL, NOFSUB, IFLAG
      INTEGER           PERI2E(NEQFUL)        , XSUPER(NSUPER+1)      ,
     $                  XLINDX(NSUPER+1)      , XLNZ(NSUPER+1)        ,
     $                  LINDX(NOFSUB)
      DOUBLE PRECISION  LNZ(*)                ,
     $                  B(NEQ)                , C(NEQ)

C     --------------------------------------------------
C --- TO PRE-MULTIPLY A COMPRESSED MATRIX WITH A VECTOR.
C     --------------------------------------------------
C     Same as SPRTMQ, but account for the internal permutation PERI2E.

      INTEGER           FSTVAR, I     , II    , IPLNZ , ISTRT ,
     $                  ISTOP , J     , JJ    , JS    , LSTVAR

      DOUBLE PRECISION  ZERO
      PARAMETER       ( ZERO = 0.0D0 )

C     -------------
C     INITIALIZE C.
C     -------------
      DO 100 I = 1, NEQ
         C(I) = ZERO
  100 CONTINUE
C
      IF ( IFLAG.EQ.1 .OR. IFLAG.EQ.11 ) THEN
C        ---------------------------------------
C                                 T
C        COMPUTE THE PRODUCT C = L B, I.E. C=UB.
C        ---------------------------------------
         DO 400 JS = 1, NSUPER
            FSTVAR = XSUPER(JS)
            LSTVAR = XSUPER(JS+1) - 1
            ISTRT  = XLINDX(JS)
            ISTOP  = XLINDX(JS+1) - 1
            IPLNZ  = XLNZ(JS) - 1
            DO 300 J = FSTVAR, LSTVAR
               IPLNZ = IPLNZ + 1
               JJ    = PERI2E(J)
               IF ( JJ.LE.NEQ ) THEN
                  C(JJ) = C(JJ) + LNZ(IPLNZ)*B(JJ)
                  DO 200 I = ISTRT+1, ISTOP
                     IPLNZ = IPLNZ + 1
                     II    = PERI2E(LINDX(I))
                     IF ( II.LE.NEQ ) THEN
                        C(JJ) = C(JJ) + LNZ(IPLNZ)*B(II)
                     ENDIF
  200             CONTINUE
               ELSEIF ( ISTOP.GT.ISTRT ) THEN
                  IPLNZ = IPLNZ + ISTOP-ISTRT
               ENDIF
               ISTRT = ISTRT + 1
  300       CONTINUE
  400    CONTINUE

      ELSEIF ( IFLAG.EQ.2 .OR. IFLAG.EQ.22 ) THEN
C        ---------------------------------------
C                                            T
C        COMPUTE THE PRODUCT C = LB, I.E. C=U B.
C        ---------------------------------------
         DO 700 JS = 1, NSUPER
            FSTVAR = XSUPER(JS)
            LSTVAR = XSUPER(JS+1) - 1
            ISTRT  = XLINDX(JS)
            ISTOP  = XLINDX(JS+1) - 1
            IPLNZ  = XLNZ(JS) - 1
            DO 600 J = FSTVAR, LSTVAR
               IPLNZ = IPLNZ + 1
               JJ    = PERI2E(J)
               IF ( JJ.LE.NEQ ) THEN
                  C(JJ) = C(JJ) + LNZ(IPLNZ)*B(JJ)
                  DO 500 I = ISTRT+1, ISTOP
                     IPLNZ = IPLNZ + 1
                     II    = PERI2E(LINDX(I))
                     IF ( II.LE.NEQ ) THEN
                        C(II) = C(II) + LNZ(IPLNZ)*B(JJ)
                     ENDIF
  500             CONTINUE
               ELSEIF ( ISTOP.GT.ISTRT ) THEN
                  IPLNZ = IPLNZ + ISTOP-ISTRT
               ENDIF
               ISTRT = ISTRT + 1
  600       CONTINUE
  700    CONTINUE

      ENDIF

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRTMR
      END
