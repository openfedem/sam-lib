C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRPRM   ( A     , B     , C     , WA    , MSPAR ,
     $                        MTREES, MSIFA , M     , N     , KSA   ,
     $                        IFLAG , LPU   , IERR                    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           M     , N     , IFLAG , KSA   , LPU   , IERR
      INTEGER           MSPAR(*)             ,
     $                  MTREES(*)             , MSIFA(*)
      DOUBLE PRECISION  A(*)                  , WA(M)                 ,
     $                  B(M,N)                , C(M,N)

C ======================================================================
C  S A M  library routine :  SPRPRM                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To compute C = aAB or C = C + aAB, where A is a compressed
C     matrix, B and C are M,N matrices and a=+-1.
C     N O T E: Since WA is used B and C may be the same array.
*              If N=1 and IFLAG=88 C is not referenced and the result
*              is returned in WA.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
*                 Jun. 01, 1999 (acd)
*                 Modified in SPRPRQ in order to account for evaluation
*                 of internal stiffness only.
C                 Feb. 25, 2003 (kmo)
C                 Get pointers from MSPAR instead of recalculating them.
C                 Added internal permutation array PERI2E.
C                 Added SPRPRR and use that instead of SPRPRQ when
C                 an internal permutation is given.
C                 Replaced the out-zeroing of external components in
C                 SPRPRQ (and now SPRPRR) by equivalent calls to
C                 SPRZR0 and SPRZR1 in SPRPRM instead.
C                 Aug. 27, 2005 (kmo)
C                 Removed calls to SPRZR0 and SPRZR1 again.
C                 Instead it is assumed that the external equations
C                 are ordered last, and that M = # internal equations
C                 when MSPAR(57) = 1.
C                 Dec. 08, 2005 (kmo)
C                 Corrected handling of IFLAG = +/-2
C                 This error was discovered by Tore Holmaas in Nov. 2005
C                 but the correction here is slightly different and more
C                 straigth forward (also involving SPRPRQ/SPRPRR) than
C                 suggested by Holmaas.
C
C     A - DOUBLE PRECISION(*)
C      Entry : A symmetric matrix where the lower triangle is stored
C              in compressed format.
C      Exit  : Unchanged, except when MSPAR(57)=1.
C     B - DOUBLE PRECISION(M,N)
C      Entry : The rectangular matrix to be premultipied by A.
C      Exit  : Unchanged.
C     C - DOUBLE PRECISION(M,N)
C      Entry : Irrelevant if IFLAG=ABS(1).
C              If IFLAG=ABS(2) it must contain appropriate coefficients.
C      Exit  : Updated with or assigned to the matrix product AB.
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
C      Entry : The number of columns in B and C.
C              It is required that N .GT. 0.
C      Exit  : Unchanged.
C     KSA - INTEGER
C      Entry : .GT. 0 : A holds a compressed matrix on sparse format.
C              .LE. 0 : A is diagonal and held as a vector.
C      Exit  : Unchanged.
C     IFLAG - INTEGER
C      Entry : +-1: C <- +-AB.
C              +-2: C <- C +- AB.
C              +88: WA <- AB (only for N=1).
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
C      Exit  : Not defined except that last column of matrix vector
C              is held. Also used to avoid C when N=1 and IFLAG=88.
C
C     Procedures
C     ----------
C
C     SPRPRQ
C     SPRPRR
C     SPRER1
C
C     Intrinsic
C     ---------
C
C     ABS
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
     $                  PERI2E, XLINDX, XLNZ  , XSUPER
      EXTERNAL          SPRPRQ, SPRPRR, SPRER1
      INTRINSIC         ABS

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRPRM, IFLAG =',IFLAG
#endif

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1).LT.0 ) THEN
         MSPAR(1) = -9
         IERR = -1
         CALL SPRER1 ( 30, 'SPRPRM', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         IERR = 0
      ENDIF

      AVOIDC = IFLAG.EQ.88 .AND. N.EQ.1

      IF (.NOT.AVOIDC .AND. ABS(IFLAG).NE.1 .AND. ABS(IFLAG).NE.2) THEN
         MSPAR(1) = -9
         IERR = -2
         CALL SPRER1 ( 25, 'SPRPRM', IFLAG, N, 0, LPU, IERR )
         RETURN
      ELSEIF ( KSA.LE.0 .AND. MSPAR(58).GE.1 ) THEN
         MSPAR(1) = -9
         IERR = -3
         CALL SPRER1 ( 26, 'SPRPRM', KSA, 0, 0, LPU, IERR )
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
      PERI2E = 1
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

C     ==================================================================

C     Use only the first M (internal) equations when M < NEQ
      IF ( M.LT.NEQ ) THEN
         IF ( MSPAR(57).EQ.1 ) THEN
            NEQ = M
         ELSE
            MSPAR(1) = -9
            IERR = -4
            CALL SPRER1 ( 24, 'SPRPRM', 0, NEQ, M, LPU, IERR )
            RETURN
         ENDIF
      ENDIF

C     --------------------------------------
C     COMPUTE C BY N CALLS TO SPRPRQ/SPRPRR.
C     --------------------------------------
      DO 200 IN = 1, N

         IF ( MSPAR(58).GE.1 ) THEN
            CALL SPRPRR ( NSUPER, NEQ, NOFSUB, IFLAG, MSIFA(PERI2E),
     $                    MTREES(XSUPER), MTREES(XLINDX), MSIFA(XLNZ),
     $                    MSIFA(LINDX), A, B(1,IN), WA )
         ELSE
            CALL SPRPRQ ( NSUPER, NEQ, NOFSUB, KSA, IFLAG,
     $                    MTREES(XSUPER), MTREES(XLINDX), MSIFA(XLNZ),
     $                    MSIFA(LINDX), A, B(1,IN), WA )
         ENDIF

         IF ( AVOIDC ) GOTO 200

         IF ( ABS(IFLAG).EQ.2 ) THEN
            DO 20 JN = 1, M
               C(JN,IN) = C(JN,IN) + WA(JN)
   20       CONTINUE
         ELSE
            DO 100 JN = 1, M
               C(JN,IN) = WA(JN)
  100       CONTINUE
         ENDIF

  200 CONTINUE

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRPRM'
#endif
      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRPRM
      END


      SUBROUTINE SPRPRQ ( NSUPER, NEQ   , NOFSUB, KSA   , IFLAG ,
     $                    XSUPER, XLINDX, XLNZ  , LINDX , LNZ   ,
     $                    B     , C       )

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , NOFSUB, KSA   , IFLAG
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  XLNZ(NSUPER+1)        , LINDX(NOFSUB)
      DOUBLE PRECISION  LNZ(*)                ,
     $                  B(NEQ)                , C(NEQ)

C     --------------------------------------------------
C --- TO PRE-MULTIPLY A COMPRESSED MATRIX WITH A VECTOR.
C     --------------------------------------------------

      INTEGER           FSTVAR, I     , II    , IPLNZ , ISTRT ,
     $                  ISTOP , J     , JS    , LSTVAR

      DOUBLE PRECISION  ALFA  , ONE   , ZERO
      PARAMETER       ( ONE = 1.0D0, ZERO = 0.0D0 )

      IF ( IFLAG.LT.0 ) THEN
C        ---------
C        SUBTRACT.
C        ---------
         ALFA = -ONE
      ELSE
C        ----
C        ADD.
C        ----
         ALFA = ONE
      ENDIF

      IF ( KSA .GT. 0 ) THEN

         DO 100 I = 1, NEQ
            C(I) = ZERO
  100    CONTINUE

C        ------------------------------
C        SYMMETRIC NON-DIAGONAL MATRIX.
C        ------------------------------
         DO 400 JS = 1, NSUPER
            FSTVAR = XSUPER(JS)
            LSTVAR = XSUPER(JS+1) - 1
            ISTRT  = XLINDX(JS)
            ISTOP  = XLINDX(JS+1) - 1
            IPLNZ  = XLNZ(JS) - 1
            DO 300 J = FSTVAR, LSTVAR
               IPLNZ = IPLNZ + 1
               IF ( J.LE.NEQ ) THEN
                  C(J) = C(J) + ALFA*(LNZ(IPLNZ)*B(J))
                  DO 200 I = ISTRT+1, ISTOP
                     IPLNZ = IPLNZ + 1
                     II    = LINDX(I)
                     IF ( II.LE.NEQ ) THEN
                        C(J)  = C(J)  + ALFA*(LNZ(IPLNZ)*B(II))
                        C(II) = C(II) + ALFA*(LNZ(IPLNZ)*B(J))
                     ENDIF
  200             CONTINUE
               ENDIF
               ISTRT = ISTRT + 1
  300       CONTINUE
  400    CONTINUE

      ELSE

C        ----------------
C        DIAGONAL MATRIX.
C        ----------------
         DO 500 J = 1, NEQ
            C(J)  = ALFA*(LNZ(J)*B(J))
  500    CONTINUE

      ENDIF

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRPRQ
      END


      SUBROUTINE SPRPRR ( NSUPER, NEQ   , NOFSUB, IFLAG , PERI2E,
     $                    XSUPER, XLINDX, XLNZ  , LINDX , LNZ   ,
     $                    B     , C       )

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , NOFSUB, IFLAG
      INTEGER           PERI2E(NEQ)           , XSUPER(NSUPER+1)      ,
     $                  XLINDX(NSUPER+1)      , XLNZ(NSUPER+1)        ,
     $                  LINDX(NOFSUB)
      DOUBLE PRECISION  LNZ(*)                ,
     $                  B(NEQ)                , C(NEQ)

C     --------------------------------------------------
C --- TO PRE-MULTIPLY A COMPRESSED MATRIX WITH A VECTOR.
C     --------------------------------------------------
C     Same as SPRPRQ, but account for the internal permutation PERI2E.

      INTEGER           FSTVAR, I     , II    , IPLNZ , ISTRT ,
     $                  ISTOP , J     , JJ    , JS    , LSTVAR

      DOUBLE PRECISION  ALFA  , ONE   , ZERO
      PARAMETER       ( ONE = 1.0D0, ZERO = 0.0D0 )

      IF ( IFLAG.LT.0 ) THEN
C        ---------
C        SUBTRACT.
C        ---------
         ALFA = -ONE
      ELSE
C        ----
C        ADD.
C        ----
         ALFA = ONE
      ENDIF

      DO 100 I = 1, NEQ
         C(I) = ZERO
  100 CONTINUE

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
               C(JJ) = C(JJ) + ALFA*(LNZ(IPLNZ)*B(JJ))
               DO 200 I = ISTRT+1, ISTOP
                  IPLNZ = IPLNZ + 1
                  II    = PERI2E(LINDX(I))
                  IF ( II.LE.NEQ ) THEN
                     C(JJ) = C(JJ) + ALFA*(LNZ(IPLNZ)*B(II))
                     C(II) = C(II) + ALFA*(LNZ(IPLNZ)*B(JJ))
                  ENDIF
  200          CONTINUE
            ELSEIF ( ISTOP.GT.ISTRT ) THEN
               IPLNZ = IPLNZ + ISTOP-ISTRT
            ENDIF
            ISTRT = ISTRT + 1
  300    CONTINUE
  400 CONTINUE

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRPRR
      END
