C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRDG1   ( MSPAR , MTREES, MSIFA , IPDIAG, LPU   ,
     $                        IERR                                    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           LPU , IERR
      INTEGER           MSPAR(*)              , MTREES(*)             ,
     $                  MSIFA(*)              , IPDIAG(*)


C ======================================================================
C  S A M  library routine :  SPRDG1                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To extract pointers to diagonal entries of L. The pointers are
C     incremented by NEQ to be used directly on SM.
C
C     Created   : Dec. 27, 1996 (acd)
C     Revisions : Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 Feb. 24, 2003 (kmo)
C                 Removed some superflouos error messages, let SPRER1
C                 handle all error message output.
C                 Get pointers from MSPAR instead of recalculating them.
C                 Added permutation array PERI2E.
C                 Replaced SPRDG2 by SPRDG2 and SPRDG3.
C                 Sep. 22, 2004 (kmo)
C                 MSPAR(1) is not changed if successful exit.
C
C     MSPAR - INTEGER(*)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed, except for MSPAR(1) that is set to -66
C              in case of unsuccessful exit.
C     MTREES - INTEGER(NTREES)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     IPDIAG - INTEGER(NEQ)
C      Entry : Not defined
C      Exit  : Holds the pointers into L, incremented by NEQ.
C     LPU - INTEGER
C      Entry : Output unit.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call,
C              = -1 if error from SPRTRS is not cleared.
C              = -2 if MTREES, MSIFA, and/ or IPDIAG
C                   are not great enough.
C
C     Working arrays
C     --------------
C
C     None
C
C     Procedures
C     ----------
C
C     SPRDG2
C     SPRDG3
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

      INTEGER           LINDX , MLUSED, NEQ   , NMSIFA, NSUPER, NTREES,
     $                  PERI2E, XLINDX, XLNZ  , XSUPER

      EXTERNAL          SPRDG2, SPRDG3, SPRER1

C     ==================================================================
C
#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRDG1'
#endif
C
C     -------------------
C     GET INFO FROM MSPAR
C     -------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -66
         IERR = -1
         CALL SPRER1 ( 30, 'SPRDG1', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         IERR = 0
      ENDIF

      NMSIFA = MSPAR( 3)
      NTREES = MSPAR( 4)
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)

C     -------------------------------------------
C     SET POINTERS TO MTREES AND CHECK ITS LENGTH
C     -------------------------------------------
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      MLUSED = XLINDX + NSUPER

      IF ( MLUSED .GT. NTREES ) THEN
         MSPAR(1) = -2
         IERR = -2
         CALL SPRER1 ( 38, 'SPRDG1', MLUSED, NTREES, 0, LPU, IERR )
         RETURN
      ENDIF

C     ------------------------------------------
C     SET POINTERS TO MSIFA AND CHECK ITS LENGTH
C     ------------------------------------------
      PERI2E = 1
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)
      MLUSED = XLNZ   + NSUPER + 1

      IF ( MLUSED .GT. NMSIFA ) THEN
         MSPAR(1) = -2
         IERR = -2
         CALL SPRER1 ( 39, 'SPRDG1', MLUSED, NMSIFA, 0, LPU, IERR )
         RETURN
      ENDIF

C     ==================================================================

C     --------------------------------------
C     Get the pointers to the diagonal of L.
C     --------------------------------------
      IF ( MSPAR(58) .LT. 1 ) THEN
         CALL SPRDG2 ( NSUPER, NEQ, MTREES(XSUPER), MTREES(XLINDX),
     $                 MSIFA(XLNZ), IPDIAG )
      ELSE
         CALL SPRDG3 ( NSUPER, NEQ, MTREES(XSUPER), MTREES(XLINDX),
     $                 MSIFA(XLNZ), MSIFA(PERI2E), IPDIAG )
      ENDIF
C     ==================================================================
C
#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRDG1'
#endif
      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRDG1
      END


      SUBROUTINE   SPRDG2   ( NSUPER, NEQ   , XSUPER, XLINDX, XLNZ  ,
     $                        IPDIAG  )

      IMPLICIT NONE

C     ------------------------------------------------------------------

      INTEGER           NSUPER, NEQ
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  XLNZ(NSUPER+1)        , IPDIAG(NEQ)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRDG2
C
C --- To extract pointers to diagonal entries of L.
C
C     CREATED   : Dec. 27, 1996 (ACD)
C     REVISIONS : Feb. 24, 2003 (KMO)
C                 Removed unused error flag.
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS, Length of IPDIAG.
C     XSUPER : INTEGER(NSUPER+1)
C              THE FUNDAMENTAL SUPERNODE PARTITION.
C     (XLINDX,
C     LINDX) : INTEGER(NSUPER+1) (Only XLINDX is used here.)
C              THE NODE BLOCK NON-ZERO STRUCTURE OF L.
C     XLNZ   : INTEGER(NSUPER+1)
C              POINTER TO START OF SUBMATRICES IN LNZ FOR EACH
C              SUPERNODE.
C
C     ON EXIT
C     -------
C
C     IPDIAG : INTEGER(NEQ)
C              The pointers to diagonal elements of L.
C
C     PROCEDURES
C     ----------
C
C     None
C
C     INTRINSIC
C     ---------
C
C     None
C
C     ------------------------------------------------------------------

      INTEGER           CHDLEN, IPJS  , J     , JS

C     ==================================================================

C     --------------------------------
C     Get the diagonal positions of L.
C     --------------------------------
      DO 200 JS = 1, NSUPER
         CHDLEN = XLINDX(JS+1) - XLINDX(JS)
         IPJS = XLNZ(JS)
         DO 100 J = XSUPER(JS), XSUPER(JS+1)-1
C                             Note: INCREMENTED BY NEQ.
            IPDIAG(J) = IPJS + NEQ
            IPJS   = IPJS + CHDLEN
            CHDLEN = CHDLEN - 1
  100    CONTINUE
  200 CONTINUE

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRDG2
      END


      SUBROUTINE   SPRDG3   ( NSUPER, NEQ   , XSUPER, XLINDX, XLNZ  ,
     $                        PERI2E, IPDIAG  )

      IMPLICIT NONE

C     ------------------------------------------------------------------

      INTEGER           NSUPER, NEQ
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  XLNZ(NSUPER+1)        , PERI2E(NEQ)           ,
     $                  IPDIAG(NEQ)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRDG3
C
C --- To extract pointers to diagonal entries of L.
C
C     CREATED   : Dec. 27, 1996 (ACD)
C     REVISIONS : Feb. 24, 2003 (KMO)
C                 Made SPRDG3 from SPRDG2 by adding PERI2E.
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS, Length of IPDIAG.
C     XSUPER : INTEGER(NSUPER+1)
C              THE FUNDAMENTAL SUPERNODE PARTITION.
C     (XLINDX,
C     LINDX) : INTEGER(NSUPER+1) (Only XLINDX is used here.)
C              THE NODE BLOCK NON-ZERO STRUCTURE OF L.
C     XLNZ   : INTEGER(NSUPER+1)
C              POINTER TO START OF SUBMATRICES IN LNZ FOR EACH
C              SUPERNODE.
C     PERI2E : Permutation from internal to external equation ordering
C
C     ON EXIT
C     -------
C
C     IPDIAG : INTEGER(NEQ)
C              The pointers to diagonal elements of L.
C
C     PROCEDURES
C     ----------
C
C     None
C
C     INTRINSIC
C     ---------
C
C     None
C
C     ------------------------------------------------------------------

      INTEGER           CHDLEN, IPJS  , J     , JS

C     ==================================================================

C     --------------------------------
C     Get the diagonal positions of L.
C     --------------------------------
      DO 200 JS = 1, NSUPER
         CHDLEN = XLINDX(JS+1) - XLINDX(JS)
         IPJS = XLNZ(JS)
         DO 100 J = XSUPER(JS), XSUPER(JS+1)-1
            IPDIAG(PERI2E(J)) = IPJS + NEQ
            IPJS   = IPJS + CHDLEN
            CHDLEN = CHDLEN - 1
  100    CONTINUE
  200 CONTINUE

C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRDG3
      END
