C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRFS1   ( MSPAR , MTREES, MSIFA , SM    , B     ,
     $                        LDB   , NRHS  , IWORK , RWORK , LPU   ,
     $                        IERR                                    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           LDB   , NRHS  , LPU   , IERR
      INTEGER           MSPAR(*)              , MTREES(*)             ,
     $                  MSIFA(*)              , IWORK(*)
      DOUBLE PRECISION  SM(*)                 , B(LDB,NRHS)           ,
     $                  RWORK(*)

C ======================================================================
C  S A M  library routine :  SPRFS1                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To forward and diagonal solve the NRHS right-hand sides stored
C     column by column in B.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
*                 Jun. 01, 1999 (acd)
*                 Modification in order to account for evaluation of
*                 internal stiffness only.
C                 Feb. 25, 2003 (kmo)
C                 Suppress error messages if LPU < 0.
C                 Get pointers from MSPAR instead of recalculating them.
C                 Added internal permutation array PERI2E.
C                 Replaced SPRFS5 by SPRFS3 and SPRFS4.
C                 Replaced SPRBS9 by SPRZR0 and SPRZR1.
C                 Aug. 14, 2003 (kmo)
C                 Added MSPAR(57) in calls to SPRFS3 and SPRFS4.
C                 Aug. 27, 2005 (kmo)
C                 Removed calls to SPRZR0 and SPRZR1 again.
C                 Instead NEQ is reset to LDB when MSPAR(57) = 1
C                 assuming the internal equations are ordered first.
C
C     MSPAR - INTEGER(*)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed, except for MSPAR(1) that is set to 7 in
C              case of sucessful exit, and to -7 else.
C     MTREES - INTEGER(NTREES)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     SM - DOUBLE PRECISION(NSM)
C      Entry : As output from SPRFC1.
C      Exit  : Not changed.
C     B - DOUBLE PRECISION(LDB,NRHS)
C      Entry : The right-hand sides stored column by column.
C      Exit  : The forward and diagonal solved right-hand sides
C              stored column by column.
C     LDB - INTEGER
C      Entry : The leading dimension of B.
C      Exit  : Not changed.
C     NRHS - INTEGER
C      Entry : The number of right-hand sides stored in B.
C      Exit  : Not changed.
C     LPU - INTEGER
C      Entry : Output unit.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call,
C              = -1 if error from SPRFC1 is not cleared.
C              = -2 if MTREES, MSIFA, SM, IWORK and/or RWORK
C                   are not great enough.
C
C     Working arrays
C     --------------
C
C     IWORK - INTEGER(NIWORK)
C      Entry : Not defined. Partitioned to internal work array,
C              see below.
C      Exit  : Need not be saved.
C     RWORK - DOUBLE PRECISION(NRWORK)
C      Entry : Not defined.
C              See below.
C      Exit  : Need not be saved.
C
C     Procedures
C     ----------
C
C     SPRFS3
C     SPRFS4
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

      INTEGER
     $  NIWORK, NMSIFA, NRWORK, NSM   , NTREES

      INTEGER
     $  DIAG  , ELSEQ , FIRST , IWUSED, IRHS  , LINDX , LNZ   , MAXSUP,
     $  MFRONT, MLUSED, MTUSED, NB    , NELACT, NEQ   , NOFSUB, NSUPER,
     $  NZEROL, PERI2E, RWUSED, SMUSED, STACK , SUPSUP, XBLOCK, XELSEQ,
     $  XLINDX, XSUPER

      EXTERNAL          SPRFS3, SPRFS4, SPRER1

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRFS1'
#endif

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -7
         IERR = -1
         CALL SPRER1 ( 30, 'SPRFS1', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         MSPAR(1) = 7
         IERR = 0
      ENDIF

      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NB     = MSPAR(12)
      MAXSUP = MSPAR(13)
      MFRONT = MSPAR(14)
      NOFSUB = MSPAR(15)
      NZEROL = MSPAR(16)
      NELACT = MSPAR(19)

      NMSIFA = MSPAR( 3)
      NTREES = MSPAR( 4)
      NSM    = NEQ + NZEROL
      NIWORK = MAXSUP
      NRWORK = MFRONT
C     ----------------------------------
C     DIVIDE IWORK AND CHECK ITS LENGTH.
C     ----------------------------------
      FIRST  = 1
      IWUSED = FIRST  + MAXSUP - 1

C     --------------------------------------------
C     SET POINTERS TO MTREES AND CHECK ITS LENGTH.
C     --------------------------------------------
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)

* === SUPSUP-acd
      SUPSUP = XLINDX + NSUPER + 1
      ELSEQ  = SUPSUP + NSUPER
      XELSEQ = ELSEQ  + NELACT
      XBLOCK = XELSEQ + NB     + 1
      MTUSED = XBLOCK + NB

C     -------------------------------------------
C     SET POINTERS TO MSIFA AND CHECK ITS LENGTH.
C     -------------------------------------------
      PERI2E = 1
      LINDX  = MSPAR(47)
      MLUSED = LINDX  + NOFSUB - 1

C     ----------------------------------------
C     SET POINTERS TO SM AND CHECK ITS LENGTH.
C     ----------------------------------------
      DIAG   = 1
      LNZ    = DIAG   + NEQ
      SMUSED = LNZ    + NZEROL - 1

C     -------------------------------------------
C     SET POINTERS TO RWORK AND CHECK ITS LENGTH.
C     -------------------------------------------
      STACK  = 1
      RWUSED = STACK  + MFRONT - 1

C     Solve for only the first LDB (internal) equations when LDB < NEQ
      IF ( MSPAR(57).EQ.1 .AND. LDB.LT.NEQ ) NEQ = LDB

      IF ( IWUSED .GT. NIWORK .OR. MTUSED .GT. NTREES .OR.
     $     MLUSED .GT. NMSIFA .OR. SMUSED .GT. NSM    .OR.
     $     RWUSED .GT. NRWORK .OR. LDB    .LT. NEQ         ) THEN
         MSPAR(1) = -2
         IF ( LPU .GE. 0 ) THEN
            WRITE(LPU,'(/A)') ' *** ERROR from SPRFS1'
            IF ( IWUSED .GT. NIWORK ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'IWORK too small',NIWORK,IWUSED
            ENDIF
            IF ( MTUSED .GT. NTREES ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'MTREES too small',NTREES,MTUSED
            ENDIF
            IF ( MLUSED .GT. NMSIFA ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'MSIFA too small',NMSIFA,MLUSED
            ENDIF
            IF ( SMUSED .GT. NSM ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'SM too small',NSM,SMUSED
            ENDIF
            IF ( RWUSED .GT. NRWORK ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'RWORK too small',NRWORK,RWUSED
            ENDIF
            IF ( LDB .LT. NEQ ) THEN
               WRITE(LPU,'(/A,I8,I8/)') 'LDB < NEQ', LDB, NEQ
            ENDIF
         ENDIF
         IERR = -2
         CALL SPRER1 ( 37, 'SPRFS1', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

C     ==================================================================

      DO 100 IRHS = 1, NRHS

C        ----------------------------------------------------
C        FORWARD AND DIAGONAL SOLVE IRHS BY A CALL TO SPRFS3.
C        ----------------------------------------------------
         IF ( MSPAR(58).GE.1 ) THEN
            CALL SPRFS4 (NSUPER, NEQ, MAXSUP, MFRONT, NOFSUB, NZEROL,
     $                   MSPAR(57), MTREES(XSUPER), MTREES(XLINDX),
     $                   MTREES(SUPSUP), MSIFA(LINDX), MSIFA(PERI2E),
     $                   B(1,IRHS), SM(LNZ), IWORK(FIRST), RWORK(STACK))
         ELSE
            CALL SPRFS3 (NSUPER, NEQ, MAXSUP, MFRONT, NOFSUB, NZEROL,
     $                   MSPAR(57), MTREES(XSUPER), MTREES(XLINDX),
     $                   MTREES(SUPSUP), MSIFA(LINDX),
     $                   B(1,IRHS), SM(LNZ), IWORK(FIRST), RWORK(STACK))
         ENDIF

  100 CONTINUE

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRFS1'
#endif
      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS1
      END


      SUBROUTINE   SPRFS3   ( NSUPER, NEQ   , MAXSUP, MFRONT, NOFSUB,
     $                        NZEROL, LPAR57, XSUPER, XLINDX, SUPSUP,
     $                        LINDX , GLBRHS, LNZ   , FIRST , LOCRHS  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , MAXSUP, MFRONT, NOFSUB, NZEROL,
     $                  LPAR57
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  SUPSUP(NSUPER)        , LINDX(NOFSUB)         ,
     $                  FIRST(MAXSUP)

      DOUBLE PRECISION  GLBRHS(NEQ)           , LNZ(NZEROL)           ,
     $                  LOCRHS(MFRONT)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFS3
C
C --- TO PERFORM AN IN-CORE MULTIFRONT FORWARD SOLUTION STEP, THIS IS
C     THE TOP LEVEL DRIVER.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Renamed from SPRFS5 to SPRFS3.
C                 Removed unused arguments IRHS and INFO.
C                 Replaced SPRFS6 by SPRFS5.
C                 AUG. 14, 2003 (KMO)
C                 Added input parameter LPAR57.
C                 AUG. 27, 2005 (KMO)
C                 Added NEQ in calls to SPRFB1 and SPRFB2.
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     MAXSUP : INTEGER
C              LARGEST SUPERNODE.
C     MFRONT : INTEGER
C              LARGEST FRONT MATRIX.
C     NOFSUB : INTEGER
C              LENGTH OF LINDX.
C     NZEROL : INTEGER
C              LENGTH OF LNZ.
C     LPAR57 : INTEGER
C              IF (LPAR57.LE.0) DO A FULL SOLVE EVEN IF THE EQUATION
C              SYSTEM IS PARTITIONED INTO INTERNAL AND RETAINED DOFS.
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION.
C     (XLINDX,
C     LINDX) : INTEGER(NSUPER+1), INTEGER(NOFSUB)
C              THE NON-ZERO STRUCTURE OF L.
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              LDL' FACTOR MATRICES.
C
C     INPUT/OUTPUT
C     ------------
C
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              ON ENTRY : RIGHT HAND SIDE.
C              ON EXIT  : FORWARD SOLVED VECTOR AND DIAGONAL SOLVED
C                         VECTOR.
C
C     WORKING ARRAYS
C     --------------
C
C     FIRST  : INTEGER(MAXSUP)
C              STORES POINTERS INTO LOCAL SUBMATRICES FOR EACH SOLVE.
C     LOCRHS : DOUBLE PRECISION(MFRONT)
C              STORES ON A DENSE FORMAT ALL "SUBMATRIX" RHS'S.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     SPRFB1 : GATHER GLOBAL CONTRIBUTIONS TO LOCAL STORAGE.
C     SPRFS5 : SUBMATRIX FORWARD SOLVE.
C     SPRFB2 : SCATTER LOCAL CONTRIBUTIONS TO GLOBAL STORAGE.
C
C     ------------------------------------------------------------------

      INTEGER           IPLNZ , JS    , JSLEN , JSSTRT, SUPSZE, SUPVAR
      EXTERNAL          SPRFB1, SPRFS5, SPRFB2

C     ==================================================================

      IPLNZ = 1

C     ==================================================================

      DO 100 JS = 1, NSUPER

         SUPVAR = XSUPER(JS+1) - XSUPER(JS)
         JSSTRT = XLINDX(JS)
         JSLEN  = XLINDX(JS+1) - JSSTRT
         SUPSZE = (JSLEN*(JSLEN+1) - (JSLEN-SUPVAR)*(JSLEN-SUPVAR+1))/2

         IF ( SUPSUP(JS).LE.0 .OR. LPAR57.LE.0 ) THEN

C           -----------------------------------------------------
C           PERFORM SUPVAR STEPS OF ELIMINATION since this is not
*           a retained supernode.
C           -----------------------------------------------------
            CALL SPRFB1 ( NEQ, JSLEN, LINDX(JSSTRT), GLBRHS, LOCRHS )
            CALL SPRFS5 ( SUPVAR, JSLEN, LNZ(IPLNZ), LOCRHS, FIRST )
            CALL SPRFB2 ( NEQ, JSLEN, LINDX(JSSTRT), GLBRHS, LOCRHS )

         ENDIF

         IPLNZ  = IPLNZ  + SUPSZE

  100 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS3
      END


      SUBROUTINE   SPRFS4   ( NSUPER, NEQ   , MAXSUP, MFRONT, NOFSUB,
     $                        NZEROL, LPAR57, XSUPER, XLINDX, SUPSUP,
     $                        LINDX , PERI2E, GLBRHS, LNZ   , FIRST ,
     $                        LOCRHS  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSUPER, NEQ   , MAXSUP, MFRONT, NOFSUB, NZEROL,
     $                  LPAR57
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  SUPSUP(NSUPER)        , LINDX(NOFSUB)         ,
     $                  FIRST(MAXSUP)         , PERI2E(NEQ)

      DOUBLE PRECISION  GLBRHS(NEQ)           , LNZ(NZEROL)           ,
     $                  LOCRHS(MFRONT)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFS4
C
C --- TO PERFORM AN IN-CORE MULTIFRONT FORWARD SOLUTION STEP, THIS IS
C     THE TOP LEVEL DRIVER. Identical to SPRFS3 except for the internal
C     permutation array PERI2E.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Made SPRFS4 from SPRFS3 by adding PERI2E, and
C                 using SPRFB3 and SPRFB4 instead of SPRFB1 and SPRFB2.
C                 AUG. 14, 2003 (KMO)
C                 Added input parameter LPAR57.
C                 AUG. 27, 2005 (KMO)
C                 Added NEQ in calls to SPRFB3 and SPRFB4.
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     MAXSUP : INTEGER
C              LARGEST SUPERNODE.
C     MFRONT : INTEGER
C              LARGEST FRONT MATRIX.
C     NOFSUB : INTEGER
C              LENGTH OF LINDX.
C     NZEROL : INTEGER
C              LENGTH OF LNZ.
C     LPAR57 : INTEGER
C              IF (LPAR57.LE.0) DO A FULL SOLVE EVEN IF THE EQUATION
C              SYSTEM IS PARTITIONED INTO INTERNAL AND RETAINED DOFS.
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION.
C     (XLINDX,
C     LINDX) : INTEGER(NSUPER+1), INTEGER(NOFSUB)
C              THE NON-ZERO STRUCTURE OF L.
C     PERI2E : INTEGER(NEQ)
C              Permutation from internal to external equation order.
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              LDL' FACTOR MATRICES.
C
C     INPUT/OUTPUT
C     ------------
C
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              ON ENTRY : RIGHT HAND SIDE.
C              ON EXIT  : FORWARD SOLVED VECTOR AND DIAGONAL SOLVED
C                         VECTOR.
C
C     WORKING ARRAYS
C     --------------
C
C     FIRST  : INTEGER(MAXSUP)
C              STORES POINTERS INTO LOCAL SUBMATRICES FOR EACH SOLVE.
C     LOCRHS : DOUBLE PRECISION(MFRONT)
C              STORES ON A DENSE FORMAT ALL "SUBMATRIX" RHS'S.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     SPRFB3 : GATHER GLOBAL CONTRIBUTIONS TO LOCAL STORAGE.
C     SPRFS5 : SUBMATRIX FORWARD SOLVE.
C     SPRFB4 : SCATTER LOCAL CONTRIBUTIONS TO GLOBAL STORAGE.
C
C     ------------------------------------------------------------------

      INTEGER           IPLNZ , JS    , JSLEN , JSSTRT, SUPSZE, SUPVAR
      EXTERNAL          SPRFB3, SPRFS5, SPRFB4

C     ==================================================================

      IPLNZ = 1

C     ==================================================================

      DO 100 JS = 1, NSUPER

         SUPVAR = XSUPER(JS+1) - XSUPER(JS)
         JSSTRT = XLINDX(JS)
         JSLEN  = XLINDX(JS+1) - JSSTRT
         SUPSZE = (JSLEN*(JSLEN+1) - (JSLEN-SUPVAR)*(JSLEN-SUPVAR+1))/2

         IF ( SUPSUP(JS).LE.0 .OR. LPAR57.LE.0 ) THEN

C           -----------------------------------------------------
C           PERFORM SUPVAR STEPS OF ELIMINATION since this is not
*           a retained supernode.
C           -----------------------------------------------------
            CALL SPRFB3 ( NEQ, JSLEN, LINDX(JSSTRT), PERI2E,
     $                    GLBRHS, LOCRHS )
            CALL SPRFS5 ( SUPVAR, JSLEN, LNZ(IPLNZ), LOCRHS, FIRST )
            CALL SPRFB4 ( NEQ, JSLEN, LINDX(JSSTRT), PERI2E,
     $                    GLBRHS, LOCRHS )

         ENDIF

         IPLNZ  = IPLNZ  + SUPSZE

  100 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS4
      END


      SUBROUTINE   SPRFS5   ( SUPVAR, LENVAR, LNZ   , LOCRHS, FIRST   )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           SUPVAR, LENVAR
      INTEGER           FIRST(SUPVAR)
      DOUBLE PRECISION  LNZ(*)                , LOCRHS(LENVAR)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFS5
C
C --- TO PERFORM A PARTIAL FORWARD SOLVE FOR THE MULTIFRONTAL METHOD.
C     AND TO PERFORM A DIAGONAL SOLVE ON THE FULLY SUMMED PART.
C
C     CREATED   : MAY  25, 1992 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Renamed from SPRFS6 to SPRFS5.
C                 Replaced SPRFS4 by SPRFS6.
C
C
C     ON ENTRY
C     --------
C
C     SUPVAR : INTEGER
C              NUMBER OF VARIABLES TO SOLVE FOR IN CURRENT STEP.
C     LENVAR : INTEGER
C              ORDER OF THE ASSOCIATED FRONT MATRIX IN CURRENT STEP.
C     LNZ    : DOUBLE PRECISION(LENVAR*SUPVAR-SUPVAR*(SUPVAR-1)/2)
C              THE CURRENT FACTOR SUBMATRIX.
C
C     ON EXIT
C     -------
C
C     LOCRHS : DOUBLE PRECISION(LENVAR)
C              FORWARD SOLVED LOCAL VECTOR, AND DIAGONAL SOLVED
C              FOR THE FULLY SUMMED VARIABLES.
C
C     WORKING ARRAYS
C     --------------
C
C     FIRST  : INTEGER(SUPVAR)
C              POINTERS INTO CHOLESKY FACTOR SUBMATRIX.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     SPRFS6
C
C     ------------------------------------------------------------------

      INTEGER           I     , IPLNZ , J     , NROWS , OFFSET
      DOUBLE PRECISION  RHSJ
      EXTERNAL          SPRFS6

C     ==================================================================

      OFFSET(I,J) = (I*(J-1))-((J-1)*(J-2)/2)

C     ==================================================================

C     --------------------
C     SOLVE TRIANGLE PART.
C     --------------------
      DO 200 J = 1, SUPVAR
         IPLNZ  = 1 + OFFSET(LENVAR,J)
         RHSJ   = LOCRHS(J)
         DO 100 I = J + 1, SUPVAR
            IPLNZ  = IPLNZ  + 1
            LOCRHS(I) = LOCRHS(I) - RHSJ*LNZ(IPLNZ)
  100    CONTINUE
         FIRST(J) = IPLNZ  + 1
  200 CONTINUE
      NROWS  = LENVAR - SUPVAR
      IF ( NROWS .GT. 0 ) THEN

C        ------------------------
C        REDUCE RHS BY RECTANGLE.
C        ------------------------
         CALL SPRFS6 ( NROWS, SUPVAR, LNZ, FIRST, LOCRHS,
     $                 LOCRHS(SUPVAR+1) )
      ENDIF

C     -------------------------------------
C     DIAGONAL SOLVE THE FULLY SUMMED PART.
C     -------------------------------------
      DO 300 J = 1, SUPVAR
         IPLNZ  = 1 + OFFSET(LENVAR,J)
         LOCRHS(J) = LOCRHS(J) / LNZ(IPLNZ)
  300 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS5
      END


      SUBROUTINE   SPRFB2   ( NEQ   , JSLEN , LINDX , GLBRHS, LOCRHS  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEQ   , JSLEN
      INTEGER           LINDX(JSLEN)
      DOUBLE PRECISION  GLBRHS(NEQ)           , LOCRHS(JSLEN)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFB2
C
C --- TO SCATTER LOCAL RIGHT HAND SIDE CONTRIBUTIONS TO GLOBAL VECTOR.
C
C     CREATED   : SEP. 14, 1993 (ACD)
C     REVISIONS : AUG. 27, 2005 (KMO)
C                 Added NEQ as input parameter.
C
C
C     ON ENTRY
C     --------
C
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     JSLEN  : INTEGER
C              THE CURRENT SUPERNODE TO SCATTER CONTRIBUTIONS TO.
C     LINDX  : INTEGER(JSLEN)
C              THE NODE BLOCK NON-ZERO STRUCTURE OF THE CURRENT COLUMN.
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              GLOBAL VECTOR.
C
C     ON EXIT
C     -------
C
C     LOCRHS : DOUBLE PRECISION(*)
C              LOCAL VECTOR
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     NONE
C
C     ------------------------------------------------------------------
      INTEGER           JS    , ROWIND
C     ==================================================================
      DO 100 JS = 1, JSLEN
         ROWIND = LINDX(JS)
         IF ( ROWIND.LE.NEQ ) THEN
            GLBRHS(ROWIND) = LOCRHS(JS)
         ENDIF
  100 CONTINUE
C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRFB2
      END


      SUBROUTINE   SPRFB4   ( NEQ   , JSLEN , LINDX , PERI2E,
     $                        GLBRHS, LOCRHS  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEQ   , JSLEN
      INTEGER           LINDX(JSLEN)          , PERI2E(NEQ)
      DOUBLE PRECISION  GLBRHS(NEQ)           , LOCRHS(JSLEN)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFB4
C
C --- TO SCATTER LOCAL RIGHT HAND SIDE CONTRIBUTIONS TO GLOBAL VECTOR.
C     Same as SPRFB2, but account for permutation array PERI2E.
C
C     CREATED   : SEP. 14, 1993 (ACD)
C     REVISIONS : FEB. 24, 2003 (KMO)
C                 Made SPRFB4 from SPRFB2 by adding PERI2E.
C                 AUG. 27, 2005 (KMO)
C                 Added NEQ as input parameter.
C
C     ON ENTRY
C     --------
C
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     JSLEN  : INTEGER
C              THE CURRENT SUPERNODE TO SCATTER CONTRIBUTIONS TO.
C     LINDX  : INTEGER(JSLEN)
C              THE NODE BLOCK NON-ZERO STRUCTURE OF THE CURRENT COLUMN.
C     PERI2E : INTEGER(NEQ)
C              Permutation from internal to external equation ordering.
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              GLOBAL VECTOR.
C
C     ON EXIT
C     -------
C
C     LOCRHS : DOUBLE PRECISION(*)
C              LOCAL VECTOR
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     NONE
C
C     ------------------------------------------------------------------
      INTEGER           JS    , ROWIND
C     ==================================================================
      DO 100 JS = 1, JSLEN
         ROWIND = PERI2E(LINDX(JS))
         IF ( ROWIND.LE.NEQ ) THEN
            GLBRHS(ROWIND) = LOCRHS(JS)
         ENDIF
  100 CONTINUE
C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRFB4
      END


      SUBROUTINE   SPRFB1   ( NEQ   , JSLEN , LINDX , GLBRHS, LOCRHS  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEQ   , JSLEN
      INTEGER           LINDX(JSLEN)
      DOUBLE PRECISION  GLBRHS(NEQ)           , LOCRHS(JSLEN)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFB1
C
C --- TO GATHER GLOBAL RIGHT HAND SIDE CONTRIBUTIONS FROM GLOBAL VECTOR
C     AND LOAD IT INTO A LOCAL VECTOR.
C
C     CREATED   : SEP. 14, 1993 (ACD)
C     REVISIONS : AUG. 27, 2005 (KMO)
C                 Added NEQ as input parameter.
C
C
C     ON ENTRY
C     --------
C
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     JSLEN  : INTEGER
C              THE LENGTH OF THE CURRENT SUPERNODE.
C     LINDX  : INTEGER(JSLEN)
C              THE NON-ZERO STRUCTURE OF THE CURRENT COLUMN.
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              GLOBAL VECTOR.
C
C     ON EXIT
C     -------
C
C     LOCRHS : DOUBLE PRECISION(*)
C              LOCAL VECTOR
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     NONE
C
C     ------------------------------------------------------------------
      INTEGER           JS    , ROWIND
      DOUBLE PRECISION  ZERO
      PARAMETER       ( ZERO = 0.0D0 )
C     ==================================================================
      DO 100 JS = 1, JSLEN
         ROWIND = LINDX(JS)
         IF ( ROWIND.LE.NEQ ) THEN
            LOCRHS(JS) = GLBRHS(ROWIND)
         ELSE
            LOCRHS(JS) = ZERO
         ENDIF
 100  CONTINUE
C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRFB1
      END


      SUBROUTINE   SPRFB3   ( NEQ   , JSLEN , LINDX , PERI2E,
     $                        GLBRHS, LOCRHS  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEQ   , JSLEN
      INTEGER           LINDX(JSLEN)          , PERI2E(NEQ)
      DOUBLE PRECISION  GLBRHS(NEQ)           , LOCRHS(JSLEN)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRFB3
C
C --- TO GATHER GLOBAL RIGHT HAND SIDE CONTRIBUTIONS FROM GLOBAL VECTOR
C     AND LOAD IT INTO A LOCAL VECTOR.
C     Same as SPRFB2, but account for permutation array PERI2E.
C
C     CREATED   : SEP. 14, 1993 (ACD)
C     REVISIONS : FEB. 24, 2003 (KMO)
C                 Made SPRFB3 from SPRFB1 by adding PERI2E.
C                 AUG. 27, 2005 (KMO)
C                 Added NEQ as input parameter.
C
C
C     ON ENTRY
C     --------
C
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     JSLEN  : INTEGER
C              THE LENGTH OF THE CURRENT SUPERNODE.
C     LINDX  : INTEGER(JSLEN)
C              THE NON-ZERO STRUCTURE OF THE CURRENT COLUMN.
C     PERI2E : INTEGER(NEQ)
C              Permutation from internal to external equation ordering.
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              GLOBAL VECTOR.
C
C     ON EXIT
C     -------
C
C     LOCRHS : DOUBLE PRECISION(*)
C              LOCAL VECTOR
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     NONE
C
C     ------------------------------------------------------------------
      INTEGER           JS    , ROWIND
      DOUBLE PRECISION  ZERO
      PARAMETER       ( ZERO = 0.0D0 )
C     ==================================================================
      DO 100 JS = 1, JSLEN
         ROWIND = PERI2E(LINDX(JS))
         IF ( ROWIND.LE.NEQ ) THEN
            LOCRHS(JS) = GLBRHS(ROWIND)
         ELSE
            LOCRHS(JS) = ZERO
         ENDIF
 100  CONTINUE
C     ==================================================================
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRFB3
      END


      SUBROUTINE   SPRFS6   ( ROWSZE, COLSZE, LNZ   , FIRST , X     ,
     $                        Y                                       )
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER           ROWSZE, COLSZE,
     $                  FIRST(COLSZE)
      DOUBLE PRECISION  LNZ(*)                , X(COLSZE)             ,
     $                  Y(ROWSZE)
C     ------------------------------------------------------------------
C
C     SPRFS6  : GENERALIZED SAXPY, UNROLLED TO DEPTH 8.
C
C     PURPOSE
C     -------
C
C     SPRFS6  : TO PERFORM MULTIPLE SAXPY'S, THAT IS, TO ADD TO A REAL
C               VECTOR Y THE MATRIX-VECTOR PRODUCT -AX, WHERE A AND
C               X ARE ALSO REAL AND THE COLUMNS OF A BEGIN IN UNRELATED
C               LOCATIONS IN MAIN MEMORY, IE IN LNZ.
C
C
C     CREATED   : DEC. 02, 1990 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Renamed from SPRFS4 to SPRFS6.
C
C
C     ON ENTRY
C     --------
C
C     ROWSZE : INTEGER
C              NUMBER OF ROWS IN MATRIX A
C     COLSZE : INTEGER
C              NUMBER OF COLUMNS IN MATRIX A
C     LNZ    : REAL(*)
C              A REAL MATRIX WHICH CONTAINS THE ROWSZE BY COLSZE
C              MATRIX A AS A SUBMATRIX
C     FIRST  : INTEGER(COLSZE)
C              FIRST(J) POINTS TO THE BEGINNING OF THE JTH
C              COLUMN OF A
C     X      : REAL(COLSZE)
C              A REAL VECTOR CONTAINING THE MULTIPLIERS
C
C     ON EXIT
C     -------
C
C     Y      : REAL(ROWSZE)
C              A REAL VECTOR WHOSE VALUE ON  INPUT IS RE-
C              PLACED BY Y + AX.
C
C     ------------------------------------------------------------------
      INTEGER           LFTOVR, K0    , K1    , K2    , K3    , K4    ,
     $                  K5    , K6    , K7    , I     , JMIN  , JMAX  ,
     $                  J
C     ------------------------------------------------------------------
      IF ( (ROWSZE .LE. 0) .OR. (COLSZE .LE. 0) )  RETURN
      LFTOVR = COLSZE - ( COLSZE / 8 ) * 8
      IF ( LFTOVR. GT. 0 ) THEN
      IF ( LFTOVR .LE. 3 ) THEN
         IF ( LFTOVR .LE. 1 ) THEN
            IF ( LFTOVR .EQ. 1 ) THEN
               K1 = FIRST(1)
CCDIR$ IVDEP
               DO 10  I = 1, ROWSZE
                  Y(I) = Y(I)
     $                 - X(1) * LNZ(K1)
                  K1 = K1 + 1
  10           CONTINUE
            ENDIF
         ELSE
            IF ( LFTOVR .EQ. 2 ) THEN
               K1 = FIRST(1)
               K2 = FIRST(2)
CCDIR$ IVDEP
               DO 20  I = 1, ROWSZE
                  Y(I) = Y(I)
     $                 - X(1) * LNZ(K1) - X(2) * LNZ(K2)
                  K1 = K1 + 1
                  K2 = K2 + 1
  20           CONTINUE
            ELSE
               K1 = FIRST(1)
               K2 = FIRST(2)
               K3 = FIRST(3)
CCDIR$ IVDEP
               DO 30  I = 1, ROWSZE
                  Y(I) = Y(I)
     $                 - X(1) * LNZ(K1) - X(2) * LNZ(K2)
     $                 - X(3) * LNZ(K3)
                  K1 = K1 + 1
                  K2 = K2 + 1
                  K3 = K3 + 1
  30           CONTINUE
            ENDIF
         ENDIF
      ELSE
         IF ( LFTOVR .LE. 5 )  THEN
            IF ( LFTOVR .EQ. 4 )  THEN
               K1 = FIRST(1)
               K2 = FIRST(2)
               K3 = FIRST(3)
               K4 = FIRST(4)
CCDIR$ IVDEP
               DO 40  I = 1, ROWSZE
                  Y(I) = Y(I)
     $                 - X(1) * LNZ(K1) - X(2) * LNZ(K2)
     $                 - X(3) * LNZ(K3) - X(4) * LNZ(K4)
                  K1 = K1 + 1
                  K2 = K2 + 1
                  K3 = K3 + 1
                  K4 = K4 + 1
  40           CONTINUE
            ELSE
               K1 = FIRST(1)
               K2 = FIRST(2)
               K3 = FIRST(3)
               K4 = FIRST(4)
               K5 = FIRST(5)
CCDIR$ IVDEP
               DO  50  I = 1, ROWSZE
                  Y(I) = Y(I)
     $                 - X(1) * LNZ(K1) - X(2) * LNZ(K2)
     $                 - X(3) * LNZ(K3) - X(4) * LNZ(K4)
     $                 - X(5) * LNZ(K5)
                  K1 = K1 + 1
                  K2 = K2 + 1
                  K3 = K3 + 1
                  K4 = K4 + 1
                  K5 = K5 + 1
  50           CONTINUE
            ENDIF
         ELSE
            IF ( LFTOVR .EQ. 6 )  THEN
               K1 = FIRST(1)
               K2 = FIRST(2)
               K3 = FIRST(3)
               K4 = FIRST(4)
               K5 = FIRST(5)
               K6 = FIRST(6)
CCDIR$ IVDEP
               DO 60  I = 1, ROWSZE
                  Y(I) = Y(I)
     $                 - X(1) * LNZ(K1) - X(2) * LNZ(K2)
     $                 - X(3) * LNZ(K3) - X(4) * LNZ(K4)
     $                 - X(5) * LNZ(K5) - X(6) * LNZ(K6)
                  K1 = K1 + 1
                  K2 = K2 + 1
                  K3 = K3 + 1
                  K4 = K4 + 1
                  K5 = K5 + 1
                  K6 = K6 + 1
  60           CONTINUE
            ELSE
               K1 = FIRST(1)
               K2 = FIRST(2)
               K3 = FIRST(3)
               K4 = FIRST(4)
               K5 = FIRST(5)
               K6 = FIRST(6)
               K7 = FIRST(7)
CCDIR$ IVDEP
               DO 70  I = 1, ROWSZE
                   Y(I) = Y(I)
     $                  - X(1) * LNZ(K1) - X(2) * LNZ(K2)
     $                  - X(3) * LNZ(K3) - X(4) * LNZ(K4)
     $                  - X(5) * LNZ(K5) - X(6) * LNZ(K6)
     $                  - X(7) * LNZ(K7)
                   K1 = K1 + 1
                   K2 = K2 + 1
                   K3 = K3 + 1
                   K4 = K4 + 1
                   K5 = K5 + 1
                   K6 = K6 + 1
                   K7 = K7 + 1
  70           CONTINUE
            ENDIF
         ENDIF
      ENDIF
      ENDIF
      JMIN = LFTOVR + 8
      JMAX = COLSZE
      IF ( JMAX .GE. JMIN )  THEN
         DO 100  J = JMIN, JMAX, 8
            K0 = FIRST(J)
            K1 = FIRST(J-1)
            K2 = FIRST(J-2)
            K3 = FIRST(J-3)
            K4 = FIRST(J-4)
            K5 = FIRST(J-5)
            K6 = FIRST(J-6)
            K7 = FIRST(J-7)
CCDIR$ IVDEP
            DO 90  I = 1, ROWSZE
               Y(I) = Y(I)
     $              - X(J  ) * LNZ(K0) - X(J-1) * LNZ(K1)
     $              - X(J-2) * LNZ(K2) - X(J-3) * LNZ(K3)
     $              - X(J-4) * LNZ(K4) - X(J-5) * LNZ(K5)
     $              - X(J-6) * LNZ(K6) - X(J-7) * LNZ(K7)
               K0 = K0 + 1
               K1 = K1 + 1
               K2 = K2 + 1
               K3 = K3 + 1
               K4 = K4 + 1
               K5 = K5 + 1
               K6 = K6 + 1
               K7 = K7 + 1
  90        CONTINUE
 100     CONTINUE
      ENDIF
      RETURN
C     ------------------------------------------------------------------
C     END OF MODULE SPRFS6
      END
