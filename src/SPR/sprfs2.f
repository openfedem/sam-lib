C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRFS2   ( MSPAR , MTREES, MSIFA , SM    , B     ,
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
C  S A M  library routine :  SPRFS2                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To forward solve the NRHS right-hand sides stored
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
C                 Replaced SPRFS7 by SPRFS7 and SPRFS8.
C                 Replaced SPRBS9 by SPRZR0 and SPRZR1.
C                 Aug. 14, 2003 (kmo)
C                 Added MSPAR(57) in calls to SPRFS7 and SPRFS8.
C                 Sep. 22, 2004 (kmo)
C                 Changed the MSPAR(1) value from 7 to 8.
C                 Aug. 27, 2005 (kmo)
C                 Removed calls to SPRZR0 and SPRZR1 again.
C                 Instead NEQ is reset to LDB when MSPAR(57) = 1
C                 assuming the internal equations are ordered first.
C
C     MSPAR - INTEGER(*)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed, except for MSPAR(1) that is set to 8
C              in case of successful exit, and to -8 otherwise.
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
C      Exit  : The forward solved right-hand sides
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
C              = -1 if error from SPRFC2 is not cleared.
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
C     SPRFS7
C     SPRFS8
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

      EXTERNAL          SPRFS7, SPRFS8, SPRER1

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'ENTERING SPRFS2'
#endif

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -8
         IERR = -1
         CALL SPRER1 ( 30, 'SPRFS2', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         MSPAR(1) = 8
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
            WRITE(LPU,'(/A)') ' *** ERROR from SPRFS2'
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
         CALL SPRER1 ( 37, 'SPRFS2', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

C     ==================================================================

      DO 100 IRHS = 1, NRHS

C        ---------------------------------------
C        FORWARD SOLVE IRHS BY A CALL TO SPRFS7.
C        ---------------------------------------
         IF ( MSPAR(58).GE.1 ) THEN
            CALL SPRFS8 (NSUPER, NEQ, MAXSUP, MFRONT, NOFSUB, NZEROL,
     $                   MSPAR(57), MTREES(XSUPER), MTREES(XLINDX),
     $                   MTREES(SUPSUP), MSIFA(LINDX), MSIFA(PERI2E),
     $                   B(1,IRHS), SM(LNZ), IWORK(FIRST), RWORK(STACK))
         ELSE
            CALL SPRFS7 (NSUPER, NEQ, MAXSUP, MFRONT, NOFSUB, NZEROL,
     $                   MSPAR(57), MTREES(XSUPER), MTREES(XLINDX),
     $                   MTREES(SUPSUP), MSIFA(LINDX),
     $                   B(1,IRHS), SM(LNZ), IWORK(FIRST), RWORK(STACK))
         ENDIF

  100 CONTINUE

C     ==================================================================

#if FT_DEBUG > 1
      IF (LPU.GE.0) WRITE(LPU,*) 'LEAVING  SPRFS2'
#endif
      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS2
      END


      SUBROUTINE   SPRFS7   ( NSUPER, NEQ   , MAXSUP, MFRONT, NOFSUB,
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
C     SPRFS7
C
C --- TO PERFORM AN IN-CORE MULTIFRONT FORWARD SOLUTION STEP, THIS IS
C     THE TOP LEVEL DRIVER.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Removed unused arguments IRHS and INFO.
C                 Replaced SPRFS8 by SPRFS9.
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
C              THE NON-ZERO STRUCTURE OF C.
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              CC' FACTOR MATRICES.
C
C     INPUT/OUTPUT
C     ------------
C
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              ON ENTRY : RIGHT HAND SIDE.
C              ON EXIT  : FORWARD SOLVED VECTOR.
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
C     SPRFS9 : SUBMATRIX FORWARD SOLVE.
C     SPRFB2 : SCATTER LOCAL CONTRIBUTIONS TO GLOBAL STORAGE.
C
C     ------------------------------------------------------------------

      INTEGER           IPLNZ , JS    , JSLEN , JSSTRT, SUPSZE, SUPVAR
      EXTERNAL          SPRFB1, SPRFS9, SPRFB2

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
            CALL SPRFS9 ( SUPVAR, JSLEN, LNZ(IPLNZ), LOCRHS, FIRST )
            CALL SPRFB2 ( NEQ, JSLEN, LINDX(JSSTRT), GLBRHS, LOCRHS )

         ENDIF

         IPLNZ  = IPLNZ  + SUPSZE

  100 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS7
      END


      SUBROUTINE   SPRFS8   ( NSUPER, NEQ   , MAXSUP, MFRONT, NOFSUB,
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
C     SPRFS8
C
C --- TO PERFORM AN IN-CORE MULTIFRONT FORWARD SOLUTION STEP, THIS IS
C     THE TOP LEVEL DRIVER. Identical to SPRFS7 except for the internal
C     permutation array PERI2E.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Made SPRFS8 from SPRFS7 by adding PERI2E, and
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
C              THE NON-ZERO STRUCTURE OF C.
C     PERI2E : INTEGER(NEQ)
C              Permutation from internal to external equation order.
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              CC' FACTOR MATRICES.
C
C     INPUT/OUTPUT
C     ------------
C
C     GLBRHS : DOUBLE PRECISION(NEQ)
C              ON ENTRY : RIGHT HAND SIDE.
C              ON EXIT  : FORWARD SOLVED VECTOR.
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
C     SPRFS9 : SUBMATRIX FORWARD SOLVE.
C     SPRFB4 : SCATTER LOCAL CONTRIBUTIONS TO GLOBAL STORAGE.
C
C     ------------------------------------------------------------------

      INTEGER           IPLNZ , JS    , JSLEN , JSSTRT, SUPSZE, SUPVAR
      EXTERNAL          SPRFB1, SPRFS9, SPRFB2

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
            CALL SPRFS9 ( SUPVAR, JSLEN, LNZ(IPLNZ), LOCRHS, FIRST )
            CALL SPRFB4 ( NEQ, JSLEN, LINDX(JSSTRT), PERI2E,
     $                    GLBRHS, LOCRHS )

         ENDIF

         IPLNZ  = IPLNZ  + SUPSZE

  100 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS8
      END


      SUBROUTINE   SPRFS9   ( SUPVAR, LENVAR, LNZ   , LOCRHS, FIRST   )

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
C     SPRFS9
C
C --- TO PERFORM A PARTIAL FORWARD SOLVE FOR THE MULTIFRONTAL METHOD.
C
C     CREATED   : MAY  25, 1992 (ACD)
C     REVISIONS : FEB. 25, 2003 (KMO)
C                 Renamed from SPRFS8 to SPRFS9.
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
C              FORWARD SOLVED LOCAL VECTOR
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
         LOCRHS(J) = LOCRHS(J) / LNZ(IPLNZ)
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

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRFS9
      END
