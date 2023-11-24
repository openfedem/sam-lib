C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRTRS   ( MSPAR , MSICA , MTREES, MEQN  ,
     $                        NSPAR , NMSICA, NTREES, NDOF  ,
     $                        NIWORK, IWORK , RINFO , LPU   , IERR )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSPAR , NMSICA, NTREES, NDOF  , NIWORK,
     $                  LPU   , IERR
      INTEGER           MSPAR(NSPAR)          , MSICA(NMSICA)         ,
     $                  MTREES(NTREES)        , MEQN(NDOF)            ,
     $                  IWORK(NIWORK)
      DOUBLE PRECISION  RINFO(4)

C ======================================================================
C  S A M  library routine :  SPRTRS                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To analyse the sparsity pattern of a matrix that is represented
C     as finite-element cliques. It uses the information stored in
C     MSICA to generate the tree structures associated with the current
C     node/variable order, and most likely it updates MEQN with a mathe-
C     matically equivalent order. The tree structures are stored in
C     MTREES, of length NTREES.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1999 (acd)
C                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 Feb. 27, 2003 (kmo)
C                 Increased MSPAR(3) by 2*NEQ if MSPAR(51) is set.
C                 May. 16, 2003 (kmo)
C                 Reduced the size of NODEL by SPRNOD.
C                 Jan. 13, 2004 (kmo)
C                 Added check for integer overflow in NZEROL.
C                 Aug. 14, 2021 (kmo)
C                 MSPAR(*) -> MSPAR(NSPAR), NSPAR is new input argument
C                 MSICA(*) -> MSICA(NMSICA), NMSICA is new argument
C                 MTREES(*) -> MTREES(NTREES), NTREES is new argument
C                 IWORK(*) -> IWORK(NIWORK), NIWORK is new argument
C
C
C     MSPAR - Integer(NSPAR)
C      Entry : As on exit from SPRSAS/SPRRNM.
C      Exit  : Finalized, see below.
C     MSICA - INTEGER(NMSICA)
C      Entry : As on exit from SPRSAS.
C      Exit  : Not changed.
C     MTREES - INTEGER(NTREES)
C      Entry : Not defined.
C      Exit  : The tree structures.
C     MEQN - INTEGER(NDOF)
C      Entry : As on exit from SPRPRP/SPRRNM.
C      Exit  : Most likely, positive integers are changed to conform
C              with the a mathematically equivalent order of the
C              SPR-nodes that is computed during the tree search.
C     NSPAR - INTEGER
C      Entry : Length of MSPAR.
C      Exit  : Not changed.
C     NMSICA - INTEGER
C      Entry : Length of MSICA.
C      Exit  : Not changed.
C     NTREES - INTEGER
C      Entry : Length of MTREES.
C      Exit  : Not changed.
C     NDOF - INTEGER
C      Entry : The number of variables in the finite-element mesh.
C      Exit  : Not changed.
C     NIWORK - INTEGER
C      Entry : Length of IWORK.
C      Exit  : Not changed.
C     LPU - INTEGER
C      Entry : Output unit for error messages.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call,
C              = -1 if error from SPRSAS/SPRRNM is not cleared.
C              = -2 if MSICA and/or MTREES and/or IWORK are not great
C                   enough.
C              = -3 if error from SPRTR4.
C              = -4 if integer overflow due to NZEROL > MAXINT_p.
*              = -6 Error from SPRTR3 concerning error in the separator
*                   nodes for superelement technique.
C
C     Working arrays
C     --------------
C
C     IWORK - INTEGER(NIWORK)
C      Entry : Not defined.
C              The array is partitioned by the routine for use in SPRMMD
C              as follows:
C              INVP   : On exit from SPRMMD, the inverse permutation.
C              XFLNOD : Pointer to FLNOD.
C              FLNOD  : List of nodes in element headed by its length.
C              XFNODL : Pointer to FNODL.
C              FNODL  : List of elements at nodes headed by its length.
C              HEAD   : Head of doubly linked degree structure.
C              QSIZE  : Size of supernodes.
C              ELMLNK : List of new elements in each elimination step.
C              FLAG   : Marker vector.
C              STOREL : List of nodes that may be indist or outmatched.
C              MARKER : Marker vector.
C      Exit  : Need not be saved.
C
C     Procedures
C     ----------
C
C     SPRTR1
C     SPRTR2
C     SPRTR3
C     SPRTR4
C     SPRTR5
C     SPRCPY
C     SPREQN
C     SPRCN1
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
C     FINAL VERSION OF MSPAR:
C
C     ------------------------------------------------------------------

      INTEGER
     $  CHILD , COLCNT, DOFRET, DPSTCK, ELNOD , ELSEQ , ELSIZE, ERROR ,
     $  FDESC , INVP  , IWUSED, LCOEFA, LEVEL , MARKER, MASK  , MAXFRT,
     $  MAXSUP, MSUSED, MTUSED, NB    , NCHILD, NEL   , NELACT, NEQ   ,
     $  NEWINV, NODEL , NODES , NOFSUB, NSIBL , NSUPER, NWPERM, NXSIBL,
     $  NMSIFA, NXSUPR, NZSTCK, PARENT, PERM  , PRVLF , PRVNBR, I     ,
     $  ROWCNT, SEPARR, SEPSZE, SET   , SIBL  , SPRCON, SPRNOD, SUINVP,
     $  SUPERM, SUPLEN, SUPRET, SUPSUP, SUPSZE, SZEREP, SZERET, VARCNT,
     $  VARWGH, WEIGHT, XBLOCK, XELNOD, XELSEQ, XNODEL, XNOD  , XNODES,
     $  XSIBL , XSUPER

      include 'integer_kind.inc'

      INTEGER(kind=i8)  NZEROL

      INTEGER           OPTION
      PARAMETER       ( OPTION = 1 )
      EXTERNAL          SPRTR1, SPRTR2, SPRTR3, SPRTR4, SPRTR5, SPRCPY,
     $                  SPREQN, SPRCN1, SPRER1

C     ==================================================================

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -3
         IERR = -1
         CALL SPRER1 ( 30, 'SPRTRS', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         IERR = 0
      ENDIF

      NEL    = MSPAR( 5)
      SPRNOD = MSPAR( 6)
      SPRCON = MSPAR( 7)
      NEQ    = MSPAR( 8)
      SEPSZE = MSPAR(40)

C     -------------------------------------------
C     SET POINTERS TO MSICA AND CHECK ITS LENGTH.
C     -------------------------------------------
      XELNOD = 1
      XNODEL = XELNOD + NEL    + 1
      XNODES = XNODEL + SPRNOD + 1
      NODES  = XNODES + SPRNOD + 1
      PERM   = NODES  + NEQ
      ELNOD  = PERM   + SPRNOD
      NODEL  = ELNOD  + SPRCON
      SEPARR = NODEL  + SPRCON
      MSUSED = SEPARR + SPRNOD - 1

      IF ( MSUSED .GT. NMSICA ) THEN
         MSPAR(1) = -2
         IERR = -2
         CALL SPRER1 ( 31, 'SPRTRS', MSUSED,  NMSICA, 0, LPU, IERR )
         RETURN
      ENDIF

C     -------------------------------------------
C     SET POINTERS TO IWORK AND CHECK ITS LENGTH.
C     -------------------------------------------
      INVP   = 1
      PARENT = INVP   + SPRNOD
      CHILD  = PARENT + SPRNOD
      SIBL   = CHILD  + SPRNOD
      NEWINV = SIBL   + SPRNOD
      IWUSED = NEWINV + SPRNOD - 1
      IF ( IWUSED .GT. NIWORK ) THEN
         MSPAR(1) = -3
         IERR = -2
         CALL SPRER1 ( 31, 'SPRTRS', IWUSED, NIWORK, 0, LPU, IERR )
         RETURN
      ENDIF

C     ==================================================================

      IF ( MSPAR(1) .EQ. 1 ) THEN

C        -----------------------------------------
C        SPRRNM IS NOT USED, COMPUTE XNODEL,NODEL.
C        -----------------------------------------
         CALL SPRCN1 ( MSICA(XELNOD), MSICA(ELNOD), MSICA(XNODEL),
     $                 MSICA(NODEL), NEL, SPRNOD, SPRCON, IWORK(INVP))
      ENDIF

C     ----------------------------------------------
C     SET MSPAR( 1) = 3 TO FLAG THAT SPRTRS IS USED.
C     ----------------------------------------------
      MSPAR( 1) = 3

C     ==================================================================

C     -------------------------------------------
C     COMPUTE AND POSTORDER THE ELIMINATION TREE.
C     -------------------------------------------

      CALL SPRTR1 ( SPRNOD, NEL, SPRCON, MSICA(XELNOD), MSICA(ELNOD),
     $              MSICA(XNODEL), MSICA(NODEL), MSICA(PERM),
     $              IWORK(INVP), IWORK(PARENT), IWORK(CHILD),
     $              IWORK(SIBL), IWORK(NEWINV) )

C     ==================================================================

C     -------------------------------------------------------------
C     COMPUTE THE COLUMN LENGTHS REFERENCED TO NODES AND VARIABLES.
C     -------------------------------------------------------------
      IWUSED = CHILD
      XNOD   = IWUSED
      COLCNT = XNOD   + SPRNOD + 1
      VARCNT = COLCNT + SPRNOD
      ROWCNT = VARCNT + SPRNOD
      PRVLF  = ROWCNT + SPRNOD
      PRVNBR = PRVLF  + SPRNOD
      SET    = PRVNBR + SPRNOD
      IWUSED = SET    + SPRNOD - 1

C     ------------------------------------
C     USE MTREES TO HOLD SOME WORK ARRAYS.
C     ------------------------------------
      MARKER = 1
      NCHILD = MARKER + SPRNOD
      WEIGHT = NCHILD + SPRNOD + 1
      VARWGH = WEIGHT + SPRNOD + 1
      FDESC  = VARWGH + SPRNOD + 1
      LEVEL  = FDESC  + SPRNOD + 1
      MTUSED = LEVEL  + SPRNOD + NEL
C                                ^--- USED TO CHECK THE ASKED LENGTH.
C                                ------------------------------------

      IF ( IWUSED .GT. NIWORK .OR. MTUSED .GT. NTREES ) THEN
         MSPAR(1) = -3
         IERR = -2
         IF ( MTUSED .GT. NTREES ) THEN
            CALL SPRER1 ( 33, 'SPRTRS', MTUSED, NTREES, 0, LPU, IERR )
         ENDIF
         IF ( IWUSED .GT. NIWORK ) THEN
            CALL SPRER1 ( 34, 'SPRTRS', IWUSED, NIWORK, 0, LPU, IERR )
         ENDIF
         RETURN
      ENDIF

      CALL SPRTR2 ( SPRNOD, NEL, SPRCON, MSICA(XELNOD), MSICA(ELNOD),
     $              MSICA(XNODEL), MSICA(NODEL), MSICA(PERM),
     $              MSICA(XNODES), IWORK(INVP), IWORK(PARENT),
     $              IWORK(XNOD), IWORK(COLCNT), IWORK(VARCNT),
     $              IWORK(ROWCNT), IWORK(PRVLF), IWORK(PRVNBR),
     $              IWORK(SET), MTREES(MARKER), MTREES(NCHILD),
     $              MTREES(WEIGHT), MTREES(VARWGH), MTREES(FDESC),
     $              MTREES(LEVEL) )

C     ==================================================================

C     --------------------------------------------
C     COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION.
C     --------------------------------------------
      IWUSED = ROWCNT
      XSUPER = IWUSED
      SUPLEN = XSUPER + SPRNOD + 1
      SUPSZE = SUPLEN + SPRNOD

* === SUPSUP-acd
      SUPSUP = SUPSZE + SPRNOD
      IWUSED = SUPSUP + SPRNOD - 1

      IF ( IWUSED .GT. NIWORK ) THEN
         MSPAR(1) = -3
         IERR = -2
         CALL SPRER1 ( 34, 'SPRTRS', IWUSED, NIWORK, 0, LPU, IERR )
         RETURN
      ENDIF

C     --------------------------------------
C     PUT (XSIBL,SIBL) AT THE END OF MTREES.
C     --------------------------------------
      XSIBL  = FDESC
      SIBL   = LEVEL

* === SUPSUP-acd
      CALL SPRTR3 ( SEPSZE, SPRNOD, IWORK(PARENT), IWORK(XNOD),
     $              IWORK(COLCNT), IWORK(VARCNT), MTREES(NCHILD),
     $              NSUPER, IWORK(XSUPER), IWORK(SUPLEN),
     $              IWORK(SUPSZE), IWORK(SUPSUP),MTREES(XSIBL),
     $              MTREES(SIBL), MSICA(PERM), SUPRET, MSICA(SEPARR),
     $              LPU, ERROR )
      IF ( ERROR.NE.0 ) THEN
         MSPAR(1) = -3
         IERR = -6
         WRITE ( LPU, '(//A,I10)' ) '*** Error from SPRTRS', ERROR
         RETURN
      ENDIF

C     -----------------------------------
C     COPY XSUPER, XSIBL, SIBL, SUPLEN TO
C     THEIR FINAL LOCATIONS IN MTREES.
C     -----------------------------------
      MTUSED = 1
      CALL SPRCPY ( NSUPER+1, IWORK(XSUPER), 1, MTREES(MTUSED), 1 )
      XSUPER = MTUSED
      MTUSED = XSUPER + NSUPER + 1
      CALL SPRCPY ( NSUPER+1, MTREES(XSIBL), 1, MTREES(MTUSED), 1 )
      XSIBL  = MTUSED
      MTUSED = XSIBL  + NSUPER + 1
      CALL SPRCPY ( NSUPER+1, MTREES(SIBL), 1, MTREES(MTUSED), 1 )
      SIBL   = MTUSED
      MTUSED = SIBL   + NSUPER + 1
      CALL SPRCPY ( NSUPER, IWORK(SUPLEN), 1, MTREES(MTUSED), 1 )
      SUPLEN = MTUSED

* --- This is important - SUPLEN has length NSUPER+1.
      MTUSED = SUPLEN + NSUPER + 1

* === SUPSUP-acd
      CALL SPRCPY ( NSUPER, IWORK(SUPSUP), 1, MTREES(MTUSED), 1 )
      SUPSUP = MTUSED
      MTUSED = SUPSUP + NSUPER

C     ---------------------------------------------------
C     MOVE SUPSZE TO THE FRONT OF IWORK, JUST AFTER INVP.
C     ---------------------------------------------------
      IWUSED = PARENT
      CALL SPRCPY ( NSUPER, IWORK(SUPSZE), 1, IWORK(IWUSED), 1 )
      SUPSZE = IWUSED
      IWUSED = SUPSZE + NSUPER

C     ------------------------------------
C     ALLOCATE SPACE IN IWORK FOR OPTIMAL
C     POSTORDER OF THE SUPERNODE CHILDREN.
C     ------------------------------------
      SUPERM = IWUSED
      SUINVP = SUPERM + NSUPER
      NWPERM = SUINVP + NSUPER
      NXSUPR = NWPERM + SPRNOD
      IWUSED = NXSUPR + NSUPER

C     ------------------------------------------------
C     USE MTREES TO STORE SOME WORK ARRAYS AT THE END.
C     ------------------------------------------------
      NXSIBL = MTUSED
      NSIBL  = NXSIBL + NSUPER + 1
      MTUSED = NSIBL  + NSUPER

      IF ( IWUSED .GT. NIWORK .OR. MTUSED .GT. NTREES ) THEN
         MSPAR(1) = -3
         IERR = -2
         IF ( MTUSED .GT. NTREES ) THEN
            CALL SPRER1 ( 33, 'SPRTRS', MTUSED, NTREES, 0, LPU, IERR )
         ENDIF
         IF ( IWUSED .GT. NIWORK ) THEN
            CALL SPRER1 ( 34, 'SPRTRS', IWUSED, NIWORK, 0, LPU, IERR )
         ENDIF
         RETURN
      ENDIF

C     ==================================================================

C     -------------------------------
C     COMPUTE THE OPTIMAL POSTORDER,
C     AND UPDATE THE TREE STRUCTURES.
C     -------------------------------

* === SUPSUP-acd
      CALL SPRTR4
     $( SEPSZE,NSUPER,SPRNOD,OPTION,MTREES(XSIBL),MTREES(SIBL),
     $  MTREES(XSUPER),MTREES(SUPLEN),MTREES(SUPSUP),MSICA(PERM),
     $  IWORK(INVP),IWORK(SUPSZE),IWORK(SUPERM),IWORK(SUINVP),
     $  IWORK(NWPERM),IWORK(NXSUPR),MTREES(NXSIBL),MTREES(NSIBL),
     $  ERROR,LPU )
      IF ( ERROR .NE. 0 ) THEN
         IERR = -3
         CALL SPRER1 ( 23, 'SPRTRS', 0, 0, 0, LPU, ERROR )
         RETURN
      ENDIF

C     ==================================================================

C     -----------------------------------
C     ALLOCATE THE FINAL SPACE IN MTREES.
C     -----------------------------------
*     MTUSED = SUPLEN + NSUPER + 1
* === SUPSUP-acd
      MTUSED = SUPSUP + NSUPER
      ELSEQ  = MTUSED
      XELSEQ = ELSEQ  + NEL
      XBLOCK = XELSEQ + NSUPER + 1
      MTUSED = XBLOCK + NSUPER

C     -------------------
C     WORKSPACE IN IWORK.
C     -------------------
      IWUSED = SUPERM
      MASK   = IWUSED
      IWUSED = MASK   + NEL - 1

      IF ( IWUSED .GT. NIWORK .OR. MTUSED .GT. NTREES ) THEN
         MSPAR(1) = -3
         IERR = -2
         IF ( MTUSED .GT. NTREES ) THEN
            CALL SPRER1 ( 33, 'SPRTRS', MTUSED, NTREES, 0, LPU, IERR )
         ENDIF
         IF ( IWUSED .GT. NIWORK ) THEN
            CALL SPRER1 ( 34, 'SPRTRS', IWUSED, NIWORK, 0, LPU, IERR )
         ENDIF
         RETURN
      ENDIF

C     ---------------------------------
C     COMPUTE THE ELEMENT SEQUENCY,
C     THE NUMBER OF BLOCK FACTORIZATION
C     STEPS AND SOME FACTOR STATISTICS.
C     ---------------------------------
      CALL SPRTR5
     $( NSUPER,SPRNOD,NEQ,NEL,SPRCON,MSICA(XELNOD),MSICA(ELNOD),
     $  MSICA(XNODEL),MSICA(NODEL),MSICA(XNODES),
     $  MTREES(XSIBL),MTREES(SIBL),IWORK(SUPSZE),MTREES(SUPSUP),
     $  MSICA(PERM),NB,MAXSUP,MAXFRT,NOFSUB,NZEROL,
     $  NZSTCK,DPSTCK,NELACT,LCOEFA,ELSIZE,DOFRET,SZERET,SZEREP,
     $  MTREES(XELSEQ),MTREES(ELSEQ),MTREES(XBLOCK),
     $  MTREES(XSUPER),MTREES(SUPLEN),RINFO,IWORK(MASK))

      IF ( NZEROL .GT. MAXINT_p ) THEN
         MSPAR(1) = -3
         IERR = -4
         WRITE(LPU,'(//A)') '*** Error from SPRTRS'
         WRITE(LPU,'(4X,A)') 'The model is too large (integer overflow)'
         WRITE(LPU,'(4X,A,I16)') 'NEQ    =',MSPAR(8)
         WRITE(LPU,'(4X,A,I16)') 'NZEROL =',NZEROL
         RETURN
      ENDIF

C     -----------------------------------------------
C     IF NELACT < NEL OR NB < NSUPER COMPRESS MTREES.
C     -----------------------------------------------
      IF ( NELACT .LT. NEL ) THEN
         MTUSED = ELSEQ + NELACT
         DO 100 I = 0, NB
            MTREES(MTUSED) = MTREES(XELSEQ+I)
            MTUSED = MTUSED + 1
  100    CONTINUE
         XELSEQ = ELSEQ + NELACT
      ENDIF
      IF ( NELACT .LT. NEL .OR. NB .LT. NSUPER ) THEN
         MTUSED = XELSEQ + NB + 1
         DO 200 I = 0, NB
            MTREES(MTUSED) = MTREES(XBLOCK+I)
            MTUSED = MTUSED + 1
  200    CONTINUE
         XBLOCK = XELSEQ + NB + 1
      ENDIF

C     ==================================================================

C     ----------------------------
C     UPDATE MEQN ACCORDING TO THE
C     SPR-NODE ORDER HELD IN PERM.
C     ----------------------------
      CALL SPREQN ( SPRNOD, NEQ, NDOF, MSICA(PERM), MSICA(XNODES),
     $              MSICA(NODES), MEQN )

C     ==================================================================

      NMSIFA = 2*NSUPER + NOFSUB + NEQ + 2

C     ---------------
C     FINALIZE MSPAR.
C     ---------------
      MSPAR(10) = 0
      MSPAR(11) = NSUPER
      MSPAR(12) = NB
      MSPAR(13) = MAXSUP
      MSPAR(14) = MAXFRT
      MSPAR(15) = NOFSUB
      MSPAR(16) = NZEROL
      MSPAR(17) = NZSTCK
      MSPAR(18) = DPSTCK
      MSPAR(19) = NELACT
      MSPAR(20) = LCOEFA
      MSPAR(25) = ELSIZE

* --- New quantities to control superelement technique.
      IF ( NSPAR .GT. 50 ) THEN
         MSPAR(53) = SUPRET
         MSPAR(54) = DOFRET
         MSPAR(55) = SZERET
         MSPAR(56) = SZEREP
         IF (MSPAR(51) .EQ. 1) NMSIFA = NMSIFA + 2*NEQ
      ENDIF

C*DAM MSPAR( 3) SOM SATT HER GIR MAXIMUM BEHOV FOR MSIFA
      MSPAR( 3) = NMSIFA
*ACD  MSPAR( 4) = 4*NSUPER + 2*NB + NELACT + 6
      MSPAR( 4) = 5*NSUPER + 2*NB + NELACT + 6
* ---------------- +NSUPER to hold SUPSUP.

C     ==================================================================


      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTRS
      END


      SUBROUTINE   SPRTR1   ( NNODS , NELS  , NELNOD, XELNOD, ELNOD ,
     $                        XNODEL, NODEL , PERM  , INVP  , PARENT,
     $                        CHILD , SIBL  , NEWINV                  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NNODS , NELS  , NELNOD
      INTEGER           XELNOD(NELS+1)        , ELNOD(NELNOD)         ,
     $                  XNODEL(NNODS+1)       , NODEL(NELNOD)         ,
     $                  PERM(NNODS)           , INVP(NNODS)           ,
     $                  PARENT(NNODS)         , CHILD(NNODS)          ,
     $                  SIBL(NNODS)           , NEWINV(NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRTR1
C
C --- TO SET UP INVP AND TO COMPUTE THE ELIMINATION TREE AND TO
C     POSTORDER IT.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF NODES IN THE MESH.
C              NOT ALTERED BY THE ROUTINE.
C     NELS   : INTEGER
C              NUMBER OF ELEMENTS IN THE MESH.
C              NOT ALTERED BY THE ROUTINE.
C     NELNOD : INTEGER
C              SIZE OF ELNOD$NODEL.
C              NOT ALTERED BY THE ROUTINE.
C     XELNOD : INTEGER(NELS+1)
C              POINTER INTO LISTS OF NODES FOR EACH ELEMENT.
C              NOT ALTERED BY THE ROUTINE.
C     ELNOD  : INTEGER(NELNOD)
C              ARRAY THAT HOLDS LISTS OF NODES FOR THE ELEMENTS.
C              NOT ALTERED BY THE ROUTINE.
C     XNODEL : INTEGER(NNODS+1)
C              POINTER INTO LISTS OF ELEMENTS FOR EACH NODE.
C              NOT ALTERED BY THE ROUTINE.
C     NODEL  : INTEGER(NELNOD)
C              ARRAY THAT HOLDS LISTS OF ELEMENTS FOR THE NODES.
C              NOT ALTERED BY THE ROUTINE.
C
C     UPDATED PARAMETERS
C     ------------------
C
C     PERM   : INTEGER(NNODS)
C              ON ENTRY:
C              THE PERMUTATION AFTER THE ORDERING STEP.
C              ON EXIT:
C              THE EQUIVALENT PERMUTATION ORDERING.
C     INVP   : INTEGER(NNODS)
C              ON ENTRY:
C              NOT DEFINED.
C              ON EXIT:
C              THE EQUIVALENT INVERSE PERMUTATION ORDERING.
C
C     ON EXIT
C     -------
C
C     PARENT : INTEGER(NNODS)
C              ON ENTRY:
C              NEEDS NOT BE SET.
C              ON EXIT:
C              THE ELIMINATION TREE REPRESENTATION CONSISTENT WITH
C              THE POSTORDERING.
C
C     WORKING ARRAYS
C     --------------
C
C     CHILD  : INTEGER(NNODS)
C              ON ENTRY:
C              NEEDS NOT BE SET.
C              STORES THE FIRST CHILD, IF ANY, FOR A NODE.
C              MARKS NODES USED IN A STEP OF ELMTR3 TO AVOID USING A
C              NODE MORE THAN ONCE.
C     SIBL   : INTEGER(NNODS)
C              ON ENTRY:
C              NEEDS NOT BE SET.
C              STORES THE REST OF THE, IF ANY,  CHILDREN FOR A NODE.
C     NEWINV : INTEGER(NNODS)
C              ON ENTRY:
C              NEEDS NOT BE SET.
C              STORES THE NEW INVERSE PERMUTATION VECTOR WHICH IS
C              USED TO UPDATE PERM, INVP.
C
C     SUBPROGRAMS
C     -----------
C
C     SPRT11 - COMPUTE INVP.
C     SPRT12 - COMPUTE THE PARENT REPRESENTATION OF THE EL. TREE.
C     SPRT13 - COMPUTE THE FIRST CHILD, SIBLING VECTORS.
C     SPRT14 - POSTORDER THE ELIMINATION TREE.
C     SPRT15 - UPDATE (PERM,INVP) WITH THE POSTORDERING.
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      EXTERNAL          SPRT11, SPRT12, SPRT13, SPRT14 , SPRT15

C     ==================================================================

      CALL  SPRT11 ( NNODS, PERM, INVP )
      CALL  SPRT12 ( NNODS, NELS, NELNOD, XELNOD, ELNOD, XNODEL, NODEL,
     $               PERM, INVP, PARENT, NEWINV, CHILD )
      CALL  SPRT13 ( NNODS, PARENT, CHILD, SIBL )
      CALL  SPRT14 ( NNODS, NNODS, CHILD, SIBL, PARENT, NEWINV, PERM )
      CALL  SPRT15 ( NNODS, NEWINV, PERM, INVP )

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTR1
      END


      SUBROUTINE   SPRTR2   ( NNODS , NELS  , NELNOD, XELNOD, ELNOD ,
     $                        XNODEL, NODEL , PERM  , XNODES, INVP  ,
     $                        PARENT, XNOD  , COLCNT, VARCNT, ROWCNT,
     $                        PRVLF , PRVNBR, SET   , MARKER, NCHILD,
     $                        WEIGHT, VARWGH, FDESC , LEVEL           )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NNODS , NELS  , NELNOD
      INTEGER           XELNOD(NELS+1)        , ELNOD(NELNOD)         ,
     $                  XNODEL(NNODS+1)       , NODEL(NELNOD)         ,
     $                  PERM(NNODS)           , XNODES(NNODS+1)       ,
     $                  INVP(NNODS)           , PARENT(NNODS)         ,
     $                  XNOD(NNODS+1)         , COLCNT(NNODS)         ,
     $                  VARCNT(NNODS)
      INTEGER           ROWCNT(NNODS)         , PRVLF(NNODS)          ,
     $                  PRVNBR(NNODS)         , SET(NNODS)            ,
     $                  MARKER(NNODS)         , NCHILD(0:NNODS)       ,
     $                  WEIGHT(0:NNODS)       , VARWGH(0:NNODS)       ,
     $                  FDESC(0:NNODS)        , LEVEL(0:NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRTR2
C
C --- TO COMPUTE THE LENGTH OF EACH NODE COLUMN IN THE FACTOR MATRIX
C     FOR THE NODE AND VARIABLE VALUE. IN ADDITION IT COMPUTES THE ROW
C     COUNT REFERENCED TO NODES. THIS IS A TRANSLATION OF AN ALGORITHM
C     DUE TO GILBERT, NG AND PEYTON. IT IS MODIFIED TO WORK ON THE
C     IMPLICIT REPRESENTATION OF THE MATRIX, AND TO COMPUTE THE COLUMN
C     LENGTHS WITH REFERENCE TO VARIABLES.
C     SEVERAL ERRORS HAVE BEEN FOUND IN THE ORIGINAL ROUTINE, THEY ARE
C     CORRECTED.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF NODES IN THE MESH.
C              NOT ALTERED BY THE ROUTINE.
C     NELS   : INTEGER
C              NUMBER OF ELEMENTS IN THE MESH.
C              NOT ALTERED BY THE ROUTINE.
C     NELNOD : INTEGER
C              SIZE OF ELNOD$NODEL.
C              NOT ALTERED BY THE ROUTINE.
C     XELNOD : INTEGER(NELS+1)
C              POINTER INTO LISTS OF NODES FOR EACH ELEMENT.
C              NOT ALTERED BY THE ROUTINE.
C     ELNOD  : INTEGER(NELNOD)
C              ARRAY THAT HOLDS LISTS OF NODES FOR THE ELEMENTS.
C              NOT ALTERED BY THE ROUTINE.
C     XNODEL : INTEGER(NNODS+1)
C              POINTER INTO LISTS OF ELEMENTS FOR EACH NODE.
C              NOT ALTERED BY THE ROUTINE.
C     NODEL  : INTEGER(NELNOD)
C              ARRAY THAT HOLDS LISTS OF ELEMENTS FOR THE NODES.
C              NOT ALTERED BY THE ROUTINE.
C     PERM   : INTEGER(NNODS)
C              THE PERMUTATION AFTER THE ORDERING AND POSTORDERING STEP.
C              NOT ALTERED BY THE ROUTINE.
C     XNODES : INTEGER(NNODS+1)
C              THE NODE PARTITION REFERENCED TO INPUT ORDER.
C              NOT ALTERED BY THE ROUTINE.
C     INVP   : INTEGER(NNODS)
C              THE INVERSE PERMUTATION AFTER THE ORDERING AND POST-
C              ORDERING STEP.
C              NOT ALTERED BY THE ROUTINE.
C     PARENT : INTEGER(NNODS)
C              THE PARENT REPRESENTATION OF THE POSTORDERED ELIMINATION
C              TREE.
C              NOT ALTERED BY THE ROUTINE.
C
C     ON EXIT
C     -------
C
C     XNOD   : INTEGER(NNODS+1)
C              THE FUNDAMENTAL NODE PARTITION.
C     COLCNT : INTEGER(NNODS)
C              THE COLUMN LENGTH COUNTED BY NODES, INCLUDING THE
C              DIAGONAL NODE BLOCK.
C     VARCNT : INTEGER(NNODS)
C              THE COLUMN LENGTH COUNTED BY VARIABLES, INCLUDING THE
C              DIAGONAL VARIABLE ELEMENT.
C     ROWCNT : INTEGER(NNODS)
C              THE LENGTH OF EACH ROW COUNTED BY NODES, INCLUDING THE
C              DIAGONAL ELEMENT.
C     NCHILD : INTEGER(0:NNODS)
C              CONTAINS THE NUMBER OF CHILDREN.
C
C     WORKING ARRAYS
C     --------------
C
C     PRVLF  : INTEGER(NNODS)
C              RECORDS THE PREVIOUS OF EACH ROW SUBTREE.
C     PRVNBR : INTEGER(NNODS)
C              RECORDS THE PREVIOUS LOWER NEIGHBOUR OF EACH NODE.
C     SET    : INTEGER(NNODS)
C              MAINTAINS THE DISJOINT SETS ( I.E. THE SUBTREES. )
C     MARKER : INTEGER(NNODS)
C              MASKS NODES ALREADY CONSIDERED IN A STEP.
C     WEIGTH : INTEGER(0:NNODS)
C              CONTAINS THE WEIGHTS USED TO COMPUTE THE NODE COLUMN
C              LENGTHS.
C     VARWGH : INTEGER(0:NNODS)
C              CONTAINS THE WEIGHTS USED TO COMPUTE THE VARIABLE COLUMN
C              LENGTHS.
C     FDESC  : INTEGER(0:NNODS)
C              CONTAINS THE FIRST ( I.E. LOWEST-NUMBERED ) DESCENDANT.
C     LEVEL  : INTEGER(0:NNODS)
C              CONTAINS THE LEVEL, DISTANCE FROM THE ROOT FOR EACH NODE.
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           ELEMNT, FATHER, HINBR , I     , IFDESC, ISTOP ,
     $                  ISTRT , J     , JSTOP , JSTRT , K     , LAST1 ,
     $                  LAST2 , LCA   , LFLAG , LOWNBR, OLDNBR, PLEAF ,
     $                  TEMP  , TEMP1 , XSUP

C     ==================================================================

      XSUP = 0

      WEIGHT(0) = 0
      VARWGH(0) = 0
      LEVEL(0) = 0

      XNOD(1) = 1
      DO 10 K = 1, NNODS
         OLDNBR = PERM(K)
         XNOD(K+1) = XNOD(K) + XNODES(OLDNBR+1) - XNODES(OLDNBR)
   10 CONTINUE

      DO 100 K = NNODS, 1, -1
         ROWCNT(K) = 1
         COLCNT(K) = 0
         VARCNT(K) = 0
         SET(K) = K
         PRVLF(K) = 0
         LEVEL(K) = LEVEL(PARENT(K)) + 1
         WEIGHT(K) = 1
         VARWGH(K) = XNOD(K+1)-XNOD(K)
         FDESC(K) = K
         NCHILD(K) = 0
         PRVNBR(K) = 0
         MARKER(K) = 0
  100 CONTINUE

      NCHILD(0) = 0
      FDESC(0) = 0

      DO 200 K = 1, NNODS
         FATHER = PARENT(K)
         WEIGHT(FATHER) = 0
         VARWGH(FATHER) = 0
         NCHILD(FATHER) = NCHILD(FATHER) + 1
         IFDESC = FDESC(K)
         IF ( IFDESC .LT. FDESC(FATHER) )
     $   THEN
            FDESC(FATHER) = IFDESC
         ENDIF
  200 CONTINUE

C     ==================================================================

      DO 700 LOWNBR = 1, NNODS
         LFLAG = 0
         IFDESC = FDESC(LOWNBR)
         OLDNBR = PERM(LOWNBR)
         MARKER(LOWNBR) = LOWNBR
         ISTRT = XNODEL(OLDNBR)
         ISTOP = XNODEL(OLDNBR+1) - 1
         DO 600 I = ISTRT, ISTOP
            ELEMNT = NODEL(I)
            JSTRT = XELNOD(ELEMNT)
            JSTOP = XELNOD(ELEMNT+1) - 1
            DO 500 J = JSTRT, JSTOP
               HINBR = INVP(ELNOD(J))
               IF ( HINBR .GT. LOWNBR .AND.
     $              MARKER(HINBR) .NE. LOWNBR )
     $         THEN
                  MARKER(HINBR) = LOWNBR
                  IF ( IFDESC .GT. PRVNBR(HINBR) )
     $            THEN
                     WEIGHT(LOWNBR) = WEIGHT(LOWNBR) + 1
                     VARWGH(LOWNBR) = VARWGH(LOWNBR)
     $                              + XNOD(HINBR+1) - XNOD(HINBR)
                     PLEAF = PRVLF(HINBR)
                     IF ( PLEAF .EQ. 0 )
     $               THEN
                        ROWCNT(HINBR) = ROWCNT(HINBR) +
     $                              LEVEL(LOWNBR) - LEVEL(HINBR)
                     ELSE
                        LAST1 = PLEAF
                        LAST2 = SET(LAST1)
                        LCA = SET(LAST2)
 300                    IF ( LCA .NE. LAST2 )
     $                  THEN
                           SET(LAST1) = LCA
                           LAST1 = LCA
                           LAST2 = SET(LAST1)
                           LCA = SET(LAST2)
                           GO TO 300
                        ENDIF
                        ROWCNT(HINBR) = ROWCNT(HINBR)
     $                                + LEVEL(LOWNBR) - LEVEL(LCA)
                        WEIGHT(LCA) = WEIGHT(LCA) - 1
                        VARWGH(LCA) = VARWGH(LCA)
     $                              - XNOD(HINBR+1) + XNOD(HINBR)
                     ENDIF
                     PRVLF(HINBR) = LOWNBR
                     LFLAG = 1
                  ENDIF
                  PRVNBR(HINBR) = LOWNBR
               ENDIF
 500        CONTINUE
 600     CONTINUE
         FATHER = PARENT(LOWNBR)
         WEIGHT(FATHER) = WEIGHT(FATHER) - 1
         VARWGH(FATHER) = VARWGH(FATHER)
     $                  - XNOD(LOWNBR+1) + XNOD(LOWNBR)
         IF ( LFLAG .EQ. 1     .OR.
     $        NCHILD(LOWNBR) .GE. 2 )
     $   THEN
            XSUP = LOWNBR
         ENDIF
         IF ( XSUP .GT. 0 )
     $      SET(XSUP) = FATHER
 700  CONTINUE
      DO 800 K = 1, NNODS
         TEMP = COLCNT(K) + WEIGHT(K)
         TEMP1 = VARCNT(K) + VARWGH(K)
         COLCNT(K) = TEMP
         VARCNT(K) = TEMP1
         FATHER = PARENT(K)
         IF ( FATHER .NE. 0 )
     $   THEN
            COLCNT(FATHER) = COLCNT(FATHER) + TEMP
            VARCNT(FATHER) = VARCNT(FATHER) + TEMP1
         ENDIF
  800 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTR2
      END


      SUBROUTINE   SPRTR3
     $( SEPSZE, NNODS , PARENT, XNOD  , COLCNT, VARCNT, NCHILD, NSUPER,
     $  XSUPER, SUPLEN, SUPSZE, SUPSUP, XSIBL , SIBL  , PERM  , SUPRET,
     $  SEPARR, LPU   , ERROR  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           SEPSZE, NNODS , NSUPER, SUPRET, LPU   , ERROR
      INTEGER           PARENT(NNODS)         , XNOD(NNODS+1)         ,
     $                  COLCNT(NNODS)         , VARCNT(NNODS)         ,
     $                  NCHILD(0:NNODS)       , XSUPER(NNODS+1)       ,
     $                  SUPLEN(NNODS)         , SUPSZE(NNODS)         ,
     $                  SUPSUP(NNODS)         , XSIBL(NNODS+1)        ,
     $                  SIBL(NNODS+1)         , PERM(NNODS)           ,
     $                  SEPARR(NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRTR3
C
C --- TO COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION AND THE ADJACENCY
C     SET REPRESENTATION OF THE FUNDAMENTAL SUPERNODE ELIMINATION TREE.
C     AS A BYPRODUCT, THE ROUTINE ALSO COMPUTES THE SUPERNODE PARENT
C     INFORMATION OF THE SUP. EL. TREE.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : Oct. 03, 2007 (KMO)
C                 Return if NNODS is zero to avoid crash.
C
C
C     ON ENTRY
C     --------
C
*     SEPSZE : Integer
*              The size of the separator. Is always zero if this is not
*              a superelement or if superelement info. is suppressed.
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     PARENT : INTEGER(NNODS)
C              THE PARENT REPRESENTAION OF THE NODE ELIMINATION TREE
C              IN POSTORDER.
C     XNOD   : INTEGER(NNODS+1)
C              THE NODE PARTITION IN NEW ORDER OF THE NODES.
C     COLCNT : INTEGER(NNODS)
C              FOR EACH NODE, THIS ARRAY GIVES THE LENGTH OF THE COLUMN,
C              WITH RESPECT TO NODES.
C     VARCNT : INTEGER(NNODS)
C              FOR EACH NODE, THIS ARRAY GIVES THE LENGTH OF THE COLUMN,
C              WITH RESPECT TO VARIABLES.
C     NCHILD : INTEGER(0:NNODS)
C              FOR EACH NODE, THIS ARRAY GIVES ITS NUMBER OF CHILDREN.
*     PERM   : Integer(NNODS)
*              The permutation vector related to the current postorder.
*     SEPARR : Integer(NNODS)
*              The nodes i in separator set s has SEPARR(i) = s > 0.
*     LPU    : Integer
*              Unit for error messages.
C
C     ON EXIT
C     -------
C
C     COLCNT : INTEGER(NNODS)
C              DESTROYED.
C     NCHILD : INTEGER(0:NNODS)
C              DETROYED.
C     NSUPER : INTEGER
C              NUMBER OF FUNDAMENTAL SUPERNODES.
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION
C     SUPLEN : INTEGER(NSUPER)
C              PARENT INFO OF THE SUPERNODE EL. TREE.
C     SUPSUP : INTEGER(NSUPER)
C              Array to control the processing in numerical steps.
C     XSIBL  : INTEGER(NSUPER+1)
C              POINTER INTO SIBL FOR THE CHILDREN OF EACH NODE IN THE
C              SUPERNODE ELIMINATION TREE.
C     SIBL   : INTEGER(NSUPER+1)
C              LISTS OF CHILDREN OF THE SUPERNODES.
C              THE SUPERNODE TREE ADJACENCY SET LISTS.
C              SIBL(NSUPER+1) STORES THE NUMBER OF ROOTS IN THE
C                             SUPERNODE ELIMINATION TREE.
C              IF IP=SIBL(NSUPER+1) THEN THE ROOTS ARE FOUND IN THE
C              LOCATIONS SIBL(NSUPER-IP+1), ... SIBL(NSUPER).
*     ERROR  : Integer
*              =  0: Normal return.
*              = -1: Error in the separators (only checked for SEPSZE
*                    greater than zero ).
C
C     WORKING ARRAYS
C     --------------
C
C     COLCNT : INTEGER(NNODS)
C              HOLDS THE SUPERNODE PARENT VECTOR.
C     NCHILD : INTEGER(NNODS)
C              USED AS WORK ARRAY IN ROUTINE SPRT31.
C     SIBL   : INTEGER(NNODS)
C              USED TO HOLD THE NODE SUPERNODE MAP VECTOR.
C
C     SUBPROGRAMS
C     -----------
C
C     SPRT31
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER
     $  CHILDS, COLLEN, FATHER, J     , JS    , JSTOP , JSTRT , NEXT  ,
     $  NEXTS , PREV  , PREVS , PRVLEN

      EXTERNAL          SPRT31

C     ==================================================================

* --- Initialize.
      ERROR = 0
      SUPRET = 0
      NSUPER = 0
      XSUPER(1) = 1
      IF ( NNODS .LE. 0 ) RETURN

C     COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION.
      NSUPER = 1
      PRVLEN = COLCNT(1)
      FATHER = PARENT(1)

      PREV = PERM(1)

* --- Compute the array SUPSUP that controls the actions on the
*     supernodes in the numerical steps. While doing this also
*     check that the separator is valid for the superelement.

* --- SUPSUP(JS) > 0 : The nodes/variables that are contained in JS
*                      are retained nodes. We compute also SUPRET,
*                      the number of supernodes that are retained.
      IF ( SEPSZE.GT.0 ) THEN
         PREVS = SEPARR(PREV)
         IF ( PREVS.GT.0 ) THEN
              SUPSUP(NSUPER) = NSUPER
              SUPRET = SUPRET + 1
         ELSE
            SUPSUP(NSUPER) = 0
         ENDIF
      ELSE
         SUPSUP(NSUPER) = 0
      ENDIF

      DO 100 J = 2, NNODS

         CHILDS = NCHILD(J)
         COLLEN = COLCNT(J) + 1

         NEXT = PERM(J)

        IF ( SEPSZE.GT.0 ) THEN

            PREVS = SEPARR(PREV)
            NEXTS = SEPARR(NEXT)

            IF ( ( FATHER.EQ.J      ).AND.
     $           ( CHILDS.EQ.1      ).AND.
     $           ( COLLEN.EQ.PRVLEN ).AND.
     $           ( NEXTS.EQ.PREVS   )      ) THEN

C - NO ACTION, NODE J IS IN THE SAME FUNDAMENTAL SUPERNODE AS ITS CHILD.

            ELSE

C - NODE J IS THE FIRST NODE IN A NEW SUPERNODE.
               NSUPER = NSUPER + 1
               XSUPER(NSUPER) = J

               IF ( NEXTS.GT.0 ) THEN
                  SUPSUP(NSUPER) = NSUPER
                  SUPRET = SUPRET + 1
               ELSE
                  SUPSUP(NSUPER) = 0
               ENDIF

            ENDIF

         ELSE

            IF ( ( FATHER.EQ.J      ).AND.
     $           ( CHILDS.EQ.1      ).AND.
     $           ( COLLEN.EQ.PRVLEN )      ) THEN

C - NO ACTION, NODE J IS IN THE SAME FUNDAMENTAL SUPERNODE AS ITS CHILD.

            ELSE

C - NODE J IS THE FIRST NODE IN A NEW SUPERNODE.
               NSUPER = NSUPER + 1
               XSUPER(NSUPER) = J
               SUPSUP(NSUPER) = 0

            ENDIF

         ENDIF

         PRVLEN = COLLEN - 1
         FATHER = PARENT(J)
         SIBL(J) = NSUPER

         PREV = NEXT

  100 CONTINUE

      XSUPER(NSUPER+1) = NNODS + 1

      IF ( SEPSZE.GT.0 ) THEN

* ------ Run through the supernodes and check the separators.
*        This test is based on the fact that each supernode is a clique
*        and thus should be in the same separator.
         DO JS = 1, NSUPER

* --------- Get the previous node.
            PREV = PERM(XSUPER(JS))

            DO J = XSUPER(JS)+1, XSUPER(JS+1)-1

* ------------ Get the previous separator.
               PREVS = SEPARR(PREV)

* ------------ Get the current node and its separator.
               NEXT = PERM(J)
               NEXTS = SEPARR(NEXT)

* ------------ The rule is that all nodes within a supernode
*              must belong to the same separator.
               IF ( PREVS.NE.NEXTS ) THEN
                  ERROR = -1
                  WRITE ( LPU, '(//A)' )
     $            '*** Error from SPRTRS::SPRTR3'
                  WRITE ( LPU, '(4X,A,2I10/)' )
     $            'Both nodes should have been in the same separator',
     $            PREVS, NEXTS
                  RETURN
               ENDIF

               PREV = NEXT

            END DO

         END DO
      ENDIF

C     ==================================================================

C     ------------------------------------------------------
C     COMPUTE THE SUPERNODE ELIMINATION TREE IN PARENT FORM,
C     AND THE SUPERNODE LENGTH AND SIZE FOR EACH SUPERNODE
C     WITH REFERENCE TO VARIABLES.
C     ------------------------------------------------------

      DO 300 JS = 1, NSUPER

         JSTRT  = XSUPER(JS)
         JSTOP  = XSUPER(JS+1)-1
         FATHER = PARENT(JSTOP)
         IF ( FATHER .GT. 0 ) THEN
            COLCNT(JS) = SIBL(FATHER)
         ELSE
            COLCNT(JS) = 0
         ENDIF
         SUPSZE(JS) = 0
         SUPLEN(JS) = VARCNT(JSTRT)
         DO 200 J = JSTRT, JSTOP
            SUPSZE(JS) = SUPSZE(JS) + XNOD(J+1) - XNOD(J)
  200    CONTINUE
  300 CONTINUE

C     ----------------------------------------
C     COMPUTE THE ADJACENCY SET REPRESENTATION
C     OF THE SUPERNODE ELIMINATION TREE.
C     ----------------------------------------

      CALL SPRT31 ( NSUPER, COLCNT, XSIBL, SIBL, NCHILD )

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTR3
      END


      SUBROUTINE   SPRTR4   ( SEPSZE, NSUPER, NNODS , OPTION, XSIBL ,
     $                        SIBL  , XSUPER, SUPLEN, SUPSUP, PERM  ,
     $                        INVP  , SUPSZE, SUPERM, SUINVP, NWPERM,
     $                        NXSUPR, NXSIBL, NSIBL , IERR  , LPU     )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           SEPSZE, NSUPER, NNODS , OPTION, IERR  , LPU
      INTEGER           XSIBL(NSUPER+1)       , SIBL(NSUPER+1)        ,
     $                  XSUPER(NSUPER+1)      , SUPLEN(NSUPER)        ,
     $                  SUPSUP(NSUPER)        , PERM(NNODS)           ,
     $                  INVP(NNODS)           , SUPSZE(NSUPER)        ,
     $                  SUPERM(NSUPER)        , SUINVP(NSUPER)        ,
     $                  NWPERM(NNODS)         , NXSUPR(NSUPER+1)      ,
     $                  NXSIBL(NSUPER+1)      , NSIBL(NSUPER+1)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRTR4
C
C --- TO PERFORM A DEPTH FIRST POSTORDER OF THE SUPERNODE ELIMINATION
C     TREE, WHERE THE TREE IS OPTIMALLY SORTED, (OPTION=0), OR IS AS IS
C     (OPTION=1). FINALLY IT UPDATES THE TREE STRUCTURES BY A CALL TO
C     SPRT44.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS :
C
C     ON ENTRY
C     --------
C
*     SEPSZE : Integer
*              The size of the separator. Is always zero if this is not
*              a superelement or if superelement info. is suppressed.
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NNODS  : INTEGER
C              NUMBER OF FINITE ELEMENT NODES.
C     OPTION : INTEGER
C              OPTIMAL OR DFS POSTORDER.
C              = 0 : OPTIMAL.
C              = 1 : DFS.
C     (XSIBL,
C     SIBL)  : INTEGER(NSUPER+1), INTEGER(NSUPER)
C              THE SUPERNODE ELIMINATION TREE ON ADJACENCY FORMAT.
C              THE SUPERNODE TREE ADJACENCY SET LISTS.
C              SIBL(NSUPER+1) STORES THE NUMBER OF ROOTS IN THE
C                             SUPERNODE ELIMINATION TREE.
C              IF IP=SIBL(NSUPER+1) THEN THE ROOTS ARE FOUND IN THE
C              LOCATIONS SIBL(NSUPER-IP+1), ... SIBL(NSUPER).
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION RELATED TO NODES.
C     SUPLEN : INTEGER(NSUPER)
C              THE LENGTH OF EACH SUPERNODE RELATED TO VARIABLES.
C     SUPSUP : INTEGER(NSUPER)
C              Array to control the processing in numerical steps.
C     PERM   : INTEGER(NNODS)
C              THE PERMUTATION, I.E. THE NEW TO OLD MAP.
C     INVP   : INTEGER(NNODS)
C              THE INVERSE PERMUTATION, I.E. THE OLD TO NEW MAP.
C     SUPSZE : INTEGER(NSUPER)
C              THE SIZE OF EACH SUPERNODE.
C
C     ON EXIT
C     -------
C
C     (XSIBL,
C     SIBL)  : INTEGER(NSUPER+1), INTEGER(NSUPER)
C              THE SUPERNODE ELIMINATION TREE ON ADJACENCY FORMAT.
C              THE CHILDREN ARE ORDERED ACCORDING TO THE OPTIMAL
C              SEQUENCE.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPLEN : INTEGER(NSUPER)
C              THE LENGTH RELATED TO VARIABLES FOR EACH SUPERNODE.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPSUP : INTEGER(NSUPER)
C              Array to control the processing in numerical steps.
C              UPDATED WITH POSTORDER INFORMATION.
C     PERM   : INTEGER(NNODS)
C              THE PERMUTATION, I.E. THE NEW TO OLD MAP.
C              UPDATED WITH POSTORDER INFORMATION.
C     INVP   : INTEGER(NNODS)
C              THE INVERSE PERMUTATION, I.E. THE OLD TO NEW MAP.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPSZE : INTEGER(NSUPER)
C              THE SIZE RELATED TO VARIABLES FOR EACH SUPERNODE.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPERM : INTEGER(NSUPER)
C              THE OPTIMAL SEQUENCE OF THE SUPERNODES.
C     IERR   : INTEGER
C              ERROR FLAG.
C              =  0 : NORMAL RETURN.
C              = -1 : ERROR FROM SPRT44.
C
C     WORKING ARRAYS
C     --------------
C
C     SUPERM : INTEGER(NSUPER)
C              STORES THE WORKSPACE NEEDED FOR EACH SUPERNODE.
C     SUINVP : INTEGER(NSUPER)
C              HOLDS THE INVERSE SUPERNODE PERMUTATION.
C     NWPERM : INTEGER(NNODS)
C              PERMUTATION INCREMENT DUE TO SUPERNODE POSTORDER.
C     NXSUPR : INTEGER(NSUPER+1)
C              FIRST, SIZE OF UPDATE MATRICES FOR EACH SUPERNODE,
C              THEN THE UPDATED SUPERNODE PARTITION IN SPRT44.
C     NXSIBL : INTEGER(NSUPER+1)
C              FIRST, INDEX FOR OPTIMAL CHILDREN SEQUENCE, THEN
C              THE UPDATED POINTER TO SUPERNODE CHILDREN IN SPRT44.
C     NSIBL  : INTEGER(NSUPER+1)
C              FIRST, QUANTITY THAT DETERMINES THE SEQUENCE OF THE
C              CHILDREN, THEN THE UPDATED CHILDREN LISTS IN SPRT44.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     SPRT41
C     SPRT42
C     SPRT43
C     SPRT44
C
C     ------------------------------------------------------------------
C
      INTEGER
     $  CHILD , ERROR , F     , I     , IP    , ISTRT , ISTOP , JS    ,
     $  S

      EXTERNAL          SPRT41, SPRT42, SPRT43, SPRT44

C     ==================================================================

      IF ( OPTION .EQ. 0 ) THEN

C        ------------------------------------------------
C        SORT THE CHILDREN ACCORDING TO THE YANG FORMULA.
C        ------------------------------------------------

         DO 300 JS = 1, NSUPER

            F = SUPLEN(JS)
            S = SUPSZE(JS)

C           -------------------------------------------------
C           NXSUPR(JS) IS SET TO THE SIZE OF THE FRONT MATRIX.
C           -------------------------------------------------

            NXSUPR(JS) = ((F*(F+1))/2)
            ISTRT = XSIBL(JS)
            ISTOP = XSIBL(JS+1) - 1
            IF ( ISTRT .LE. ISTOP )
     $      THEN
               IP = 0
               DO 100 I = ISTRT, ISTOP
                  IP = IP + 1
                  CHILD = SIBL(I)
                  NSIBL(IP) = SUPERM(CHILD) - NXSUPR(CHILD)
  100          CONTINUE
               IF ( IP .GT. 1 )
     $         THEN
                  DO 200 I = 1, IP
                     NXSIBL(I) = I
  200             CONTINUE
                  CALL SPRT41 ( IP, NXSIBL, IP, NSIBL )
               ELSE
                  NXSIBL(1) = 1
               ENDIF
               CALL SPRT42 ( IP, JS, SUPERM, SIBL(ISTRT), NXSUPR,
     $                       NXSIBL, NSIBL )
            ELSE
               SUPERM(JS) = NXSUPR(JS)
            ENDIF

C           ---------------------------------------------------------
C           FINALLY, SET NXSUPR(JS) TO THE SIZE OF THE UPDATE MATRIX.
C           ---------------------------------------------------------

            NXSUPR(JS) = (((F-S)*(F-S+1))/2)

  300    CONTINUE

      ENDIF

C     ---------------------------------------
C     COMPUTE SUPERM BY A DEPTH FIRST SEARCH.
C     ---------------------------------------

      CALL SPRT43 ( NSUPER, XSIBL, SIBL, SUPERM, NXSUPR, NXSIBL )

C     -----------------------------------------------
C     UPDATE THE TREE STRUCTURES BY A CALL TO SPRT44.
C     -----------------------------------------------
      CALL SPRT44
     $( SEPSZE,NNODS,NSUPER,XSUPER,XSIBL,SIBL,PERM,INVP,SUPSZE,SUPLEN,
     $  SUPSUP,SUPERM,NWPERM,NXSUPR,NXSIBL,NSIBL,SUINVP,ERROR,LPU )

      IF ( ERROR.NE.0 ) THEN
         IERR = -1
         CALL SPRER1 ( 100, 'SPRTR4', 0, 0, 0, LPU, ERROR )
         RETURN
      ENDIF

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTR4
      END


      SUBROUTINE   SPRTR5
     $( NSUPER, SPRNOD, NEQ   , NEL   , SPRCON, XELNOD, ELNOD , XNODEL,
     $  NODEL , XNODES, XSIBL , SIBL  , SUPSZE, SUPSUP, PERM  ,
     $  NB    , MAXSUP, MAXFRT, NOFSUB, NZEROL, NZSTCK, DPSTCK,
     $  NELACT, LCOEFA, ELSIZE, DOFRET, SZERET, SZEREP, XELSEQ, ELSEQ ,
     $  XBLOCK, XSUPER, SUPLEN, RINFO , MASK                )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      include 'integer_kind.inc'

      INTEGER           NSUPER, SPRNOD, NEQ   , NEL   , SPRCON, NB    ,
     $                  MAXSUP, MAXFRT, NOFSUB, NZSTCK, DPSTCK,
     $                  NELACT, LCOEFA, ELSIZE, DOFRET, SZERET, SZEREP
      INTEGER(kind=i8)  NZEROL

      INTEGER           XELNOD(NEL+1)         , ELNOD(SPRCON)         ,
     $                  XNODEL(SPRNOD+1)      , NODEL(SPRCON)         ,
     $                  XNODES(SPRNOD+1)      ,
     $                  XSIBL(NSUPER+1)       , SIBL(NSUPER+1)        ,
     $                  SUPSZE(NSUPER)        , SUPSUP(NSUPER)

      INTEGER           PERM(SPRNOD)          ,
     $                  XELSEQ(SPRNOD+1)      , ELSEQ(NEL)            ,
     $                  XBLOCK(NSUPER+1)      , XSUPER(NSUPER+1)      ,
     $                  SUPLEN(NSUPER)        , MASK(NEL)

      DOUBLE PRECISION  RINFO(4)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRTR5
C
C --- TO COMPUTE THE NUMBER OF EXTERNAL ASSEMBLY STEPS NB, THE ELEMENT
C     SEQUENCE VECTORS (XELSEQ,ELSEQ), THE ELEMENT ASSEMBLY BLOCK
C     STRUCTURE XBLOCK, TO UPDATE XSUPER TO RELATE TO VARIABLES INSTEAD
C     OF TO NODES, AND FINALLY TO COMPUTE STORAGE REQUIREMENT AND
C     FACTORIZATION STATISTICS.
C
C --- N O T E: IT IS ASSUMED THAT ALL DATA STRUCTURES ARE UPDATED
C              ACCORDING TO THE FINAL OPTIMAL/DFS POSTORDER.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : JAN. 13, 2004 (KMO)
C                 Changed to kind=i8 for NZEROL to capture overflow.
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     SPRNOD : INTEGER
C              NUMBER OF NODES WITH AT LEAST ONE FREE VARIABLE.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     NEL    : INTEGER
C              NUMBER OF ORIGINAL FINITE-ELEMENST IN THE MODEL.
C     SPRCON : INTEGER
C              SIZE OF THE ACTIVE CONNECTIVITY ARRAYS.
C     XELNOD : INTEGER(NEL+1)
C              POINTER INTO ELNOD.
C     ELNOD  : INTEGER(SPRCON)
C              LIST OF NODES CONNECTED TO EACH ELEMENT.
C     XNODEL : INTEGER(SPRNOD+1)
C              POINTER INTO THE LIST OF ELEMENTS CONNECTED TO
C              EACH ACTIVE NODE.
C     NODEL  : INTEGER(SPRCON)
C              LIST OF THE ELMENTS CONNECTED TO EACH ACTIVE NODE.
C     XNODES : INTEGER(SPRNOD+1)
C              THE SPR-NODE PARTITION, AND POINTER TO NODES.
C     XSIBL  : INTEGER(NSUPER+1)
C              POINTER TO THE SUPERNODE TREE ADJACENCY SET.
C     SIBL   : INTEGER(NSUPER+1)
C              THE SUPERNODE TREE ADJACENCY SET LISTS.
C              SIBL(NSUPER+1) STORES THE NUMBER OF ROOTS IN THE
C                             SUPERNODE ELIMINATION TREE.
C              IF IP=SIBL(NSUPER+1) THEN THE ROOTS ARE FOUND IN THE
C              LOCATIONS SIBL(NSUPER-IP+1), ... SIBL(NSUPER).
C     SUPSZE : INTEGER(NSUPER)
C              THE SIZE OF EACH SUPERNODE RELATED TO VARIABLES
C              IN THE FINAL OPTIMAL POSTORDER.
C     PERM   : INTEGER(SPRNOD)
C              THE FINAL POSTORDER OF THE NODES.
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION VECTOR RELATED TO NODES.
C
C     ON EXIT
C     -------
C
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION VECTOR RELATED VARIABLES.
C     NB     : INTEGER
C              NUMBER OF EXTERNAL ASSEMBLY STEP, IE. THE NUMBER OF
C              CALLS TO MA46BD IN ORDER TO FACTORIZE THE MATRIX.
C     MAXSUP : INTEGER
C              THE LARGEST SUPERNODE.
C     MAXFRT : INTEGER
C              THE LARGEST FRONT MATRIX.
C     NOFSUB : INTEGER
C              THE NUMBER OF ROW/COLUMN INDICES NEEDED TO REPRESENT
C              L/U IF THE TENTATIVE PIVOT ORDER IS USED BY MA46BD.
C     NZEROL : INTEGER
C              THE SIZE OF L.
C     NZSTCK : INTEGER
C              SIZE OF THE WORKING STACK IF TENTATIVE ORDER IS USED.
C     DPSTCK : INTEGER
C              THE DEPTH OF THE STACK.
C     NELACT : INTEGER
C              NUMBER OF ORIGINAL FINITE ELEMENTS THAT ARE NEEDED IN THE
C              ASSEMBLY STEP.
C     LCOEFA : INTEGER
C              THE LENGTH OF THE COEFFICIENT MATRIX IF IT IS STORED
C              AS A SEQUENCE OF FRONT MATRICES IN EACH ASSEMBLY STEP.
C     ELSIZE : INTEGER
C              MAXIMUM ORDER OF AN EMBEDDED ELEMENT MATRIX WITH RESPECT
C              TO VARIABLES.
C     XBLOCK : INTEGER(NB+1)
C              THE ASSEMBLY PARTITION OF THE SUPERNODES.
C     XELSEQ : INTEGER(NB+1)
C              POINTER TO THE LIST OF ORIGINAL FINITE ELEMENTS
C              NEEDED IN EACH ASSEMBLY STEP.
C     ELSEQ  : INTEGER(NEL)
C              THE LISTS OF ORIGINAL FINITE ELEMENTS THAT ARE NEEDED
C              IN EACH ASSEMBLY STEP.
C     RINFO  : DOUBLE PRECISION(4), OR REAL(4).
C              ASSEMBLY AND FACTORIZATION COUNTS.
C              RINFO(1) : THE NUMBER OF FLOPS USED TO ASSEMBLE THE
C                         ORIGINAL FINITE ELEMENT COEFFICIENTS.
C              RINFO(2) : THE NUMBER OF FLOPS USED TO ASSEMBLE THE
C                         GENERATED ELEMENT COEFFICIENTS.
C              RINFO(3) : THE NUMBER OF FLOPS USED TO FACTORIZE THE
C                         MATRIX INTO ITS TRIANGULAR FACTORS.
C              RINFO(4) : THE NUMBER OF FLOPS TO SOLVE ONE SYSTEM.
C
C     WORKING ARRAYS
C     --------------
C
C     MASK   : INTEGER(NEL)
C              MASKS OFF ELEMENTS THAT ARE ALREADY ASSEMBLED.
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     MAX
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           CHILD , COLLEN, DEPTH , EF    , EFNODE, ELEMNT,
     $                  F     , FCH   , I     , IP    , J     , JS    ,
     $                  K     , KSTRT , KSTOP , NODE  ,
     $                  S     , SCH   , STACK
      LOGICAL           LFIRST

      DOUBLE PRECISION  ONE         , TWO         , SIX         ,
     *                  ZERO
      PARAMETER       ( ONE  = 1.0D0, TWO  = 2.0D0, SIX  = 6.0D0,
     *                  ZERO = 0.0D0 )
      DOUBLE PRECISION  SUMS1, SUMS2, SS, FF, SSCH, FFCH

      INTRINSIC         MAX

C     ==================================================================

C     -----------------------------------------------
C     INITIALIZE MASK TO ZERO IN ORDER TO FACILITATE
C     MASK OF EACH ELEMENTS WHICH HAS BEEN ASSEMBLED.
C     -----------------------------------------------
      DO 100 I = 1, NEL
         MASK(I) = 0
  100 CONTINUE

      RINFO(1) = ZERO

C     ------------------------------------
C     NB NUMBER OF EXTERNAL ASSEMBLY STEP.
C     IP POINTER INTO ELSEQ.
C     NELACT THE NUMBER OF ELEMENTS THAT
C            ARE ACCESSED IN THE ASSEMBLY.
C     ------------------------------------
      NB = 0
      IP = 1
      NELACT = 0
      LCOEFA = 0
      ELSIZE = 0

* --- The number of retained degrees of freedom.
      DOFRET = 0

* --- Size of the generated element.
      SZERET = 0

* --- Size of the index representation for the generated element.
      SZEREP = 0

      DO 400 JS = 1, NSUPER

C        ------------------------------
C        SET THE FIRST FLAG TO TRUE FOR
C        THE NEW SUPERNODE.
C        ------------------------------
         LFIRST = .TRUE.
         DO 300 I = XSUPER(JS), XSUPER(JS+1) - 1

C           ------------------------------------------
C           RUN THROUGH THE NODES OF THE SUPERNODE ...
C           ------------------------------------------
            NODE = PERM(I)
            DO 200 J = XNODEL(NODE), XNODEL(NODE+1) - 1

C           ---------------------------------------------
C           ... AND LOOK AT THE ELEMENTS CONNECTED TO THE
C               NODES TO SEE IF THEY NEED TO BE MERGED.
C           ---------------------------------------------
               ELEMNT = NODEL(J)
               IF ( MASK(ELEMNT) .LE. 0 )
     $         THEN

C                 -------------------------------
C                 THIS ELEMENT IS NOT YET MERGED.
C                 -------------------------------
                  IF ( LFIRST )
     $            THEN

C                    -----------------------------------
C                    THE ELEMENT IS THE FIRST IN A NEW
C                    EXTERNAL ASSEMBLY STEP, MARK AS SO.
C                    -----------------------------------
                     NB = NB + 1
                     XELSEQ(NB) = IP
                     XBLOCK(NB) = JS
                     COLLEN = SUPLEN(JS)
                     LCOEFA = LCOEFA + (((COLLEN+1)*COLLEN)/2)
                     LFIRST = .FALSE.
                  ENDIF

C                 -----------------------------------------
C                 COLLECT THE NUMBER OF ASSEMBLY OPERATIONS
C                 OF THE ORIGINAL FINITE-ELEMENTS.
C                 -----------------------------------------
                  KSTRT = XELNOD(ELEMNT)
                  KSTOP = XELNOD(ELEMNT+1) - 1

C                 -------------------------------------
C                 F IS THE ORDER OF THE ELEMENT MATRIX.
C                 EF IS THE ORDER OF THE NODE.
C                 -------------------------------------
                  F = 0
                  DO 150 K = KSTRT, KSTOP
                     EFNODE = ELNOD(K)
                     EF = XNODES(EFNODE+1)-XNODES(EFNODE)
                     IF ( EF .GT. 0 )
     $                  F = F + EF
  150             CONTINUE
C --------------- SYMMETRIC MATRIX
                  FF = DBLE(F)
                  RINFO(1) = RINFO(1) + (((FF+ONE)*FF)/TWO)
C*DAM * --------------- UNSYMMETRIC MATRIX
C*DAM                   RINFO(1) = RINFO(1) + (F*F)

C                 ---------------------------
C                 KEEP LARGEST ELEMENT ORDER.
C                 ---------------------------
                  ELSIZE = MAX(ELSIZE,F)

C                 -------------------------------
C                 STORE THE ELEMENT IN THE LIST
C                 AND INCREMENT THE LIST POINTER,
C                 AND THE NUMBER OF ELEMENTS USED.
C                 -------------------------------
                  ELSEQ(IP) = ELEMNT
                  IP = IP + 1
                  NELACT = NELACT + 1

C                 ---------------------------------------
C                 IF ALL THE ORIGINAL ELEMENTS ARE MERGED
C                 WE CAN JUMP BECAUSE THE REST OF THE
C                 SUPERNODES MAY BE MERGED TO THIS STEP.
C                 ---------------------------------------
                  IF ( IP .GT. NEL )
     $               GOTO 500

                  MASK(ELEMNT) = ELEMNT
               ENDIF
  200       CONTINUE
  300    CONTINUE
  400 CONTINUE

C     --------------------
C     SET THE LAST LIMITS.
C     --------------------
  500 XELSEQ(NB+1) = NELACT + 1
      XBLOCK(NB+1) = NSUPER + 1

C     ==================================================================

      MAXSUP = 0
      MAXFRT = 0
      NOFSUB = 0
      NZEROL = 0
      NZSTCK = 0
      DPSTCK = 0

      DEPTH  = 0
      STACK  = 0

      RINFO(2) = ZERO
      RINFO(3) = ZERO
      RINFO(4) = ZERO

      XSUPER(1) = 1

      DO 700 JS = 1, NSUPER

C        -------------------------
C        S SIZE OF THE SUPERNODE.
C        F ORDER OF THE SUPERNODE.
C        -------------------------
         S = SUPSZE(JS)
         F = SUPLEN(JS)

C        ---------------------------------------
C        MAXSUP LARGEST SUPERNODE.
C        MAXFRT MAXIMUM ORDER OF A FRONT MATRIX.
C        ---------------------------------------
         MAXSUP = MAX(MAXSUP,S)
         MAXFRT = MAX(MAXFRT,F)

C        -------------------------------------------
C        FZNONZ FACTOR NONZEROS.
C        ELBOWR SPACE NECESSARY TO DO FACTORIZATION.
C        -------------------------------------------
C ------ SYMMETRIC MATRIX.
C*KMO          ELBOWR = (((F+1)*F)/2) + FZNONZ + STACK

C*DAM * ------ UNSYMMETRIC MATRIX.
C*DAM          ELBOWR = ((F-S)*S) + (F*F) + FZNONZ + STACK

C        -------------------------------------------------
C        LA MINIMUM SIZE OF A FOR SUCCESSFUL FACTORIZATION
C           OF A MATRIX THAT DOES NOT SUFFER FROM DELAYED
C           ELIMIANTIONS.
C        -------------------------------------------------
C*KMO          LA = MAX(LA,ELBOWR)

C        -------------------------------------------------
C        NOFSUB NUMBER OF COLUMN INDICES TO REPRESENT L/U.
C        NZEROL SIZE OF L.
C        -------------------------------------------------
         NOFSUB = NOFSUB + F
         NZEROL = NZEROL + int((F*S) - ((S*(S-1))/2),i8)

C ------ SYMMETRIC MATRIX.
C*KMO          FZNONZ = FZNONZ + ((F*S) - ((S*(S-1))/2))

C*DAM * ------ UNSYMMETRIC MATRIX.
C*DAM          FZNONZ = FZNONZ + 2*((F*S) - ((S*(S-1))/2)) - S

C        ----------------------------------------
C        STACK CURRENT SIZE OF THE WORKING STACK.
C        ----------------------------------------
C ------ SYMMETRIC MATRIX.
         STACK  = STACK  + (((F+1)*F)/2)

C*DAM * ------ UNSYMMETRIC MATRIX.
C*DAM          STACK  = STACK  + (F*F)

C        --------------------------------------------------------------
C        DEPTH NUMBER OF UPDATE MATRICES CURRENTLY STORED ON THE STACK.
C        --------------------------------------------------------------
         DEPTH  = DEPTH  + 1


C        ------------------------------------------
C        NZSTCK MAXIMUM REQUIRED SIZE OF THE STACK.
C        DPSTCK MAXIMUM DEPTH OF THE STACK.
C        ------------------------------------------
         NZSTCK = MAX(NZSTCK,STACK)
         DPSTCK = MAX(DPSTCK,DEPTH)

C        ----------------------------------------------
C        XSUPER IS NOW UPDATED TO HOLD THE SIZE RELATED
C               TO VARIABLES FOR EACH SUPERNODE.
C        ----------------------------------------------
         XSUPER(JS+1) = XSUPER(JS) + S

C        -----------------------------------------
C        RUN THROUGH THE CHILDREN OF THE SUPERNODE
C        AND REDUCE THE SIZE OF STACK AND DEPTH.
C        -----------------------------------------
         DO 600 I = XSIBL(JS), XSIBL(JS+1) - 1
            CHILD = SIBL(I)
            DEPTH = DEPTH - 1
            SCH = SUPSZE(CHILD)
            FCH = SUPLEN(CHILD)

C --------- SYMMETRIC MATRIX.
            STACK = STACK - (((FCH-SCH+1)*(FCH-SCH))/2)

C*DAM * --------- UNSYMMETRIC MATRIX.
C*DAM             STACK = STACK - ((FCH-SCH)*(FCH-SCH))

C           ----------------------------------
C           RINFO(2) NUMBER OF ASSEMBLY
C           OPERATIONS FOR GENERATED ELEMENTS.
C           ----------------------------------
C --------- SYMMETRIC MATRIX
            SSCH = DBLE(SCH)
            FFCH = DBLE(FCH)
            RINFO(2) = RINFO(2)
     $               + (((FFCH-SSCH+ONE)*(FFCH-SSCH))/TWO)
C*DAM * --------- UNSYMMETRIC MATRIX
C*DAM             RINFO(2) = RINFO(2)
C*DAM      $               + ((FCH-SCH)*(FCH-SCH))

  600    CONTINUE

         IF ( SUPSUP(JS).LE.0 ) THEN

C           --------------------------------------------
C           RINFO(3) NUMBER OF FACTORIZATION OPERATIONS.
C           --------------------------------------------
            SS = DBLE(S)
            FF = DBLE(F)
            SUMS1 = (SS*(SS+ONE))/TWO
            SUMS2 = (SS*(SS+ONE)*(TWO*SS+ONE))/SIX

C --------- SYMMETRIC MATRIX
            RINFO(3) = RINFO(3)
     $            + (((FF*FF)+(TWO*FF))*SS) - (TWO*(ONE+FF)*SUMS1)
     $            + SUMS2
C*DAM * ------ UNSYMMETRIC MATRIX
C*DAM          RINFO(3) = RINFO(3)
C*DAM      $            + (((F*F)+F)*S) - ((1+(2*F))*SUMS1)
C*DAM      $            + SUMS2

C           ---------------------------------------
C           RINFO(4) NUMBER OF SOLUTION OPERATIONS.
C                    ACCUMULATE FORWARD SOLVE OPS.
C           ---------------------------------------
            RINFO(4) = RINFO(4) + (FF*SS) - (SS*(SS-ONE)/TWO) - SS

         ELSE

            DOFRET = DOFRET + S
            SZEREP = SZEREP + F
            SZERET = SZERET + ( S*F ) - ( ( (S-1)*S )/2 )

         ENDIF

C        ----------------------------------------
C        STORE THE TRIANGULAR FACTORS IN THEIR
C        PERMANENT STORAGE, I.E REDUCE THE STACK.
C        ----------------------------------------
C ------ SYMMETRIC MATRIX.
         STACK = STACK
     $         - ((F*S) - ((S*(S-1))/2))

C*DAM * ------ UNSYMMETRIC MATRIX.
C*DAM          STACK = STACK
C*DAM      $         - (((F*S) - ((S*(S-1))/2)) - S )
C*DAM      $         - ((F*S) - ((S*(S-1))/2))

C        ------------------------------------------
C        IF SIZE S = ORDER F, THEN THE STACK SPACE
C        IS RELEASED IMMEDIATELY AFTER ELIMINATION.
C        ------------------------------------------
         IF ( F .EQ. S )
     $      DEPTH = DEPTH - 1

  700 CONTINUE

C*DAM *     -----------
C*DAM *     SET LKEEPB.
C*DAM *     -----------
C*DAM       LKEEPB = (4*NSUPER) + NEQ + NOFSUB

C     ------------------
C     FINALIZE RINFO(4).
C     ------------------
      RINFO(4) = TWO*RINFO(4) + DBLE(NEQ)

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRTR5
      END


      SUBROUTINE   SPRT11   ( N     , PERM  , INVP   )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           N
      INTEGER           PERM(N)       , INVP(N)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT11 : TO COMPUTE THE ARRAY INVP, GIVEN THE ARRAY PERM.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     N     : INTEGER
C             DIMENSION OF THE ARRAYS PERM AND INVP
C     PERM  : INTEGER(N)
C             PERM(I) = K, MEANS THAT THE ORIGINAL NODE NUMBER K IS THE
C             I'TH NODE IN THE NEW ORDERING. THAT IS THE ARRAY PERM DE-
C             FINES A MAP FROM NEW NODE NUMBERS TO OLD NODE NUMBERS.
C
C     ON EXIT
C     -------
C
C     INVP  : INTEGER(N)
C             INVP(PERM(I)) = I, IE. INVP(K) = I, MEANS THAT THE ORI-
C             GINAL NODE K IS THE I'TH NODE IN THE NEW ORDERING. IE.
C             INVP DEFINES A MAP FROM OLD NODE NUMBERS TO NEW NODE
C             NUMBERS.
C
C     ------------------------------------------------------------------

      INTEGER           I

C     ==================================================================

      DO 100 I = 1, N
         INVP(PERM(I)) = I
 100  CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRT11
      END


      SUBROUTINE   SPRT12   ( NNODS , NELS  , NELNOD, XELNOD, ELNOD ,
     $                        XNODEL, NODEL , PERM  , INVP  , PARENT,
     $                        ANCSTR, MARKER                          )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER           NNODS , NELS  , NELNOD
      INTEGER           XELNOD(NELS+1)        , ELNOD(NELNOD)         ,
     $                  XNODEL(NNODS+1)       , NODEL(NELNOD)         ,
     $                  PERM(NNODS)           , INVP(NNODS)           ,
     $                  PARENT(NNODS)         , ANCSTR(NNODS)         ,
     $                  MARKER(NNODS)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT12
C
C --- TO COMPUTE THE ELIMINATION TREE.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : Jun. 01, 2020 (KMO)
C                 Removed uneccessary initialization of MARKER.
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF NODES, NOT ALTERED.
C     NELS   : INTEGER
C              NUMBER OF ELEMENTS IN THE MESH.
C              NOT ALTERED BY THE ROUTINE.
C     NELNOD : INTEGER
C              SIZE OF ELNOD$NODEL.
C              NOT ALTERED BY THE ROUTINE.
C     XELNOD : INTEGER(NELS+1)
C              POINTER INTO LISTS OF NODES FOR EACH ELEMENT.
C              NOT ALTERED BY THE ROUTINE.
C     ELNOD  : INTEGER(NELNOD)
C              ARRAY THAT HOLDS LISTS OF NODES FOR THE ELEMENTS.
C              NOT ALTERED BY THE ROUTINE.
C     XNODEL : INTEGER(NNODS+1)
C              POINTER INTO LISTS OF ELEMENTS FOR EACH NODE.
C              NOT ALTERED BY THE ROUTINE.
C     NODEL  : INTEGER(NELNOD)
C              ARRAY THAT HOLDS LISTS OF ELEMENTS FOR THE NODES.
C              NOT ALTERED BY THE ROUTINE.
C     (PERM,
C      INVP) : 2*INTEGER(NNODS)
C              THE PERMUTATION.
C              NOT ALTERED BY THE ROUTINE.
C
C     ON EXIT
C     -------
C
C     PARENT : INTEGER(NNODS)
C              THE PARENT REPRESENTATION OF THE ELIMINATION TREE.
C
C     WORKING ARRAYS
C     --------------
C
C     ANCSTR : INTEGER(NNODS)
C              A VECTOR WHICH STORES ANCESTOR RELATIONS SO THAT
C              PATH COMPRESSION IS IN EFFECT.
C     MARKER : INTEGER(NNODS)
C              IN EACH STEP, MARKS NODES WHICH ARE USED TO AVOID
C              ACCESSING A NODE MORE THAN ONCE.
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------
C
      INTEGER           ELEMNT, I     , J     , JSTOP , JSTRT , K     ,
     $                  KSTRT , KSTOP , NABOR , NEXT  , NODE
C
C     ==================================================================
C
      IF ( NNODS .LE. 0 )
     $   RETURN
C
C     ==================================================================
C
      DO 400 I = 1, NNODS
         PARENT(I) = 0
         ANCSTR(I) = 0
         MARKER(I) = I
         NODE = PERM(I)
         JSTRT = XNODEL(NODE)
         JSTOP = XNODEL(NODE+1) - 1
         DO 300 J = JSTRT, JSTOP
            ELEMNT = NODEL(J)
            KSTRT = XELNOD(ELEMNT)
            KSTOP = XELNOD(ELEMNT+1) - 1
            DO 200 K = KSTRT, KSTOP
               NABOR = INVP(ELNOD(K))
               IF ( NABOR .LT. I         .AND.
     $              MARKER(NABOR) .NE. I       )
     $         THEN
                  MARKER(NABOR) = I
 100              IF ( ANCSTR(NABOR) .EQ. I )
     $               GOTO 200
                  IF ( ANCSTR(NABOR) .GT. 0 )
     $            THEN
                     NEXT = ANCSTR(NABOR)
                     ANCSTR(NABOR) = I
                     NABOR = NEXT
                     GOTO 100
                  ENDIF
                  PARENT(NABOR) = I
                  ANCSTR(NABOR) = I
               ENDIF
 200        CONTINUE
 300     CONTINUE
 400  CONTINUE
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT12
      END


      SUBROUTINE   SPRT13   ( NNODS , PARENT, CHILD , SIBL   )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER           NNODS
      INTEGER           PARENT(NNODS)         ,
     $                  CHILD(NNODS)          , SIBL(NNODS)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT13
C
C --- TO COMPUTE THE FIRST CHILD-SIBLINGS VECTORS OF THE ELIMINATION
C     TO FACLITATE FAST POSTORDERING OF THE ELIMINATION TREE.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     PARENT : INTEGER(NNODS)
C              THE ELIMINATION TREE PARENT VECTOR.
C
C     ON EXIT
C     -------
C
C     CHILD  : INTEGER(NNODS)
C              THE FIRST CHILD VECTOR.
C     SIBL   : INTEGER(NNODS)
C              THE SIBLINGS OF THE CHILD.
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------
C
      INTEGER           LROOT , NODE  , NDPAR
C
C     ==================================================================
C
      IF ( NNODS .LE. 0 )
     $   RETURN
C
      DO 100 NODE = 1, NNODS
         CHILD(NODE) = 0
         SIBL(NODE) = 0
  100 CONTINUE
C
      IF ( NNODS .LE. 1 )
     $   RETURN
C
      LROOT = NNODS
C
C     ==================================================================
C
C     ------------------------------------------------------------
C     FOR EACH NODE := NNODS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING.
C     ------------------------------------------------------------
C
      DO 300 NODE = NNODS-1, 1, -1
         NDPAR = PARENT(NODE)
         IF ( NDPAR .LE. 0  .OR.  NDPAR .EQ. NODE )
     $   THEN
C
C           -------------------------------------------------
C           NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
C           SET NODE TO BE ONE OF THE ROOTS OF THE TREES.
C           -------------------------------------------------
C
            SIBL(LROOT) = NODE
            LROOT = NODE
         ELSE
C
C           -------------------------------------------
C           OTHERWISE, BECOMES FIRST SON OF ITS PARENT.
C           -------------------------------------------
C
            SIBL(NODE) = CHILD(NDPAR)
            CHILD(NDPAR) = NODE
         ENDIF
  300 CONTINUE
C
      SIBL(LROOT) = 0
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT13
      END


      SUBROUTINE   SPRT14   ( ROOT  , NNODS , CHILD , SIBL  , PARENT,
     $                        NEWINV, STACK                           )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER           ROOT  , NNODS
      INTEGER           CHILD(NNODS)          , SIBL(NNODS)           ,
     $                  PARENT(NNODS)         , NEWINV(NNODS)         ,
     $                  STACK(NNODS)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT14
C
C --- TO POSTORDER THE ELIMINATION TREE AND TO UPDATE THE PARENT
C     VECTOR ACCORDING TO THE EQUIVALENT POSTORDER. THE POSTORDERING
C     IS RETURNED IN NEWINV.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : Oct. 03, 2007 (KMO)
C                 Return if NNODS is zero to avoid crash.
C
C
C     ON ENTRY
C     --------
C
C     ROOT   : INTEGER
C              THE FIRST ROOT OF THE ELIMINATION TREE.
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     CHILD  : INTEGER(NNODS)
C              THE FIRST CHILD VECTOR.
C     SIBL   : INTEGER(NNODS)
C              THE SIBLING VECTOR.
C
C     UPDATED PARAMETERS
C     ------------------
C
C     PARENT : INTEGER(NNODS)
C              THE PARENT REPRESENTATION OF THE ELIMINATION TREE.
C
C     ON EXIT
C     -------
C
C     NEWINV : INTEGER(NNODS)
C              THE NEW INVERSE PERMUTATION VECTOR FOR THE POSTORDERING.
C
C     WORKING ARRAYS
C     --------------
C
C     STACK  : INTEGER(NNODS)
C              THE STACK FOR POSTORDER TRAVERSAL OF THE TREE.
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------
C
      INTEGER           ITOP  , NDPAR , NODE  , NUM   , NUNODE
C
C     ==================================================================
C
      IF ( NNODS .LE. 0 ) RETURN

      NUM = 0
      ITOP = 0
      NODE = ROOT
C
C     -------------------------------------------------------------
C     TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES
C     ALONG THE TRAVERSAL INTO THE STACK.
C     -------------------------------------------------------------
C
  100 CONTINUE
         ITOP = ITOP + 1
         STACK(ITOP) = NODE
         NODE = CHILD(NODE)
         IF ( NODE .GT. 0 )
     $      GOTO 100
C
C        ----------------------------------------------------------
C        IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT.
C        ----------------------------------------------------------
C
  200    CONTINUE
            IF ( ITOP .LE. 0 )
     $         GOTO 300
            NODE = STACK(ITOP)
            ITOP = ITOP - 1
            NUM = NUM + 1
            NEWINV(NODE) = NUM
C
C           ----------------------------------------------------
C           THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE.
C           ----------------------------------------------------
C
            NODE = SIBL(NODE)
            IF ( NODE .LE. 0 )
     $         GOTO 200
         GOTO 100
C
  300 CONTINUE
C
C     ------------------------------------------------------------
C     DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  SIBL
C     IS USED TEMPORARILY FOR THE NEW PARENT VECTOR.
C     ------------------------------------------------------------
C
      DO 400 NODE = 1, NUM
         NUNODE = NEWINV(NODE)
         NDPAR = PARENT(NODE)
         IF ( NDPAR .GT. 0 )
     $      NDPAR = NEWINV(NDPAR)
         SIBL(NUNODE) = NDPAR
  400 CONTINUE
C
      DO 500 NUNODE = 1, NUM
         PARENT(NUNODE) = SIBL(NUNODE)
  500 CONTINUE
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT14
      END


      SUBROUTINE   SPRT15   ( NNODS , NEWINV, PERM  , INVP   )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER           NNODS
      INTEGER           NEWINV(NNODS)         ,
     $                  PERM(NNODS)           , INVP(NNODS)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT15
C
C --- TO UPDATE PERM,INVP WITH THE PERMUTATION INCREMENT FROM THE
C     EQUIVALENT POSTORDERING STORED IN NEWINV.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     NEWINV : INTEGER(NNODS)
C              THE INVERSE PERMUTATION INCREMENT.
C
C     ON EXIT
C     -------
C
C     PERM   : INTEGER
C              THE UPDATED PERMUTATION VECTOR.
C     INVP   : INTEGER
C              THE UPDATED INVERSE PERMUTATION VECTOR.
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------
C
      INTEGER           I     , INTERM, NODE
C
C     ==================================================================
C
      DO 100 I = 1, NNODS
         INTERM = INVP(I)
         INVP(I) = NEWINV(INTERM)
  100 CONTINUE
C
      DO 200 I = 1, NNODS
         NODE = INVP(I)
         PERM(NODE) = I
  200 CONTINUE
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT15
      END


      SUBROUTINE   SPRT31   ( NNODS , PARENT, XSIBL , SIBL  , MASK    )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER           NNODS
      INTEGER           PARENT(NNODS)         , MASK(NNODS)           ,
     $                  XSIBL(NNODS+1)        , SIBL(NNODS+1)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT31
C
* --- TO COMPUTE THE ADJACENCY SET REPRESENTATION OF THE ELIMINATION
C     TREE REPRESENTED BY PARENT. THE ROUTINE HANDLES MORE THAN ONE
C     COMPONENT IF PARENT HAS ZERO OR NEGATIVE VALUES FOR TREE ROOTS
C     IN THE FORREST. THE ROOTS ARE STORED LAST IN SIBL WITH NEGATIVE
C     SIGN, AND THEY MAY EASILY BE COUNTED IN ROUTINES WHICH USES THE
C     RESULTS.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF DISTINCT VERTICES IN THE ELIMINATION TREE
C     PARENT : INTEGER(NNODS)
C              REPRESENTATION OF THE ELIMINATION TREE
C
C     ON EXIT
C     -------
C
C     (XSIBL,
C      SIBL) : INTEGER(NNODS+1), INTEGER(NNODS)
C              ADJACENCY REPRESENTATION OF THE ELIMINATION TREE.
C              THE SUPERNODE TREE ADJACENCY SET LISTS.
C              SIBL(NNODS+1) STORES THE NUMBER OF ROOTS IN THE
C                             SUPERNODE ELIMINATION TREE.
C              IF IP=SIBL(NNODS+1) THEN THE ROOTS ARE FOUND IN THE
C              LOCATIONS SIBL(NNODS-IP+1), ... SIBL(NNODS).
C
C     WORKING ARRAYS
C     --------------
C
C     MASK   : INTEGER(NNODS)
C              NUMBER OF CHILDREN OF EACH VERTEX AND MASK OF ORDERED
C              VERTICES
C
C     SUBROUTINES USED
C     ----------------
C
C     NONE
C
C     ------------------------------------------------------------------
C
      INTEGER           FATHER, I     , IPSIBL, NROOTS, ROOT
C
C     ==================================================================
C
      DO 10 I = 1, NNODS
         MASK(I) = 0
         SIBL(I) = 0
   10 CONTINUE
C
C     ==================================================================
C
C     ---------------------------
C     COMPUTE THE TREE ADJACENCY.
C     ---------------------------
C
      DO 100 I = 1, NNODS
         FATHER = PARENT(I)
         IF ( FATHER .GT. 0 ) THEN
            MASK(FATHER) = MASK(FATHER) + 1
         ENDIF
 100  CONTINUE
C
      XSIBL(1) = 1
C
      DO 200 I = 1, NNODS
         XSIBL(I+1) = XSIBL(I) + MASK(I)
 200  CONTINUE
C
      DO 210 I = 1, NNODS
         MASK(I) = 0
  210 CONTINUE
C
      DO 300 I = 1, NNODS
         FATHER = PARENT(I)
         IF ( FATHER .GT. 0 ) THEN
            SIBL(XSIBL(FATHER) + MASK(FATHER)) = I
            MASK(FATHER) = MASK(FATHER) + 1
         ENDIF
 300  CONTINUE
C
C     ----------------
C     STORE THE ROOTS.
C     ----------------
C
      DO 310 I = 1, NNODS
         MASK(I) = 0
  310 CONTINUE
C
      DO 400 I = 1, NNODS
         ROOT = SIBL(I)
         IF ( ROOT .GT. 0 ) THEN
            MASK(ROOT) = ROOT
         ENDIF
  400 CONTINUE
C
      IPSIBL = NNODS
      NROOTS = 0
      DO 500 I = NNODS, 1, -1
         ROOT = MASK(I)
         IF ( ROOT .EQ. 0 ) THEN
            SIBL(IPSIBL) = -I
            IPSIBL = IPSIBL - 1
            NROOTS = NROOTS + 1
         ENDIF
  500 CONTINUE
      SIBL(NNODS+1) = NROOTS
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT31
      END


      SUBROUTINE   SPRT41   ( NL    , LIST  , NK    , KEY    )
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT41 : ORDER A LIST OF INTEGERS IN DECREASING SEQUENCE OF THEIR
C              KEYS USING INSERTION SORT
C
C     WRITTEN   : AUG. 12, 1990 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C             - NL     : INTEGER
C                        LENGTH OF LIST
C             - LIST   : INTEGER(NL)
C                        A LIST OF INTEGERS
C             - NK     : INTEGER
C                        LENGTH OF KEY, NOTE!! NK GE NL
C             - KEY    : INTEGER(NK)
C                        A LIST OF INTEGER KEYS
C
C     ON EXIT
C     -------
C
C             - LIST   : INTEGER(NL)
C                        A LIST OF INTEGERS SORTED IN DECREASING
C                        SEQUENCE OF KEY.
C
C     ------------------------------------------------------------------
C
C     ----------
C     PARAMETERS
C     ----------
C
      IMPLICIT NONE

      INTEGER           NL    , NK    , I     , J     , T     , VALUE,
     $                  LIST(NL)      , KEY(NK)
C
      DO 20 I = 2, NL
         T     = LIST(I)
         VALUE = KEY(T)
         DO 10 J = I - 1, 1, -1
            IF ( VALUE .LE. KEY(LIST(J)) ) THEN
               LIST(J+1) = T
               GOTO 20
            ENDIF
            LIST(J+1) = LIST(J)
 10      CONTINUE
         LIST(1) = T
 20   CONTINUE
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT41
      END


      SUBROUTINE   SPRT42   ( IP    , JS    , WSTORE, SIBL  , DELTA ,
     $                        INDEX , WORK                            )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER           IP    , JS
      INTEGER           WSTORE(*)             , SIBL(*)               ,
     $                  DELTA(*)              , INDEX(*)              ,
     $                  WORK(*)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT42
C
* --- TO COMPUTE WORKING STACK STORAGE FOR A NON-TRIVIAL SUPERNODE.
C     METHOD OF YANG.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     IP     : INTEGER
C              NUMBER OF CHILDREN OF THE SUPERNODE.
C     JS     : INTEGER
C              CURRENT PARENT SUPERNODE.
C     WSTORE : INTEGER(NSUPER)
C              STORAGE FOR SUPERNODES NUMBERED LOWER THAN JS.
C     SIBL   : INTEGER(IP)
C              LIST OF CHILDREN NODES OF JS.
C     DELTA  : INTEGER(NSUPER)
C              SIZE OF UPDATE MATRICES.
C     INDEX  : INTEGER(IP)
C              OPTIMAL CHILDREN SEQUENCE INDEX VECTOR.
C
C     ON EXIT
C     -------
C
C     WSTORE(JS) : INTEGER
C                  COMPUTED RESULTAT STACK STORAGE FOR SUPERNODE JS.
C
C     WORKING ARRAYS
C     --------------
C
C     WORK   : INTEGER(IP)
C              TEMPORARY STORAGE OF OPTIMAL SEQUENCE.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     SPRCPY
C
C     INTRINSIC
C     ---------
C
C     MAX
C
C     ------------------------------------------------------------------
C
      INTEGER           CHILD , I     , SUMDLT, WORKST
      INTRINSIC         MAX
      EXTERNAL          SPRCPY
C
C     ==================================================================
C
      IF ( IP .LE. 0 ) THEN
         RETURN
      ELSE
         WORKST = 0
         SUMDLT = 0
         DO 100 I = 1, IP
            CHILD  = SIBL(INDEX(I))
            WORK(I) = CHILD
            WORKST = MAX( WORKST, WSTORE(CHILD) + SUMDLT )
            SUMDLT = SUMDLT + DELTA(CHILD)
 100     CONTINUE
         WSTORE(JS) = MAX( DELTA(JS) + SUMDLT, WORKST )
      ENDIF
C
      IF ( IP .GT. 1 ) THEN
         CALL SPRCPY ( IP, WORK, 1, SIBL, 1 )
      ENDIF
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT42
      END


      SUBROUTINE   SPRT43   ( NNODS , XSIBL , SIBL  , PERM  , STACK,
     $                        MASK                                    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NNODS
      INTEGER           PERM(NNODS)           ,
     $                  STACK(NNODS)          , MASK(NNODS)           ,
     $                  XSIBL(NNODS+1)        , SIBL(NNODS+1)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT43
C
* --- TO FIND A POSTORDERING OF THE SUPERNODE ELIMINATION TREE REPRESEN-
C     TED BY (XSIBL,SIBL) THE POSTORDER OF THE SUPERNODE ELIMINATION
C     TREE IS RETURNED IN VECTOR PERM. THE ROUTINE ALSO HANDLES THE CASE
C     OF MORE THAN ONE CONNECTED COMPONENT. NOTE THAT (XSIBL,SIBL) MUST
C     BE LEFT UNALTERED AFTER ITS COMPUTATION IN SPRT31.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF DISTINCT VERTICES IN THE ELIMINATION TREE
C     (XSIBL,
C      SIBL) : INTEGER(NNODS+1), INTEGER(NNODS+1)
C              ADJACENCY REPRESENTATION OF THE SUPERNODE ELIMINATION
C              TREE.
C              SIBL(NNODS+1) STORES THE NUMBER OF ROOTS IN THE
C                             SUPERNODE ELIMINATION TREE.
C              IF IP=SIBL(NNODS+1) THEN THE ROOTS ARE FOUND IN THE
C              LOCATIONS SIBL(NNODS-IP+1), ... SIBL(NNODS).
C
C     ON EXIT
C     -------
C
C     PERM   : INTEGER(NNODS)
C              POSTORDER OF THE SUPERNODE ELIMINATION TREE VERTICES
C
C     WORKING ARRAYS
C     --------------
C
C     MASK   : INTEGER(NNODS)
C              NUMBER OF CHILDREN OF EACH VERTEX AND MASK OF ORDERED
C              VERTICES
C     STACK  : INTEGER(NNODS)
C              STORAGE OF PATH FROM ROOT TO CURRENT NODE
C
C     SUBROUTINES USED
C     ----------------
C
C     SPRT45 : FIND RIGHT SIBLING OF A VERTEX
C     SPRT46 : FIND LEFTMOST UNEXPLORED CHILD OF A VERTEX
C
C     ------------------------------------------------------------------

      INTEGER           BROTHR, FATHER, IPSIBL, LCHILD, NUMBER, RSIBL ,
     $                  SPOINT

      EXTERNAL          SPRT45, SPRT46
      INTRINSIC         ABS

C     ==================================================================

      DO 10 NUMBER = 1, NNODS
         MASK(NUMBER) = 0
   10 CONTINUE

      NUMBER = 1
      SPOINT = 0
      IPSIBL = NNODS

C     ==================================================================

C     -------------------------------------------------------
C     UTILIZE THAT SPRT31 STORES THE ROOTS AT THE END OF SIBL
C     WITH NEGATIVE SIGN AND TRAVERSE EACH OF THE TREES.
C     -------------------------------------------------------

 100  IF ( IPSIBL .GT. 0 ) THEN

         FATHER = SIBL(IPSIBL)

C        ---------------------
C        FATHER IS INITIALIZED
C        ---------------------

         IF ( FATHER .LT. 0 ) THEN

C           -----------------------------------
C           FATHER IS A ROOT, DECREMENT POINTER
C           INTO SIBL AND SET FATHER POSITIVE
C           -----------------------------------

            IPSIBL = IPSIBL - 1
            FATHER = ABS(FATHER)

  200       IF ( FATHER .GT. 0 ) THEN

C              ---------------------------------
C              CURRENT TREE IS NOT YET EXHAUSTED
C              ---------------------------------

               SPOINT = SPOINT + 1
               STACK(SPOINT) = FATHER
               MASK(FATHER) = FATHER

C              --------------------------------------------
C              FIND THE LEFTMOST UNEXPLORED CHILD OF FATHER
C              --------------------------------------------

               CALL SPRT46 ( FATHER, XSIBL, SIBL, MASK, LCHILD )

C              ------------------------------------------
C              SET BROTHR TO BE FATHER AND SET FATHER TO
C              BE THE LEFTMOST CHILD IN ORDER TO CONTINUE
C              DOWN IN THE CURRENT SUBTREE OF THE TREE.
C              ------------------------------------------

               BROTHR = FATHER
               FATHER = LCHILD

            ENDIF
            IF ( FATHER .EQ. 0 ) THEN

C              ---------------------------------------
C              WE HAVE REACHED THE BOTTOM OF A SUBTREE
C              ---------------------------------------

  300          IF ( SPOINT .EQ. 0 ) THEN

C                 -------------------------------------
C                 NO MORE UNEXPLORED NODES IN THIS TREE
C                 -------------------------------------

                  GOTO 100
               ENDIF

C              -------------------------------------------------
C              FIND THE THE BROTHER OF BROTHR TO THE RIGHT OF IT
C              -------------------------------------------------

               CALL SPRT45 ( BROTHR, SPOINT, STACK, XSIBL,
     $                       SIBL, MASK, RSIBL )

C              ------------------------
C              POP A NODE OFF THE STACK
C              SINCE IT IS NECESSARY TO
C              GO UPWARDS AGAIN
C              ------------------------

               SPOINT = SPOINT - 1
               IF ( RSIBL .GT. 0 ) THEN
                  FATHER = RSIBL
                  PERM(NUMBER) = STACK(SPOINT+1)
                  NUMBER = NUMBER + 1
               ELSE
                  IF ( SPOINT .GT. 0 ) THEN
                     BROTHR = STACK(SPOINT)
                  ENDIF
                  PERM(NUMBER) = STACK(SPOINT+1)
                  NUMBER = NUMBER + 1
                  GOTO 300
               ENDIF
            ENDIF
            GOTO 200
         ENDIF
      ENDIF

C     ==================================================================

C     ------------------------------------------------------------------
C     END OF MODULE SPRT43
      END


      SUBROUTINE   SPRT44   ( SEPSZE, NNODS , NSUPER, XSUPER, XSIBL ,
     $                        SIBL  , PERM  , INVP  , SUPSZE, SUPLEN,
     $                        SUPSUP, SUPERM, NWPERM, NXSUP , NXSIBL,
     $                        NSIBL , SUINVP, ERROR , LPU             )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           SEPSZE, NNODS , NSUPER, ERROR , LPU
      INTEGER           XSUPER(NSUPER+1)      , XSIBL(NSUPER+1)       ,
     $                  SIBL(NSUPER+1)        , PERM(NNODS)           ,
     $                  INVP(NNODS)           , SUPSZE(NSUPER)        ,
     $                  SUPLEN(NSUPER)        , SUPSUP(NSUPER)        ,
     $                  SUPERM(NSUPER)        , NWPERM(NNODS)         ,
     $                  NXSUP(NSUPER+1)       , NXSIBL(NSUPER+1)      ,
     $                  NSIBL(NSUPER+1)       , SUINVP(NSUPER)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT44
C
* --- TO UPDATE THE PERMUTATION VECTORS AND THE SUPERNODE ELIMINATION
C     TREE AFTER A POSTORDER OF THE SUPERNODE ELIMINATION TREE IS DONE.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
*     SEPSZE : Integer
*              The size of the separator. Is always zero if this is not
*              a superelement or if superelement info. is suppressed.
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     XSUPER : INTEGER(NSUPER+1)
C              THE FUNDAMENTAL SUPERNODE PARTITION.
C     XSIBL  : INTEGER(NSUPER+1)
C              POINTER INTO THE SUPERNODE ELIMINATION TREE
C              ADJACENCY SET.
C     SIBL   : INTEGER(NSUPER)
C              LIST OF CHILDREN FOR EACH SUPERNODE.
C     PERM   : INTEGER(NNODS)
C              NEW TO OLD PERMUTATION.
C     INVP   : INTEGER(NNODS)
C              OLD TO NEW PERMUTATION.
C     SUPSZE : INTEGER(NSUPER)
C              THE SIZE OF EACH SUPERNODE.
C     SUPLEN : INTEGER(NSUPER)
C              THE LENGTH OF EACH SUPERNODE RELATED TO VARIABLES.
C     SUPSUP : INTEGER(NSUPER)
C              Array to control the processing in numerical steps.
C     SUPERM : INTEGER(NSUPER)
C              NEW TO OLD PERMUTATION OF THE SUPERNODES.
C
C     ON EXIT
C     -------
C
C     ERROR  : INTEGER
C              =  0 : NORMAL RETURN.
C              = -1 : ERROR FROM INTERNAL SORT ROUTINE.
C     XSUPER : INTEGER(NSUPER+1)
C              THE FUNDAMENTAL SUPERNODE PARTITION.
C              UPDATED WITH POSTORDER INFORMATION.
C     XSIBL  : INTEGER(NSUPER+1)
C              POINTER INTO THE SUPERNODE ELIMINATION TREE
C              ADJACENCY SET.
C              UPDATED WITH POSTORDER INFORMATION.
C     SIBL   : INTEGER(NSUPER)
C              LIST OF CHILDREN FOR EACH SUPERNODE.
C              UPDATED WITH POSTORDER INFORMATION.
C     PERM   : INTEGER(NNODS)
C              NEW TO OLD PERMUTATION.
C              UPDATED WITH POSTORDER INFORMATION.
C     INVP   : INTEGER(NNODS)
C              OLD TO NEW PERMUTATION.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPSZE : INTEGER(NSUPER)
C              THE SIZE RELATED TO VARIABLES FOR EACH SUPERNODE.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPLEN : INTEGER(NSUPER)
C              THE LENGTH RELATED TO VARIABLES FOR EACH SUPERNODE.
C              UPDATED WITH POSTORDER INFORMATION.
C     SUPSUP : INTEGER(NSUPER)
C              Array to control the processing in numerical steps.
C              UPDATED WITH POSTORDER INFORMATION.
C
C     WORKING ARRAYS
C     --------------
C
C     NWPERM : INTEGER(NNODS)
C              PERMUTATION INCREMENT DUE TO SUPERNODE POSTORDER.
C     NXSUP  : INTEGER(NSUPER+1)
C              HOLDS THE UPDATED SUPERNODE PARTITION.
C     NXSIBL : INTEGER(NSUPER+1)
C              HOLDS THE UPDATED TREE ADJACENCY POINTERS
C     NSIBL  : INTEGER(NSUPER+1)
C              HOLDS THE UPDATED CHILDREN LISTS FOR EACH SUPERNODE.
C     SUINVP : INTEGER(NSUPER)
C              OLD TO NEW PERMUTATION OF THE SUPERNODES.
C
C     SUBPROGRAMS
C     -----------
C
C     SPRT11
C     SPRCPY
C     SPRT48
C     SPRT47
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     ABS
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER
     $  FATHER, I     , IP1   , IP2   , ISTRT , ISTOP , JS    , OLDSUP,
     $  NROOTS, RETNOD, ROOT

      LOGICAL           LERROR

      INTRINSIC         ABS

      EXTERNAL          SPRT11, SPRCPY, SPRT48, SPRT47

C     ==================================================================

      ERROR  = 0
      LERROR = .FALSE.

C     ==================================================================

C     COMPUTE THE INVERSE SUPERNODE PERMUTATION.
      CALL SPRT11 ( NSUPER, SUPERM, SUINVP )

C     ==================================================================

C     UPDATE XSUPER, XSIBL, SIBL AND CONSTRUCT NWPERM.
      IP1 = 1
      IP2 = 1
      DO 300 JS = 1, NSUPER
         NXSUP(JS) = IP1
         NXSIBL(JS) = IP2
         OLDSUP = SUPERM(JS)

         DO 100 I = XSUPER(OLDSUP), XSUPER(OLDSUP+1) - 1
            NWPERM(IP1) = I
            IP1 = IP1 + 1
  100    CONTINUE

         DO 200 I = XSIBL(OLDSUP), XSIBL(OLDSUP+1) - 1
            NSIBL(IP2) = SUINVP(SIBL(I))
            IP2 = IP2 + 1
  200    CONTINUE

  300 CONTINUE

      NXSUP(NSUPER+1) = IP1
      NXSIBL(NSUPER+1) = IP2

C     UPDATE THE ROOTS.
      NROOTS = SIBL(NSUPER+1)
      DO 400 JS = NSUPER, NSUPER-NROOTS+1, -1
         ROOT = ABS(SIBL(JS))
         NSIBL(JS) = SUINVP(ROOT)
  400 CONTINUE
      NSIBL(NSUPER+1) = NROOTS
      IF ( NROOTS .GT. 1 ) THEN
         CALL SPRT48 ( NROOTS, NSIBL(NSUPER-NROOTS+1), LERROR )
         IF ( LERROR ) THEN
            ERROR = -1
            WRITE ( LPU, '(//A)' )
     $      '*** Error from SPRTRS::SPRTR4::SPRT44'
            WRITE ( LPU, '(4X,A/)' )
     $      'Error detected in SPRT48 when sorting roots'
            RETURN
         ENDIF
      ENDIF
      DO 500 JS = NSUPER, NSUPER-NROOTS+1, -1
         NSIBL(JS) = -NSIBL(JS)
  500 CONTINUE

C     COPY UPDATED ARRAYS BACK TO ORIGINAL POSITIONS.
      CALL SPRCPY ( NSUPER+1, NXSUP, 1, XSUPER, 1 )
      CALL SPRCPY ( NSUPER+1, NXSIBL, 1, XSIBL, 1 )
      CALL SPRCPY ( NSUPER+1, NSIBL, 1, SIBL, 1 )

C     UPDATE THE PERMUTATION VECTORS.
      CALL SPRT47 ( NNODS, PERM, INVP, NWPERM )

C     UPDATE SUPLEN AND SUPSZE.
      DO 600 JS = 1, NSUPER
         OLDSUP = SUPERM(JS)
         NXSUP(JS) = SUPLEN(OLDSUP)
         NXSIBL(JS) = SUPSZE(OLDSUP)
  600 CONTINUE
      CALL SPRCPY ( NSUPER, NXSUP, 1, SUPLEN, 1 )
      CALL SPRCPY ( NSUPER, NXSIBL, 1, SUPSZE, 1 )

* --- Update SUPSUP.
      IF ( SEPSZE.GT.0 ) THEN

         DO JS = 1, NSUPER
            OLDSUP = SUPERM(JS)
            NXSUP(JS) = SUPSUP(OLDSUP)
         END DO
         CALL SPRCPY ( NSUPER, NXSUP, 1, SUPSUP, 1 )

      ENDIF

* --- Do final control of the separators.
      IF ( SEPSZE.GT.0 ) THEN

         DO JS = 1, NSUPER
* --------- For control, NXSIBL will be used as Parent vector.
            NXSIBL(JS) = 0
         END DO

* ------ Set up Parent vector.
         DO JS = 1, NSUPER
            ISTRT = XSIBL(JS)
            ISTOP = XSIBL(JS+1)-1
            DO I = ISTRT, ISTOP
               NXSIBL(SIBL(I)) = JS
            END DO
         END DO

* ------ Control the final supernodal sequence.
*        NOTE - that the supernodes need to be in postorder for the
*               test to be correct.
         RETNOD = 0
         DO JS = 1, NSUPER

            FATHER = NXSIBL(JS)

* --------- If we have a retained supernode, turn on RETNOD.
            IF ( SUPSUP(JS).GT.0 ) THEN
               RETNOD = 1
            ENDIF

            IF ( FATHER.LE.0 ) THEN

* ------------ Since the supernodes are ordered in postorder this node
*              is the last in a component. If at least one node previous
*              in the component is a retained node, then also the root
*              must be a retained supernode.
               IF ( RETNOD.EQ.1 ) THEN
                  IF ( SUPSUP(JS).LE.0 ) THEN
                     ERROR = -2
                     WRITE ( LPU, '(//A)' )
     $               '*** Error from SPRTRS::SPRTR4::SPRT44'
                     WRITE ( LPU, '(4X,A,I10/)' )
     $               'Supernode JS is a root that is not a separator',JS
                     RETURN
                  ENDIF
               ENDIF

* ------------ Reset RETNOD since this is the end of a component.
               RETNOD = 0

            ELSEIF ( FATHER.GT.0 ) THEN

* ------------ This node has a parent in the elimination tree. If it is
*              a retained node then also the parent need be retained.
               IF ( SUPSUP(JS).GT.0 ) THEN
                  IF ( SUPSUP(FATHER).LE.0 ) THEN
                     ERROR = -3
                     WRITE ( LPU, '(//A)' )
     $               '*** Error from SPRTRS::SPRTR4::SPRT44'
                     WRITE ( LPU, '(4X,A,I10/)' )
     $               'Separator supernodes does not define a chain',JS
                     RETURN
                  ENDIF
               ENDIF

            ENDIF

         END DO

      ENDIF

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRT44
      END


      SUBROUTINE SPRT45 ( BROTHR, SPOINT, STACK, XSIBL, SIBL, MASK,
     $                    RSIBL )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER           BROTHR, SPOINT, RSIBL
      INTEGER           STACK(*), XSIBL(*), SIBL(*), MASK(*)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT45
C
* --- TO FIND RIGHT SIBLING RSIBL OF BROTHR. USED IN SUBROUTINE SPRT43.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     BROTHR : INTEGER
C              VERTEX OF WHICH RIGTH SIBLING IS TO BE FOUND
C     SPOINT : INTEGER
C              POINTER TO TOP OF STACK
C     STACK  : INTEGER(NNODS)
C              STACK WHICH STORES THE PATH FROM ROOT TO THE CURRENT
C              VERTEX
C     (XSIBL,
C      SIBL) : INTEGER(NNODS+1), INTEGER(NNODS-1)
C              ADJACENCY OF THE TREE STRUCTURE
C     MASK   : INTEGER(NNODS)
C              MASK OF EXPLORED VERTICES
C
C     ON EXIT
C     -------
C
C     RSIBL  : INTEGER
C              RIGTH SIBLING OF BROTHR
C
C     ------------------------------------------------------------------
C
      INTEGER           FATHER, I     , ISTART, ISTOP , CHILD , RCHILD,
     $                  J     , JSTART, JSTOP
C
      RSIBL = 0
      IF ( SPOINT .LE. 1 ) RETURN
C
      FATHER = STACK(SPOINT-1)
      ISTART = XSIBL(FATHER)
      ISTOP  = XSIBL(FATHER+1) - 1
C
      IF ( ISTART .GT. ISTOP ) RETURN
C
      DO 200 I = ISTART, ISTOP
C
         CHILD = SIBL(I)
C
         IF ( CHILD .EQ. BROTHR ) THEN
C
            IF ( I .EQ. ISTOP ) RETURN
C
               JSTART = I + 1
               JSTOP  = ISTOP
C
               DO 100 J = JSTART, JSTOP
C
                  RCHILD = SIBL(J)
C
                  IF ( MASK(RCHILD) .EQ. 0 ) THEN
C
                     RSIBL = RCHILD
                     RETURN
C
                  ENDIF
C
 100           CONTINUE
C
               RETURN
C
         ENDIF
C
 200  CONTINUE
C
      RETURN
      END


      SUBROUTINE SPRT46 ( FATHER, XSIBL, SIBL, MASK, LCHILD )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER           FATHER, LCHILD
      INTEGER           XSIBL(*), SIBL(*), MASK(*)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT46
C
* --- TO FIND LEFTMOST UNEXPLORED CHILD LCHILD OF FATHER. USED IN SUB-
C     ROUTINE SPRT43.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     FATHER : INTEGER
C              VERTIX FOR WHICH LEFTMOST CHILD IS TO BE FOUND
C     (XSIBL,
C      SIBL) : INTEGER(NNODS+1), INTEGER(NNODS-1)
C              ADJACENCY OF THE TREE STRUCTURE
C     MASK   : INTEGER(NNODS)
C              MASK OF EXPLORED VERTICES
C
C     ON EXIT
C     -------
C
C     LCHILD : INTEGER
C              LEFTMOST UNEXPLORED CHILD OF FATHER
C
C     ------------------------------------------------------------------
C
      INTEGER           CHILD , I     ,ISTART, ISTOP
C
      LCHILD = 0
      ISTART = XSIBL(FATHER)
      ISTOP  = XSIBL(FATHER+1) - 1
C
      IF ( ISTART .GT. ISTOP ) RETURN
C
      DO 100 I = ISTART, ISTOP
C
         CHILD = SIBL(I)
C
         IF ( MASK(CHILD) .EQ. 0 ) THEN
C
            LCHILD = CHILD
            RETURN
C
         ENDIF
C
 100  CONTINUE
C
      RETURN
      END


      SUBROUTINE   SPRT47   ( NNODS , PERM  , INVP  , NEWPRM )
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE

      INTEGER           NNODS
      INTEGER           PERM(NNODS)           , INVP(NNODS)           ,
     $                  NEWPRM(NNODS)
C
C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRT47
C
* --- TO UPDATE THE PERMUTATION VECTOR, WHICH GIVES THE CONTACT FROM
C     LAST PERMUTATION ORDER TO THE INITIAL ORDER OF THE NODES, WITH THE
C     NEW PERMUTATION DEFINED FROM PERM TO NEWPRM. NEWPRM IS A PERMU-
C     TATION WHERE PERM DEFINES THE INITIAL STATE.
C
C     *** N O T E ***
C     ON INPUT THE ROUTINE REQUIRES THAT PERM CONTAINS THE INITIAL
C     PERMUTATION VECTOR AND THAT NEWPRM CONTAINS THE PERMUTATION
C     INCREMENT.
C     BUT THE VECTOR INVP NEEDS NOT CONTAIN THE INITIAL INVERSE PERMU-
C     TATION SINCE IT IS RESTATET BY MEANS OF THE UPDATED PERMUTATION
C     VECTOR.
C
C     CREATED   : May  01, 1997 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF VERTICES
C     PERM   : INTEGER(NNODS)
C              REPRESENTS THE INITIAL PERMUTATION VECTOR
C     NEWPRM : INTEGER(NNODS)
C              DEFINES THE PERMUTATION INCREMENT FROM OLD
C              REPRESENTATION TO NEW REPRESENTATION.
C
C     ON EXIT
C     -------
C
C     (PERM,
C     INVP)  :  INTEGER(NNODS), INTEGER(NNODS)
C               NEW REPRESENTATION OF THE PERMUTATION.
C
C     ------------------------------------------------------------------
C
      INTEGER           I
C
C     ==================================================================
C
      DO 100 I = 1, NNODS
         INVP(I) = PERM(NEWPRM(I))
 100  CONTINUE
C
      DO 200 I = 1, NNODS
         PERM(I) = INVP(I)
 200  CONTINUE
C
      DO 300 I = 1, NNODS
         INVP(PERM(I)) = I
 300  CONTINUE
C
C     ==================================================================
C
      RETURN
C
C     ------------------------------------------------------------------
C     END OF MODULE SPRT47
      END


      SUBROUTINE   SPRT48   ( N, KEY, ERROR )
C                  SPRT48 -- SORT INTEGERS INTO ASCENDING ORDER

      IMPLICIT NONE

      INTEGER                 N
      LOGICAL                 ERROR
      INTEGER                 KEY (N)
      INTEGER                 STKLEN, TINY
      PARAMETER               ( STKLEN = 50,
     $                          TINY   =  9 )
      INTEGER                 I     , IP1   , J     , JM1   , K     ,
     $                        LEFT  , LLEN  , RIGHT , RLEN  , TOP   ,
     $                        V
      LOGICAL                 DONE
      INTEGER                 STACK (STKLEN)
      ERROR = .FALSE.
      IF  (N .EQ. 1)  RETURN
      IF  (N .LT. 1)  GO TO 6000
      TOP = 1
      LEFT = 1
      RIGHT = N
      DONE = (N .LE. TINY)
  100 IF  (DONE)  GO TO 2000
          K = KEY ((LEFT+RIGHT)/2)
          KEY ((LEFT+RIGHT)/2) = KEY (LEFT)
          KEY (LEFT) = K
          IF  ( KEY(LEFT+1) .LE. KEY(RIGHT) ) GO TO 200
              K = KEY (LEFT+1)
              KEY (LEFT+1) = KEY (RIGHT)
              KEY (RIGHT) = K
  200     IF  ( KEY(LEFT) .LE. KEY(RIGHT) )  GO TO 300
              K = KEY (LEFT)
              KEY (LEFT) = KEY (RIGHT)
              KEY (RIGHT) = K
  300     IF  ( KEY (LEFT+1) .LE. KEY (LEFT) )  GO TO 400
              K = KEY (LEFT+1)
              KEY (LEFT+1) = KEY (LEFT)
              KEY (LEFT) = K
  400     V = KEY (LEFT)
          I = LEFT+1
          J = RIGHT
  500     CONTINUE
  600         I  = I + 1
              IF  ( KEY(I) .LT. V )  GO TO 600
  700         J = J - 1
              IF  ( KEY(J) .GT. V )  GO TO 700
          IF  (J .LT. I)  GO TO 800
              K = KEY (I)
              KEY (I) = KEY (J)
              KEY (J) = K
          GO TO 500
  800     K = KEY (LEFT)
          KEY (LEFT) = KEY (J)
          KEY (J) = K
          LLEN = J-LEFT
          RLEN = RIGHT - I + 1
          IF  ( MAX0 (LLEN, RLEN) .GT. TINY )  GO TO 1100
              IF  (TOP .EQ. 1)  GO TO 900
                  TOP = TOP - 2
                  LEFT = STACK (TOP)
                  RIGHT = STACK (TOP+1)
                  GO TO 1000
  900             DONE = .TRUE.
 1000             GO TO 1700
 1100     IF  (MIN0 (LLEN, RLEN) .GT. TINY)  GO TO 1400
              IF  ( LLEN .GT. RLEN )  GO TO 1200
                  LEFT = I
                  GO TO 1300
 1200             RIGHT = J - 1
 1300         GO TO 1700
 1400     IF  ( TOP+1 .GT. STKLEN )  GO TO 6000
          IF  ( LLEN .GT. RLEN )  GO TO 1500
              STACK (TOP) = I
              STACK (TOP+1) = RIGHT
              RIGHT = J-1
              GO TO 1600
 1500         STACK (TOP) = LEFT
              STACK (TOP+1) = J-1
              LEFT = I
 1600     TOP = TOP + 2
 1700 GO TO 100
 2000 I = N - 1
      IP1 = N
 2100     IF  ( KEY (I) .LE. KEY (IP1) )  GO TO 2400
              K = KEY (I)
              J = IP1
              JM1 = I
 2200             KEY (JM1) = KEY (J)
                  JM1 = J
                  J = J + 1
                  IF  ( J .GT. N )  GO TO 2300
                  IF  (KEY (J) .LT. K)  GO TO 2200
 2300         KEY (JM1) = K
 2400     IP1 = I
          I = I - 1
          IF  ( I .GT. 0 )  GO TO 2100
      RETURN
 6000 ERROR = .TRUE.
      RETURN
      END


      SUBROUTINE   SPREQN   ( SPRNOD, NEQ   , NDOF  , PERM  , XNODES,
     $                        NODES , MEQN                            )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           SPRNOD, NEQ   , NDOF
      INTEGER           PERM(SPRNOD)          , XNODES(SPRNOD+1)      ,
     $                  NODES(NEQ)            , MEQN(NDOF)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPREQN
C
C --- TO UPDATE MEQN GIVEN A PERMUTATION OF THE NODES HELD IN PERM.
C     NOTE THAT THERE ARE ONLY THE POSITIVE INTEGERS IN MEQN THAT
C     CORRESPOND TO FREE VARIABLES IN THE FINITE-ELEMENT MESH THAT ARE
C     UPDATED.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     SPRNOD - INTEGER
C      ENTRY : NUMBER OF SPR-NODES IN THE FINITE-ELEMENT MESH.
C      EXIT  : NOT CHANGED.
C     NEQ - INTEGER
C      ENTRY : NUMBER OF FREE VARIABLES, I.E. EQUATIONS, IN THE
C              FINITE-ELEMENT MESH.
C      EXIT  : NOT CHANGED.
C     NDOF - INTEGER
C      ENTRY : NUMBER OF VARIABLES IN THE FINITE-ELEMENT MESH.
C      EXIT  : NOT CHANGED.
C     PERM - INTEGER(SPRNOD)
C      ENTRY : A PERMUTATION OF THE SPR-NODES, MAY BE IDENTITY.
C      EXIT  : NOT CHANGED.
C     XNODES - INTEGER(SPRNOD+1)
C      ENTRY : THE SPR-NODE PARTITION, I.E. MADOF FOR SPR.
C      EXIT  : NOT CHANGED.
C     NODES - INTEGER(NEQ)
C      ENTRY : THE DOF NUMBER, I.E. POINTER INTO MEQN, FOR EACH FREE
C              VARIABLE IN THE FINITE-ELEMENT MESH.
C      EXIT  : NOT CHANGED.
C     MEQN - INTEGER(NDOF)
C      ENTRY : FOR ZERO OR NEGATIVE VALUES, AS DEFINED IN SAM-LIBRARY,
C              THE POSITIVE VALUES CORRESPOND TO EQUATION NUMBERS.
C      EXIT  : THE POSITIVE VALUES ARE UPDATED ACCORDING TO PERM.
C
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C      ENTRY : NOT DEFINED.
C      EXIT  : NEED NOT BE SAVED.
C
C     PROCEDURES
C     ----------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           I     , IDOF  , J     , NODE  , VAR

C     ==================================================================

      VAR = 0
      DO 200 I = 1, SPRNOD
         NODE = PERM(I)
         DO 100 J = XNODES(NODE), XNODES(NODE+1)-1
            IDOF = NODES(J)
            VAR = VAR + 1
            MEQN(IDOF) = VAR
  100    CONTINUE
  200 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPREQN
      END


      SUBROUTINE   SPRCN1   ( MPMNPC, MMNPC , XNODEL, NODEL , NEL   ,
     $                        NANOD , NMMNPC, OFFSET                  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEL   , NANOD , NMMNPC
      INTEGER           MPMNPC(NEL+1)         , MMNPC(NMMNPC)         ,
     $                  XNODEL(NANOD+1)       , NODEL(NMMNPC)         ,
     $                  OFFSET(NANOD)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRCN1
C
C --- TO COMPUTE THE NODE-ELEMENT CONNECTIVITY ARRAY (XNODEL,NODEL).
C     NOTE THAT THE ROUTINE ASSUMES THAT EACH NODE IS REPRESENTED ONLY
C     ONCE FOR EACH ELEMENT IT IS CONNECTED TO.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     NEL - INTEGER
C      ENTRY : NUMBER OF ORIGINAL FINITE-ELEMENTS.
C      EXIT  : NOT CHANGED.
C     NANOD - INTEGER
C      ENTRY : NUMBER OF ACTIVE FINITE-ELEMENT NODES.
C      EXIT  : NOT CHANGED.
C     NMMNPC - INTEGER
C      ENTRY : LENGTH OF MMNPC AND NODEL.
C      EXIT  : NOT CHANGED.
C     MPMNPC - INTEGER(NEL+1)
C      ENTRY : POINTER TO MMNPC.
C      EXIT  : NOT CHANGED.
C     MMNPC - INTEGER(NMMNPC)
C      ENTRY : NODES CONNECTED TO EACH ELEMENT.
C      EXIT  : NOT CHANGED.
C     XNODEL - INTEGER(NANOD+1)
C      ENTRY : NOT DEFINED.
C      EXIT  : POINTER TO LISTS OF ELEMENTS FOR EACH NODE IN NODEL.
C     NODEL - INTEGER(NMMNPC)
C      ENTRY : NOT DEFINED.
C      EXIT  : LISTS OF ELEMENTS CONNECTED FOR EACH NODE.
C
C     WORKING ARRAYS : NOT DEFINED ON ENTRY/EXIT.
C     -------------------------------------------
C
C     OFFSET : INTEGER(NANOD)
C              OFFSET INTO NODEL FOR EACH NODE WHEN THE ARRAY IS BUILT.
C
C     SUBPROGRAMS
C     -----------
C
C     NONE
C
C     FUNCTIONS
C     ---------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     INCLUDE BLOCKS
C     --------------
C
C     NONE
C
C     COMMON BLOCKS
C     -------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           I     , IP    , J     , NODE

C     ==================================================================

      DO 100 I = 1, NANOD
         OFFSET(I) = 0
         XNODEL(I) = 0
  100 CONTINUE
      XNODEL(NANOD+1) = 0

C     ==================================================================

C     -------------------------------------------
C     RUN THROUGH MMNPC AND COUNT THE APPEARANCES
C     OF EACH NODE, IT IS ASSUMED THAT A NODE IS
C     REPRESENTED ONLY ONE TIME FOR EACH ELEMENT
C     IT IS CONNECTED TO. THIS LOOP IS WRONG IF
C     THIS IS NOT THE CASE.
C     -------------------------------------------
      DO 300 I = 1, NMMNPC
         NODE = MMNPC(I)
         XNODEL(NODE+1) = XNODEL(NODE+1) + 1
  300 CONTINUE

C     -----------------------------------
C     BUILD THE NODE-ELEMENT CONNECTIVITY
C     AND STORE IT IN (XNODEL,NODEL).
C     -----------------------------------
      XNODEL(1) = 1
      DO 400 I = 1, NANOD
         XNODEL(I+1) = XNODEL(I+1) + XNODEL(I)
  400 CONTINUE
      DO 600 J = 1, NEL
         DO 500 I = MPMNPC(J), MPMNPC(J+1)-1
            NODE = MMNPC(I)
            IP = XNODEL(NODE) + OFFSET(NODE)
            NODEL(IP) = J
            OFFSET(NODE) = OFFSET(NODE) + 1
  500    CONTINUE
  600 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRCN1
      END
