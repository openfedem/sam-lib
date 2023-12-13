C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRSMB   ( MSPAR , MSICA , MTREES, MSIFA , MEQN  ,
     $                        IWORK , LPU   , IERR                    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           LPU   ,  IERR
      INTEGER           MSPAR(*)              , MSICA(*)              ,
     $                  MTREES(*)             , MSIFA(*)              ,
     $                  MEQN(*)               , IWORK(*)

C ======================================================================
C  S A M  library routine :  SPRSMB                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To compute the structure of the Cholesky factor L. The index info
C     is returned in MSIFA of length NMSIFA.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
*                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 Feb. 27, 2003 (kmo)
C                 Increased the size of MSIFA by MSPAR(58)*NEQ.
C
C
C     MSPAR - INTEGER(*)
C      Entry : As on exit from SPRTRS.
C      Exit  : Not changed, except for:
C              MSPAR(1) that is set to 4 in case of sucessful exit, and
C                       to -4 else.
C              MSPAR(22) Number of panels in L.
C              Also:  MSPAR(3) = NMSIFA  is updatad
C     MSICA - INTEGER(NMSICA)
C      Entry : As on exit from SPRSAS.
C      Exit  : Not changed.
C     MEQN - INTEGER(NDOF)
C      Entry : As on exit from SPRTRS.
C      Exit  : Not changed.
C     MTREES - INTEGER(NTREES)
C      Entry : The tree structures of the matrix corresponding to the
C              final order of the spr-nodes, i.e. the variables.
C      Exit  : The pointer array SUPLEN->XLINDX, i.e. it is converted
C              to a pointer array into MSIFA=LINDX.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : Not defined.
C      Exit  : The index lists for each supernode in order to represent
C              the structure of the Cholesky factor.
C     LPU - INTEGER
C      Entry : Output unit.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call,
C              = -1 if error from SPRTRS is not cleared.
C              = -2 if MSICA, MTREES, MSIFA and/or IWORK
C                   are not great enough.
C              = -3 Error from SPRSM1.
C
C     Working arrays
C     --------------
C
C     IWORK - INTEGER(NIWORK)
C      Entry : Not defined. Partitioned to internal work array,
C              see below.
C      Exit  : Need not be saved.
C
C     Procedures
C     ----------
C
C     SPRSM1
C     SPRSM2
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
     $  NMSICA, NTREES, NMSIFA, NIWORK

      INTEGER
     $  CACHE , ELNOD , ELSEQ , INVP  , IWUSED, LINK  , LINDX , MASK  ,
     $  MLUSED, MSUSED, MTUSED, NB    , NBLOCK, NEL   , NELACT, NEQ   ,
     $  NODES , NODMAP, NOFSUB, NSUPER, PERM  , SIBL  , SPLIT , SPRCON,
     $  SPRNOD, SUPSUP, XBLOCK, XELNOD, XELSEQ, XLINDX, XLNZ  , XNOD  ,
     $  XNODES, XSIBL , XSPLIT, XSUPER

C*DAM       INTEGER           I     , ILENGT, IP    , J

      EXTERNAL          SPRSM1, SPRSM2, SPRER1

C     ==================================================================

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      NMSICA = MSPAR(2)
      NTREES = MSPAR(4)
      NMSIFA = MSPAR(3)
      NIWORK = 3*MSPAR(6) + 2*MSPAR(8) + 1
C
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -4
         IERR = -1
         CALL SPRER1 ( 30, 'SPRSMB', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         MSPAR(1) = 4
         IERR = 0
      ENDIF

      NEL    = MSPAR( 5)
      SPRNOD = MSPAR( 6)
      SPRCON = MSPAR( 7)
      NEQ    = MSPAR( 8)
      CACHE  = MSPAR( 9)
      NSUPER = MSPAR(11)
      NB     = MSPAR(12)
      NOFSUB = MSPAR(15)
      NELACT = MSPAR(19)

C     ----------------------------------
C     DIVIDE IWORK AND CHECK ITS LENGTH.
C     ----------------------------------
      INVP   = 1
      LINK   = INVP   + SPRNOD
      MASK   = LINK   + NEQ
      NODMAP = MASK   + SPRNOD
      XNOD   = NODMAP + NEQ
      IWUSED = XNOD   + SPRNOD

C     -------------------------------------------
C     SET POINTERS TO MSICA AND CHECK ITS LENGTH.
C     -------------------------------------------
      XELNOD = 1
      XNODES = XELNOD + NEL    + SPRNOD + 2
      NODES  = XNODES + SPRNOD + 1
      PERM   = NODES  + NEQ
      ELNOD  = PERM   + SPRNOD
      MSUSED = ELNOD  + SPRCON - 1

C     --------------------------------------------
C     SET POINTERS TO MTREES AND CHECK ITS LENGTH.
C     --------------------------------------------
      XSUPER = 1
      XSIBL  = XSUPER + NSUPER + 1
      SIBL   = XSIBL  + NSUPER + 1
      XLINDX = SIBL   + NSUPER + 1

* === SUPSUP-acd
      SUPSUP = XLINDX + NSUPER + 1
      ELSEQ  = SUPSUP + NSUPER
      XELSEQ = ELSEQ  + NELACT
      XBLOCK = XELSEQ + NB     + 1
      MTUSED = XBLOCK + NB

C     -------------------------------------------
C     SET POINTERS TO MSIFA AND CHECK ITS LENGTH.
C     -------------------------------------------
      LINDX  = MSPAR(58)*NEQ + 1
      XLNZ   = LINDX  + NOFSUB
      XSPLIT = XLNZ   + NSUPER + 1
      SPLIT  = XSPLIT + NSUPER
      MLUSED = SPLIT  + NEQ

      IF ( IWUSED .GT. NIWORK .OR. MSUSED .GT. NMSICA .OR.
     $     MTUSED .GT. NTREES .OR. MLUSED .GT. NMSIFA      ) THEN
         MSPAR(1) = -2
         IF ( IWUSED .GT. NIWORK ) THEN
            WRITE ( LPU, '(/A,I8,I8/)' )'IWORK too small',NIWORK,IWUSED
            NIWORK = IWUSED
         ENDIF
         IF ( MSUSED .GT. NMSICA ) THEN
            WRITE ( LPU, '(/A,I8,I8/)' )'MSICA too small',NMSICA,MSUSED
            NMSICA = MSUSED
         ENDIF
         IF ( MTUSED .GT. NTREES ) THEN
            WRITE ( LPU, '(/A,I8,I8/)' )'MTREES too small',NTREES,MTUSED
            NTREES = MTUSED
         ENDIF
         IF ( MLUSED .GT. NMSIFA ) THEN
            WRITE ( LPU, '(/A,I8,I8/)' )'MSIFA too small',NMSIFA,MLUSED
            NMSIFA = MLUSED
         ENDIF
         IERR = -2
         CALL SPRER1 ( 35, 'SPRSMB', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

C     ==================================================================

C     --------------------------------------------------
C     COMPUTE XLINDX,LINDX AND XLNZ BY A CALL TO SPRSM1.
C     --------------------------------------------------
      CALL SPRSM1 ( NB, NSUPER, SPRNOD, NEL, NELACT,  NEQ, SPRCON,
     $              NOFSUB, MTREES(XBLOCK), MEQN, MTREES(XELSEQ),
     $              MTREES(ELSEQ), MSICA(XELNOD), MSICA(ELNOD),
     $              MTREES(XSIBL), MTREES(SIBL), MSICA(PERM),
     $              MTREES(XSUPER), MSICA(XNODES), MSICA(NODES),
     $              MTREES(XLINDX), MSIFA(LINDX), MSIFA(XLNZ),
     $              IWORK(INVP), IWORK(LINK), IWORK(MASK),
     $              IWORK(NODMAP), IWORK(XNOD), IERR )

      IF ( IERR .NE. 0 ) THEN
         CALL SPRER1 ( 23, 'SPRSMB', 0, 0, 0, LPU, IERR )
         IERR = -3
         RETURN
      ENDIF

C     --------------------------------------------------------------
C     COMPUTE THE PANELS OF L FOR OPTIMIZED USE OF THE CACHE MEMORY.
C     --------------------------------------------------------------
      CALL SPRSM2 ( NSUPER, CACHE, MTREES(XSUPER), MTREES(XLINDX),
     $              NBLOCK, MSIFA(XSPLIT), MSIFA(SPLIT) )

C     -----------------------------------------------------
C     FINALLY, SET SIZE OF MSIFA AND NUMBER OF PANELS IN L.
C     -----------------------------------------------------
      MSPAR(3) = MLUSED - NEQ + NBLOCK
C*DAM This is the final size of MSIFA.
      MSPAR(22) = NBLOCK

* --- SWAP XSIBL <-> XLINDX to simplify access to matrix structures.

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSMB
      END


      SUBROUTINE   SPRSM1   ( NB    , NSUPER, NNODS , NEL   , NELACT,
     $                        NEQ   , NELNOD, NOFSUB, XBLOCK,
     $                        MEQN  , XELSEQ, ELSEQ , XELNOD, ELNOD ,
     $                        XSIBL , SIBL  , PERM  , XSUPER, XNODES,
     $                        NODES , XLINDX, LINDX , XLNZ  , INVP  ,
     $                        LINK  , MASK  , NODMAP, XNOD  , IERR    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NB    , NSUPER, NNODS , NEL   , NELACT,
     $                  NEQ   , NELNOD, NOFSUB, IERR
      INTEGER           XBLOCK(NB+1)          , MEQN(*)               ,
     $                  XELSEQ(NB+1)          , ELSEQ(NELACT)         ,
     $                  XELNOD(NEL+1)         , ELNOD(NELNOD)         ,
     $                  XSIBL(NSUPER+1)       , SIBL(NSUPER+1)        ,
     $                  PERM(NNODS)           , XSUPER(NSUPER+1)      ,
     $                  XNODES(NNODS+1)       , NODES(NEQ)            ,
     $                  XLINDX(NSUPER+1)      , LINDX(NOFSUB)         ,
     $                  XLNZ(NSUPER+1)        , INVP(NNODS)           ,
     $                  LINK(NEQ)             , MASK(NNODS)           ,
     $                  NODMAP(NEQ)           , XNOD(NNODS+1)

C ======================================================================
C  S A M  library routine :  SPRSM1                  GROUP 9 / PRIVATE
C ======================================================================
C
C     PURPOSE
C     -------
C
C --- TO COMPUTE THE STRUCTURE OF THE CHOLESKY FACTOR, I.E. TO DO THE
C     SYMBOLIC FACTORIZATION STEP. THE ROUTINE GENERATES THE ARRAY PAIR
C     XLINDX,LINDX. IT ALSO COMPUTES THE POINTER ARRAY XLNZ THAT POINTS
C     TO THE START IN LNZ FOR EACH SUPERNODE.
C
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     NB - INTEGER
C      ENTRY : NUMBER EXTERNAL ASSEMBLY STEPS.
C      EXIT  : NOT CHANGED.
C     NSUPER - INTEGER
C      ENTRY : NUMBER OF SUPERNODES.
C      EXIT  : NOT CHANGED.
C     NNODS - INTEGER
C      ENTRY : NUMBER OF NODES.
C      EXIT  : NOT CHANGED.
C     NEL - INTEGER
C      ENTRY : NUMBER OF ELEMENTS.
C      EXIT  : NOT CHANGED.
C     NELACT - INTEGER
C      ENTRY : NUMBER OF ELEMENTS THAT IS ACTIVE IN THE ASSEMBLY.
C      EXIT  : NOT CHANGED.
C     NEQ - INTEGER
C      ENTRY : NUMBER OF FREE VARIABLES
C      EXIT  : NOT CHANGED.
C     NOFSUB - INTEGER
C      ENTRY : LENGTH OF LINDX.
C      EXIT  : NOT CHANGED.
C     XBLOCK - INTEGER(NB+1)
C      ENTRY : THE ASSEMBLY STEP PARTITION.
C      EXIT  : NOT CHANGED.
C     MEQN - INTEGER(NDOF)
C      ENTRY : FOR EACH DOF, IF MEQN(IDOF)>0 IT GIVES THE EQUATION
C              NUMBER FOR A DEGREE OF FREEDOM.
C      EXIT  : NOT CHANGED.
C     XELSEQ - INTEGER(NB+1)
C      ENTRY : THE ELEMENT PARTITON ACCORDING TO THE ASSEMBLY STEPS.
C      EXIT  : NOT CHANGED.
C     ELSEQ - INTEGER(NEL)
C      ENTRY : THE ORIGINAL FINITE-ELEMENTS THAT ARE ACTIVE IN EACH
C              ASSEMBLY STEP.
C      EXIT  : NOT CHANGED.
C     XSIBL - INTEGER(NSUPER+1)
C      ENTRY : POINTER TO CHILDREN STORED IN SIBL FOR CHILDREN LIST.
C      EXIT  : NOT CHANGED.
C     SIBL - INTEGER(NSUPER+1)
C      ENTRY : LISTS OF CHILDREN FOR EACH SUPERNODE.
C      EXIT  : NOT CHANGED.
C     PERM - INTEGER(NNODS)
C      ENTRY : NEW TO OLD MAP OF THE NODES.
C      EXIT  : NOT CHANGED.
C     XSUPER - INTEGER(NSUPER)
C      ENTRY : THE SUPERNODE PARTITION WITH REFERENCE TO VARIABLES.
C      EXIT  : NOT CHANGED.
C     XNODES - INTEGER(NNODS+1)
C      ENTRY : THE NODE PARTITION ACCORDING TO INITIAL ORDER.
C      EXIT  : NOT CHANGED.
C     NODES - INTEGER(NEQ)
C      ENTRY : THE DEGREE OF FREEDOM IDENTIFICATORS FOR EACH ACTIVE
C              VARIABLE IN A NODE.
C     XLINDX - INTEGER(NSUPER+1)
C      ENTRY : THE LENGTH OF EACH SUPERNODE.
C      EXIT  : POINTERS TO THE START OF THE VARIABLE REFERENCED INDEX
C              LIST FOR EACH SUPERNODE THAT IS STORED IN LINDX.
C     LINDX - INTEGER(NOFSUB)
C      ENTRY : NOT DEFINED.
C      EXIT  : FOR EACH SUPERNODE, THE LIST OF VARIABLES IN THE SUPER-
C              NODE INCLUDING EXPLICIT INDEX FOR THE DIAGONAL ENTRY.
C     XLNZ - INTEGER(NSUPER+1)
C      ENTRY : NOT DEFINED.
C      EXIT  : POINTERS TO THE FIRST LOCATION IN LNZ FOR EACH SUPERNODE.
C     IERR - INTEGER
C      ENTRY : NOT DEFINED.
C      EXIT  : ERROR FLAG.
C              =  0 : NORMAL RETURN.
C              = -1 : LINDX IS NOT GREAT ENOUGH, SIGNALS INPUT DATA
C                     ERROR.
C              = -2 : TOO FEW INDICES ARE FOUND IN A SUPERNODE,SIGNALS
C                     INPUT DATA ERROR.
C      ENTRY :
C      EXIT  : NOT CHANGED.
C
C
C     WORKING ARRAYS
C     --------------
C
C     INVP - INTEGER(NNODS)
C      ENTRY : NOT DEFINED.
C              OLD TO NEW MAP OF THE NODES.
C      EXIT  : NEED NOT BE SAVED.
C     LINK - INTEGER(NEQ)
C      ENTRY : NOT DEFINED.
C              LINKS TOGETHER INDICES TO VARIABLES IN THE FRONT.
C      EXIT  : NEED NOT BE SAVED.
C     MASK - INTEGER(NNODS)
C      ENTRY : NOT DEFINED.
C              MASKS NODES FOR VARIABLE INDICES THAT ARE MERGED TO
C              REDUCE SEARCH.
C      EXIT  : NEED NOT BE SAVED.
C     NODMAP - INTEGER(NEQ)
C      ENTRY : NOT DEFINED.
C              MAPS EACH VARIABLE TO ITS NODE, THIS IS THE POSITIVE
C              INTEGERS IN MEQN, NOT THE DOF IDENTIFICATORS.
C      EXIT  : NEED NOT BE SAVED.
C     XNOD - INTEGER(NNODS+1)
C      ENTRY : NOT DEFINED.
C              THE NODE PARTITION.
C      EXIT  : NEED NOT BE SAVED.
C
C     PROCEDURES
C     ----------
C
C     SPRS11
C     SPRS12
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

      INTEGER           COLLEN, FSTSUP, FSTVAR, I     , IBL   , IDOF  ,
     $                  INDEX , IP    , IPELM , ISTOP , ISTRT , IVAR  ,
     $                  J     , JS    , JSLEN , JSTOP , JSTRT , LSTSUP,
     $                  LSTVAR, NCHILD, NFOUND, NUMELM, NXTIND, OLDN
      EXTERNAL          SPRS11, SPRS12

C     ==================================================================

      IERR = 0

C     -----------------------------------------------
C     INITIALIZE FOR THE SYMBOLIC FACTORIZATION STEP.
C     -----------------------------------------------
      XNOD(1) = 1
      DO 200 I = 1, NNODS
         OLDN = PERM(I)
         INVP(OLDN) = I
         JSTRT = XNODES(OLDN)
         JSTOP = XNODES(OLDN+1)
         XNOD(I+1) = XNOD(I) + JSTOP - JSTRT
         MASK(I) = 0
         DO 100 J = JSTRT, JSTOP-1
            IDOF = NODES(J)
            IVAR = MEQN(IDOF)
            NODMAP(IVAR) = I
  100    CONTINUE
  200 CONTINUE

C     -----------------------------
C     CONVERT XLINDX TO FINAL FORM.
C     -----------------------------
      JSLEN = XLINDX(1)
      XLINDX(1) = 1
      DO 300 JS = 1, NSUPER
         COLLEN = XLINDX(JS+1)
         XLINDX(JS+1) = XLINDX(JS) + JSLEN
         JSLEN = COLLEN
  300 CONTINUE

C     ==================================================================

C     -----------
C     START XLNZ.
C     -----------
      XLNZ(1) = 1

C     --------------------------
C     FOR EACH ASSEMBLY STEP ...
C     --------------------------
      DO 700 IBL = 1, NB

C        --------------------------------
C        SET UP FIRST AND LAST SUPERNODE,
C        POINTER TO LIST OF ORIGINAL
C        FINITE-ELEMENTS AND THE NUMBER.
C        --------------------------------
         FSTSUP = XBLOCK(IBL)
         LSTSUP = XBLOCK(IBL+1) - 1
         IPELM  = XELSEQ(IBL)
         NUMELM = XELSEQ(IBL+1) - IPELM

C        ----------------------
C        FOR EACH SUPERNODE ...
C        ----------------------
         DO 600 JS = FSTSUP, LSTSUP

C           -------------------------------------------------
C           SET UP POINTER TO LINDX, THE LENGTH OF COLUMN JS,
C           POINTER TO CHILDREN LIST IN SIBL, FIRST AND LAST
C           VARIABLE AND INITIATE LINK FOR THIS SUPERNODE.
C           -------------------------------------------------
            IP = XLINDX(JS)
            JSLEN = XLINDX(JS+1) - IP
            ISTRT = XSIBL(JS)
            ISTOP = XSIBL(JS+1) - 1
            FSTVAR = XSUPER(JS)
            LSTVAR = XSUPER(JS+1) - 1
            LINK(LSTVAR) = NEQ + 1

C           --------------------------------
C           CHECK THAT LINDX IS LONG ENOUGH.
C           --------------------------------
            IF ( (IP+JSLEN-1) .GT. NOFSUB ) THEN
               IERR = -1
               RETURN
            ENDIF

C           -------------------------------------------------
C           STORE THE "IN-SUPER" INDICES DIRECTLY INTO LINDX,
C           AND INITIALIZE THE NUMBER OF INDICES FOUND.
C           -------------------------------------------------
            DO 400 I = FSTVAR, LSTVAR
               LINDX(IP) = I
               IP = IP + 1
  400       CONTINUE
            NFOUND = LSTVAR - FSTVAR + 1

C           --------------------------------
C           COMPUTE XLNZ FOR NEXT SUPERNODE.
C           --------------------------------
            XLNZ(JS+1) = XLNZ(JS) + JSLEN*NFOUND - (NFOUND*(NFOUND-1)/2)

            IF ( ISTRT .LE. ISTOP .AND. NFOUND .LT. JSLEN ) THEN

C              -----------------------------------------------
C              ASSEMBLE THE INDEX LISTS OF GENERATED ELEMENTS.
C              -----------------------------------------------
               NCHILD = ISTOP - ISTRT + 1
               CALL SPRS11 ( JS, NSUPER, NFOUND, NCHILD, NNODS, NEQ,
     $                       NOFSUB, SIBL(ISTRT), XNOD, XSUPER, XLINDX,
     $                       LINDX, LINK, NODMAP, MASK )
            ENDIF

            IF ( JS .EQ. FSTSUP .AND. NFOUND .LT. JSLEN ) THEN

C              ---------------------------------------------------
C              ASSEMBLE INDEX LISTS FROM ORIGINAL FINITE-ELEMENTS.
C              ---------------------------------------------------
               CALL SPRS12 ( JS, LSTVAR, NNODS, NEL, NELNOD, NUMELM,
     $                       NEQ, NFOUND, INVP, XNOD, XELNOD, ELNOD,
     $                       ELSEQ(IPELM), LINK, MASK, NODMAP )
            ENDIF

            IF ( NFOUND .LT. JSLEN ) THEN

C              ---------------------------
C              INDEX LIST IS NOT COMPLETE,
C              EXIT WITH AN ERROR.
C              ---------------------------
               IERR = -2
               RETURN
            ENDIF

C           ----------------------------------------
C           FINALLY, COPY CONTENTS OF LINK TO LINDX.
C           ----------------------------------------
            JSTRT = XLINDX(JS) + LSTVAR - FSTVAR + 1
            JSTOP = XLINDX(JS+1) - 1

            IF ( JSTRT .LE. JSTOP ) THEN

C              -----------------------------------------------
C              THE FIRST INDEX BELOW THE DIAGONAL BLOCK OF THE
C              SUPERNODE IS FOUND IN LINK(LSTVAR). STORE IT IN
C              THE FIRST OFF-DIAGONAL POSITION OF LINDX AND
C              CONTINUE DOWN THE COLUMN WHILE FETCHING INDICES
C              FROM LINK.
C              -----------------------------------------------
               NXTIND = LINK(LSTVAR)
               LINDX(JSTRT) = NXTIND
               DO 500 J = JSTRT+1, JSTOP
                  INDEX = LINK(NXTIND)
                  LINDX(J) = INDEX
                  NXTIND = INDEX
  500          CONTINUE
            ENDIF

  600    CONTINUE
  700 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSM1
      END


      SUBROUTINE   SPRSM2   ( NSUPER, CACHSZ, XSUPER, XLINDX, NBLOCK,
     $                        XSPLIT, SPLIT                           )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSUPER, CACHSZ, NBLOCK
      INTEGER           XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  XSPLIT(NSUPER)        , SPLIT(*)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRSM2
C
C --- TO BUILD THE CACHE PANELS FOR MORE EFFICIENT USE OF CACHE.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES
C     CACHSZ : INTEGER
C              CACHE SIZE IN KBYTES, .LE. 0 => UNBLOCKED OUTPUT, AND
C              (XSPLIT,SPLIT) ARE NOT REFERENCED IN FACTOR ROUTINES.
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION.
C     XLINDX : INTEGER(NSUPER+1)
C              GIVES THE LENGTH OF EACH SUPERNODE.
C
C     ON EXIT
C     -------
C
C     NBLOCK : INTEGER
C              TOTAL NUMBER OF CACHE BLOCKS, USED TO DIMENSION THE
C              ARRAY SPLIT.
C     (XSPLIT,
C      SPLIT): INTEGER(NSUPER), INTEGER(NBLOCK+1)
C              INFORMATION ABOUT THE CACHE BLOCKS. WORST CASE REQ.
C              FOR SPLIT IS NEQNS+1, I.E. PASS THE ARRAY WITH SIZE
C              NEQNS+1 TO THE ROUTINE. THIS IMPOSES THE REQUIREMENT
C              THAT TWO COLUMNS OF THE MATRIX ARE ASSUMED TO FIT
C              IN CACHE.
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

      INTEGER           CACHE , FSTVAR, ICACHE, J     , JS    , LENVAR,
     $                  LSTVAR, UCACHE

C     ==================================================================

C     ----------------------------------------------
C     QUICK RETURN IF NO CACHE BLOCKING IS REQUIRED.
C     ----------------------------------------------
      NBLOCK = 0
      IF ( CACHSZ .LE. 0 ) RETURN

C     --------------------------------
C     COMPUTE THE CACHE SIZE IN WORDS.
C     --------------------------------
      CACHE  = CACHSZ * 128 ! 1024 / 8

C     ==================================================================

      DO 200 JS = 1, NSUPER
C        -------------------
C        GET SUPERNODE INFO.
C        -------------------
         LENVAR = XLINDX(JS+1) - XLINDX(JS)
         FSTVAR = XSUPER(JS)
         LSTVAR = XSUPER(JS+1) - 1

         UCACHE = CACHE  - ( 2*LENVAR )

C        -------------------
C        SET START BLOCK FOR
C        CURRENT SUPERNODE.
C        -------------------
         NBLOCK = NBLOCK + 1
         XSPLIT(JS) = NBLOCK
         SPLIT(NBLOCK) = FSTVAR

C        ----------------
C        PUT FIRST COLUMN
C        INTO THE CACHE.
C        ----------------
         ICACHE = LENVAR
         DO 100 J = FSTVAR+1, LSTVAR
C           --------------------
C           FOR THE NEXT COLUMNS
C           BLOCK THE SUPERNODE.
C           --------------------
            LENVAR = LENVAR - 1
            IF ( ICACHE + LENVAR .GT. UCACHE ) THEN
C              -----------------
C              CREATE NEW BLOCK.
C              -----------------
               NBLOCK = NBLOCK + 1
               SPLIT(NBLOCK) = J
               ICACHE = 0
            ELSE
C              ------------------------
C              ADD IT TO CURRENT BLOCK.
C              ------------------------
               ICACHE = ICACHE + LENVAR
            ENDIF
  100    CONTINUE
  200 CONTINUE

C     -------------------
C     SET THE LAST LIMIT.
C     -------------------
      SPLIT(NBLOCK+1) = XSUPER(NSUPER+1)

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSM2
      END


      SUBROUTINE   SPRS11   ( JS    , NSUPER, NFOUND, NCHILD, NNODS ,
     $                        NEQNS , NOFSUB, SIBL  , XNOD  , XSUPER,
     $                        XLINDX, LINDX , LINK  , NODMAP, MASK    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           JS    , NSUPER, NFOUND, NCHILD, NNODS , NEQNS ,
     $                  NOFSUB
      INTEGER           SIBL(NCHILD)          , XNOD(NNODS+1)         ,
     $                  XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  LINDX(NOFSUB)         , LINK(NEQNS)           ,
     $                  MASK(NNODS)           , NODMAP(NEQNS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRS11
C
C --- TO PERFORM THE SYMBOLIC ASSEMBLY STEP OF GENERATED ELEMENTS. THE
C     MAIN TASK OF THE ROUTINE IS THUS TO START THE SET UP OF THE ARRAYS
C     (XLINDX,LINDX) THAT STORE THE INDEX INFORMATION FOR EACH GENERATED
C     ELEMENT. THIS SUBROUTINE DOES AN INCOMPLETE SET UP OF THE ARRAYS.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     JS     : INTEGER
C              THE CURRENT SUPERNODE.
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES AS COMPUTED IN MA46AD.
C     NFOUND : INTEGER
C              THE NUMBER OF INDICES IN THE FRONT FOUND SO FAR.
C     NCHILD : INTEGER
C              NUMBER OF CHILDREN TO THE CURRENT SUPERNODE.
C     NNODS  : INTEGER
C              THE TOTAL NUMBER OF NODES IN THE MODEL.
C     NEQNS  : INTEGER
C              NUMBER OF VARIABLES.
C     NOFSUB : INTEGER
C              NUMBER OF INDICES ALLOACTED FOR LINDX. SHOULD BE SOMEWHAT
C              GREATER THAN WHAT WAS COMPUTED IN MA46AD.
C     SIBL   : INTEGER(NCHILD)
C              LIST OF CHILDREN TO THE CURRENT SUPERNODE.
C     XNOD   : INTEGER(NNODS+1)
C              THE NODE PARTITION OF THE VARIABLES.
C     XSUPER : INTEGER(NSUPER+1)
C              THE SUPERNODE PARTITION.
C     XLINDX : INTEGER(NSUPER)
C              POINTER ARAY INTO LINDX FOR THE LIST OF INDICES TO THE
C              GENERATED ELEMENTS THAT ARE COMPLETED AT THIS STAGE.
C              NOTE THAT IT POINTS TO THE FIRST LOCATION IN LINDX WHERE
C              THERE IS INFORMATION ON THE SUPERNODE IN QUESTION.
C     LINDX  : INTEGER(NOFSUB)
C              FOR SUPERNODES/GENERATED ELEMENTS THAT ARE COMPLETED THIS
C              ARRAY STORES THE INDEX LISTS.
C     LINK   : INTEGER(NEQNS)
C              AN ARRAY THAT RECEIVES THE OFF-SUPER ROW INDICES, ON ENTR
C              IT IS ASSUMED THAT LINK(LSTVAR) = 'NEXT IN LIST' IN ORDER
C              TO START THE LINK PROCESS.
C     MASK   : INTEGER(NACTIV)
C              A VECTOR USED TO MASK OFF NODES THAT ARE MERGED IN THE
C              CURRENT STEP.
C     NODMAP : INTEGER(NEQNS)
C              A VECTOR THAT MAPS EACH VARIABLE TO ITS NODE.
C
C     ON EXIT
C     -------
C
C     NFOUND : INTEGER
C              UPDATED WITH THE INDICES FOUND.
C     LINDX  : INTEGER(NOFSUB)
C     LINK   : INTEGER(NEQNS)
C              STORES THE INDICES THAT BELONGS TO THE CURRENT
C              GENERATED ELEMENT AS A LINKED LIST. THE FIRST
C              INDEX AFTER THE SUPERNODE IS FOUND IN LINK(LSTVAR).
C              THE NEXT INDICES FOLLOWS IN LINKED ORDER.
C     MASK   : INTEGER(NACTIV)
C              UPDATED WITH THE NODES MERGED IN THE CURRENT STEP.
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
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

      INTEGER           CHDSZE, CHILD , FSTIND, I     , INDEX , J     ,
     $                  JSTRT , JSTOP , K     , LSTVAR, NEWNOD, NODE  ,
     $                  NXTIND, RINDEX, SIZE

C     ==================================================================

      LSTVAR = XSUPER(JS+1) - 1

      DO 400 I = 1, NCHILD
         CHILD  = SIBL(I)
         FSTIND = LSTVAR
         CHDSZE = XSUPER(CHILD+1) - XSUPER(CHILD)
         JSTRT  = XLINDX(CHILD) + CHDSZE
         JSTOP  = XLINDX(CHILD+1) - 1

C        -------------------------------------
C        TAKE THE ENTRIES FOR THIS GENERATED
C        ELEMENT AND STORE IN THE LINKED LIST.
C        -------------------------------------
         DO 300 J = JSTRT, JSTOP
            INDEX = LINDX(J)
            NODE  = NODMAP(INDEX)
            IF ( MASK(NODE) .NE. JS .AND. INDEX .GT. LSTVAR ) THEN

C              ----------------------------------------------
C              INDEX IS THE FIRST VARIABLE IN A NODE THAT HAS
C              NOT YET BEEN MERGED INTO THE INDEX SET OF THE
C              NEW GENERATED ELEMENT. NOTE THAT THE VARIABLES
C              AT THE NODE ARE STORED CONSEQUTIVELY WHICH
C              MAKES IT POSSIBLE TO STORE SETS OF VARIABLES
C              AT THE SAME TIME.
C              ----------------------------------------------
               MASK(NODE) = JS

C              -------------------------------
C              GET THE NEXT INDEX IN THE LIST.
C              -------------------------------
  100          NXTIND = LINK(FSTIND)
               IF ( INDEX .GT. NXTIND ) THEN
C                 --------------------------------------
C                 INDEX SHOULD BE PLACED ABOVE THIS ONE,
C                 SET FIRST INDEX TO NEXT INDEX AND TRY
C                 AGAIN.
C                 --------------------------------------
                  NEWNOD = NODMAP(NXTIND)
                  SIZE = XNOD(NEWNOD+1) - XNOD(NEWNOD)
                  FSTIND = NXTIND + SIZE - 1
                  GOTO 100
               ENDIF

C              -------------------------------------------
C              THE NEXT INDEX IS GREATER THAN THIS INDEX.
C              FIRST, SET INDEX TO NXTIND FOR FSTIND, THEN
C              PUSH INDEX, AND THE OTHER INDICES IN ITS
C              NODE INTO THE LINKED LIST. SECOND, SET THE
C              LAST INDEX IN THE NODE TO FSTIND, AND PUT
C              NXTIND INTO LINK(FSTIND).
C              -------------------------------------------
               SIZE = XNOD(NODE+1) - XNOD(NODE)
               LINK(FSTIND) = INDEX
               RINDEX = INDEX
               DO 200 K = 1, SIZE-1
                  LINK(RINDEX) = RINDEX + 1
                  RINDEX = RINDEX + 1
  200          CONTINUE
               FSTIND = RINDEX
               LINK(FSTIND) = NXTIND
               NFOUND = NFOUND + SIZE
            ENDIF
  300    CONTINUE
  400 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRS11
      END


      SUBROUTINE   SPRS12   ( JS    , LSTVAR, NNODS , NELS  , NELNOD,
     $                        NUMEL , NEQNS , NFOUND, INVP  , XNOD  ,
     $                        XELNOD, ELNOD , ELSEQ , LINK  , MASK  ,
     $                        NODMAP                                  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           JS    , LSTVAR, NNODS , NELS  , NELNOD, NUMEL ,
     $                  NEQNS , NFOUND
      INTEGER           INVP(NNODS)           , XNOD(NNODS+1)         ,
     $                  XELNOD(NELS+1)        , ELNOD(NELNOD)         ,
     $                  ELSEQ(NUMEL)          , LINK(NEQNS)           ,
     $                  MASK(NNODS)           , NODMAP(NEQNS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRS12
C
C --- TO PERFORM THE SYMBOLIC ASSEMBLY STEP OF ORIGINAL FINITE-ELEMENTS.
C
C     CREATED   : NOV. 24, 1993 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     JS     : INTEGER
C              THE CURRENT SUPERNODE.
C     LSTVAR : INTEGER
C              THE LAST VARIABLE IN THE SUPERNODE. WE KNOW THAT FSTVAR
C              TO LSTVAR DEFINES A FULL MATRIX ON THE DIAGONAL, THEREFOR
C              INDICES IN THE AREA FSTVAR=<INDEX=<LSTVAR WILL NOT BE
C              GATHERED INTO THE LINKED LIST.
C     NNODS  : INTEGER
C              THE TOTAL NUMBER OF NODES IN THE MODEL.
C     NELS   : INTEGER
C              NUMBER OF FINITE-ELEMENTS IN THE MODEL.
C     NELNOD : INTEGER
C              LENGTH OF ELNOD.
C     NUMEL  : INTEGER
C              NUMBER OF ELEMENTS TO BE ASSEMBLED AT THIS STAGE.
C     NEQNS  : INTEGER
C              NUMBER OF VARIABLES.
C     NFOUND : INTEGER
C              NUMBER OF TENTATIVE INDICES IN THE FRONT FOUND SO FAR.
C     INVP   : INTEGER(NNODS)
C              OLD TO NEW MAP.
C     XNOD   : INTEGER(NNODS+1)
C              THE NODE PARTITION OF THE VARIABLES.
C     XELNOD : INTEGER(NELS+1)
C              POINTER INTO THE ELEMENT-NODE CONNECTIVITY ARRAY ELNOD.
C     ELNOD  : INTEGER(NELNOD)
C              THE LIST OF NODES CONNECTED TO EACH ELEMENT.
C     ELSEQ  : INTEGER(NUMEL)
C              LIST OF ELEMENTS TO BE ASSEMBLED AT THIS STAGE.
C     LINK   : INTEGER(NEQNS)
C              ON ENTRY IT IS ASSUMED THAT LINK(LSTVAR) = 'NEXT IN LIST'
C     NODMAP : INTEGER(NEQNS)
C              A VECTOR THAT MAPS EACH VARIABLE TO ITS NODE.
C     MASK   : INTEGER(NNODS)
C              A VECTOR USED TO MASK OFF NODES THAT ARE MERGED IN THE
C              CURRENT STEP.
C
C     ON EXIT
C     -------
C
C     NFOUND : INTEGER
C              INCREMENTED WITH THE NUMBER OF INDICES FOUND.
C     LINK   : INTEGER(NEQNS)
C              STORES THE INDICES THAT BELONGS TO THE CURRENT ORIGINAL
C              FINITE-ELEMENT AS A LINKED LIST. THE FIRST INDEX AFTER
C              THE SUPERNODE IS FOUND IN LINK(LSTVAR).
C              THE NEXT INDICES FOLLOWS IN LINKED ORDER.
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     PROCEDURES
C     ----------
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

      INTEGER           ELEMNT, FSTIND, I     , INDEX , J     , JSTRT ,
     $                  JSTOP , K     , NEWNOD, NODE  , NNODE , NSIZE ,
     $                  NXTIND, RINDEX, SIZE

C     ==================================================================

      DO 400 I = 1, NUMEL
         ELEMNT = ELSEQ(I)
         FSTIND = LSTVAR
         JSTRT = XELNOD(ELEMNT)
         JSTOP = XELNOD(ELEMNT+1) - 1
         DO 300 J = JSTRT, JSTOP
            NODE = ELNOD(J)
            FSTIND = LSTVAR
            NNODE = INVP(NODE)
            NSIZE = XNOD(NNODE+1) - XNOD(NNODE)
            IF ( MASK(NNODE) .NE. JS .AND. NSIZE .GT. 0 ) THEN
C              -------------------------------------------
C              NNODE IS NOT YET MERGED INTO THE STRUCTURE,
C              MERGE IT IN AND DO THE PROPER UPDATES.
C              -------------------------------------------
               MASK(NNODE) = JS
               INDEX = XNOD(NNODE)

               IF ( INDEX .GT. LSTVAR ) THEN

C                 -------------------------------
C                 GET THE NEXT INDEX IN THE LIST.
C                 -------------------------------
  100             NXTIND = LINK(FSTIND)
                  IF ( INDEX .GT. NXTIND ) THEN
C                    --------------------------------------
C                    INDEX SHOULD BE PLACED ABOVE THIS ONE,
C                    SET FIRST INDEX TO NEXT INDEX AND TRY
C                    AGAIN.
C                    --------------------------------------
                     NEWNOD = NODMAP(NXTIND)
                     SIZE = XNOD(NEWNOD+1) - XNOD(NEWNOD)
                     FSTIND = NXTIND + SIZE - 1
                     GOTO 100
                  ENDIF

C                 -------------------------------------------
C                 THE NEXT INDEX IS GREATER THAN THIS INDEX.
C                 FIRST, SET INDEX TO NXTIND FOR FSTIND, THEN
C                 PUSH INDEX, AND THE OTHER INDICES IN ITS
C                 NODE INTO THE LINKED LIST. SECOND, SET THE
C                 LAST INDEX IN THE NODE TO FSTIND, AND PUT
C                 NXTIND INTO LINK(FSTIND).
C                 -------------------------------------------
                  LINK(FSTIND) = INDEX
                  RINDEX = INDEX
                  DO 200 K = 1, NSIZE-1
                     LINK(RINDEX) = RINDEX + 1
                     RINDEX = RINDEX + 1
  200             CONTINUE
                  FSTIND = RINDEX
                  LINK(FSTIND) = NXTIND
                  NFOUND = NFOUND + NSIZE
               ENDIF
            ENDIF
  300    CONTINUE
  400 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRS12
      END
