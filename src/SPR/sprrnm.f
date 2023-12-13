C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRRNM   ( MSPAR , MSICA , IWORK , LPU   , IERR    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           LPU   , IERR
      INTEGER           MSPAR(*)      , MSICA(*) , IWORK(*)

C ======================================================================
C  S A M  library routine :  SPRRNM                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To order the finite-element variables by means of the SPR-node
C     partition.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
C                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 May. 16, 2003 (kmo)
C                 Added checks that MSICA and IWORK are great enough.
C                 Reduced the size of NODEL by SPRNOD.
C
C
C     MSPAR - INTEGER(*)
C      Entry : As on exit from SPRSAS.
C      Exit  : Not changed, except for MSPAR(1) that is set to 2 in case
C              of sucessful exit, and to -2 else.
C     MSICA - INTEGER(NMSICA)
C      Entry : As on exit from SPRSAS.
C      Exit  : Not changed.
C     NMSICA - INTEGER
C      Entry : The length of MSICA.
C      Exit  : Not changed.
C     NIWORK - INTEGER
C      Entry : The length of IWORK.
C      Exit  : Not changed.
C     LPU - INTEGER
C      Entry : Output unit for error messages.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call,
C              = -1 if error from SPRSAS is not cleared.
C              = -2 if MSICA and/or IWORK are not great enough.
C
C     Working arrays
C     --------------
C
C     IWORK - INTEGER(MSPAR(37))
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
C     SPRCN1
C     SPRMMD
C     SPRCMD
C     SPRER1
C
C     Intrinsic
C     ---------
C
C     MAX
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

      INTEGER           DELTA , ELMLNK, ELNOD , FLAG  , FLNOD , FNODL ,
     $                  HEAD  , INVP  , IWUSED, MARKER, MAXMRK, MAXNEL,
     $                  NCOMPR, NEL   , NEQ   , NODEL , NODES ,
     $                  PERM  , QSIZE , SPRCON, SPRNOD, STOREL, XELNOD,
     $                  XFLNOD, XFNODL, XNODEL, XNODES, SEPSZE, SEPNUM,
     $                  SEPCNT, SEPARR, VWGHT

      LOGICAL           LWEIGH

      INTRINSIC         MAX
      EXTERNAL          SPRCN1, SPRMMD, SPRCMD, SPRER1

C     ==================================================================

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -2
         IERR = -1
         CALL SPRER1 ( 30, 'SPRRNM', 0, 0, 0, LPU, IERR )
         RETURN
      ELSE
         MSPAR(1) = 2
         IERR = 0
      ENDIF

      NEL    = MSPAR( 5)
      SPRNOD = MSPAR( 6)
      SPRCON = MSPAR( 7)
      NEQ    = MSPAR( 8)
      SEPSZE = MSPAR(40)

C     ----------------------------------
C     DIVIDE IWORK AND CHECK ITS LENGTH.
C     ----------------------------------
      INVP   = 1
      XFLNOD = INVP   + SPRNOD
      FLNOD  = XFLNOD + NEL
      XFNODL = FLNOD  + SPRCON + NEL
      FNODL  = XFNODL + SPRNOD
      HEAD   = FNODL  + SPRCON + SPRNOD
      QSIZE  = HEAD   + NEQ    + 1
      VWGHT  = QSIZE  + SPRNOD
      ELMLNK = VWGHT  + SPRNOD
      FLAG   = ELMLNK + NEL
      STOREL = FLAG   + SPRNOD
      MARKER = STOREL + SPRNOD
      IWUSED = MARKER + SPRNOD - 1
      IF ( IWUSED .GT. MSPAR(37) ) THEN
         MSPAR(1) = -2
         IERR = -2
         CALL SPRER1 ( 34, 'SPRRNM', IWUSED, MSPAR(37), 0, LPU, IERR)
         RETURN
      ENDIF

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
      IWUSED = SEPARR + SPRNOD - 1
      IF ( IWUSED .GT. MSPAR(2) ) THEN
         MSPAR(1) = -2
         IERR = -2
         CALL SPRER1 ( 32, 'SPRRNM', IWUSED, MSPAR(35), 0, LPU, IERR)
         RETURN
      ENDIF

C     ==================================================================

C     ---------------------
C     COMPUTE XNODEL,NODEL.
C     ---------------------
      CALL SPRCN1 ( MSICA(XELNOD), MSICA(ELNOD), MSICA(XNODEL),
     $              MSICA(NODEL), NEL, SPRNOD, SPRCON, IWORK(MARKER) )

C     ==================================================================

C     ------------------------------------
C     ORDER THE NODES BY A CALL TO SPRMMD.
C     ------------------------------------
      MAXNEL = SPRCON + NEL
      MAXMRK = MAX(NEL,SPRNOD) + 1
      DELTA  = 0

      IF ( SEPSZE.LE.0 ) THEN

* ------ Standard ordering
         CALL SPRMMD ( NEL, SPRNOD, SPRCON, MAXNEL, MAXMRK, DELTA,
     $                 MSICA(XELNOD), MSICA(ELNOD), MSICA(XNODEL),
     $                 MSICA(NODEL), NCOMPR, MSICA(PERM), IWORK(INVP),
     $                 IWORK(XFLNOD), IWORK(FLNOD), IWORK(XFNODL),
     $                 IWORK(FNODL), IWORK(HEAD), IWORK(QSIZE),
     $                 IWORK(ELMLNK), IWORK(FLAG), IWORK(STOREL),
     $                 IWORK(MARKER) )

      ELSE IF ( SEPSZE.LE.SPRNOD ) THEN

* ------ Constrained ordering for superelement technique.
*        Reset SEPNUM and SEPCNT to the only legal values per now.
*        Use standard minimum degree condition (LWEIGH=.FALSE.).
         SEPNUM = 0
         SEPCNT = 1
         LWEIGH = .FALSE.

         CALL SPRCMD
     $    ( NEL,NEQ,SPRNOD,SPRCON,MAXNEL,MAXMRK,DELTA,LWEIGH,
     $      NCOMPR,MSICA(XELNOD),MSICA(ELNOD),MSICA(XNODEL),
     $      MSICA(NODEL),MSICA(XNODES),MSICA(PERM),IWORK(INVP),
     $      IWORK(XFLNOD),IWORK(FLNOD),IWORK(XFNODL),IWORK(FNODL),
     $      IWORK(HEAD),IWORK(QSIZE),IWORK(VWGHT),IWORK(ELMLNK),
     $      IWORK(FLAG),IWORK(STOREL),IWORK(MARKER),SEPNUM,
     $      SEPSZE,SEPCNT,MSICA(SEPARR),LPU )

      ELSE

* ------ There are no subgraphs to order. It is assumed that
*        PERM is already initialized according to the initial
*        order of the nodes.

      ENDIF

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRRNM
      END


      SUBROUTINE   SPRMMD   ( NELS  , NNODS , NELNOD, MAXNEL, MAXMRK,
     $                        DELTA , XELNOD, ELNOD , XNODEL, NODEL ,
     $                        NCOMPR, PERM  , INVP  , XFLNOD, FLNOD ,
     $                        XFNODL, FNODL , HEAD  , QSIZE , ELMLNK,
     $                        FLAG  , STOREL, MARKER                  )

C     ------------------------------------------------------------------

      INTEGER           NELS  , NNODS , NELNOD, MAXNEL, MAXMRK, DELTA ,
     $                  NCOMPR
      INTEGER           XELNOD(NELS+1)        , ELNOD(NELNOD)         ,
     $                  XNODEL(NNODS+1)       , NODEL(NELNOD)         ,
     $                  PERM(NNODS)           , INVP(NNODS)           ,
     $                  XFLNOD(NELS)          , FLNOD(MAXNEL)         ,
     $                  XFNODL(NNODS)         , FNODL(NNODS+NELNOD)   ,
     $                  HEAD(NNODS)           , QSIZE(NNODS)          ,
     $                  ELMLNK(NELS)          , FLAG(NNODS)           ,
     $                  STOREL(NNODS)         , MARKER(NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRMMD
C
C --- TO PERFORM A MINIMUM DEGREE PERMUTATION. THE ROUTINE UTILIZES THE
C     INHERENT CLIQUE STRUCTURE OF THE FINITE ELEMENT MESH.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NELS   : INTEGER
C              NUMBER OF ELEMENTS.
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     NELNOD : INTEGER
C              SIZE OF CONNECTIVITY ARRAYS.
C     MAXNEL : INTEGER
C              SIZE OF ARRAY FLNOD. THE SIZE MUST BE AT LEAST
C              NELS+NELNOD, BUT THE NUMBER OF WORKSPACE COMPRESSIONS
C              WILL DECREASE IF THE SIZE IS INCREASED.
C     MAXMRK : INTEGER
C              USED AS UPPER LIMIT FOR MARKING NODES AND ELEMENTS, ANY
C              NUMBER LARGER THAN MAX(NNODS,NELS)+1 IS OK.
C     DELTA  : INTEGER
C              TOLERANCE FOR MULTIPLE ELIMINATION.
C              .LT. 0 : WITHOUT MULTIPLE ELIMINATION.
C              .EQ. 0 : WITH MULTIPLE ELIMINATION.
C              .GT. 0 : WITH MULTIPLE ELIMINATION AND RELAXED MINIMUM
C                       DEGREE CONDITION.
C     (XELNOD,
C     ELNOD) : INTEGER(NELS+1), INTEGER(NELNOD)
C              ELEMENT-NODE CONNECTIVITY.
C     (XNODEL,
C     NODEL) : INTEGER(NNODS+1), INTEGER(NELNOD)
C              NODE-ELEMENT CONNECTIVITY.
C
C     ON EXIT
C     -------
C
C     NCOMPR : INTEGER
C              NUMBER OF WORKSPACE COMPRESSIONS NEEDED IN SPRMM2.
C     (PERM,
C     INVP)  : 2*INTEGER(NNODS)
C              THE PERMUTATION.
C              DURING EXECUTION THEY ARE USED AS WORKSPACE AS EXPLAINED
C              BELOW.
C
C     WORKING ARRAYS
C     --------------
C
C     XFLNOD : INTEGER(NELS)
C              STORES THE POINTERS INTO THE LISTS OF NODES WHICH ARE
C              CONNECTED TO THE ELEMENT. IF ELEMENTS ARE ABORBED, THEN
C              XFLNOD=0 FOR THESE ELEMENTS.
C     FLNOD  : INTEGER(MAXNEL)
C              STORES FOR EACH ELEMENT A LIST OF NODES CONNECTED TO THE
C              ELEMENT. EACH LIST OF NODES IS HEADED BY ITS LENGTH.
C     XFNODL : INTEGER(NNODS)
C              STORES THE POINTERS INTO THE LISTS OF ELEMENTS WHICH ARE
C              CONNECTED TO THE NODE. IF NODES ARE ELIMINATED, THEN
C              XFNODL=0 FOR THESE NODES.
C     FNODL  : INTEGER(NNODS+NELNOD)
C              STORES FOR EACH NODE A LIST OF ELEMENTS CONNECTED TO THE
C              NODE. EACH LIST OF ELEMENTS IS HEADED BY ITS LENGTH.
C     HEAD   : INTEGER(NNODS)
C              HEAD OF DEGREE LIST.
C     QSIZE  : INTEGER(NNODS)
C              SUPERNODE SIZES.
C     ELMLNK : INTEGER(NNODS)
C              GENERALIZED ELEMENTS WHICH NEED DEGREE UPDATES.
C     FLAG   : INTEGER(NNODS)
C              FLAGS ACTIVE NODES IN THE DEGREE UPDATE STEP.
C     STOREL : INTEGER(NNODS)
C              LIST OF NODES IN NEW GENERATED ELEMENTS AND LISTS
C              OF NODES OF DIFFERENT TYPES IN DEGREE UPDATE STEP.
C     MARKER : INTEGER(NNODS)
C              MARK OF NODES.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     SPRMM1
C     SPRMM2
C     SPRMM3
C     SPRMM4
C
C     ------------------------------------------------------------------

      INTEGER           CMINDG, CURMDN, I     , LIMIT , MARK  , NUM   ,
     $                  NUMELM, NXTMDN, NXTUSE

      EXTERNAL          SPRMM1, SPRMM2, SPRMM3, SPRMM4

C     ==================================================================

C     ----------------------------------
C     SET NUM TO # ELIMINATED NODES + 1.
C     NXTUSE IS SET TO NELS+NELNOD+1.
C     ----------------------------------

      NCOMPR = 0
      NUM    = 1
      NXTUSE = NELS + NELNOD + 1

C     ----------------------------------------------
C     INITIALIZE FOR THE MINIMUM DEGREE ELIMINATION.
C     ----------------------------------------------

      CALL SPRMM1 ( NELS, NNODS, NELNOD, MAXNEL, MAXMRK, XELNOD, ELNOD,
     $              XNODEL, NODEL, XFLNOD, FLNOD, XFNODL, FNODL, HEAD,
     $              INVP, PERM, QSIZE, ELMLNK, MARKER, FLAG, STOREL )

C     ==================================================================

      IF ( HEAD(1) .GT. 0 ) THEN

C        ------------------------------------
C        ELIMINATE THE ISOLATED NODAL BLOCKS.
C        ------------------------------------

         NXTMDN = HEAD(1)
 100     IF ( NXTMDN .GT. 0 ) THEN
            CURMDN = NXTMDN
            NXTMDN = INVP(CURMDN)
            MARKER(CURMDN) = MAXMRK
            INVP(CURMDN) = -NUM
            NUM = NUM + QSIZE(CURMDN)
            GOTO 100
         ENDIF

C        -----------------------------------------------
C        GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
C        -----------------------------------------------

         IF ( NUM .GT. NNODS )
     $      GOTO 600

      ENDIF

C     --------------------------------------------
C     INITIALIZE FOR THE NON-TRIVIAL ELIMINATIONS.
C     --------------------------------------------

      MARK   = 1
      CMINDG = 2
      HEAD(1) = 0

C     ================================
C     START THE MAIN ELIMINATION LOOP.
C     ================================

C     -----------------------------
C     FIND THE NEXT MINIMUM DEGREE.
C     -----------------------------

 200  IF ( HEAD(CMINDG) .LE. 0 ) THEN
            CMINDG = CMINDG + 1
            GOTO 200
         ENDIF

      LIMIT  = CMINDG + DELTA
      NUMELM = 0

C        ====================================
C        START THE MULTIPLE ELIMINATION STEP.
C        ====================================

C        -----------------------------------------
C        RETRIEVE THE CURRENT MINIMUM DEGREE NODE.
C        -----------------------------------------

 300     CURMDN = HEAD(CMINDG)

         IF ( CURMDN .LE. 0 ) THEN
            CMINDG = CMINDG + 1
            IF ( CMINDG .GT. LIMIT )
     $         GOTO 500
            GOTO 300
         ENDIF

C        -----------------------------------
C        REMOVE THE MINIMUM DEGREE NODE FROM
C        THE DOUBLY LINKED DEGREE STRUCTURE.
C        AND SET ITS INVERSE PERMUTATION
C        POSITION IN INVP.
C        -----------------------------------

         NXTMDN = INVP(CURMDN)
         HEAD(CMINDG) = NXTMDN
         IF ( NXTMDN .GT. 0 ) PERM(NXTMDN) = -CMINDG
         INVP(CURMDN) = -NUM

C        -----------------------------------------------
C        GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
C        -----------------------------------------------

         IF ( NUM + QSIZE(CURMDN) .GT. NNODS )
     $      GOTO 600

C        -----------------------
C        UPDATE THE NODE MARKER.
C        -----------------------

         MARK = MARK + 1

         IF ( MARK .GE. MAXMRK ) THEN

C           ------------------------
C           RESET THE MARKER VECTOR.
C           ------------------------

            MARK = 1
            DO 400 I = 1, NNODS
               IF ( MARKER(I) .LT. MAXMRK )
     $            MARKER(I) = 0
 400        CONTINUE
         ENDIF

C        ---------------------------------
C        UPDATE THE GRAPH REPRESENTATION
C        DUE TO THE ELIMINATION OF CURMDN.
C        ---------------------------------

         CALL SPRMM2 ( CURMDN, NELS, NNODS, NELNOD, MAXNEL, MARK,
     $                 MAXMRK, NUMELM, NXTUSE, NCOMPR, XFNODL, FNODL,
     $                 XFLNOD, FLNOD, HEAD, PERM, INVP, QSIZE, ELMLNK,
     $                 MARKER, STOREL )
         NUM = NUM + QSIZE(CURMDN)

C        -----------------------------------------------
C        GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
C        -----------------------------------------------

         IF ( NUM .GT. NNODS )
     $      GOTO 600

C        -------------------------------
C        GOTO A NEW NODE ELIMINATION IN
C        THIS MULTIPLE ELIMINATION STEP.
C        -------------------------------

         IF ( DELTA .GE. 0 )
     $      GOTO 300

C        -----------------------------------------------
C        GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
C        -----------------------------------------------

 500     IF ( NUM .GT. NNODS )
     $      GOTO 600

C        ----------------------------------
C        UPDATE THE DEGREE OF ALL THE NODES
C        WHICH ARE AFFECTED IN THE CURRENT
C        MULTIPLE ELIMINATION STEP, EXCEPT
C        THOSE WHICH ARE OUTMATCHED.
C        ----------------------------------

         CALL SPRMM3 ( NUMELM, CMINDG, MARK, NELS, NNODS,
     $                 NELNOD, MAXNEL, MAXMRK, XFNODL, FNODL,
     $                 XFLNOD, FLNOD, HEAD, INVP, PERM, QSIZE,
     $                 ELMLNK, MARKER, FLAG, STOREL )

C     ----------------------------------------
C     GOTO THE NEXT MULTIPLE ELIMINATION STEP.
C     ----------------------------------------

      GOTO 200

C     =============================
C     END OF MAIN ELIMINATION LOOP.
C     =============================

 600  CONTINUE

C     -----------------------------------
C     DETERMINE THE FINAL NODAL ORDERING.
C     -----------------------------------

      CALL SPRMM4 ( NNODS, PERM, INVP, QSIZE )

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRMMD
      END


      SUBROUTINE   SPRMM1   ( NELS  , NNODS , NELNOD, MAXNEL, MAXMRK,
     $                        XELNOD, ELNOD , XNODEL, NODEL , XFLNOD,
     $                        FLNOD , XFNODL, FNODL , HEAD  , INVP  ,
     $                        PERM  , QSIZE , ELMLNK, MARKER, FLAG  ,
     $                        MASK    )

C     ------------------------------------------------------------------

      INTEGER           NELS  , NNODS , NELNOD, MAXNEL, MAXMRK
      INTEGER           XELNOD(NELS+1)        , ELNOD(NELNOD)         ,
     $                  XNODEL(NNODS+1)       , NODEL(NELNOD)         ,
     $                  XFLNOD(NELS)          , FLNOD(MAXNEL)         ,
     $                  XFNODL(NNODS)         , FNODL(NNODS+NELNOD)   ,
     $                  HEAD(NNODS)           , INVP(NNODS)           ,
     $                  PERM(NNODS)           , QSIZE(NNODS)          ,
     $                  ELMLNK(NELS)          , MARKER(NNODS)         ,
     $                  FLAG(NNODS)           , MASK(NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRMM1
C
C --- TO INITIALIZE FOR SPRMMD.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NELS   : INTEGER
C              NUMBER OF ELEMENTS.
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     NELNOD : INTEGER
C              SIZE OF CONNECTIVITY ARRAYS.
C     MAXNEL : INTEGER
C              WORKSPACE FOR FLNOD. AT LEAST NELS+NELNOD.
C     (XELNOD,
C     ELNOD) : INTEGER(NELS), INTEGER(NELNOD)
C              ELEMENT-NODE CONNECTIVITY.
C     (XNODEL,
C     NODEL) : INTEGER(NNODS+1), INTEGER(NELNOD)
C              NODE-ELEMENT CONNECTIVITY.
C
C     ON EXIT
C     -------
C
C     XFLNOD : INTEGER(NELS)
C              POINTER INTO FLNOD FOR LISTS OF NODES CONNECTED TO THE
C              ELEMENTS.
C     FLNOD  : INTEGER(NELNOD+NELS)
C              LISTS OF NODES, WHERE EACH LIST IS HEADED BY ITS LENGTH.
C     XFNODL : INTEGER(NNODS)
C              POINTER INTO FNODL FOR LISTS OF ELEMENTS CONNECTED TO
C              THE NODES.
C     FNODL  : INTEGER(NELNOD+NNODS)
C              LISTS OF ELEMENTS, WHERE EACH LIST IS HEADED BY ITS
C              LENGTH.
C     (HEAD,
C     INVP,
C     PERM)  : 3*INTEGER(NNODS)
C              DOUBLY LINKED DEGREE STRUCTURE WITH INITIAL DEGREE OF
C              THE NODES.
C     QSIZE  : INTEGER(NNODS)
C              SUPERNODE SIZES INITIALIZED TO 1.
C     ELMLNK : INTEGER(NELS)
C              LIST OF ELEMENTS IS NOW THE ZERO LIST.
C     MARKER : INTEGER(NNODS)
C              INITIALIZED TO ZERO.
C     FLAG   : INTEGER(NNODS)
C              FLAGS ACTIVE NODES IN THE DEGREE UPDATE STEP.
C
C     WORKING ARRAYS
C     --------------
C
C     MASK   : INTEGER(NNODS)
C              Used to hold the clique sizes for nodes that are
C              candidates for merge.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           DEG   , ELEMNT, I     , IDUMMY, IEL   , INOD  ,
     $                  IP    , ISTRT , ISTOP , J     , JSTRT , JSTOP ,
     $                  K     , KSTRT , KSTOP , MERGED, NEWNOD, NEXT  ,
     $                  NODE  , NXTNOD, PREV

      INTRINSIC         ABS

C     ==================================================================

C     -----------
C     INITIALIZE.
C     -----------
      DO 10 I = 1, NNODS
         FLAG(I) = 0
         HEAD(I) = 0
         MARKER(I) = 0
         QSIZE(I) = 1
 10   CONTINUE
      DO 20 I = 1, NELS
         ELMLNK(I) = 0
 20   CONTINUE

C     ==================================================================

C     -------------------------------------
C     COMPUTE THE INITIAL DEGREE STRUCTURE.
C     -------------------------------------
      DO 300 NODE = 1, NNODS
         MARKER(NODE) = NODE
         ISTRT = XNODEL(NODE)
         ISTOP = XNODEL(NODE+1) - 1
         DEG = 1
         DO 200 I = ISTRT, ISTOP
            ELEMNT = NODEL(I)
            JSTRT = XELNOD(ELEMNT)
            JSTOP = XELNOD(ELEMNT+1) - 1
            DO 100 J = JSTRT, JSTOP
               NEWNOD = ELNOD(J)
               IF ( MARKER(NEWNOD) .NE. NODE )
     $         THEN
                  MARKER(NEWNOD) = NODE
                  DEG = DEG + 1
               ENDIF
 100        CONTINUE
 200     CONTINUE
         NEWNOD = HEAD(DEG)
         INVP(NODE) = NEWNOD
         HEAD(DEG)  = NODE
         IF ( NEWNOD .GT. 0 )
     $      PERM(NEWNOD) = NODE
         PERM(NODE) = -DEG
 300  CONTINUE

C     ---------------------
C     RESET MARKER TO ZERO.
C     ---------------------
      DO 400 I = 1, NNODS
         MARKER(I) = 0
 400  CONTINUE

C     ==================================================================

C     -----------------------------------------------------
C     BUILD THE GRAPH REPRESENTATION NEEDED BY THE ROUTINE.
C     -----------------------------------------------------
      IP = NELS + NELNOD
      DO 600 I = NELS, 1, -1
         JSTRT = XELNOD(I)
         JSTOP = XELNOD(I+1)-1
         DO 500 J = JSTOP, JSTRT, -1
            FLNOD(IP) = ELNOD(J)
            IP = IP - 1
 500     CONTINUE
         XFLNOD(I) = IP
         FLNOD(IP) = JSTOP + 1 - JSTRT
         IP = IP - 1
 600  CONTINUE

      IP = NNODS + NELNOD
      DO 800 I = NNODS, 1, -1
         JSTRT = XNODEL(I)
         JSTOP = XNODEL(I+1)-1
         DO 700 J = JSTOP, JSTRT, -1
            FNODL(IP) = NODEL(J)
            IP = IP - 1
 700     CONTINUE
         XFNODL(I) = IP
         FNODL(IP) = JSTOP + 1 - JSTRT
         IP = IP - 1
 800  CONTINUE

C     ==================================================================

C --- FIND THE STRUCTURALLY INDISTINGUISHABLE NODES ACCORDING TO
C     ALGORITHM 2.1 IN MY THESIS.
      MERGED = 0
      DO 1300 IEL = 1, NELS
         ISTRT = XFLNOD(IEL)
         ISTOP = ISTRT + FLNOD(ISTRT)
         DO 1200 I = ISTRT+1, ISTOP
            INOD = FLNOD(I)
            IF ( FLAG(INOD) .EQ. 0 ) THEN
               FLAG(INOD) = INOD
C ------------ SET NXTNOD =-INOD IN ORDER TO INITIATE THE LINKED
C              LIST OF NODES TO BE SCANNED.
               NXTNOD = -INOD
               JSTRT = XFNODL(INOD)
               JSTOP = JSTRT + FNODL(JSTRT)
               DO 1000 J = JSTRT+1, JSTOP
                  ELEMNT = FNODL(J)
                  KSTRT = XFLNOD(ELEMNT)
                  KSTOP = KSTRT + FLNOD(KSTRT)
                  DO 900 K = KSTRT+1, KSTOP
                     NODE = FLNOD(K)
                     IF ( FLAG(NODE) .EQ. 0 ) THEN
C --------------------- THIS NODE IS NOT YET VISITED, INITIALIZE THE
C                       GENERATED ELEMENT COUNTER, SET FLAG IN ORDER
C                       TO UPDATE LINKED LIST AND SET -NODE AS NEXT
C                       NODE.
                        MASK(NODE) = 1
                        FLAG(NODE) = NXTNOD
                        NXTNOD = -NODE
                     ELSEIF ( FLAG(NODE) .LT. 0 ) THEN
C --------------------- THIS NODE IS VISITED, INCREMENT ITS GENERATED
C                       ELEMENT COUNTER BY ONE.
                        MASK(NODE) = MASK(NODE) + 1
                     ENDIF
 900              CONTINUE
 1000          CONTINUE
C ------------ JSTOP IS USED AS CLIQUE NUMBER COUNTER FOR INOD.
               JSTOP = FNODL(JSTRT)
               DO 1100 IDUMMY = 1, NNODS
C --------------- DUMMY LOOP, WILL ONLY INCREMENT ACCORDING TO
C                 THE NUMBER OF MERGED NODES.
                  NODE = ABS(NXTNOD)
                  IF ( FLAG(NODE) .LT. 0 ) THEN
C ------------------ HAS NOT YET REACHED THE BOTTOM OF LIST. SET NXTNOD
C                    TO BE NEXT NODE IN LINKED LIST STORED IN FLAG.
                     NXTNOD = FLAG(NODE)
C ------------------ KSTOP IS USED AS CLIQUE NUMBER COUNTER FOR NODE.
                     KSTOP = FNODL(XFNODL(NODE))
                     IF ( ( MASK(NODE) .LT. KSTOP ).OR.
     $                    (      JSTOP .NE. KSTOP )      ) THEN
C --------------------- NODE CANNOT BE MERGED WITH INOD.
                        FLAG(NODE) = 0
                     ELSE
C --------------------- NODE CAN BE MERGED WITH INOD
C                       FIRST REMOVE NODE FROM THE DEGREE STRUCTURE.
                        PREV = PERM(NODE)
                        IF ( PREV .NE. 0 .AND. PREV .NE. -MAXMRK ) THEN
C ------------------------ IT IS IN THE DEGREE STRUCTURE, REMOVE IT.
                           NEXT = INVP(NODE)
                           IF ( NEXT .GT. 0 ) PERM(NEXT) = PREV
                           IF ( PREV .GT. 0 ) INVP(PREV) = NEXT
                           IF ( PREV .LT. 0 ) HEAD(-PREV) = NEXT
                        ENDIF
C --------------------- MERGE NODE WITH INOD.
                        QSIZE(INOD) = QSIZE(INOD) + QSIZE(NODE)
                        QSIZE(NODE) = 0
                        MARKER(NODE) = MAXMRK
                        INVP(NODE) = -INOD
                        XFNODL(NODE) = 0
C --------------------- ACCUMULATE THE NUMBER OF MERGED NODES
C                       AND MASK IT OFF IN FLAG.
                        MERGED = MERGED + 1
                        FLAG(NODE) = NODE
                     ENDIF
                  ELSE
C ------------------ END OF LIST, JUMP TO NEXT NODE.
                     GOTO 1150
                  ENDIF
 1100          CONTINUE
 1150          CONTINUE
            ENDIF
 1200    CONTINUE
 1300 CONTINUE

      DO 1400 I = 1, NNODS
         FLAG(I) = 0
 1400 CONTINUE

**DAM       WRITE ( 6, '(/A,I6)' ) ' Number of nodes     : ', NNODS
**DAM       WRITE ( 6, '(A,I6/)' ) ' New number of nodes : ', NNODS-MERGED

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRMM1
      END


      SUBROUTINE   SPRMM2   ( NODE  , NELS  , NNODS , NELNOD, MAXNEL,
     $                        MARK  , MAXMRK, NUMELM, NXTUSE, NCOMPR,
     $                        XFNODL, FNODL , XFLNOD, FLNOD , HEAD  ,
     $                        PERM  , INVP  , QSIZE , ELMLNK, MARKER,
     $                        STOREL                                  )

C     ------------------------------------------------------------------

      INTEGER           NODE  , NELS  , NNODS , NELNOD, MAXNEL, MARK  ,
     $                  MAXMRK, NUMELM, NXTUSE, NCOMPR
      INTEGER           XFNODL(NNODS)         , FNODL(NNODS+NELNOD)   ,
     $                  XFLNOD(NELS)          , FLNOD(MAXNEL)         ,
     $                  HEAD(NNODS)           , PERM(NNODS)           ,
     $                  INVP(NNODS)           , MARKER(NNODS)         ,
     $                  QSIZE(NNODS)          , ELMLNK(NELS)          ,
     $                  STOREL(NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRMM2
C
C --- TO UPDATE THE REPRESENTATION OF THE GRAPH DUE TO THE ELIMINATION
C     OF NODE. THE ROUTINE WORKS ON THE IMPLIVIT REPRESENTATION OF THE
C     GRAPH. THE IMPLICIT GRAPH STRUCTURE IS REPRESENTED BY THE ARRAYS
C     (XFNODL,FNODL) $ (XFLNOD,FLNOD).
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NODE   : INTEGER
C              NODE TO BE ELIMINATED.
C     NELS   : INTEGER
C              NUMBER OF FINITE ELEMENTS IN THE MESH.
C     NNODS  : INTEGER
C              NUMBER OF NODES IN THE MESH.
C     NELNOD : INTEGER
C              SIZE OF INDICES STORED IN FNODL AND FLNOD AT FIRST ENTRY.
C     MAXNEL : INTEGER
C              SIZE OF FLNOD, WHERE MAXNEL .GE. NELS+NELNOD.
C     MARK   : INTEGER
C              STAGE COUNTER, FACILITATES MARKING NODES.
C     MAXMRK : INTEGER
C              A LARGEST POSSIBLE INTEGER, HOWEVER,
C              MAXMRK .GT. MAX(NELS,NNODS) IS ENOUGH.
C
C     UPDATED PARAMETERS
C     ------------------
C
C     NUMELM : INTEGER
C              NUMBER OF ACTIVE NEW GENERATED ELEMENTS AT THIS STAGE.
C     NXTUSE : INTEGER
C              THE NEXT FREE LOCATION OF FLNOD.
C              ON THE FIRST ENTRY, IT MUST BE SET TO NELNOD+NELS+1.
C              UPDATED BY THE ROUTINE TO ALWAYS POINT TO THE FIRST
C              FREE LOCATION IN FLNOD.
C     NCOMPR : INTEGER
C              UPDATED WITH THE NUMBER OF WORKSPACE COMPRESSIONS NEEDED
C              IN THIS GRAPH UPDATE.
C     XFNODL : INTEGER(NNODS)
C              POINTER INTO THE GENERATED ELEMENTS FOR EACH NODE.
C              ABSORBED AND ELIMINATED NODES HAVE XFNODL(I)=0.
C              ON EXIT FOR AN ELIMINATED NODE I, XFNODL(I)=0.
C     FNODL  : INTEGER(NELNOD+NNODS)
C              THE LIST OF GENERATED ELEMENTS FOR EACH NODE HEADED BY
C              THE LENGTH OF THE LIST.
C              ON EXIT REFERENCE TO ABSORBED ELEMENTS IS REMOVED.
C     XFLNOD : INTEGER(NELS)
C              POINTER INTO THE LIST OF NODES FOR EACH ELEMENT.
C              IF AN ELEMENT I IS ABSORBED, THEN XFLNOD(I)=0.
C     FLNOD  : INTEGER(MAXNEL)
C              THE NODES IN A GENERATED ELEMENT HEADED BY THE LENGTH
C              OF THE LIST.
C     HEAD   : INTEGER(NNODS)
C              HEAD OF DEGREE LIST.
C     PERM   : INTEGER(NNODS)
C              PREVIOUS NODE.
C     INVP   : INTEGER(NNODS)
C              NEXT NODE.
C     QSIZE  : INTEGER(NNODS)
C              SIZE OF THE SUPERNODES.
C     ELMLNK : INTEGER(NELS)
C              LIST OF GENERATED ELEMENTS WHERE SOME NODES NEED
C              DEGREE UPDATE.
C     MARKER : INTEGER(NNODS)
C              MARKER VECTOR FOR EFFICIENT MASK OF TOUCHED NODES.
C
C     ON EXIT
C     -------
C
C     NO PURE EXIT PARAMETERS.
C
C     WORKING ARRAYS
C     --------------
C
C     STOREL : INTEGER(NNODS)
C              TEMPORARY STORAGE OF THE NODES BELONGING TO THE NEW
C              GENERATED ELEMENT.
C
C     SUBPROGRAMS
C     -----------
C
C     SPRMM5
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

      INTEGER           COUNT , ELEMNT, ELMSZE, FSTELM, GEFREE, I     ,
     $                  IP    , IPFSTL, ISTRT , ISTOP , J     , JSTRT ,
     $                  JSTOP , MFLNOD, NEXT  , PREV  , RCHLOC, RCHNOD,
     $                  RCHSZE
      EXTERNAL          SPRMM5

C     ------------------------------------------------------------------
C     COUNT COUNTS THE NUMBER OF GENERATED ELEMENTS CONNECTED TO A NODE.
C     ELEMNT ELEMENT INDEX.
C     FSTELM THE FIRST ELEMENT IN THE LIST OF GENERATED ELEMENTS FOR
C            NODE. USED AS REPRESENTATIVE FOR THE NEW ELEMENT.
C     FSTIP FIRST POSITION IN FLNOD FOR AN ELEMENT WHEN LIST IS COM-
C           PRESSED.
C     GEFREE STORES THE SIZE OF FREE SPACE AT THE END OF FLNOD.
C     IP USED AS ARRAY POINTER.
C     I, J ARE LOOP COUNTERS
C     ISTRT, ISTOP ARE ARRAY/LIST BOUNDS FOR LOOP INDEX I.
C     JSTRT, JSTOP ARE ARRAY/LIST BOUNDS FOR LOOP INDEX J.
C     MFLNOD IN CALLS TO SPRMM5, LAST LOCATION OF FLNOD IN USE.
C     NEXT NEXT NODE IN DEGREE LIST.
C     PREV PREVIOUS NODE IN THE DEGREE LIST.
C     RCHLOC POINTER INTO THE GENERATED ELEMENT.
C     RCHNOD NODE IN THE NEW GENERATED ELEMENT.
C     RCHSZE SIZE OF THE NEW GENERATED ELEMENT.
C     ------------------------------------------------------------------

C     ==================================================================

C     ---------------------------------------------------
C     FETCH THE GENERATED ELEMENT AND STORE IT IN STOREL.
C     ---------------------------------------------------

      ISTRT = XFNODL(NODE)
      ISTOP = ISTRT + FNODL(ISTRT)
      RCHSZE = 0
      ELMSZE = 0
      FSTELM = 0
      IPFSTL = 0

C     ------------------------
C     MARK NODE AS ELIMINATED.
C     ------------------------

      MARKER(NODE) = MAXMRK
      PERM(NODE) = -MAXMRK
      XFNODL(NODE) = 0
C
      DO 200 I = ISTRT+1, ISTOP
         ELEMNT = FNODL(I)
         JSTRT = XFLNOD(ELEMNT)
         JSTOP = JSTRT + FLNOD(JSTRT)
         IF ( FLNOD(JSTRT) .GT. ELMSZE )
     $   THEN
            ELMSZE = FLNOD(JSTRT)
            FSTELM = ELEMNT
            IPFSTL = JSTRT
         ENDIF
         XFLNOD(ELEMNT) = 0
         DO 100 J = JSTRT+1, JSTOP
            RCHNOD = FLNOD(J)
            IF ( MARKER(RCHNOD) .LT. MARK )
     $      THEN

               MARKER(RCHNOD) = MARK

C              ---------------------------------
C              THE NODE HAS NOT YET BEEN MERGED.
C              ---------------------------------

               PREV = PERM(RCHNOD)
               IF ( PREV .NE. 0 .AND. PREV .NE. -MAXMRK )
     $         THEN

C                 -------------------------------
C                 THE NODE IS STILL IN THE DEGREE
C                 STRUCTURE, REMOVE IT.
C                 -------------------------------

                  NEXT = INVP(RCHNOD)
                  IF ( NEXT .GT. 0 ) PERM(NEXT) = PREV
                  IF ( PREV .GT. 0 ) INVP(PREV) = NEXT
                  IF ( PREV .LT. 0 ) HEAD(-PREV) = NEXT

C                 ------------------------------
C                 AND MARK IT FOR DEGREE UPDATE.
C                 ------------------------------

                  PERM(RCHNOD) = 0
               ENDIF
               IF ( XFNODL(RCHNOD) .GT. 0 )
     $         THEN

C                 ---------------------
C                 NODE IS STILL ACTIVE.
C                 ---------------------

                  RCHSZE = RCHSZE + 1
                  STOREL(RCHSZE) = RCHNOD
                  PERM(RCHNOD) = 0
               ENDIF
            ENDIF
 100     CONTINUE
 200  CONTINUE

C     ----------------------------------------------------------
C     FSTELM IS THE REPRESENTATIVE OF THE NEW GENERATED ELEMENT.
C     STOREL NOW CONTAINS THE NODES OF THE GENERATED ELEMENT
C            WITHOUT DUPLICATE REFERENCES.
C     RCHSZE NOW STORES THE SIZE OF THE NEW GENERATED ELEMENT.
C     XFLNOD IS ZEROED FOR EACH ABSORBED ELEMENT.
C     ----------------------------------------------------------

C     -------------------------
C     UPDATE NUMELM AND ELMLNK.
C     -------------------------

      NUMELM = NUMELM + 1
      ELMLNK(NUMELM) = FSTELM

C     ==================================================================

      IF ( RCHSZE .LE. ELMSZE )
     $THEN

C        -------------------------------------------------------
C        STORE THE NEW ELEMENT IN THE LOCATIONS USED FOR FSTELM.
C        -------------------------------------------------------

         XFLNOD(FSTELM) = IPFSTL
         FLNOD(IPFSTL) = RCHSZE
         ISTRT = IPFSTL
         ISTOP = ISTRT + RCHSZE
         IP = 0
         DO 300 I = ISTRT+1, ISTOP
            IP = IP + 1
            FLNOD(I) = STOREL(IP)
 300     CONTINUE
      ELSE

C        --------------------------------------
C        STORE THE ELEMENT AT THE END OF FLNOD.
C        --------------------------------------

         GEFREE = MAXNEL - NXTUSE
         IF ( GEFREE .LT. RCHSZE )
     $   THEN

C           ---------------
C           COMPRESS FLNOD.
C           ---------------

            MFLNOD = NXTUSE-1
            CALL SPRMM5 ( NELS, XFLNOD, FLNOD, MFLNOD, NXTUSE, NCOMPR )
         ENDIF

C        ------------------
C        STORE THE ELEMENT.
C        ------------------

         XFLNOD(FSTELM) = NXTUSE
         FLNOD(NXTUSE) = RCHSZE
         ISTRT = NXTUSE
         ISTOP = ISTRT + RCHSZE
         IP = 0
         DO 400 I = ISTRT+1, ISTOP
            IP = IP + 1
            FLNOD(I) = STOREL(IP)
 400     CONTINUE
         NXTUSE = ISTOP + 1
      ENDIF

C     --------------------------------------------------------------
C     NXTUSE NOW POINTS TO THE NEXT FREE LOCATION AFTER NEW ELEMENT.
C     XFLNOD HAS GOT THE POINTER.
C     FLNOD STORES THE NODES.
C     --------------------------------------------------------------

C     ==================================================================

C     -------------------------------------------------------------
C     GO THROUGH THE NEW ELEMENT, MARK THE NODES FOR DEGREE UPDATE,
C     AND LOOK FOR NODES THAT ARE INTERNAL TO THE NEW ELEMENT.
C     -------------------------------------------------------------

      ISTRT = XFLNOD(FSTELM)
      ISTOP = ISTRT + FLNOD(ISTRT)
      RCHLOC = ISTRT
      RCHSZE = 0
      DO 600 I = ISTRT+1, ISTOP
         RCHNOD = FLNOD(I)

C           ---------------------------
C           THE NODE IS NOT ELIMINATED.
C           ---------------------------

            JSTRT = XFNODL(RCHNOD)
            JSTOP = JSTRT + FNODL(JSTRT)
            COUNT = 0
            IP    = JSTRT
            DO 500 J = JSTRT+1, JSTOP
               ELEMNT = FNODL(J)
               IF ( XFLNOD(ELEMNT) .GT. 0 .AND.
     $              ELEMNT .NE. FSTELM          )
     $         THEN
                  COUNT = COUNT + 1
                  IP = IP + 1
                  FNODL(IP) = ELEMNT
               ENDIF
 500        CONTINUE

C           ------------------------------------------
C           THE NODE IS BY DEFAULT ATTACHED TO FSTELM.
C           ------------------------------------------

            COUNT = COUNT + 1
            IP = IP + 1
            FNODL(IP) = FSTELM
            IF ( COUNT .EQ. 1 )
     $      THEN

C              ------------------------------------------
C              THE NODE RCHNOD IS INTERIOR TO THE ELEMENT
C              AND MAY BE ELIMINATED TOGETHER WITH NODE.
C              ------------------------------------------

               QSIZE(NODE) = QSIZE(NODE) + QSIZE(RCHNOD)
               QSIZE(RCHNOD) = 0
               INVP(RCHNOD) = -NODE
               PERM(RCHNOD) = -MAXMRK
               MARKER(RCHNOD) = MAXMRK
               XFNODL(RCHNOD) = 0
            ELSE

C              -------------------------------------
C              THE NODE IS AT LEAST ON TWO ELEMENTS.
C              -------------------------------------
               RCHLOC = RCHLOC + 1
               RCHSZE = RCHSZE + 1
               FLNOD(RCHLOC) = RCHNOD
               FNODL(JSTRT) = COUNT
            ENDIF

 600  CONTINUE

C     -------------------------------------
C     UPDATE THE SIZE OF THE REACHABLE SET.
C     -------------------------------------

      FLNOD(ISTRT) = RCHSZE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRMM2
      END


      SUBROUTINE   SPRMM3   ( NUMELM, CMINDG, MARK  , NELS  , NNODS ,
     $                        NELNOD, MAXNEL, MAXMRK, XFNODL, FNODL ,
     $                        XFLNOD, FLNOD , HEAD  , INVP  , PERM  ,
     $                        QSIZE , ELMLNK, MARKER, FLAG  , STOREL  )

C     ------------------------------------------------------------------

      INTEGER           NUMELM, CMINDG, MARK  , NELS  , NNODS , NELNOD,
     $                  MAXNEL, MAXMRK
      INTEGER           XFNODL(NNODS)         , FNODL(NNODS+NELNOD)   ,
     $                  XFLNOD(NELS)          , FLNOD(MAXNEL)         ,
     $                  HEAD(NNODS)           , INVP(NNODS)           ,
     $                  PERM(NNODS)           , QSIZE(NNODS)          ,
     $                  ELMLNK(NELS)          , MARKER(NNODS)         ,
     $                  FLAG(NNODS)           , STOREL(NNODS)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRMM3
C
C --- TO PERFORM DEGREE UPDATE. IT ALSO SEARCH FOR INDISTINGUISHABLE AND
C     OUTMATCHED NODES.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NUMELM : INTEGER
C              HEAD TO THE LIST OF GENERALIZED ELEMENTS OF WHICH THE
C              NODES NEED DEGREE UPDATE.
C     CMINDG : INTEGER
C              CURRENT MINIMUM DEGREE.
C     MARK   : INTEGER
C              USED TO PREVENT A NODE TO BE TREATED MORE THAN ONCE IN
C              EACH STEP OF THE ROUTINE.
C     NELS   : INTEGER
C              NUMBER OF ELEMENTS.
C     NNODS  : INTEGER
C              NUMBER OF NODES.
C     NELNOD : INTEGER
C              SIZE OF ELIMINATION GRAPH REPRESENTATION.
C     MAXNEL : INTEGER
C              SIZE OF FLNOD.
C     MAXMRK : INTEGER
C              USED AS UPPER LIMIT FOR MARKING NODES, ANY NUMBER LARGER
C              THAN NNODS IS OK.
C     (XFNODL,
C     FNODL,
C     XFLNOD,
C     FLNOD) : INTEGER(NNODS), INTEGER(NNODS+NELNOD), INTEGER(NELS),
C              INTEGER(MAXNEL)
C              REPRESENTATION OF THE ELIMINATION GRAPH.
C     QSIZE  : INTEGER(NNODS)
C              SIZE OF THE SUPERNODES.
C     ELMLNK : INTEGER(NNODS)
C              LIST OF GENERALIZED ELEMENTS OF WHICH THE NODES NEED
C              DEGREE UPDATE.
C     MARKER : INTEGER(NNODS)
C              FACILITATES MARKING NODES.
C     FLAG   : INTEGER(NNODS)
C              MARKS ACTIVE NODES IN THE DEGREE UPDATE STEP. THE ROUTINE
C              ASSUMES THAT THE ARRAY IS ZERO ON ENTRY. THE ARRAY IS
C              RESET TO ZERO ON EXIT.
C
C     ON EXIT
C     -------
C
C     (HEAD,
C     INVP,
C     PERM)  : 3*INTEGER(NNODS)
C              UPDATED DEGREE STRUCTURE.
C     XFNODL : INTEGER(NNODS)
C              FOR AN ABSORBED NODE I, XFNODL(I)=0.
C     QSIZE  : INTEGER(NNODS)
C              UPDATED SIZE OF THE SUPERNODES.
C
C     WORKING ARRAYS
C     --------------
C
C     STOREL : INTEGER(NNODS)
C              LIST OF NODES WHICH NEED DEGREE UPDATE.
C
C     SUBPROGRAMS/FUNCTIONS
C     ---------------------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER           CLQSZE, DEG0  , DEG1  , ELM   , I     , IEL   ,
     $                  IP    , ISTOP , ISTRT , J     , JSTOP , JSTRT ,
     $                  K     , KSTOP , KSTRT , NELM  , NEWELM, NEWNOD,
     $                  NEXT  , NODE  , PREV  , QLX   , STORED

C     ==================================================================

      DO 800 IEL = 1, NUMELM

C        ---------------------------------------------
C        START THE DEGREE UPDATE FOR THE NODES ON ELM.
C        ---------------------------------------------

         ELM = ELMLNK(IEL)
         DEG0 = 0
         ISTRT = XFLNOD(ELM)
         ISTOP = ISTRT + FLNOD(ISTRT)
         NELM = 0
         IP = ISTRT

C        -----------------------------------------
C        FIRST INITIALIZE FLAG AND COMPRESS FLNOD.
C        -----------------------------------------

         DO 100 I = ISTRT+1, ISTOP
            NODE = FLNOD(I)
            IF ( XFNODL(NODE) .GT. 0 )
     $      THEN
               DEG0 = DEG0 + QSIZE(NODE)
               NELM = NELM + 1
               IP = IP + 1
               FLNOD(IP) = NODE
               FLAG(NODE) = 1
            ENDIF
 100     CONTINUE

C        -----------------------------------------------------
C        RUN THROUGH THE ELEMENT, SEARCH FOR INDISTINGUISHABLE
C        AND OUTMATCHED NODES, AND FINALLY UPDATE THE DEGREE.
C        -----------------------------------------------------

         FLNOD(ISTRT) = NELM
         ISTOP = ISTRT + NELM
         DO 600 I = ISTRT+1, ISTOP
            NODE = FLNOD(I)

            IF ( PERM(NODE) .EQ. 0 )
     $      THEN

C              ------------------------------
C              THIS NODE NEEDS DEGREE UPDATE.
C              ------------------------------

               DEG1 = DEG0
               MARK = MARK + 1
               IF ( MARK .GE. MAXMRK )
     $         THEN

C                 ------------------------
C                 RESET THE MARKER VECTOR.
C                 ------------------------

                  MARK = 1
                  DO 200 J = 1, NNODS
                     IF ( MARKER(J) .LT. MAXMRK )
     $                  MARKER(J) = 0
 200              CONTINUE
               ENDIF
               MARKER(NODE) = MARK

               JSTRT = XFNODL(NODE)
               JSTOP = JSTRT + FNODL(JSTRT)
               STORED = 0
               DO 400 J = JSTRT+1, JSTOP

C                 --------------------------------------
C                 RUN THROUGH THE REACHABLE SET OF NODE.
C                 --------------------------------------

                  NEWELM = FNODL(J)
                  IF ( NEWELM .NE. ELM )
     $            THEN

                     KSTRT = XFLNOD(NEWELM)
                     IF ( KSTRT .GT. 0 )
     $               THEN
                        KSTOP = KSTRT + FLNOD(KSTRT)
                        DO 300 K = KSTRT + 1, KSTOP
                           NEWNOD = FLNOD(K)
                           IF ( XFNODL(NEWNOD) .GT. 0 .AND.
     $                          NEWNOD .NE. NODE )
     $                     THEN

                              IF ( FLAG(NEWNOD) .GE. 1 )
     $                        THEN

C                                ---------------------------
C                                NEWNOD IS ON ELM. COUNT ITS
C                                APPEARANCES TO SEE IF ITS
C                                INDISTINGUISHABLE OR OUT-
C                                MATCHED BY NODE.
C                                ---------------------------

                                 IF ( FLAG(NEWNOD) .EQ. 1 )
     $                           THEN
                                    FLAG(NEWNOD) = 2
                                    STORED = STORED + 1
                                    STOREL(STORED) = NEWNOD
                                 ELSE
                                    FLAG(NEWNOD) = FLAG(NEWNOD) + 1
                                 ENDIF
                              ELSEIF ( MARKER(NEWNOD) .LT. MARK )
     $                        THEN

C                                --------------------------
C                                THIS NODE IS NOT ON ELM
C                                ACCUMULATE ITS SIZE IN THE
C                                DEGREE COUNT AND MARK IT.
C                                --------------------------
C
                                 DEG1 = DEG1 + QSIZE(NEWNOD)
                                 MARKER(NEWNOD) = MARK
                              ENDIF

                           ENDIF
 300                    CONTINUE
                     ENDIF
                  ENDIF
 400           CONTINUE
               IF ( STORED .GT. 0 )
     $         THEN

C                 ------------------------------
C                 QLX IS NOW SET TO THE SIZE OF
C                 THE CLIQUE MEMBERSHIP OF NODE.
C                 ------------------------------

                  QLX = FNODL(JSTRT)

                  DO 500 J = 1, STORED
                     NEWNOD = STOREL(J)

                     IF ( FLAG(NEWNOD) .EQ. QLX )
     $               THEN

C                       ----------------------------------
C                       THIS NODE HAS APPEARED QLX TIMES
C                       WHEN RUNNING THROUGH THE REACHABLE
C                       SET OF NODE. THIS MEANS THAT IT IS
C                       EITHER INDISTINGUISHABLE FROM NODE
C                       OR OUTMATCHED BY NODE.
C                       ----------------------------------

                        PREV = PERM(NEWNOD)
                        IF ( PREV .NE. 0 .AND. PREV .NE. -MAXMRK ) THEN

C                          -----------------------------------------
C                          IT IS IN THE DEGREE STRUCTURE, REMOVE IT.
C                          -----------------------------------------

                           NEXT = INVP(NEWNOD)
                           IF ( NEXT .GT. 0 ) PERM(NEXT) = PREV
                           IF ( PREV .GT. 0 ) INVP(PREV) = NEXT
                           IF ( PREV .LT. 0 ) HEAD(-PREV) = NEXT
                        ENDIF

                        CLQSZE = FNODL(XFNODL(NEWNOD))
                        IF ( CLQSZE .EQ. QLX )
     $                  THEN

C                          ----------------------------------------
C                          THE NODE IS INDISTINGUISHABLE FROM NODE.
C                          MERGE IT WITH NODE.
C                          ----------------------------------------

                           QSIZE(NODE) = QSIZE(NODE) + QSIZE(NEWNOD)
                           QSIZE(NEWNOD) = 0
                           MARKER(NEWNOD) = MAXMRK
                           INVP(NEWNOD) = -NODE
                           XFNODL(NEWNOD) = 0
                           FLAG(NEWNOD) = 0
                        ELSE

C                          -------------------------------
C                          THE NODE IS OUTMATCHED BY NODE.
C                          -------------------------------

                           FLAG(NEWNOD) = 1
                        ENDIF

C                       ------------------------------
C                       FOR BOTH INDISTINGUISHABLE AND
C                       OUTMATCHED NODES MARK IT TO
C                       AVOID DEGREE UPDATE OF OUTM.
C                       ------------------------------

                        PERM(NEWNOD) = -MAXMRK
                     ELSE

C                       -------------------------------------
C                       THE NODE IS NEITHER INDISTINGUISHABLE
C                       NOR OUTMATCHED FROM/BY NODE.
C                       -------------------------------------

                        FLAG(NEWNOD) = 1
                     ENDIF

 500              CONTINUE
               ENDIF

C              ----------------------------------
C              DO THE FINAL DEGREE COUNT OF NODE,
C              AND PUT IT INTO THE DEGREE LIST.
C              ----------------------------------

               DEG1 = DEG1 - QSIZE(NODE) + 1
               NEXT = HEAD(DEG1)
               INVP(NODE) = NEXT
               IF ( NEXT .GT. 0 ) PERM(NEXT) = NODE
               HEAD(DEG1) = NODE
               PERM(NODE) = -DEG1
               IF ( DEG1 .LT. CMINDG ) CMINDG = DEG1

            ENDIF
 600     CONTINUE

C        ---------------------------
C        FINALLY RESET FLAG TO ZERO,
C        AND COMPRESS FLNOD.
C        ---------------------------

         NELM = 0
         IP = ISTRT
         DO 700 I = ISTRT + 1, ISTOP
            NODE = FLNOD(I)
            IF ( XFNODL(NODE) .GT. 0 )
     $      THEN
               NELM = NELM + 1
               IP = IP + 1
               FLNOD(IP) = NODE
               FLAG(NODE) = 0
            ENDIF
 700     CONTINUE
         FLNOD(ISTRT) = NELM

 800  CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRMM3
      END


      SUBROUTINE   SPRMM4   ( NNODS , PERM  , INVP  , QSIZE )

C     ==================================================================

      INTEGER           NNODS

      INTEGER           PERM  (NNODS)         , INVP  (NNODS)         ,
     $                  QSIZE (NNODS)

C     ==================================================================
C
C     PURPOSE
C     -------
C
C     SPRMM4 : TO PERFORM THE FINAL STEP IN PRODUCING THE PERMUTATION
C              AND INVERSE PERM VECTORS IN THE MULTIPLE MINIMUM DEGREE
C              ORDERING ALGORITHM.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : MNT. XX, 199X (ACD)
C
C
C     ON ENTRY
C     --------
C
C     NNODS  : INTEGER
C              NUMBER OF NODES IN THE GRAPH.
C     QSIZE  : INTEGER(NNODS)
C              SIZE OF SUPERNODES AT ELIMINATION.
C
C     ON EXIT
C     -------
C
C     INVP   : INTEGER(NNODS)
C              INVERSE PERM VECTOR.  ON INPUT, IF QSIZE(NODE)=0, THEN
C              NODE HAS BEEN MERGED INTO THE NODE -INVP(NODE);
C              OTHERWISE, -INVP(NODE) IS ITS INVERSE LABELLING.
C     PERM   : INTEGER(NNODS)
C              THE PERMUTATION VECTOR.
C
C     ==================================================================

      INTEGER           NODE  , FATHER, NUM   , ROOT  , NEXTF , NQSIZE

C     ==================================================================

      DO 100 NODE = 1, NNODS
         NQSIZE = QSIZE(NODE)
         IF ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
         IF ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
  100 CONTINUE

C     ---------------------------------------
C     FOR EACH NODE WHICH HAS BEEN MERGED, DO
C     ---------------------------------------

      DO 500 NODE = 1, NNODS
         IF ( PERM(NODE) .GT. 0 )  GOTO 500

C        -----------------------------------------
C        TRACE THE MERGED TREE UNTIL ONE WHICH HAS
C        NOT BEEN MERGED, CALL IT ROOT
C        -----------------------------------------

         FATHER = NODE

  200    IF ( PERM(FATHER) .GT. 0 )  GOTO 300
            FATHER = - PERM(FATHER)
            GOTO 200

C        -----------------------
C        NUMBER NODE AFTER ROOT.
C        -----------------------

  300    ROOT       = FATHER
         NUM        = PERM(ROOT) + 1
         INVP(NODE) = - NUM
         PERM(ROOT) = NUM

C        ------------------------
C        SHORTEN THE MERGED TREE.
C        ------------------------

         FATHER = NODE

  400    NEXTF = - PERM(FATHER)
         IF ( NEXTF .LE. 0 )  GOTO 500
            PERM(FATHER) = - ROOT
            FATHER       = NEXTF
            GOTO 400

  500 CONTINUE

C     ----------------------
C     READY TO COMPUTE PERM.
C     ----------------------

      DO 600 NODE = 1, NNODS
         NUM        = - INVP(NODE)
         INVP(NODE) = NUM
         PERM(NUM)  = NODE
  600 CONTINUE

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRMM4
      END


      SUBROUTINE SPRMM5(N, IPE, IW, LW, IWFR, ICOMP)

C COMPRESS LISTS HELD BY SPRMM2 IN IW AND ADJUST POINTERS
C     IN IPE TO CORRESPOND.
C     CREATED : JAN. 11, 1994 (ACD)

      INTEGER N, LW, IWFR, ICOMP
      INTEGER IPE(N)
      INTEGER IW(LW)

C N IS THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
C     ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
C IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
C     LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
C LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
C     LOCATION IN IW.

      INTEGER I, IR, K, K1, K2, LWFR

      ICOMP = ICOMP + 1

C PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
C     LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
C     -(LIST NUMBER).

      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE

C COMPRESS
C IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
C LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.

      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70

C SEARCH FOR THE NEXT NEGATIVE ENTRY.
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70

C PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
C     AND PREPARE TO COPY LIST.
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50

C COPY LIST TO NEW POSITION.
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END


      SUBROUTINE   SPRCMD
     $( NELS  , NEQ   , NNODS , NELNOD, MAXNEL, MAXMRK, DELTA , LWEIGH,
     $  NCOMPR, XELNOD, ELNOD , XNODEL, NODEL , XNODES, PERM  , INVP  ,
     $  XFLNOD, FLNOD , XFNODL, FNODL , DEGREE, QSIZE , VWGHT , ELMLNK,
     $  FLAG  , STOREL, MARKER, SEPNUM, SEPSZE, SEPCNT, SEPARR, LPU    )

*     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER
     $  NELS  , NEQ   , NNODS , NELNOD, MAXNEL, MAXMRK, DELTA , NCOMPR,
     $  SEPNUM, SEPSZE, SEPCNT, LPU

      LOGICAL
     $  LWEIGH

      INTEGER
     $  XELNOD(NELS+1), ELNOD(NELNOD)      , XNODEL(NNODS+1),
     $  NODEL(NELNOD) , XNODES(NNODS+1)    , PERM(NNODS)    ,
     $  INVP(NNODS)   , XFLNOD(NELS)       , FLNOD(MAXNEL)  ,
     $  XFNODL(NNODS) , FNODL(NNODS+NELNOD),
     $  DEGREE(NEQ)   , QSIZE(NNODS)       , VWGHT(NNODS)   ,
     $  ELMLNK(NELS)  , FLAG(NNODS)        , STOREL(NNODS)  ,
     $  MARKER(NNODS) , SEPARR(NNODS)

*     ------------------------------------------------------------------
*
*     PURPOSE
*     -------
*
*     SPRCMD
*
* --- TO PERFORM A MINIMUM DEGREE PERMUTATION. THE ROUTINE UTILIZES THE
*     INHERENT CLIQUE STRUCTURE OF THE FINITE ELEMENT MESH.
*     This version may do also a constrained ordering of the nodes.
*
*     CREATED   : May  01, 1997 (ACD)
*     REVISIONS : MNT. XX, 199X (ACD)
*
*     ON ENTRY
*     --------
*
*     NELS   : INTEGER
*              NUMBER OF ELEMENTS.
*     NNODS  : INTEGER
*              NUMBER OF NODES.
*     NELNOD : INTEGER
*              SIZE OF CONNECTIVITY ARRAYS.
*     MAXNEL : INTEGER
*              SIZE OF ARRAY FLNOD. THE SIZE MUST BE AT LEAST
*              NELS+NELNOD, BUT THE NUMBER OF WORKSPACE COMPRESSIONS
*              WILL DECREASE IF THE SIZE IS INCREASED.
*     MAXMRK : INTEGER
*              USED AS UPPER LIMIT FOR MARKING NODES AND ELEMENTS, ANY
*              NUMBER LARGER THAN MAX(NNODS,NELS)+1 IS OK.
*     DELTA  : INTEGER
*              TOLERANCE FOR MULTIPLE ELIMINATION.
*              .LT. 0 : WITHOUT MULTIPLE ELIMINATION.
*              .EQ. 0 : WITH MULTIPLE ELIMINATION.
*              .GT. 0 : WITH MULTIPLE ELIMINATION AND RELAXED MINIMUM
*                       DEGREE CONDITION.
*     LWEIGH : Logical
*              Option for wieghted or standard minimum degree selection.
*              = TRUE : Weighted selection.
*              = FALSE: Standard selection.
*     (XELNOD,
*     ELNOD) : INTEGER(NELS+1), INTEGER(NELNOD)
*              ELEMENT-NODE CONNECTIVITY.
*     (XNODEL,
*     NODEL) : INTEGER(NNODS+1), INTEGER(NELNOD)
*              NODE-ELEMENT CONNECTIVITY.
*     SEPNUM : Integer
*              Type of ordering to be performed.
*              =    0 : Order the nodes that are not members of the
*                       separator set.
*                       The nodes in the separator sets are numbered
*                       as follows: If there is SEPCNT separator sets
*                       the separators are ordered backwards from
*                       1, 2, ... SEPCNT. The nodes in a separator set
*                       are ordered backward initial order.
*                       The only option for SEPSZE = 0.
*              =    1 : Order the nodes of the separator set as well.
*                       The non-separator nodes are ordered first and
*                       then we continue with a mmd order also for
*                       the separator nodes.
*                   DO NOT USE THIS OPTION.
*              =    2 : Order the nodes of the separator set as well.
*                       The non-separator nodes are ordered first and
*                       the separators are then merged to supervariab-
*                       les, and these supervariables are then ordered
*                       according to mmd.
*                   DO NOT USE THIS OPTION.
*
*     SEPSZE : Integer
*              The number of nodes in the separator sets. The separator
*              sets are distinguished by means of the array SEPARR.
*     SEPCNT : Integer
*            - Second pass for hybrid nested dissection & mmd:
*              The number of separator sets that have been computed
*              during the incomplete nested dissection step.
*            - Superelement Technique:
*              The number of separator sets as defined by user, in
*              normal this should be one.
*     SEPARR : Integer(NNODS)
*            - Second pass for hybrid nested dissection & mmd:
*              If there are SEPCNT separator sets and i is a number
*              such that 0<i=<SEPCNT, then SEPARR(v) = i for all the
*              nodes v that belong to separator set i. This is the
*              only way to find the separator sets. All the nodes u
*              that are not in a separator set have SEPARR(u)=<0.
*            - Superelement Technique:
*              The same rule as above applies for the nodes. We note
*              that SEPCNT=1 is the normal situation in a user defined
*              superelement assembly.
*
*     ON EXIT
*     -------
*
*     NCOMPR : INTEGER
*              NUMBER OF WORKSPACE COMPRESSIONS NEEDED IN SPRMD2.
*     (PERM,
*     INVP)  : 2*INTEGER(NNODS)
*              THE PERMUTATION.
*              DURING EXECUTION THEY ARE USED AS WORKSPACE AS EXPLAINED
*              BELOW.
*
*     WORKING ARRAYS
*     --------------
*
*     XFLNOD : INTEGER(NELS)
*              STORES THE POINTERS INTO THE LISTS OF NODES WHICH ARE
*              CONNECTED TO THE ELEMENT. IF ELEMENTS ARE ABORBED, THEN
*              XFLNOD=0 FOR THESE ELEMENTS.
*     FLNOD  : INTEGER(MAXNEL)
*              STORES FOR EACH ELEMENT A LIST OF NODES CONNECTED TO THE
*              ELEMENT. EACH LIST OF NODES IS HEADED BY ITS LENGTH.
*     XFNODL : INTEGER(NNODS)
*              STORES THE POINTERS INTO THE LISTS OF ELEMENTS WHICH ARE
*              CONNECTED TO THE NODE. IF NODES ARE ELIMINATED, THEN
*              XFNODL=0 FOR THESE NODES.
*     FNODL  : INTEGER(NNODS+NELNOD)
*              STORES FOR EACH NODE A LIST OF ELEMENTS CONNECTED TO THE
*              NODE. EACH LIST OF ELEMENTS IS HEADED BY ITS LENGTH.
*     DEGREE : INTEGER(NEQ)
*              HEAD OF DEGREE LIST.
*     QSIZE  : INTEGER(NNODS)
*              SUPERNODE SIZES.
*     VWGHT  : INTEGER(NNODS)
*              The weights of each vertex.
*     ELMLNK : INTEGER(NELS)
*              GENERALIZED ELEMENTS WHICH NEED DEGREE UPDATES.
*     FLAG   : INTEGER(NNODS)
*              FLAGS ACTIVE NODES IN THE DEGREE UPDATE STEP.
*     STOREL : INTEGER(NNODS)
*              LIST OF NODES IN NEW GENERATED ELEMENTS AND LISTS
*              OF NODES OF DIFFERENT TYPES IN DEGREE UPDATE STEP.
*     MARKER : INTEGER(NNODS)
*              MARK OF NODES.
*
*     SUBPROGRAMS/FUNCTIONS
*     ---------------------
*     SPRMD1
*     SPRMD2
*     SPRMD3
*     SPRMD4
*
*     ------------------------------------------------------------------
*
*     DATA STRUCTURE FOR DEGREE LISTS
*
*     DEGREE, INVP AND PERM  --   DOUBLY LINKED LIST OF NODES BY DEGREE
*                                 WITH OVERLOADED MEANINGS AS BELOW
*
*         INVP   -- >  0 --        NEXT NODE IN LIST
*                   =  0 --        Node is separator and should be
*                                  deferred until all other nodes are
*                                  ordered.
*                   <  0 --       -INVERSE PERMUTATION POSITION OF NODE
*                                    (MEANS NODE HAS BEEN ELIMINATED)
*
*         PERM   -- =  MAXMRK --   Node is separator and should be
*                                  deferred until all other nodes are
*                                  ordered.
*                   >  0 --        PREVIOUS NODE
*                   =  0 --        NODE NEEDS DEGREE UPDATING
*                   <  0 --       -DEGREE OF NODE
*                                    (POINTER TO HEAD OF LIST)
*                   = -MAXINT --   NODE DOES NOT REQUIRE DEGREE UPDATE
*                                  AND IS OUTMATCHED BY ANOTHER NODE
*
*         MARK   --                USED TO NOTE STAGE AT WHICH NODE HAS
*                                  BEEN MARKED
*                                    (MULTIPLE VALUES USED TO REDUCE
*                                     AMOUNT OF RESETTING REQUIRED)
*
*         MARKER -- =  MAXINT --   NODE HAS BEEN COMBINED INTO A CLIQUE
*                                  EITHER AN ELIMINATED SUPERELEMENT OR
*                                  A CLIQUE OF INDISTINGUISHABLE NODES
*                                  THESE CAUSE FURTHER OVERLOADINGS OF
*                                  INVP AND PERM.
*
*     ------------------------------------------------------------------

      INTEGER
     $ CMINDG, CMPGRP, CURMDN, I     , LIMIT , MARK  , NODE  , NOWORD,
     $ NUM   , NUMELM, NUMORD, NXTMDN, NXTUSE, xwrong

      EXTERNAL
     $ SPRMD1, SPRMD2, SPRMD3, SPRMD4

*     ==================================================================

*     SET NUM TO # ELIMINATED NODES+1.
*     NXTUSE IS SET TO NELS+NELNOD+1.
      CMPGRP = NNODS
      NCOMPR = 0
      NUM    = 1
      NXTUSE = NELS + NELNOD + 1

* --- The number of nodes to order in the first pass.
      NOWORD = NNODS - SEPSZE

* --- Abort the ordering if there are no internal nodes.
      IF ( NOWORD.LE.0 ) RETURN

* --- INITIALIZE FOR THE MINIMUM DEGREE ELIMINATION.
      CALL SPRMD1
     $( NELS,NEQ,NNODS,NELNOD,MAXNEL,MAXMRK,LWEIGH,XELNOD,ELNOD,XNODEL,
     $  NODEL,XNODES,XFLNOD,FLNOD,XFNODL,FNODL,DEGREE,INVP,PERM,QSIZE,
     $  VWGHT,ELMLNK,MARKER,FLAG,CMPGRP,SEPARR )

*     ==================================================================

      IF ( DEGREE(1).GT.0 ) THEN

* ------ ELIMINATE THE ISOLATED NODAL BLOCKS.
         NXTMDN = DEGREE(1)
 100     IF ( NXTMDN.GT.0 ) THEN
            CURMDN = NXTMDN
            NXTMDN = INVP(CURMDN)
            MARKER(CURMDN) = MAXMRK
            INVP(CURMDN) = -NUM
            NUM = NUM + QSIZE(CURMDN)
            GOTO 100
         ENDIF

* ------ GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
         IF ( NUM.GT.NOWORD ) GOTO 600

      ENDIF

* --- INITIALIZE FOR THE NON-TRIVIAL ELIMINATIONS.
      MARK   = 1
      CMINDG = 2
      DEGREE(1) = 0

* --- START THE MAIN ELIMINATION LOOP.

* --- FIND THE NEXT MINIMUM DEGREE.
 200  IF ( DEGREE(CMINDG).LE.0 ) THEN
            CMINDG = CMINDG + 1
            GOTO 200
      ENDIF

      LIMIT  = CMINDG + DELTA
      NUMELM = 0

* --- START THE MULTIPLE ELIMINATION STEP.

* --- RETRIEVE THE CURRENT MINIMUM DEGREE NODE.
 300  CURMDN = DEGREE(CMINDG)

      IF ( CURMDN.LE.0 ) THEN
         CMINDG = CMINDG + 1
         IF ( CMINDG.GT.LIMIT ) GOTO 500
         GOTO 300
      ENDIF

* --- REMOVE THE MINIMUM DEGREE NODE FROM THE DOUBLY LINKED DEGREE
*     STRUCTURE AND SET ITS INVERSE PERMUTATION POSITION IN INVP.
      NXTMDN = INVP(CURMDN)
      DEGREE(CMINDG) = NXTMDN
      IF ( NXTMDN.GT.0 ) PERM(NXTMDN) = -CMINDG
      INVP(CURMDN) = -NUM

* --- GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
      IF ( NOWORD.EQ.NNODS ) THEN
         IF ( NUM + QSIZE(CURMDN).GT.NOWORD ) GOTO 600
      ENDIF

* --- UPDATE THE NODE MARKER.
      MARK = MARK + 1

      IF ( MARK.GE.MAXMRK ) THEN

* ------ RESET THE MARKER VECTOR.
         MARK = 1
         DO 400 I = 1, NNODS
            IF ( MARKER(I).LT.MAXMRK ) MARKER(I) = 0
 400     CONTINUE
      ENDIF

* --- UPDATE THE GRAPH REPRESENTATION DUE TO THE ELIM. OF CURMDN.
      CALL SPRMD2
     $( CURMDN,NELS,NEQ,NNODS,NELNOD,MAXNEL,MARK,MAXMRK,NUMELM,NXTUSE,
     $  NCOMPR,XFNODL,FNODL,XFLNOD,FLNOD,DEGREE,PERM,INVP,QSIZE,VWGHT,
     $  ELMLNK,MARKER,STOREL )
      NUM = NUM + QSIZE(CURMDN)

* --- GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
      IF ( NUM.GT.NOWORD ) GOTO 600

* --- GOTO A NEW NODE ELIMINATION IN THIS MULTIPLE ELIMINATION STEP.
      IF ( DELTA.GE.0 ) GOTO 300

* --- GOTO NUMBERING IF ALL THE NODES ARE ELIMINATED.
 500  IF ( NUM.GT.NOWORD ) GOTO 600

* --- UPDATE THE DEGREE OF ALL THE NODES AFFECTED IN THE CURRENT
*     MULTIPLE ELIMINATION STEP, EXCEPT THOSE THAT ARE OUTMATCHED.
      CALL SPRMD3
     $( NUMELM,CMINDG,MARK,NELS,NEQ,NNODS,NELNOD,MAXNEL,MAXMRK,XFNODL,
     $  FNODL,XFLNOD,FLNOD,DEGREE,INVP,PERM,QSIZE,VWGHT,ELMLNK,MARKER,
     $  FLAG,STOREL )

* --- GOTO THE NEXT MULTIPLE ELIMINATION STEP.
      GOTO 200

* --- END OF MAIN ELIMINATION LOOP.
 600  CONTINUE
      IF ( NOWORD.LT.NNODS ) THEN

* ------ First pass on the subgraphs is finished.

         IF ( SEPNUM.EQ.1 ) THEN

* --------- Do not merge the separator nodes prior to post-pass.
            NUMORD = 0
            DO I = 1, NNODS
               IF ( SEPARR(I).GT.0 ) THEN
                  NUMORD = NUMORD + 1
                  FLAG(NUMORD) = I
               ENDIF
            END DO

         ELSE

* --------- Merge the separator nodes together as a preparation
*           for Post-Pass or exit with relative initial order.
            NUMORD = SEPCNT
            DO I = 1, NUMORD
               FLAG(I) = 0
            END DO

            xwrong = 0

            DO I = 1, NNODS
               IF ( SEPARR(I).GT.0 ) THEN

                  IF ( QSIZE(I).LE.0 ) THEN
                     xwrong = xwrong + 1
                     write ( lpu, '(a,i10)' ) 'Error in node ',i
                  ENDIF

* --------------- Node I is in a separator
                  IF ( FLAG(SEPARR(I)).EQ.0 ) THEN

* ------------------ Node I is the first node in this separator.
*                    Set flag(this separator) = I in order to merge all
*                    other nodes in this separator woth node I.
                     FLAG(SEPARR(I)) = I

                  ELSE

* ------------------ Node I belongs to a separator that has been
*                    initialised. Merge node I with the first node
*                    in this separator.
                     NODE = FLAG(SEPARR(I))
                     QSIZE(NODE) = QSIZE(NODE) + QSIZE(I)
                     QSIZE(I) = 0
                     PERM(I) = -MAXMRK
                     MARKER(I) = MAXMRK
                     INVP(I) = -NODE
                     XFNODL(I) = 0

                  ENDIF

               ENDIF
            END DO

         ENDIF

         IF ( SEPNUM.GT.0 ) THEN

* --------- Prepare for a Post-Pass over the separators.
            CALL SPRMD6
     $      ( NUMORD,CMINDG,MARK,NELS,NEQ,NNODS,NELNOD,MAXNEL,MAXMRK,
     $        XFNODL,FNODL,XFLNOD,FLNOD,DEGREE,INVP,PERM,VWGHT,
     $        MARKER,FLAG )

* --------- Goto to Post-Pass
            NOWORD = NNODS
            GOTO 200

         ENDIF

* ------ Order the nodes according to initial order.
         DO I = NUMORD, 1, -1
            NODE = FLAG(I)
            INVP(NODE) = -NUM
            NUM = NUM + QSIZE(NODE)
         END DO

      ENDIF

* --- DETERMINE THE FINAL NODAL ORDERING.
      CALL SPRMD4 ( NNODS, PERM, INVP, QSIZE )

*     ==================================================================

      RETURN

*     ------------------------------------------------------------------
*     END OF MODULE SPRCMD
      END


      SUBROUTINE   SPRMD1
     $( NELS  , NEQ   , NNODS , NELNOD, MAXNEL, MAXMRK, LWEIGH, XELNOD,
     $  ELNOD , XNODEL, NODEL , XNODES, XFLNOD, FLNOD , XFNODL, FNODL ,
     $  DEGREE, INVP  , PERM  , QSIZE , VWGHT , ELMLNK, MARKER, FLAG  ,
     $  CMPGRP, SEPARR )

*     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER
     $  NELS  , NEQ   , NNODS , NELNOD, MAXNEL, MAXMRK, CMPGRP
      LOGICAL
     $  LWEIGH
      INTEGER
     $  XELNOD(NELS+1), ELNOD(NELNOD), XNODEL(NNODS+1), NODEL(NELNOD),
     $  XNODES(NNODS+1), XFLNOD(NELS), FLNOD(MAXNEL),XFNODL(NNODS),
     $  FNODL(NNODS+NELNOD), DEGREE(NEQ), INVP(NNODS),
     $  PERM(NNODS), QSIZE(NNODS), VWGHT(NNODS), ELMLNK(NELS),
     $  MARKER(NNODS), FLAG(NNODS), SEPARR(NNODS)

*     ------------------------------------------------------------------
*
*     PURPOSE
*     -------
*
*     SPRMD1
*
* --- TO INITIALIZE FOR SPRCMD.
*
*     CREATED   : May  01, 1997 (ACD)
*     REVISIONS : MNT. XX, 199X (ACD)
*
*     ON ENTRY
*     --------
*
*     NELS   : INTEGER
*              NUMBER OF ELEMENTS.
*     NNODS  : INTEGER
*              NUMBER OF NODES.
*     NELNOD : INTEGER
*              SIZE OF CONNECTIVITY ARRAYS.
*     MAXNEL : INTEGER
*              WORKSPACE FOR FLNOD. AT LEAST NELS+NELNOD.
*     LWEIGH : Logical
*              Option for wieghted or standard minimum degree selection.
*              = TRUE : Weighted selection.
*              = FALSE: Standard selection.
*     (XELNOD,
*     ELNOD) : INTEGER(NELS), INTEGER(NELNOD)
*              ELEMENT-NODE CONNECTIVITY.
*     (XNODEL,
*     NODEL) : INTEGER(NNODS+1), INTEGER(NELNOD)
*              NODE-ELEMENT CONNECTIVITY.
*
*     ON EXIT
*     -------
*
*     XFLNOD : INTEGER(NELS)
*              POINTER INTO FLNOD FOR LISTS OF NODES CONNECTED TO THE
*              ELEMENTS.
*     FLNOD  : INTEGER(NELNOD+NELS)
*              LISTS OF NODES, WHERE EACH LIST IS HEADED BY ITS LENGTH.
*     XFNODL : INTEGER(NNODS)
*              POINTER INTO FNODL FOR LISTS OF ELEMENTS CONNECTED TO
*              THE NODES.
*     FNODL  : INTEGER(NELNOD+NNODS)
*              LISTS OF ELEMENTS, WHERE EACH LIST IS HEADED BY ITS
*              LENGTH.
*     DEGREE : Integer(NEQ)
*              Head of degree list.
*     (INVP,
*     PERM)  : 2*INTEGER(NNODS)
*              DOUBLY LINKED DEGREE STRUCTURE WITH INITIAL DEGREE OF
*              THE NODES.
*     QSIZE  : INTEGER(NNODS)
*              SUPERNODE SIZES INITIALIZED TO 1.
*     VWGHT  : INTEGER(NNODS)
*              Vertex weights.
*     ELMLNK : INTEGER(NELS)
*              LIST OF ELEMENTS IS NOW THE ZERO LIST.
*     MARKER : INTEGER(NNODS)
*              INITIALIZED TO ZERO.
*     FLAG   : INTEGER(NNODS)
*              FLAGS ACTIVE NODES IN THE DEGREE UPDATE STEP.
*     CMPGRP : Integer
*              The number of nodes after graph compression.
*
*     SUBPROGRAMS/FUNCTIONS
*     ---------------------
*
*     NONE
*
*     ------------------------------------------------------------------

      INTEGER           DEG   , ELEMNT, I     , IP    , ISTRT , ISTOP ,
     $                  J     , JSTRT , JSTOP , NEWNOD, NODE

      INTRINSIC         ABS

*     ==================================================================

* --- INITIALIZE.
      CMPGRP = 0
      DO 10 I = 1, NNODS
         FLAG(I) = 0
         MARKER(I) = 0
         QSIZE(I) = 1
         IF ( LWEIGH ) THEN
            VWGHT(I) = XNODES(I+1) - XNODES(I)
         ELSE
            VWGHT(I) = 1
         ENDIF
 10   CONTINUE
      DO 20 I = 1, NELS
         ELMLNK(I) = 0
 20   CONTINUE
      DO I = 1, NEQ
         DEGREE(I) = 0
      END DO

*     ==================================================================

* --- COMPUTE THE INITIAL DEGREE STRUCTURE.
      DO 300 NODE = 1, NNODS

         IF ( SEPARR(NODE).GT.0 ) THEN

* --------- The node is in a separator, mark it as so.
            PERM(NODE) = MAXMRK
            INVP(NODE) = 0

         ELSE

* --------- The node is interior, compute its degree.
            MARKER(NODE) = NODE
            ISTRT = XNODEL(NODE)
            ISTOP = XNODEL(NODE+1) - 1
            DEG = 1

            DO 200 I = ISTRT, ISTOP
               ELEMNT = NODEL(I)
               JSTRT = XELNOD(ELEMNT)
               JSTOP = XELNOD(ELEMNT+1) - 1
               DO 100 J = JSTRT, JSTOP
                  NEWNOD = ELNOD(J)
                  IF ( MARKER(NEWNOD).NE.NODE ) THEN
                     MARKER(NEWNOD) = NODE
                     DEG = DEG + VWGHT(NEWNOD)
                  ENDIF
 100           CONTINUE
 200        CONTINUE

            NEWNOD = DEGREE(DEG)
            INVP(NODE) = NEWNOD
            DEGREE(DEG)  = NODE
            IF ( NEWNOD.GT.0 ) PERM(NEWNOD) = NODE
            PERM(NODE) = -DEG

         ENDIF

 300  CONTINUE

* --- RESET MARKER TO ZERO.
      DO 400 I = 1, NNODS
         MARKER(I) = 0
 400  CONTINUE

*     ==================================================================

* --- BUILD THE GRAPH REPRESENTATION NEEDED BY THE ROUTINE.
      IP = NELS + NELNOD
      DO 600 I = NELS, 1, -1
         JSTRT = XELNOD(I)
         JSTOP = XELNOD(I+1)-1
         DO 500 J = JSTOP, JSTRT, -1
            FLNOD(IP) = ELNOD(J)
            IP = IP - 1
 500     CONTINUE
         XFLNOD(I) = IP
         FLNOD(IP) = JSTOP + 1 - JSTRT
         IP = IP - 1
 600  CONTINUE

      IP = NNODS + NELNOD
      DO 800 I = NNODS, 1, -1
         JSTRT = XNODEL(I)
         JSTOP = XNODEL(I+1)-1
         DO 700 J = JSTOP, JSTRT, -1
            FNODL(IP) = NODEL(J)
            IP = IP - 1
 700     CONTINUE
         XFNODL(I) = IP
         FNODL(IP) = JSTOP + 1 - JSTRT
         IP = IP - 1
 800  CONTINUE

*     ==================================================================

      RETURN

*     ------------------------------------------------------------------
*     END OF MODULE SPRMD1
      END


      SUBROUTINE   SPRMD2
     $( NODE  , NELS  , NEQ   , NNODS , NELNOD, MAXNEL, MARK  , MAXMRK,
     $  NUMELM, NXTUSE, NCOMPR, XFNODL, FNODL , XFLNOD, FLNOD , DEGREE,
     $  PERM  , INVP  , QSIZE , VWGHT , ELMLNK, MARKER, STOREL )

*     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER
     $  NODE  , NELS  , NEQ   , NNODS , NELNOD, MAXNEL, MARK  , MAXMRK,
     $  NUMELM, NXTUSE, NCOMPR

      INTEGER
     $  XFNODL(NNODS), FNODL(NNODS+NELNOD), XFLNOD(NELS),
     $  FLNOD(MAXNEL), DEGREE(NEQ)        , PERM(NNODS) ,
     $  INVP(NNODS)  , MARKER(NNODS)      , QSIZE(NNODS) ,
     $  VWGHT(NNODS) , ELMLNK(NELS)       , STOREL(NNODS)

*     ------------------------------------------------------------------
*
*     PURPOSE
*     -------
*
*     SPRMD2
*
* --- TO UPDATE THE REPRESENTATION OF THE GRAPH DUE TO THE ELIMINATION
*     OF NODE. THE ROUTINE WORKS ON THE IMPLIVIT REPRESENTATION OF THE
*     GRAPH. THE IMPLICIT GRAPH STRUCTURE IS REPRESENTED BY THE ARRAYS
*     (XFNODL,FNODL) $ (XFLNOD,FLNOD).
*
*     CREATED   : May  01, 1997 (ACD)
*     REVISIONS : MNT. XX, 199X (ACD)
*
*     ON ENTRY
*     --------
*
*     NODE   : INTEGER
*              NODE TO BE ELIMINATED.
*     NELS   : INTEGER
*              NUMBER OF FINITE ELEMENTS IN THE MESH.
*     NNODS  : INTEGER
*              NUMBER OF NODES IN THE MESH.
*     NELNOD : INTEGER
*              SIZE OF INDICES STORED IN FNODL AND FLNOD AT FIRST ENTRY.
*     MAXNEL : INTEGER
*              SIZE OF FLNOD, WHERE MAXNEL .GE. NELS+NELNOD.
*     MARK   : INTEGER
*              STAGE COUNTER, FACILITATES MARKING NODES.
*     MAXMRK : INTEGER
*              A LARGEST POSSIBLE INTEGER, HOWEVER,
*              MAXMRK .GT. MAX(NELS,NNODS) IS ENOUGH.
*
*     UPDATED PARAMETERS
*     ------------------
*
*     NUMELM : INTEGER
*              NUMBER OF ACTIVE NEW GENERATED ELEMENTS AT THIS STAGE.
*     NXTUSE : INTEGER
*              THE NEXT FREE LOCATION OF FLNOD.
*              ON THE FIRST ENTRY, IT MUST BE SET TO NELNOD+NELS+1.
*              UPDATED BY THE ROUTINE TO ALWAYS POINT TO THE FIRST
*              FREE LOCATION IN FLNOD.
*     NCOMPR : INTEGER
*              UPDATED WITH THE NUMBER OF WORKSPACE COMPRESSIONS NEEDED
*              IN THIS GRAPH UPDATE.
*     XFNODL : INTEGER(NNODS)
*              POINTER INTO THE GENERATED ELEMENTS FOR EACH NODE.
*              ABSORBED AND ELIMINATED NODES HAVE XFNODL(I)=0.
*              ON EXIT FOR AN ELIMINATED NODE I, XFNODL(I)=0.
*     FNODL  : INTEGER(NELNOD+NNODS)
*              THE LIST OF GENERATED ELEMENTS FOR EACH NODE HEADED BY
*              THE LENGTH OF THE LIST.
*              ON EXIT REFERENCE TO ABSORBED ELEMENTS IS REMOVED.
*     XFLNOD : INTEGER(NELS)
*              POINTER INTO THE LIST OF NODES FOR EACH ELEMENT.
*              IF AN ELEMENT I IS ABSORBED, THEN XFLNOD(I)=0.
*     FLNOD  : INTEGER(MAXNEL)
*              THE NODES IN A GENERATED ELEMENT HEADED BY THE LENGTH
*              OF THE LIST.
*     DEGREE : INTEGER(NEQ)
*              HEAD OF DEGREE LIST.
*     PERM   : INTEGER(NNODS)
*              PREVIOUS NODE.
*     INVP   : INTEGER(NNODS)
*              NEXT NODE.
*     QSIZE  : INTEGER(NNODS)
*              SIZE OF THE SUPERNODES.
*     VWGHT  : INTEGER(NNODS)
*              Vertex weights.
*     ELMLNK : INTEGER(NELS)
*              LIST OF GENERATED ELEMENTS WHERE SOME NODES NEED
*              DEGREE UPDATE.
*     MARKER : INTEGER(NNODS)
*              MARKER VECTOR FOR EFFICIENT MASK OF TOUCHED NODES.
*
*     ON EXIT
*     -------
*
*     NO PURE EXIT PARAMETERS.
*
*     WORKING ARRAYS
*     --------------
*
*     STOREL : INTEGER(NNODS)
*              TEMPORARY STORAGE OF THE NODES BELONGING TO THE NEW
*              GENERATED ELEMENT.
*
*     SUBPROGRAMS
*     -----------
*
*     SPRMD5
*
*     FUNCTIONS
*     ---------
*
*     NONE
*
*     INTRINSIC
*     ---------
*
*     NONE
*
*     INCLUDE BLOCKS
*     --------------
*
*     NONE
*
*     COMMON BLOCKS
*     -------------
*
*     NONE
*
*     ------------------------------------------------------------------

      INTEGER           COUNT , ELEMNT, ELMSZE, FSTELM, GEFREE, I     ,
     $                  IP    , IPFSTL, ISTRT , ISTOP , J     , JSTRT ,
     $                  JSTOP , MFLNOD, NEXT  , PREV  , RCHLOC, RCHNOD,
     $                  RCHSZE

      EXTERNAL          SPRMD5

*     ------------------------------------------------------------------
*     COUNT COUNTS THE NUMBER OF GENERATED ELEMENTS CONNECTED TO A NODE.
*     ELEMNT ELEMENT INDEX.
*     FSTELM THE FIRST ELEMENT IN THE LIST OF GENERATED ELEMENTS FOR
*            NODE. USED AS REPRESENTATIVE FOR THE NEW ELEMENT.
*     FSTIP FIRST POSITION IN FLNOD FOR AN ELEMENT WHEN LIST IS COM-
*           PRESSED.
*     GEFREE STORES THE SIZE OF FREE SPACE AT THE END OF FLNOD.
*     IP USED AS ARRAY POINTER.
*     I, J ARE LOOP COUNTERS
*     ISTRT, ISTOP ARE ARRAY/LIST BOUNDS FOR LOOP INDEX I.
*     JSTRT, JSTOP ARE ARRAY/LIST BOUNDS FOR LOOP INDEX J.
*     MFLNOD IN CALLS TO SPRMD5, LAST LOCATION OF FLNOD IN USE.
*     NEXT NEXT NODE IN DEGREE LIST.
*     PREV PREVIOUS NODE IN THE DEGREE LIST.
*     RCHLOC POINTER INTO THE GENERATED ELEMENT.
*     RCHNOD NODE IN THE NEW GENERATED ELEMENT.
*     RCHSZE SIZE OF THE NEW GENERATED ELEMENT.
*     ------------------------------------------------------------------

* --- FETCH THE GENERATED ELEMENT AND STORE IT IN STOREL.
      ISTRT = XFNODL(NODE)
      ISTOP = ISTRT + FNODL(ISTRT)
      RCHSZE = 0
      ELMSZE = 0
      FSTELM = 0
      IPFSTL = 0

* --- MARK NODE AS ELIMINATED.
      MARKER(NODE) = MAXMRK
      PERM(NODE) = -MAXMRK
      XFNODL(NODE) = 0

      DO 200 I = ISTRT+1, ISTOP
         ELEMNT = FNODL(I)
         JSTRT = XFLNOD(ELEMNT)
         JSTOP = JSTRT + FLNOD(JSTRT)
         IF ( FLNOD(JSTRT).GT.ELMSZE ) THEN
            ELMSZE = FLNOD(JSTRT)
            FSTELM = ELEMNT
            IPFSTL = JSTRT
         ENDIF

         XFLNOD(ELEMNT) = 0

         DO 100 J = JSTRT+1, JSTOP

            RCHNOD = FLNOD(J)

            IF ( MARKER(RCHNOD).LT.MARK ) THEN

               MARKER(RCHNOD) = MARK

* ------------ THE NODE HAS NOT YET BEEN MERGED.
               PREV = PERM(RCHNOD)

               IF ( PREV.EQ.MAXMRK ) THEN

* --------------- Node is in separator.
                  IF ( XFNODL(RCHNOD).GT.0 ) THEN

* ------------------ NODE IS STILL ACTIVE.
                     RCHSZE = RCHSZE + 1
                     STOREL(RCHSZE) = RCHNOD
                  ENDIF

               ELSE

* --------------- Node is interior.
                  IF ( ( PREV.NE.0 ).AND.( PREV.NE.-MAXMRK ) ) THEN

* ------------------ THE NODE IS STILL IN THE DEGREE STRUCTURE,
*                    REMOVE IT.
                     NEXT = INVP(RCHNOD)
                     IF ( NEXT.GT.0 ) PERM(NEXT) = PREV
                     IF ( PREV.GT.0 ) INVP(PREV) = NEXT
                     IF ( PREV.LT.0 ) DEGREE(-PREV) = NEXT

* ------------------ AND MARK IT FOR DEGREE UPDATE.
                     PERM(RCHNOD) = 0
                  ENDIF

                  IF ( XFNODL(RCHNOD).GT.0 ) THEN

* ------------------ NODE IS STILL ACTIVE.
                     RCHSZE = RCHSZE + 1
                     STOREL(RCHSZE) = RCHNOD
                     PERM(RCHNOD) = 0
                  ENDIF

               ENDIF

            ENDIF

 100     CONTINUE
 200  CONTINUE

*     ----------------------------------------------------------
*     FSTELM IS THE REPRESENTATIVE OF THE NEW GENERATED ELEMENT.
*     STOREL NOW CONTAINS THE NODES OF THE GENERATED ELEMENT
*            WITHOUT DUPLICATE REFERENCES.
*     RCHSZE NOW STORES THE SIZE OF THE NEW GENERATED ELEMENT.
*     XFLNOD IS ZEROED FOR EACH ABSORBED ELEMENT.
*     ----------------------------------------------------------

* --- UPDATE NUMELM AND ELMLNK.
      NUMELM = NUMELM + 1
      ELMLNK(NUMELM) = FSTELM

*     ==================================================================

      IF ( RCHSZE.LE.ELMSZE ) THEN

* ------ STORE THE NEW ELEMENT IN THE LOCATIONS USED FOR FSTELM.
         XFLNOD(FSTELM) = IPFSTL
         FLNOD(IPFSTL) = RCHSZE
         ISTRT = IPFSTL
         ISTOP = ISTRT + RCHSZE
         IP = 0
         DO 300 I = ISTRT+1, ISTOP
            IP = IP + 1
            FLNOD(I) = STOREL(IP)
 300     CONTINUE
      ELSE

* ------ STORE THE ELEMENT AT THE END OF FLNOD.
         GEFREE = MAXNEL - NXTUSE
         IF ( GEFREE.LT.RCHSZE ) THEN

* --------- COMPRESS FLNOD.
            MFLNOD = NXTUSE-1
            CALL SPRMD5 ( NELS, XFLNOD, FLNOD, MFLNOD, NXTUSE, NCOMPR )
         ENDIF

* ------ STORE THE ELEMENT.
         XFLNOD(FSTELM) = NXTUSE
         FLNOD(NXTUSE) = RCHSZE
         ISTRT = NXTUSE
         ISTOP = ISTRT + RCHSZE
         IP = 0
         DO 400 I = ISTRT+1, ISTOP
            IP = IP + 1
            FLNOD(I) = STOREL(IP)
 400     CONTINUE
         NXTUSE = ISTOP + 1

      ENDIF

*     --------------------------------------------------------------
*     NXTUSE NOW POINTS TO THE NEXT FREE LOCATION AFTER NEW ELEMENT.
*     XFLNOD HAS GOT THE POINTER.
*     FLNOD STORES THE NODES.
*     --------------------------------------------------------------

* --- GO THROUGH THE NEW ELEMENT, MARK THE NODES FOR DEGREE UPDATE,
*     AND LOOK FOR NODES THAT ARE INTERNAL TO THE NEW ELEMENT.
      ISTRT = XFLNOD(FSTELM)
      ISTOP = ISTRT + FLNOD(ISTRT)
      RCHLOC = ISTRT
      RCHSZE = 0

      DO 600 I = ISTRT+1, ISTOP

         RCHNOD = FLNOD(I)

* ------ THE NODE IS NOT ELIMINATED.
         JSTRT = XFNODL(RCHNOD)
         JSTOP = JSTRT + FNODL(JSTRT)
         COUNT = 0
         IP    = JSTRT

         DO 500 J = JSTRT+1, JSTOP

            ELEMNT = FNODL(J)

            IF ( ( XFLNOD(ELEMNT).GT.0 ).AND.( ELEMNT.NE.FSTELM ) ) THEN
               COUNT = COUNT + 1
               IP = IP + 1
               FNODL(IP) = ELEMNT
            ENDIF

 500     CONTINUE

* ------ THE NODE IS BY DEFAULT ATTACHED TO FSTELM.
         COUNT = COUNT + 1
         IP = IP + 1
         FNODL(IP) = FSTELM

         IF ( COUNT.EQ.1 ) THEN

            IF ( PERM(RCHNOD).EQ.MAXMRK ) THEN

* ------------ RCHNOD is a separator node and cannot be merged
*              with node even though it is interior.
               RCHLOC = RCHLOC + 1
               RCHSZE = RCHSZE + 1
               FLNOD(RCHLOC) = RCHNOD
               FNODL(JSTRT) = COUNT

            ELSE

* ------------ THE NODE RCHNOD IS INTERIOR TO THE ELEMENT
*              AND MAY BE ELIMINATED TOGETHER WITH NODE.
               QSIZE(NODE) = QSIZE(NODE) + QSIZE(RCHNOD)
               QSIZE(RCHNOD) = 0

               VWGHT(NODE) = VWGHT(NODE) + VWGHT(RCHNOD)

               INVP(RCHNOD) = -NODE
               PERM(RCHNOD) = -MAXMRK
               MARKER(RCHNOD) = MAXMRK
               XFNODL(RCHNOD) = 0

            ENDIF

         ELSE

* --------- THE NODE IS AT LEAST ON TWO ELEMENTS.
            RCHLOC = RCHLOC + 1
            RCHSZE = RCHSZE + 1
            FLNOD(RCHLOC) = RCHNOD
            FNODL(JSTRT) = COUNT

         ENDIF

 600  CONTINUE

* --- UPDATE THE SIZE OF THE REACHABLE SET.
      FLNOD(ISTRT) = RCHSZE

*     ==================================================================

      RETURN

*     ------------------------------------------------------------------
*     END OF MODULE SPRMD2
      END


      SUBROUTINE   SPRMD3
     $( NUMELM, CMINDG, MARK  , NELS  , NEQ   , NNODS , NELNOD, MAXNEL,
     $  MAXMRK, XFNODL, FNODL , XFLNOD, FLNOD , DEGREE, INVP  , PERM  ,
     $  QSIZE , VWGHT , ELMLNK, MARKER, FLAG  , STOREL  )

*     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER
     $  NUMELM, CMINDG, MARK  , NELS  , NEQ   , NNODS , NELNOD, MAXNEL,
     $  MAXMRK

      INTEGER
     $  XFNODL(NNODS),FNODL(NNODS+NELNOD),XFLNOD(NELS),FLNOD(MAXNEL),
     $  DEGREE(NEQ),INVP(NNODS),PERM(NNODS),QSIZE(NNODS),VWGHT(NNODS),
     $  ELMLNK(NELS),MARKER(NNODS),FLAG(NNODS),STOREL(NNODS)

*     ------------------------------------------------------------------

*     PURPOSE
*     -------

*     SPRMD3

* --- TO PERFORM DEGREE UPDATE. IT ALSO SEARCH FOR INDISTINGUISHABLE AND
*     OUTMATCHED NODES.

*     CREATED   : May  01, 1997 (ACD)
*     REVISIONS : MNT. XX, 199X (ACD)

*     ON ENTRY
*     --------

*     NUMELM : INTEGER
*              HEAD TO THE LIST OF GENERALIZED ELEMENTS OF WHICH THE
*              NODES NEED DEGREE UPDATE.
*     CMINDG : INTEGER
*              CURRENT MINIMUM DEGREE.
*     MARK   : INTEGER
*              USED TO PREVENT A NODE TO BE TREATED MORE THAN ONCE IN
*              EACH STEP OF THE ROUTINE.
*     NELS   : INTEGER
*              NUMBER OF ELEMENTS.
*     NNODS  : INTEGER
*              NUMBER OF NODES.
*     NELNOD : INTEGER
*              SIZE OF ELIMINATION GRAPH REPRESENTATION.
*     MAXNEL : INTEGER
*              SIZE OF FLNOD.
*     MAXMRK : INTEGER
*              USED AS UPPER LIMIT FOR MARKING NODES, ANY NUMBER LARGER
*              THAN NNODS IS OK.
*     (XFNODL,
*     FNODL,
*     XFLNOD,
*     FLNOD) : INTEGER(NNODS), INTEGER(NNODS+NELNOD), INTEGER(NELS),
*              INTEGER(MAXNEL)
*              REPRESENTATION OF THE ELIMINATION GRAPH.
*     QSIZE  : INTEGER(NNODS)
*              SIZE OF THE SUPERNODES.
*     VWGHT  : INTEGER(NNODS)
*              Vertex weights.
*     ELMLNK : INTEGER(NELS)
*              LIST OF GENERALIZED ELEMENTS OF WHICH THE NODES NEED
*              DEGREE UPDATE.
*     MARKER : INTEGER(NNODS)
*              FACILITATES MARKING NODES.
*     FLAG   : INTEGER(NNODS)
*              MARKS ACTIVE NODES IN THE DEGREE UPDATE STEP. THE ROUTINE
*              ASSUMES THAT THE ARRAY IS ZERO ON ENTRY. THE ARRAY IS
*              RESET TO ZERO ON EXIT.

*     ON EXIT
*     -------

*     DEGREE : INTEGER(NEQ)
*              Head of degree list.
*     (INVP,
*     PERM)  : 2*INTEGER(NNODS)
*              UPDATED DEGREE STRUCTURE.
*     XFNODL : INTEGER(NNODS)
*              FOR AN ABSORBED NODE I, XFNODL(I)=0.
*     QSIZE  : INTEGER(NNODS)
*              UPDATED SIZE OF THE SUPERNODES.

*     WORKING ARRAYS
*     --------------

*     STOREL : INTEGER(NNODS)
*              LIST OF NODES WHICH NEED DEGREE UPDATE.

*     SUBPROGRAMS/FUNCTIONS
*     ---------------------

*     NONE

*     ------------------------------------------------------------------

      INTEGER           CLQSZE, DEG0  , DEG1  , ELM   , I     , IEL   ,
     $                  IP    , ISTOP , ISTRT , J     , JSTOP , JSTRT ,
     $                  K     , KSTOP , KSTRT , NELM  , NEWELM, NEWNOD,
     $                  NEXT  , NODE  , PREV  , QLX   , STORED

*     ==================================================================

      DO 800 IEL = 1, NUMELM

* ------ START THE DEGREE UPDATE FOR THE NODES ON ELM.
         ELM = ELMLNK(IEL)

* ------ DEG0 is the 'degree' of the new generated element. To get the
*        externel degree count of a node on this element it must be
*        subtracted.
         DEG0 = 0
         ISTRT = XFLNOD(ELM)
         ISTOP = ISTRT + FLNOD(ISTRT)
         NELM = 0
         IP = ISTRT

* ------ FIRST INITIALIZE FLAG AND COMPRESS FLNOD.
         DO 100 I = ISTRT+1, ISTOP
            NODE = FLNOD(I)

            IF ( XFNODL(NODE) .GT. 0 ) THEN
               DEG0 = DEG0 + VWGHT(NODE)
               NELM = NELM + 1
               IP = IP + 1
               FLNOD(IP) = NODE
               FLAG(NODE) = 1
            ENDIF

 100     CONTINUE

* ------ RUN THROUGH THE ELEMENT, SEARCH FOR INDISTINGUISHABLE AND OUTM.
*        NODES, AND FINALLY UPDATE THE DEGREE.
         FLNOD(ISTRT) = NELM
         ISTOP = ISTRT + NELM
         DO 600 I = ISTRT+1, ISTOP

* --------- Get the node.
            NODE = FLNOD(I)

            IF ( PERM(NODE).EQ.0 ) THEN

* ------------ THIS NODE NEEDS DEGREE UPDATE.
               DEG1 = DEG0
               MARK = MARK + 1
               IF ( MARK.GE.MAXMRK ) THEN

* --------------- RESET THE MARKER VECTOR.
                  MARK = 1
                  DO 200 J = 1, NNODS
                     IF ( MARKER(J) .LT. MAXMRK ) MARKER(J) = 0
 200              CONTINUE
               ENDIF

* ------------ Mark it.
               MARKER(NODE) = MARK

               JSTRT = XFNODL(NODE)
               JSTOP = JSTRT + FNODL(JSTRT)
               STORED = 0

               DO 400 J = JSTRT+1, JSTOP

* --------------- RUN THROUGH THE REACHABLE SET OF NODE.
                  NEWELM = FNODL(J)

                  IF ( NEWELM .NE. ELM ) THEN

                     KSTRT = XFLNOD(NEWELM)

                     IF ( KSTRT .GT. 0 ) THEN

                        KSTOP = KSTRT + FLNOD(KSTRT)

                        DO 300 K = KSTRT + 1, KSTOP
                           NEWNOD = FLNOD(K)

                           IF ( XFNODL(NEWNOD) .GT. 0 .AND.
     $                          NEWNOD .NE. NODE )          THEN

                              IF ( FLAG(NEWNOD) .GE. 1 ) THEN

* ------------------------------ NEWNOD IS ON ELM. COUNT ITS APPEARANCES
*                                TO SEE IF ITS INDISTINGUISHABLE OR OUT-
*                                MATCHED BY NODE.
                                 IF ( FLAG(NEWNOD) .EQ. 1 ) THEN
                                    FLAG(NEWNOD) = 2
                                    STORED = STORED + 1
                                    STOREL(STORED) = NEWNOD
                                 ELSE
                                    FLAG(NEWNOD) = FLAG(NEWNOD) + 1
                                 ENDIF

                              ELSEIF ( MARKER(NEWNOD) .LT. MARK ) THEN

* ------------------------------ THIS NODE IS NOT ON ELM ACCUMULATE ITS
*                                SIZE IN THE DEGREE COUNT AND MARK IT.
                                 DEG1 = DEG1 + VWGHT(NEWNOD)
                                 MARKER(NEWNOD) = MARK

                              ENDIF

                           ENDIF

 300                    CONTINUE

                     ENDIF
                  ENDIF

 400           CONTINUE

               IF ( STORED .GT. 0 ) THEN

* --------------- QLX IS NOW SET TO THE SIZE OF THE CLIQUE MEMBERSHIP OF
*                 NODE.
                  QLX = FNODL(JSTRT)

                  DO 500 J = 1, STORED

                     NEWNOD = STOREL(J)

                     IF ( PERM(NEWNOD).EQ.MAXMRK ) THEN

                        FLAG(NEWNOD) = 1

                     ELSE

                        IF ( FLAG(NEWNOD) .EQ. QLX ) THEN

* ------------------------ THIS NODE HAS APPEARED QLX TIMES WHEN RUNNING
*                          THROUGH THE REACHABLE SET OF NODE. THIS MEANS
*                          THAT IT IS EITHER INDISTINGUISHABLE FROM NODE
*                          OR OUTMATCHED BY NODE.
                           PREV = PERM(NEWNOD)
                           IF ( ( PREV.NE.0 ).AND.( PREV.NE.-MAXMRK ) )
     $                     THEN

* --------------------------- IT IS IN THE DEGREE STRUCTURE, REMOVE IT.
                              NEXT = INVP(NEWNOD)
                              IF ( NEXT .GT. 0 ) PERM(NEXT) = PREV
                              IF ( PREV .GT. 0 ) INVP(PREV) = NEXT
                              IF ( PREV .LT. 0 ) DEGREE(-PREV) = NEXT
                           ENDIF

                           CLQSZE = FNODL(XFNODL(NEWNOD))

                           IF ( CLQSZE.EQ.QLX ) THEN

* --------------------------- THE NODE IS INDISTINGUISHABLE FROM NODE.
*                             MERGE IT WITH NODE.
                              QSIZE(NODE) = QSIZE(NODE) + QSIZE(NEWNOD)
                              QSIZE(NEWNOD) = 0

                              VWGHT(NODE) = VWGHT(NODE) + VWGHT(NEWNOD)

                              MARKER(NEWNOD) = MAXMRK
                              INVP(NEWNOD) = -NODE
                              XFNODL(NEWNOD) = 0
                              FLAG(NEWNOD) = 0
                           ELSE

* --------------------------- THE NODE IS OUTMATCHED BY NODE.
                              FLAG(NEWNOD) = 1
                           ENDIF

* ------------------------ FOR BOTH INDISTINGUISHABLE AND OUTMATCHED
*                          NODES MARK IT TO AVOID DEGREE UPDATE OF
*                          OUTMATCHED NODES.
                           PERM(NEWNOD) = -MAXMRK
                        ELSE

* ------------------------ THE NODE IS NEITHER INDISTINGUISHABLE NOR
*                          OUTMATCHED FROM/BY NODE.
                           FLAG(NEWNOD) = 1
                        ENDIF

                     ENDIF

 500              CONTINUE

               ENDIF

* ------------ DO THE FINAL DEGREE COUNT OF NODE, AND PUT IT INTO THE
*              DEGREE LIST.
               DEG1 = DEG1 - VWGHT(NODE) + 1
               NEXT = DEGREE(DEG1)
               INVP(NODE) = NEXT
               IF ( NEXT.GT.0 ) PERM(NEXT) = NODE
               DEGREE(DEG1) = NODE
               PERM(NODE) = -DEG1
               IF ( DEG1.LT.CMINDG ) CMINDG = DEG1

            ENDIF

 600     CONTINUE

* ------ FINALLY RESET FLAG TO ZERO, AND COMPRESS FLNOD.
         NELM = 0
         IP = ISTRT
         DO 700 I = ISTRT + 1, ISTOP
            NODE = FLNOD(I)
            IF ( XFNODL(NODE) .GT. 0 ) THEN
               NELM = NELM + 1
               IP = IP + 1
               FLNOD(IP) = NODE
               FLAG(NODE) = 0
            ENDIF
 700     CONTINUE
         FLNOD(ISTRT) = NELM

 800  CONTINUE

*     ==================================================================

      RETURN

*     ------------------------------------------------------------------
*     END OF MODULE SPRMD3
      END


      SUBROUTINE   SPRMD4   ( NNODS , PERM  , INVP  , QSIZE )

*     ==================================================================

      IMPLICIT NONE

      INTEGER           NNODS

      INTEGER           PERM  (NNODS)         , INVP  (NNODS)         ,
     $                  QSIZE (NNODS)

*     ==================================================================

*     PURPOSE
*     -------

*     SPRMD4 : TO PERFORM THE FINAL STEP IN PRODUCING THE PERMUTATION
*              AND INVERSE PERM VECTORS IN THE MULTIPLE MINIMUM DEGREE
*              ORDERING ALGORITHM.

*     CREATED   : May  01, 1997 (ACD)
*     REVISIONS : MNT. XX, 199X (ACD)

*     ON ENTRY
*     --------

*     NNODS  : INTEGER
*              NUMBER OF NODES IN THE GRAPH.
*     QSIZE  : INTEGER(NNODS)
*              SIZE OF SUPERNODES AT ELIMINATION.

*     ON EXIT
*     -------

*     INVP   : INTEGER(NNODS)
*              INVERSE PERM VECTOR.  ON INPUT, IF QSIZE(NODE)=0, THEN
*              NODE HAS BEEN MERGED INTO THE NODE -INVP(NODE);
*              OTHERWISE, -INVP(NODE) IS ITS INVERSE LABELLING.
*     PERM   : INTEGER(NNODS)
*              THE PERMUTATION VECTOR.

*     ==================================================================

      INTEGER           NODE  , FATHER, NUM   , ROOT  , NEXTF , NQSIZE

*     ==================================================================

      DO 100 NODE = 1, NNODS
         NQSIZE = QSIZE(NODE)
         IF ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
         IF ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
  100 CONTINUE

* --- FOR EACH NODE WHICH HAS BEEN MERGED, DO
      DO 500 NODE = 1, NNODS
         IF ( PERM(NODE) .GT. 0 )  GOTO 500

* ------ TRACE THE MERGED TREE UNTIL ONE WHICH HAS
*        NOT BEEN MERGED, CALL IT ROOT
         FATHER = NODE

  200    IF ( PERM(FATHER) .GT. 0 )  GOTO 300
            FATHER = - PERM(FATHER)
            GOTO 200

* ------ NUMBER NODE AFTER ROOT.
  300    ROOT       = FATHER
         NUM        = PERM(ROOT) + 1
         INVP(NODE) = - NUM
         PERM(ROOT) = NUM

* ------ SHORTEN THE MERGED TREE.
         FATHER = NODE

  400    NEXTF = - PERM(FATHER)
         IF ( NEXTF .LE. 0 )  GOTO 500
            PERM(FATHER) = - ROOT
            FATHER       = NEXTF
            GOTO 400

  500 CONTINUE

* --- READY TO COMPUTE PERM.
      DO 600 NODE = 1, NNODS
         NUM        = - INVP(NODE)
         INVP(NODE) = NUM
         PERM(NUM)  = NODE
  600 CONTINUE

      RETURN

*     ------------------------------------------------------------------
*     END OF MODULE SPRMD4
      END


      SUBROUTINE SPRMD5(N, IPE, IW, LW, IWFR, ICOMP)

      IMPLICIT NONE

* COMPRESS LISTS HELD BY SPRMD2 IN IW AND ADJUST POINTERS
*     IN IPE TO CORRESPOND.
*     CREATED : May  01, 1997 (ACD)

      INTEGER N, LW, IWFR, ICOMP
      INTEGER IPE(N)
      INTEGER IW(LW)

* N IS THE MATRIX ORDER. IT IS NOT ALTERED.
* IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
*     ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
* IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
*     LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
* LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
* IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
*     LOCATION IN IW.

      INTEGER I, IR, K, K1, K2, LWFR

      ICOMP = ICOMP + 1

* PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
*     LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
*     -(LIST NUMBER).

      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE

* COMPRESS
* IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
* LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.

      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70

* SEARCH FOR THE NEXT NEGATIVE ENTRY.
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70

* PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
*     AND PREPARE TO COPY LIST.
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50

* COPY LIST TO NEW POSITION.
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END


      SUBROUTINE SPRMD6
     $( SEPCNT, CMINDG, MARK  , NELS  , NEQ   , NNODS , NELNOD, MAXNEL,
     $  MAXMRK, XFNODL, FNODL , XFLNOD, FLNOD , DEGREE, INVP  , PERM  ,
     $  VWGHT , MARKER, FLAG )

*     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER
     $  SEPCNT, CMINDG, MARK  , NELS  , NEQ   , NNODS , NELNOD, MAXNEL,
     $  MAXMRK
      INTEGER
     $  XFNODL(NNODS),FNODL(NNODS+NELNOD),XFLNOD(NELS),FLNOD(MAXNEL),
     $  DEGREE(NEQ),INVP(NNODS),PERM(NNODS),VWGHT(NNODS),
     $  MARKER(NNODS),FLAG(NNODS)

*     ------------------------------------------------------------------
*
*     This routine prepares for a Post-Pass on the separator nodes.
*
*     ------------------------------------------------------------------

      INTEGER
     $  DEG   , ELEMNT, I     , INOD  , ISTRT , ISTOP , J     , JSTRT ,
     $  JSTOP , NEXT  , NODE  , RCHNOD

* --- Initialize
      CMINDG = NEQ + 2

* --- The degree structure must be reset for
*     the separator nodes.
      DO I = 1, SEPCNT
         NODE = FLAG(I)
         PERM(NODE) = 0
         INVP(NODE) = 0
      END DO
      DO I = 1, NEQ
         DEGREE(I) = 0
      END DO

* --- Update the degree of each separator.
      DO INOD = 1, SEPCNT

* ------ Get the separator representative.
         NODE = FLAG(INOD)
         JSTRT = XFNODL(NODE)

         IF ( JSTRT.GT.0 ) THEN

            DEG = 1

* --------- Increase marker value.
            MARK = MARK + 1
            IF ( MARK.GE.MAXMRK ) THEN

* ------------ RESET THE MARKER VECTOR.
               MARK = 1
               DO J = 1, NNODS
                  IF ( MARKER(J) .LT. MAXMRK ) MARKER(J) = 0
               END DO
            ENDIF

* --------- Mark NODE.
            MARKER(NODE) = MARK

* --------- Run through the reachable set of node.
            JSTOP = JSTRT + XFNODL(JSTRT)

            DO J = JSTRT+1, JSTOP

* ------------ Get the element in reachable set of NODE.
               ELEMNT = FNODL(J)

               ISTRT = XFLNOD(ELEMNT)
               IF ( ISTRT.GT.0 ) THEN
                  ISTOP = ISTRT + FLNOD(ISTRT)

                  DO I = ISTRT+1, ISTOP

* ------------------ Get node in reachable set.
                     RCHNOD = FLNOD(I)

                     IF ( MARKER(RCHNOD).LT.MARK ) THEN

* --------------------- This node has not yet been visited.
                        DEG = DEG + VWGHT(RCHNOD)
                        MARKER(RCHNOD) = MARK

                     ENDIF

                  END DO
               ENDIF

            END DO

* --------- Update the degree counter.
            CMINDG = MIN( CMINDG,DEG )

* --------- Insert NODE into the degree structure.
            NEXT = DEGREE(DEG)
            INVP(NODE) = NEXT
            IF ( NEXT.GT.0 ) PERM(NEXT) = NODE
            DEGREE(DEG) = NODE
            PERM(NODE) = -DEG

         ENDIF

      END DO

* --- Set FLAG back to zero.
      DO I = 1, SEPCNT
         FLAG(I) = 0
      END DO

      RETURN
      END
