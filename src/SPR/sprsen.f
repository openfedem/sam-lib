C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRSEN   ( MPAR  , MPMNPC, MMNPC , MADOF , MSC   ,
     $                        MPMCEQ, MMCEQ , MEQN  , NEL   , NANOD ,
     $                        NDOF  , NMMNPC, NSPAR , MSPAR , MASTER,
     $                        NUMNOD, SPRMAP                          )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEL   , NANOD , NDOF  , NMMNPC, NSPAR
      INTEGER           MPAR(50)              ,
     $                  MPMNPC(NEL+1)         , MMNPC(NMMNPC)         ,
     $                  MADOF(NANOD+1)        , MSC(NDOF)             ,
     $                  MPMCEQ(*)             , MMCEQ(*)              ,
     $                  MEQN(NDOF)            , MSPAR(NSPAR)          ,
     $                  MASTER(NDOF)          , NUMNOD(NANOD)         ,
     $                  SPRMAP(NDOF)

C ======================================================================
C  S A M  library routine :  SPRSEN                  GROUP 9 / PRIVATE
C ======================================================================
C
C     Purpose
C     -------
C
C --- To compute the size of the element-node connectivity array that
C     captures all the dependencies in the finite-element mesh.
C     It is assumed by the routine that a node is represented only once
C     for each element in MMNPC and that NANOD .LE. NDOF.
C
C     The routine first scan MADOF to find the number of SPR-nodes in
C     the finite-element mesh.
C
C     For each original finite-element node we have the following
C     possibilities:
C        - Each dof with MSC=1 that is master creates a node of size 1.
C        - Each dof with MSC=2 that is master creates a node of size 1.
C        - The dofs with MSC=1 that are not masters create one node of
C          size the number of dofs with this status in the node.
C        - The dofs with MSC=2 that are not masters create one node of
C          size the number of dofs with this status in the node.
C
C     This mean that an original finite-element node may be divided
C     into at most four different types. The reason for making masters
C     singletons is to simplify the process. It should be noted that the
C     minimum-degree ordering step merge nodes that may be merged.
C
C     Next, and based on the SPR-node partition it computes the size
C     of the representation.
C
C     The total workspace needed is: 2*NDOF+NANOD.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
C                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 Jan. 16, 2003 (kmo)
C                 Initialized SPRMAP to zero. Initialized MSPAR(35:36).
C                 May. 16, 2003 (kmo)
C                 Removed duplicated code block for MSC=2 and introduced
C                 the variable ISC instead. Reduced NMSICA by SPRNOD.
C                 Sep. 21, 2004 (kmo)
C                 MSPAR(*) -> MSPAR(NSPAR), NSPAR is new input argument
C                 Oct. 10, 2005 (kmo)
C                 Check for MMCEQ(I).GT.0
C                 May. 25, 2018 (kmo)
C                 Moved some calculation of MASTER and NUMNOD to SPRSEQ.
C                 Check that free non-master DOFs are connected.
C
C     MPAR - INTEGER(50)
C      Entry : Control parameters.
C      Exit  : Not changed.
C     MPMNPC - INTEGER(NEL+1)
C      Entry : Pointers to MMNPC.
C      Exit  : Not changed.
C     MMNPC - INTEGER(NMMNPC)
C      Entry : Lists of nodes connected to each element.
C      Exit  : Not changed.
C     MADOF - INTEGER(NANOD+1)
C      Entry : Matrix of accumulated dofs.
C      Exit  : Not changed.
C     MSC - INTEGER(NDOF)
C      Entry : Matrix of status codes.
C      Exit  : Not changed.
C     MPMCEQ - INTEGER(NCEQ+1)
C      Entry : Pointer to MMCEQ, it is not referenced in NCEQ=0.
C      Exit  : Not changed.
C     MMCEQ - INTEGER(NMMCEQ)
C      Entry : Lists of slaves and corresponding master dofs.
C      Exit  : Not changed.
C     MEQN - INTEGER(NDOF)
C      Entry : Matrix of equation numbers.
C      Exit  : Not changed.
C     NEL - INTEGER
C      Entry : Number of original finite-elements.
C      Exit  : Not changed.
C     NANOD - INTEGER
C      Entry : Number of original finite-element nodes.
C      Exit  : Not changed.
C     NDOF - INTEGER
C      Entry : Number of variables in the finite-element mesh.
C      Exit  : Not changed.
C     NMMNPC - INTEGER
C      Entry : Length of MMNPC and NODEL.
C      Exit  : Not changed.
C     NSPAR - INTEGER
C      Entry : Length of MSPAR.
C      Exit  : Not changed.
C     MSPAR - INTEGER(NSPAR)
C      Entry : Not defined.
C      Exit  : Control array for SPR-use, for definition of contents,
C              see below.
C
C     Working arrays
C     --------------
C
C     MASTER - INTEGER(NDOF)
C      Entry : Contains integers > 0 for each master.
C      Exit  : In case of slaves, it marks SPR-nodes that are already
C              present in an element list, in order to remove
C              duplicate entries (need not be saved).
C     NUMNOD - INTEGER(NANOD)
C      Entry : Contains integers > 0 for each connected finite-element node.
C      Exit  : Holds the number of SPR-nodes that is generated
C              in each original finite-element node (need not be saved).
C     SPRMAP - INTEGER(NDOF)
C      Entry : Not defined.
C              Maps each dof to its SPR-node. Prescribed dofs are mapped
C              to the zero node, and slaves are mapped to -1 node.
C      Exit  : Need not be saved.
C
C     Procedures
C     ----------
C
C     None
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
C     Definition of MSPAR.
C     --------------------
C
C     ------------------------------------------------------------------

      INTEGER           ELMSZE, I     , IDOF  , IP    , ISC   , J     ,
     $                  K     , NCEQ  , NEQ   , NEWNOD, NMASTR, NODE  ,
     $                  SEPSZE, SPRCON, SPRFST, SPRNOD, NODSET(2)

C     ==================================================================
C     SPRCON THE MAXIMUM LENGTH OF THE SPR-VERSION OF THE ELEMENT-NODE
C            CONNECTIVITY ARRAY MMNPC.
C     SPRNOD THE NUMBER OF NODES PRESENT IN THE INTERNAL PARTITION
C            OF THE NODES.
C     ==================================================================

      DO 10 I = 1, NSPAR
         MSPAR(I) = 0
   10 CONTINUE

* --- The number of status 2 nodes, i.e. the size of the separator.
      SEPSZE = 0

C     -------------------------------------
C     GET INFO FROM THE CONTROL ARRAY MPAR.
C     -------------------------------------
      NCEQ = MPAR( 7)
      NEQ  = MPAR(11)

      NMASTR = 0
      SPRNOD = 0
      DO 500 J = 1, NANOD
         SPRFST = SPRNOD

C        ----------------------------------------
C        COMPUTE THE ARRAY SPRMAP THAT MAPS EACH
C        DOF TO ITS SPR-NODE AND THE ARRAY NUMNOD
C        THAT COUNTS THE NUMBER OF SPR-NODES IN
C        EACH ORIGINAL FINITE-ELEMENT NODE.
C        IN ADDITON THE NUMBER OF SPR-NODES IS
C        COMPUTED.
C        ----------------------------------------
         NODSET(1) = 0
         NODSET(2) = 0
         DO 400 IDOF = MADOF(J), MADOF(J+1)-1
            SPRMAP(IDOF) = 0
            IF ( MEQN(IDOF) .GT. 0 ) THEN
               ISC = MSC(IDOF)
               IF ( ISC .LT. 1 .OR. ISC .GT. 2 ) GOTO 400

C              -------------------
C              AN ACTIVE VARIABLE.
C              -------------------
               IF ( MASTER(IDOF) .GT. 0 ) THEN

C                 -----------------------------------
C                 A MASTER THAT CREATES ITS OWN NODE.
C                 -----------------------------------
                  NMASTR = NMASTR + 1
                  SPRNOD = SPRNOD + 1
                  SPRMAP(IDOF) = SPRNOD
                  IF ( ISC .EQ. 2 ) SEPSZE = SEPSZE + 1
               ELSE IF ( NUMNOD(J) .GT. 0 ) THEN ! kmo added 25.05.2018
C We now do this only if the node is connected to at least one finite-element
C to be in compliance with the processing in SPRSA1. Otherwise the generated
C data structures will be inconsistent in the case of master nodes that are not
C connected to finite-elements (but used as external nodes, for instance).
                  IF ( NODSET(ISC) .EQ. 0 ) THEN

C                    -----------------------------------
C                    ALL THE FREE STATUS 'ISC' DOFS THAT
C                    ARE NOT MASTERS CREATE ONE NODE.
C                    -----------------------------------
                     SPRNOD = SPRNOD + 1
                     NODSET(ISC) = SPRNOD
                     IF ( ISC .EQ. 2 ) SEPSZE = SEPSZE + 1
                  ENDIF
                  SPRMAP(IDOF) = NODSET(ISC)
               ENDIF
            ELSEIF ( MEQN(IDOF) .LT. 0 ) THEN

C              -----------------------------------
C              A PRESCRIBED OR DEPENDENT VARIABLE.
C              -----------------------------------
               IP = -MEQN(IDOF)
               IF ( MPMCEQ(IP+1) .GT. MPMCEQ(IP)+1 ) THEN

C                 -------------------------
C                 THE VARIABLE IS DEPENDENT
C                 NEGATE IN SPRMAP.
C                 -------------------------
                  SPRMAP(IDOF) = -1
               ELSE

C                 ---------------------------
C                 THE VARIABLE IS PRESCRIBED.
C                 ZERO OUT IN SPRMAP.
C                 ---------------------------
                  SPRMAP(IDOF) = 0
               ENDIF
            ENDIF
  400    CONTINUE

C        -------------------------------------
C        SET NUMBER OF SPR-NODES IN THIS NODE.
C        -------------------------------------

         NUMNOD(J) = SPRNOD-SPRFST

  500 CONTINUE

      SPRCON = 0
      IF ( NMASTR .EQ. 0 ) THEN

C        --------------------------------------
C        THERE ARE NO SLAVES PRESENT AND SPRCON
C        IS THE SUM OF THE ENTRIES IN NUMNOD.
C        NOTE THAT THIS LOOP IS WRONG IF THERE
C        ARE DUPLICATE REFERENCES TO A NODE IN
C        SOME NODE LIST FOR AN ELEMENT.
C        --------------------------------------
         DO 600 J = 1, NMMNPC
            NODE = MMNPC(J)
            SPRCON = SPRCON + NUMNOD(NODE)
  600    CONTINUE

      ELSE

C        ---------------------------------------
C        THIS LOOP GOES EXPLICITLY THROUGH EACH
C        ORIGINAL FINITE-ELEMENT AND EMBEDDS FOR
C        EACH SLAVE THE MASTERS THAT ARE NOT IN
C        ITS LIST INITIALLY.

C        FIRST, INITIALIZE MASTER <- 0 IN ORDER
C        TO FACILITATE MASK OF NODES THAT ARE
C        ALREADY ASSEMBLED INTO AN ELEMENT LIST.
C        ---------------------------------------
         DO 700 J = 1, NDOF
            MASTER(J) = 0
  700    CONTINUE

         DO 1100 J = 1, NEL
            ELMSZE = 0
            DO 1000 I = MPMNPC(J), MPMNPC(J+1)-1
               NODE = MMNPC(I)

               DO 900 IDOF = MADOF(NODE), MADOF(NODE+1)-1
                  NEWNOD = SPRMAP(IDOF)
                  IF ( NEWNOD .GT. 0 ) THEN

                     IF ( MASTER(NEWNOD) .NE. J ) THEN

C                       -------------------------------
C                       THIS NODE IS NOT YET ASSEMBLED.
C                       -------------------------------
                        MASTER(NEWNOD) = J
                        ELMSZE = ELMSZE + 1
                     ENDIF

                  ELSEIF ( NEWNOD .LT. 0 ) THEN

C                    ---------------------
C                    THIS IS A SLAVE NODE.
C                    ---------------------
                     IP = -MEQN(IDOF)
                     DO 800 K = MPMCEQ(IP)+1, MPMCEQ(IP+1)-1
                        IF ( MMCEQ(K) .LE. 0 ) GOTO 800
                        NEWNOD = SPRMAP(MMCEQ(K))
                        IF ( MASTER(NEWNOD) .NE. J ) THEN

C                          -------------------------------
C                          THIS NODE IS NOT YET ASSEMBLED.
C                          -------------------------------
                           MASTER(NEWNOD) = J
                           ELMSZE = ELMSZE + 1
                        ENDIF
  800                CONTINUE
                  ENDIF
  900          CONTINUE

 1000       CONTINUE
            SPRCON = SPRCON + ELMSZE
 1100    CONTINUE

      ENDIF

C     ---------------------------------------------------
C     UPDATE MSPAR WITH NMSICA, NEL, SPRNOD, SPRCON, NEQ.
C     ---------------------------------------------------
      MSPAR( 5) = NEL
      MSPAR( 6) = SPRNOD
      MSPAR( 7) = SPRCON
      MSPAR( 8) = NEQ

      MSPAR(35) = 4*SPRNOD + NEL + 2*SPRCON + NEQ + 3 ! NMSICA
      MSPAR(36) = 7*SPRNOD + NEL + 6                  ! NTREES
      MSPAR(40) = SEPSZE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSEN
      END


      SUBROUTINE   SPRSEL   ( NSPAR , MPAR  , MSPAR   )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSPAR
      INTEGER           MPAR(50)              , MSPAR(NSPAR)

C ======================================================================
C  S A M  library routine :  SPRSEL                  GROUP 9 / PRIVATE
C ======================================================================
C
C     Purpose
C     -------
C
C --- To compute the size of the element-node connectivity array that
C     captures all the dependencies in the finite-element mesh. It also
C     computes the size of MSICA that is needed to store the element
C     based representation of the coefficient matrix structure.
C
C     This is a lumped-matrix version of SPRSEN where we assume a
C     one-to-one correspondance between the "elements" and the degrees
C     of freedom. Thus, the only off-diagonal terms in the system matrix
C     emanate from the constraint equations, if any.
C
C     This assumption makes it possible to estimate the size of the
C     data structures without traversing the actual finite-element mesh.
C     It can be computed directly from the parameters stored in MPAR.
C
C     Created   : Mar. 14, 2003 (kmo)
C     Revisions : May. 16, 2003 (kmo)
C                 Reduced NMSICA by SPRNOD.
C                 Sep. 21, 2004 (kmo)
C                 MSPAR(*) -> MSPAR(NSPAR), NSPAR is new input argument
C
C     NSPAR - INTEGER
C      Entry : Length of MSPAR.
C      Exit  : Not changed.
C     MPAR - INTEGER(50)
C      Entry : Control parameters.
C      Exit  : Not changed.
C     MSPAR - INTEGER(NSPAR)
C      Entry : Not defined.
C      Exit  : Control array for SPR-use.
C
C     Working arrays
C     --------------
C
C     None
C
C     Procedures
C     ----------
C
C     None
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

      INTEGER           I     , NCEQ  , NEQ   , NDOF  , NDOF2 , NMMCEQ,
     $                  SEPSZE, SPRCON, SPRNOD

C     ==================================================================
C     SPRCON THE MAXIMUM LENGTH OF THE SPR-VERSION OF THE ELEMENT-NODE
C            CONNECTIVITY ARRAY MMNPC.
C     SPRNOD THE NUMBER OF NODES PRESENT IN THE INTERNAL PARTITION
C            OF THE NODES.
C     ==================================================================

      DO 10 I = 1, NSPAR
         MSPAR(I) = 0
   10 CONTINUE

C     -------------------------------------
C     GET INFO FROM THE CONTROL ARRAY MPAR.
C     -------------------------------------
      NDOF   = MPAR( 3)
      NDOF2  = MPAR( 5)
      NCEQ   = MPAR( 7)
      NEQ    = MPAR(11)
      NMMCEQ = MPAR(16)

C     -------------------------------------------------------
C     DEFINE SPRNOD, SEPSZE AND SPRCON FOR THE LUMPED MATRIX.
C     -------------------------------------------------------
      SPRNOD = NEQ
      SEPSZE = NDOF2
      SPRCON = NEQ + NMMCEQ - NCEQ

C     ----------------------------------------------------
C     UPDATE MSPAR WITH NMSICA, NDOF, SPRNOD, SPRCON, NEQ.
C     ----------------------------------------------------
      MSPAR( 5) = NDOF
      MSPAR( 6) = SPRNOD
      MSPAR( 7) = SPRCON
      MSPAR( 8) = NEQ

      MSPAR(35) = 4*SPRNOD + NDOF + 2*SPRCON + NEQ + 3 ! NMSICA
      MSPAR(36) = 7*SPRNOD + NDOF + 6                  ! NTREES
      MSPAR(40) = SEPSZE

C     ----------------------------------
C     FLAG THAT THIS IS A LUMPED MATRIX.
C     ----------------------------------
      MSPAR(51) = 1

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSEL
      END
