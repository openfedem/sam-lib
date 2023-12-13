C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRSAS   ( MPAR  , MPMNPC, MMNPC , MADOF , MSC   ,
     $                        MPMCEQ, MMCEQ , MEQN  , MSPAR , MSICA ,
     $                        IWORK , NSPAR , LPU   , IERR    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSPAR , LPU   , IERR
      INTEGER           MPAR(50)              , MPMNPC(*)             ,
     $                  MMNPC(*)              , MADOF(*)              ,
     $                  MSC(*)                , MPMCEQ(*)             ,
     $                  MMCEQ(*)              , MEQN(*)               ,
     $                  MSPAR(NSPAR)          , MSICA(*)              ,
     $                  IWORK(*)

C ======================================================================
C  S A M  library routine :  SPRSAS                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To partition the work arrays for routine SPRSA1 that computes the
C     element-based coefficient matrix representation and the SPR-node
C     partition.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Aug. 24, 1998 (acd)
C                 MSPAR(50) -> MSPAR(*), and added superelement option.
C                 Mar. 06, 2003 (kmo)
C                 Added option for lumped element matrices.
C                 May. 12, 2003 (kmo)
C                 Check whether SPRCON is less than anticipated, and
C                 move array SEPARR if that is the case.
C                 May. 16, 2003 (kmo)
C                 Reduced the size of NODEL by SPRNOD.
C                 May. 25, 2018 (kmo)
C                 Added temporary variable INTSZE.
C                 Aug. 14, 2021 (kmo)
C                 MSPAR(*) -> MSPAR(NSPAR), NSPAR is new input argument
C
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
C     MSPAR - INTEGER(NSPAR)
C      Entry : Not defined.
C      Exit  : Control array for SPR-use, for definition of contents,
C              see in routine SPRSA1.
C     MSICA - INTEGER(NMSICA)
C      Entry : Not defined.
C      Exit  : Stores the SPR-representation of the coefficient matrix
C              that is finished at this stage.
C     NSPAR - INTEGER
C      Entry : Length of MSPAR.
C      Exit  : Not changed.
C     LPU - INTEGER
C      Entry : Output unit for error messages.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined.
C      Exit  : Error flag, is set to zero in case of a normal return.
C              = -1: Number of SPR-nodes computed here does not match
C                    the number computed in SPRSEN.
C              = -2: Array NODES is not large enough. This happens
C                    if there are input errors in MADOF.
C              = -3: Array ELNOD is not large enough.
C              = -4: MPAR( 2) .ne. MSPAR(5) and/or
C                    MPAR(11) .ne. MSPAR(8).
C
C     Working arrays
C     --------------
C
C     IWORK - INTEGER(2*NDOF+NANOD)
C      Entry : Not defined.
C              For use, see definition of the work arrays in SPRSA1.
C      Exit  : Need not be saved.
C
C     Procederes
C     ----------
C
C     SPRSA1
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
C     Definition of MSPAR.
C     --------------------
C
C     ------------------------------------------------------------------

      INTEGER           NANOD , NCEQ  , NDOF  , NEL   , NEQ   ,
     $                  NMMNPC, NMMCEQ, SPRCON, SPRNOD

      INTEGER           ELNOD , MARKER, MASTER, NODEL , NODES , PERM  ,
     $                  SPRMAP, SEPARR, SEPNEW, XELNOD, XNODEL, XNODES

      INTRINSIC         MAX
      EXTERNAL          SPRSA1, SPRER1

C     ==================================================================

      IERR   = 0

C     ----------------------------------------
C     GET INFORMATION FROM CONTROL ARRAY MPAR.
C     ----------------------------------------
      NANOD  = MPAR( 1)
      IF ( NSPAR .GT. 50 .AND. MSPAR(51) .EQ. 1 ) THEN
         NEL = MPAR( 3)
      ELSE
         NEL = MPAR( 2)
      ENDIF
      NDOF   = MPAR( 3)
      NCEQ   = MPAR( 7)
      NEQ    = MPAR(11)
      NMMNPC = MPAR(15)
      NMMCEQ = MPAR(16)

C     -----------------------------------------
C     GET INFORMATION FROM CONTROL ARRAY MSPAR.
C     -----------------------------------------
      SPRNOD = MSPAR(6)
      SPRCON = MSPAR(7)

C     -----------------------------------------
C     CHECK CONSISTENCY BETWEEN MPAR AND MSPAR.
C     -----------------------------------------
      IF ( NEL .NE. MSPAR(5) .OR. NEQ .NE. MSPAR(8) ) THEN
         MSPAR(1) = -1
         IERR = -4
         CALL SPRER1 ( 21, 'SPRSAS', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

C     -----------------------------------
C     DEFINE POINTERS TO ARRAYS IN MSICA.
C     -----------------------------------

      XELNOD = 1
      XNODEL = XELNOD + NEL + 1
      XNODES = XNODEL + SPRNOD + 1
      NODES  = XNODES + SPRNOD + 1
      PERM   = NODES  + NEQ
      ELNOD  = PERM   + SPRNOD
      NODEL  = ELNOD  + SPRCON
      SEPARR = NODEL  + SPRCON

C     ----------------------------------------
C     DEFINE POINTERS TO WORK ARRAYS IN IWORK.
C     ----------------------------------------

      MASTER = 1
      SPRMAP = MASTER + NDOF
      MARKER = SPRMAP + NDOF

C     --------------------------------------
C     CALL SPRSA1 TO COMPUTE (XELNOD,ELNOD),
C     (XNODES,NODES) AND TO INITIALZE PERM.
C     --------------------------------------
      CALL SPRSA1 (MPMNPC,MMNPC,MADOF,MSC,MPMCEQ,MMCEQ,MEQN,
     $             NEL,NANOD,NDOF,NEQ,NCEQ,NMMNPC,NMMCEQ,NSPAR,
     $             SPRCON,SPRNOD,MSPAR,MSICA(XELNOD),MSICA(ELNOD),
     $             MSICA(XNODES),MSICA(NODES),MSICA(PERM),MSICA(SEPARR),
     $             IWORK(MASTER),IWORK(SPRMAP),IWORK(MARKER),LPU,IERR)
C
      IF ( MSPAR(7) .LT. SPRCON .AND. IERR .EQ. 0 ) THEN
C
CKMO --- The size of ELNOD and NODEL is less than anticipated.
C        Warn the user, and move array SEPARR accordingly.
         IERR   = 1
         NODEL  = ELNOD + MSPAR(7)
         SEPNEW = NODEL + MSPAR(7)
         CALL SPRCPY (SPRNOD,MSICA(SEPARR),1,MSICA(SEPNEW),1)
         CALL SPRER1 (61,'SPRSAS',SPRCON,MSPAR(7),0,LPU,IERR)
      ENDIF

      MSPAR( 2) = 4*SPRNOD + 2*MSPAR(7) + NEL + NEQ + 3

*ACD  MSPAR(37) = 8*MSPAR(6) + 2*MSPAR(7) + 3*MSPAR(5)
      MSPAR(37) = 8*MSPAR(6) + 2*MSPAR(7) + 3*MSPAR(5) + NEQ + 1
* --- New length due to extensions in the ordering step.
* --- HUSK at dette er NEQ+1 i MMDw og at det er MAXFIL i MMFw.
*     Dette betyr at MMFw ikke kan allokeres MAXFIL .gt. NEQ + 1
*     i MFRRNM!!!

*ACD  MSPAR(38) = 8*MSPAR(6) + MAX(MSPAR(6),MSPAR(5)) + 1
      MSPAR(38) = 9*MSPAR(6) + MAX(MSPAR(6),MSPAR(5)) + 1
* ---------------- +MSPAR(6) to hold copy of SUPSUP.

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSAS
      END


      SUBROUTINE   SPRSA1   ( MPMNPC, MMNPC , MADOF , MSC   , MPMCEQ,
     $                        MMCEQ , MEQN  , NEL   , NANOD , NDOF  ,
     $                        NEQ   , NCEQ  , NMMNPC, NMMCEQ, NSPAR ,
     $                        SPRCON, SPRNOD, MSPAR , XELNOD, ELNOD ,
     $                        XNODES, NODES , PERM  , SEPARR, MASTER,
     $                        SPRMAP, MARKER, LPU   , IERR    )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NEL   , NANOD , NDOF  , NEQ   , NCEQ  , NMMNPC,
     $                  NMMCEQ, NSPAR , SPRCON, SPRNOD, LPU   , IERR
      INTEGER           MPMNPC(NEL+1)         , MMNPC(NMMNPC)         ,
     $                  MADOF(NANOD+1)        , MSC(NDOF)             ,
     $                  MPMCEQ(NCEQ+1)        , MMCEQ(NMMCEQ)         ,
     $                  MEQN(NDOF)            , MSPAR(NSPAR)          ,
     $                  XELNOD(NEL+1)         , ELNOD(SPRCON)         ,
     $                  XNODES(SPRNOD+1)      , NODES(NEQ)            ,
     $                  PERM(SPRNOD)          , SEPARR(SPRNOD)        ,
     $                  MASTER(NDOF)          , SPRMAP(NDOF)          ,
     $                  MARKER(NANOD)

C ======================================================================
C  S A M  library routine :  SPRSA1                  GROUP 9 / PRIVATE
C ======================================================================
C
C     PURPOSE
C     -------
C
C --- THE MAIN TASK OF THIS ROUTINE IS TO COMPUTE THE SPR-VERSION OF
C     MMNPC, THIS STRUCTURE IS STORED IN THE ARRAY PAIR (XELNOD,ELNOD).
C     IT ALSO COMPUTES THE SPR-NODE PARTITION WHICH IS STORED IN THE
C     ARRAY PAIR (XNODES,NODES), AND INITIALIZE THE PERMUTATION VECTOR
C     PERM TO THE IDENTITY PERMUTATION. FINALLY THE ROUTINE INITIALIZES
C     MSPAR AND STORE THE SPR-PARAMETERS THAT ARE COMPUTED HERE.
C
C     IT IS ASSUMED BY THE ROUTINE THAT A NODE IS REPRESENTED ONLY ONCE
C     FOR EACH ELEMENT IN MMNPC AND THAT NANOD .LE. NDOF, AND THAT A DOF
C     IS REPRESENTED BY A UNIQUE INTEGER. I.E. THE INTEGERS IMPLICITLY
C     STORED IN MADOF.
C
C     THE ROUTINE FIRST SCAN MADOF TO FIND THE NUMBER OF SPR-NODES IN
C     THE FINITE ELEMENT MESH.
C
C     FOR EACH FINITE ORIGINAL FINITE-ELEMENT NODE WE HAVE THE FOLLOWING
C     POSSIBILITIES:
C        - EACH DOF WITH MSC=1 THAT IS MASTER CREATES A NODE OF SIZE 1.
C        - EACH DOF WITH MSC=2 THAT IS MASTER CREATES A NODE OF SIZE 1.
C        - THE DOFS WITH MSC=1 THAT ARE NOT MASTERS CREATE ONE NODE OF
C          SIZE THE NUMBER OF DOFS WITH THIS STATUS IN THE NODE.
C        - THE DOFS WITH MSC=2 THAT ARE NOT MASTERS CREATE ONE NODE OF
C          SIZE THE NUMBER OF DOFS WITH THIS STATUS IN THE NODE.
C
C     THIS MEAN THAT AN ORIGINAL FINITE-ELEMENT NODE MAY BE DIVIDED
C     INTO AT MOST FOUR DIFFERENT TYPES. THE REASON FOR MAKING MASTERS
C     SINGLETONS IS TO SIMPLIFY THE PROCESS. IT SHOULD BE NOTED THAT THE
C     MINIMUM-DEGREE ORDERING STEP MERGE NODES THAT MAY BE MERGED SO THI
C     IS NO BIG PROBLEM.
C
C     NEXT, AND BASED ON THE SPR-NODE PARTITION IT COMPUTES THE
C     SPR-VERSION OF MMNPC.
C
C     CREATED   : JAN. 11, 1994 (ACD)
C     REVISIONS : Mar. 06, 2003 (kmo)
C                 Added option for lumped element matrices.
C                 Moved some stuff to the new SPRSA2.
C                 May. 12, 2003 (kmo)
C                 Removed duplicated code blocks for MSC=2, introduced
C                 the variable ISC and associated arrays of length 2.
C                 Jul. 31, 2003 (kmo)
C                 Account for master DOFs that are not directly
C                 associated with elements. Minor error message changes.
C                 Oct. 10, 2005 (kmo)
C                 Check for MMCEQ(I).GT.0
C                 May. 25, 2018 (kmo)
C                 Check for SPRMAP(I).GT.0 for masters in constraint equations.
C                 Aug. 14, 2021 (kmo)
C                 MMCEQ(*) -> MMCEQ(NMMCEQ), NMMCEQ is new input argument
C                 MSPAR(*) -> MSPAR(NSPAR), NSPAR is new input argument
C
C
C     MPMNPC - INTEGER(NEL+1)
C      ENTRY : POINTERS TO MMNPC.
C      EXIT  : NOT CHANGED.
C     MMNPC - INTEGER(NMMNPC)
C      ENTRY : LISTS OF NODES CONNECTED TO EACH ELEMENT.
C      EXIT  : NOT CHANGED.
C     MADOF - INTEGER(NANOD+1)
C      ENTRY : MATRIX OF ACCUMULATED DOFS.
C      EXIT  : NOT CHANGED.
C     MSC - INTEGER(NDOF)
C      ENTRY : MATRIX OF STATUS CODES.
C      EXIT  : NOT CHANGED.
C     MPMCEQ - INTEGER(NCEQ+1)
C      ENTRY : POINTER TO MMCEQ, IT IS NOT REFERENCED IF NCEQ=0.
C      EXIT  : NOT CHANGED.
C     MMCEQ - INTEGER(NMMCEQ)
C      ENTRY : LISTS OF SLAVES AND CORRESPONDING MASTER DOFS.
C      EXIT  : NOT CHANGED.
C     MEQN - INTEGER(NDOF)
C      ENTRY : MATRIX OF EQUATION NUMBERS.
C      EXIT  : NOT CHANGED.
C     NEL - INTEGER
C      ENTRY : NUMBER OF ORIGINAL FINITE-ELEMENTS.
C      EXIT  : NOT CHANGED.
C     NANOD - INTEGER
C      ENTRY : NUMBER OF ORIGINAL FINITE-ELEMENT NODES.
C      EXIT  : NOT CHANGED.
C     NDOF - INTEGER
C      ENTRY : NUMBER OF VARIABLES IN THE FINITE-ELEMENT MESH.
C      EXIT  : NOT CHANGED.
C     NEQ - INTEGER
C      ENTRY : THE NUMBER OF EQUATIONS.
C      EXIT  : NOT CHANGED.
C     NCEQ - INTEGER
C      ENTRY : THE NUMBER OF CONSTRAINT EQUATIONS, MAY BE ZERO.
C      EXIT  : NOT CHANGED.
C     NMMNPC - INTEGER
C      ENTRY : LENGTH OF MMNPC.
C      EXIT  : NOT CHANGED.
C     NMMCEQ - INTEGER
C      Entry : Length of MMCEQ.
C      Exit  : Not changed.
C     NSPAR - INTEGER
C      Entry : Length of MSPAR.
C      Exit  : Not changed.
C     SPRCON - INTEGER
C      ENTRY : THE SIZE OF THE SPR-VERSION OF MMNPC AS COMPUTED IN
C              SPRSEN.
C      EXIT  : NOT CHANGED.
C     SPRNOD - INTEGER
C      ENTRY : THE NUMBER OF NODES PRESENT IN THE INTERNAL PARTITION
C              OF THE NODES.
C      EXIT  : NOT CHANGED.
C     MSPAR - INTEGER(NSPAR)
C      ENTRY : AS SET IN SPRSEN.
C      EXIT  : CONTROL ARRAY FOR SPR-USE, FOR DEFINITION OF CONTENTS,
C              SEE BELOW.
C     XELNOD - INTEGER(NEL+1)
C      ENTRY : NOT DEFINED.
C      EXIT  : POINTER TO LIST OF SPR-NODES FOR EACH ELEMENT IN ELNOD.
C     ELNOD - INTEGER(SPRCON)
C      ENTRY : NOT DEFINED.
C      EXIT  : IF IERR=0, THIS ARRAY STORES THE LIST OF SPR-NODES FOR
C              EACH ORIGINAL FINITE-ELEMENT.
C     XNODES - INTEGER(SPRNOD+1)
C      ENTRY : NOT DEFINED.
C      EXIT  : IF IERR=0, THIS ARRAY STORES THE SPR-NODE PARTITION, I.E.
C              THIS ARRAY IS THE SPR VERSION OF MADOF.
C     NODES - INTEGER(NEQ)
C      ENTRY : NOT DEFINED.
C      EXIT  : THE LIST OF DOFS THAT ARE CONNECTED TO EACH SPR-NODE.
C     PERM - INTEGER(SPRNOD)
C      ENTRY : NOT DEFINED.
C      EXIT  : IF IERR=0, IT HOLDS AN IDENTITY PERMUTATION OF THE SPR-
C              NODES. THIS IS NEEDED IF THE ORDERING ROUTINE SPRRNM IS
C              NOT USED.
C     LPU - INTEGER
C      Entry : Output unit for error messages.
C      Exit  : Not changed.
C     IERR - INTEGER
C      ENTRY : NOT DEFINED.
C      EXIT  : ERROR FLAG, IS SET TO ZERO IN CASE OF A NORMAL RETURN.
C              = -1: NUMBER OF SPR-NODES COMPUTED HERE DOES NOT MATCH
C                    THE NUMBER COMPUTED IN SPRSEN.
C              = -2: ARRAY NODES IS NOT LARGE ENOUGH. THIS HAPPENS
C                    IF THERE ARE INPUT ERRORS IN MADOF.
C              = -3: ARRAY ELNOD IS NOT LARGE ENOUGH.
C
C     WORKING ARRAYS
C     --------------
C
C     MASTER - INTEGER(NDOF)
C      ENTRY : NOT DEFINED.
C              FIRST, THIS ARRAY STORES INTEGERS > 0 FOR EACH MASTER.
C              SECOND, IN CASE OF SLAVES, THIS ARRAY MARKS SPR-NODES
C              THAT ARE ALREADY PRESENT IN AN ELEMENT LIST IN ORDER
C              TO REMOVE DUPLICATE ENTRIES.
C      EXIT  : NEED NOT BE SAVED.
C     SPRMAP - INTEGER(NDOF)
C      ENTRY : NOT DEFINED.
C              MAPS EACH DOF TO ITS SPR-NODE. PRESCRIBED DOFS ARE MAPPED
C              TO THE ZERO NODE, AND SLAVES ARE MAPPED TO -1 NODE. THIS
C              USE OF THE ARRAY REQUIRES THAT ALL DOFS ARE REPRESENTED
C              BY UNIQUELY DEFINED INTEGERS, BUT THEY NEED NOT BE
C              CONTIGUOUS.
C              IN THE SET UP OF (XNODES,NODES) IT IS ALSO USED AS A
C              LINKED LIST THAT HOLDS, FOR EACH NODE, THE DOFS THAT ARE
C              NOT MASTERS FOR STATUS 1 AND 2 DOFS. THIS USE OF THE
C              ARRAY REQUIRES THAT EACH DOF IN A NODE IS UNIQUELY
C              REPRESENTED WITH RESPECT TO DOF NUMBERS.
C      EXIT  : NEED NOT BE SAVED.
C     MARKER - INTEGER(NANOD)
C      Entry : Not defined.
C              The array holds MARKER(I) = 0 for each node that is not
*              connected to an element and so is not active in the
*              model. The active nodes have MARKER(I) > 0.
C
C     PROCEDURES
C     ----------
C
C     SPRSA2
C     SPRER1
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
C     DEFINITION OF MSPAR.
C     --------------------
C
C     ------------------------------------------------------------------

      INTEGER           I     , IDOF  , IEQ   , IP    , IPNODS, ISC   ,
     $                  J     , NODE  , NUMNOD, SEPSZE, INTSZE, SPRCHK,
     $                  FIRST(2), FSTNOD(2), NEXT(2), NUM(2)

      LOGICAL           LERROR, LMASTR, LUMPED

      INTRINSIC         MAX
      EXTERNAL          SPRSA2, SPRER1

C     ==================================================================

      IERR = 0

* --- The number of MFR-Nodes that are supernodes in Sesam terms.
      SEPSZE = MSPAR(40)
      INTSZE = SPRNOD - SEPSZE

* --- For use when initializing PERM.
      FSTNOD(1) = 0
      FSTNOD(2) = INTSZE
      LMASTR = .FALSE.
      LUMPED = NSPAR .GT. 50 .AND. MSPAR(51) .EQ. 1

C     ------------------
C     INITIALIZE MASTER.
C     ------------------
      DO 100 I = 1, NDOF
         MASTER(I) = 0
         SPRMAP(I) = 0
  100 CONTINUE

C     ---------------------------------
C     COMPUTE THE ARRAY MASTER THAT HAS
C     MASTER > 0 FOR VARIABLES THAT ARE
C     MASTERS IN A CONSTRAINT EQUATION.
C     ---------------------------------
      DO 300 J = 1, NCEQ
         DO 200 I = MPMCEQ(J)+1, MPMCEQ(J+1)-1
            IDOF  = MMCEQ(I)
            IF ( IDOF .GT. 0 ) THEN
               MASTER(IDOF) = MASTER(IDOF) + 1
               LMASTR = .TRUE.
            ENDIF
  200    CONTINUE
  300 CONTINUE

* --- Mask off the nodes that are not associated with elements.
      DO I = 1, NANOD
         MARKER(I) = 0
      END DO
      DO I = 1, NMMNPC
         NODE = MMNPC(I)
         IF ( NODE .GT. 0 .AND. NODE .LE. NANOD ) THEN
            MARKER(NODE) = MARKER(NODE) + 1
         ENDIF
      END DO

C     ------------------------------------------
C     SCAN MADOF AND COMPUTE THE NODE PARTITION.
C     THE NODE PARTITION FOR INTERNAL SPR USE IS
C     STORED IN THE ARRAY PAIR (XNODES,NODES).
C     N O T E: THIS LOOP WILL NOT WORK IF THE
C              DOF IDENTIFICATORS IN MADOF ARE
C              NOT UNIQUE.
C     ------------------------------------------
      SPRCHK = 0
      IPNODS = 1
      XNODES(1) = 1

      DO 700 J = 1, NANOD

* ------ An active node.

         FIRST(1) = -1
         FIRST(2) = -1
         NUM(1) = 0
         NUM(2) = 0
         DO 400 IDOF = MADOF(J), MADOF(J+1)-1
            IEQ = MEQN(IDOF)
            ISC = MSC(IDOF)
            IF ( IEQ .GT. 0 ) THEN
               IF ( ISC .LT. 1 .OR. ISC .GT. 2 ) GOTO 400

C              -------------------
C              AN ACTIVE VARIABLE.
C              -------------------
               IF ( MASTER(IDOF) .GT. 0 .OR. LUMPED ) THEN

C                 -----------------------------------
C                 A MASTER THAT CREATES ITS OWN NODE.
C                 -----------------------------------
                  SPRCHK = SPRCHK + 1
                  IF ( SPRCHK .GT. SPRNOD ) THEN
                     MSPAR(1) = -1
                     IERR = -1
                     CALL SPRER1 ( 32,'SPRSA1',0,0,0,LPU,IERR )
                     RETURN
                  ELSE IF ( IPNODS .GT. NEQ ) THEN
                     MSPAR(1) = -1
                     IERR = -2
                     CALL SPRER1 ( 32,'SPRSA1',0,0,0,LPU,IERR )
                     RETURN
                  ENDIF

* --------------- Initialize in PERM.
                  FSTNOD(ISC) = FSTNOD(ISC) + 1
                  PERM(FSTNOD(ISC)) = SPRCHK

C                 ---------------------------
C                 MAP THE DOF TO ITS SPR-NODE
C                 AND UPDATE (XNODES,NODES).
C                 ---------------------------
                  SPRMAP(IDOF) = SPRCHK
                  NODES(IPNODS) = IDOF
                  IPNODS = IPNODS + 1
                  XNODES(SPRCHK+1) = IPNODS

               ELSE IF ( MARKER(J) .GT. 0 ) THEN

C                 -----------------------------------------
C                 ALL NON-MASTERS CREATE ONE NODE TOGETHER.
C                 -----------------------------------------
                  IF ( FIRST(ISC) .EQ. -1 ) THEN

C                    ----------------------------------
C                    THIS IS THE FIRST DOF, INITIALIZE
C                    THE LINKED LIST OF FREE VARIABLES.
C                    ----------------------------------
                     FIRST(ISC) = IDOF
                     NEXT(ISC) = IDOF
                     NUM(ISC) = 1
                  ELSE

C                    ---------------------------------
C                    BUILD THE LIST OF FREE VARIABLES.
C                    ---------------------------------
                     SPRMAP(NEXT(ISC)) = IDOF
                     NEXT(ISC) = IDOF
                     NUM(ISC) = NUM(ISC) + 1
                  ENDIF
               ENDIF

            ELSEIF ( IEQ .LT. 0 .AND. MARKER(J) .GT. 0 ) THEN

C              -----------------------------------
C              A PRESCRIBED OR DEPENDENT VARIABLE.
C              -----------------------------------
               IP = -IEQ
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

C        -----------------------------------------------------
C        NOW LINK IN NON-MASTERS OF STATUS 1 AND 2 IF PRESENT.
C        -----------------------------------------------------
         DO 600 ISC = 1, 2
            IF ( NUM(ISC) .EQ. 0 ) GOTO 600

C           ----------------------------------------------
C           THERE ARE NON-MASTERS OF STATUS 'ISC' PRESENT.
C           ----------------------------------------------
            SPRCHK = SPRCHK + 1
            IF ( SPRCHK .GT. SPRNOD ) THEN
               MSPAR(1) = -1
               IERR = -1
               CALL SPRER1 ( 32, 'SPRSA1', 0, 0, 0, LPU, IERR )
               RETURN
            ELSE IF ( IPNODS+NUM(ISC)-1 .GT. NEQ ) THEN
               MSPAR(1) = -1
               IERR = -2
               CALL SPRER1 ( 32, 'SPRSA1', 0, 0, 0, LPU, IERR )
               RETURN
            ENDIF

* --------- Initialize in PERM.
            FSTNOD(ISC) = FSTNOD(ISC) + 1
            PERM(FSTNOD(ISC)) = SPRCHK

C           -------------------------------
C           LOOP THE LINKED LIST AND UPDATE
C           SPRMAP, AND (XNODES,NODES).
C           -------------------------------
            DO 500 I = 1, NUM(ISC)
               NODES(IPNODS) = FIRST(ISC)
               NEXT(ISC) = SPRMAP(FIRST(ISC))
               SPRMAP(FIRST(ISC)) = SPRCHK
               FIRST(ISC) = NEXT(ISC)
               IPNODS = IPNODS + 1
  500       CONTINUE
            XNODES(SPRCHK+1) = IPNODS

  600    CONTINUE

  700 CONTINUE

      IF ( SPRCHK .LT. SPRNOD ) THEN
         MSPAR(1) = -1
         IERR = -1
         CALL SPRER1 ( 22, 'SPRSA1', -1, SPRCHK, SPRNOD, LPU, IERR )
         RETURN
      ENDIF

C     -----------------------------
C     COMPUTE SPR-VERSION OF MMNPC.
C     -----------------------------
      IPNODS = 1
      XELNOD(1) = 1
      NUMNOD = 0

C     -------------------------------
C     INITIALIZE MASTER TO FACILITATE
C     MARKING OF THE ASSEMBLED NODES.
C     -------------------------------
      DO 800 J = 1, NDOF
         MASTER(J) = 0
  800 CONTINUE

      IF ( LUMPED ) THEN

C        ---------------------------------------------------------
C        RUN THROUGH THE LIST OF DOFS AND BUILD (XELNOD,ELNOD).
C        In the lumped element matrix case we have a one-to-one
C        correspondance between DOFs and elements in the SPR-grid.
C        ---------------------------------------------------------

         DO 900 J = 1, NDOF

            CALL SPRSA2 (J,J,NDOF,NCEQ,NMMCEQ,SPRCON,
     $                   SPRMAP,MEQN,MPMCEQ,MMCEQ,
     $                   IPNODS,ELNOD,MASTER,IERR)
            IF ( IERR .LT. 0 ) THEN
               MSPAR(1) = -1
               CALL SPRER1 ( 32,'SPRSA1',0,0,0,LPU,IERR )
               RETURN
            ENDIF

            XELNOD(J+1) = IPNODS
            NUMNOD = MAX(NUMNOD,IPNODS-XELNOD(J))

  900    CONTINUE

      ELSE

         DO 1200 J = 1, NEL

C           -----------------------------------------
C           RUN THROUGH THE LIST OF NODES/DOFS
C           IN EACH ELEMENT AND BUILD (XELNOD,ELNOD).
C           -----------------------------------------
            DO 1100 I = MPMNPC(J), MPMNPC(J+1)-1
               NODE = MMNPC(I)
               IF ( MARKER(NODE) .EQ. 0 ) GOTO 1100

* ------------ An active node.

               DO 1000 IDOF = MADOF(NODE), MADOF(NODE+1)-1

                  CALL SPRSA2 (IDOF,J,NDOF,NCEQ,NMMCEQ,SPRCON,
     $                         SPRMAP,MEQN,MPMCEQ,MMCEQ,
     $                         IPNODS,ELNOD,MASTER,IERR)
                  IF ( IERR .LT. 0 ) THEN
                     MSPAR(1) = -1
                     CALL SPRER1 ( 32,'SPRSA1',0,0,0,LPU,IERR )
                     RETURN
                  ENDIF

 1000          CONTINUE

 1100       CONTINUE

            XELNOD(J+1) = IPNODS
            NUMNOD = MAX(NUMNOD,IPNODS-XELNOD(J))

 1200    CONTINUE

      ENDIF

* *** CONTROL *** Dette er meget viktig.
      IF ( LMASTR .AND. SEPSZE .GT. 0 ) THEN

* ------ Sort the two sets of PERM to relative initial order.
         IF ( INTSZE .GT. 1 ) THEN

* --------- The first set is for internal nodes.
            CALL SPRT48 ( INTSZE, PERM, LERROR)
            IF ( LERROR ) THEN
               IERR = -41
               RETURN
            ENDIF
         ENDIF
         IF ( SEPSZE .GT. 1 ) THEN

* --------- The second set is for external nodes,
*           so-called supernodes in Sesam.
            CALL SPRT48 ( SEPSZE, PERM(INTSZE+1), LERROR )
            IF ( LERROR ) THEN
               IERR = -42
               RETURN
            ENDIF
         ENDIF

      ENDIF
* *** CONTROL ***

* --- Initialize SEPARR.
      DO I = 1, SPRNOD
         NODE = PERM(I)
         IF (I .LE. INTSZE) THEN
            SEPARR(NODE) = 0
         ELSE
            SEPARR(NODE) = 1
         ENDIF
      END DO

C     ---------------
C     STORE IN MSPAR.
C     ---------------
      MSPAR( 1) = 1
      MSPAR( 2) = 0
      MSPAR( 3) = 0
      MSPAR( 4) = 0
      MSPAR( 7) = IPNODS - 1
      MSPAR(24) = NUMNOD

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRSA1
      END


      SUBROUTINE SPRSA2 (IDOF,IEL,NDOF,NCEQ,NMMCEQ,SPRCON,
     $                   SPRMAP,MEQN,MPMCEQ,MMCEQ,
     $                   IPNODS,ELNOD,MASTER,IERR)

      IMPLICIT NONE

      INTEGER IDOF,IEL,NDOF,NCEQ,NMMCEQ,SPRCON,
     $        SPRMAP(NDOF),MEQN(NDOF),MPMCEQ(NCEQ+1),MMCEQ(NMMCEQ),
     $        IPNODS,ELNOD(SPRCON),MASTER(NDOF),IERR
      INTEGER ICEQ,IPCEQ,NEWNOD,K

      IF ( IDOF .LT. 1 .OR. IDOF .GT. NDOF ) GOTO 900

      IERR   = 0
      NEWNOD = SPRMAP(IDOF)

      IF ( NEWNOD .GT. 0 ) THEN

C        ----------------------
C        THIS IS A MASTER NODE.
C        ----------------------
         IF ( NEWNOD .GT. NDOF ) GOTO 900
         IF ( MASTER(NEWNOD) .EQ. IEL ) RETURN

C        -------------------------------
C        THIS NODE IS NOT YET ASSEMBLED.
C        -------------------------------
         IF ( IPNODS .GT. SPRCON ) GOTO 900
         MASTER(NEWNOD) = IEL
         ELNOD(IPNODS)  = NEWNOD
         IPNODS = IPNODS + 1

      ELSEIF ( NEWNOD .LT. 0 .AND. MEQN(IDOF) .LT. 0 ) THEN

C        ---------------------
C        THIS IS A SLAVE NODE.
C        ---------------------
         ICEQ = -MEQN(IDOF)
         IF ( ICEQ .GT. NCEQ ) GOTO 900

         DO 100 K = MPMCEQ(ICEQ)+1, MPMCEQ(ICEQ+1)-1
            IF ( K .LT. 1 .OR. K .GT. NMMCEQ ) GOTO 900
            IPCEQ = MMCEQ(K)
            IF ( IPCEQ .LT. 1 ) GOTO 100
            IF ( IPCEQ .GT. NDOF ) GOTO 900
            NEWNOD = SPRMAP(IPCEQ)
            IF ( NEWNOD .LT. 1 .OR. NEWNOD .GT. NDOF ) GOTO 900
            IF ( MASTER(NEWNOD) .EQ. IEL ) GOTO 100
            IF ( IPNODS .GT. SPRCON ) GOTO 900

C           -------------------------------
C           THIS NODE IS NOT YET ASSEMBLED.
C           -------------------------------
            MASTER(NEWNOD) = IEL
            ELNOD(IPNODS)  = NEWNOD
            IPNODS = IPNODS + 1
  100    CONTINUE

      ENDIF

      RETURN

  900 CONTINUE
      IERR = -3
      RETURN
      END
