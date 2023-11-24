C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRADM ( EM    , TTCC  ,
     $                    MPAR  , MSPAR ,
     $                    MADOF , MEQN  ,
     $                    MPMNPC, MMNPC ,
     $                    MPMCEQ, MMCEQ ,
     $                    MSICA , MTREES,
     $                    MSIFA , MDOFNC,
     $                    SM    , SV    ,
     $                    IWORK , IEL   , NEDOF , LPU , NSV , IERR )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRADM                  GROUP 9 / PUBLIC
C
C     T A S K :  To assemble an element matrix (EM) into the system
C                matrix (SM) and/or vector (SV).
C
C     MPAR   : INTEGER(NMPAR)
C      Entry : Matrix of parameters.
C      Exit  : Unchanged.
C     MADOF  : INTEGER(NANOD+1)
C      Entry : Matrix of accumulated degrees of freedom.
C      Exit  : Unchanged.
C     MPMNPC : INTEGER(NEL+1)
C      Entry : Matrix of pointers to MNPC.
C      Exit  : Unchanged.
C     MMNPC  : INTEGER(NMMNPC)
C      Entry : Matrix of nodal point correspondence for all elements.
C      Exit  : Unchanged.
C
C     ROUTINES CALLED/REFERENCED :  SPRASM,SPRER1,SPREEQ         (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-02-18 / 1.0
*                       98-08-24 / 1.1   A.C.Damhaug
*                                        MSPAR(50) -> MSPAR(*), and
*                                        added superelement option.
C                       03-03-03 / 2.0   K.M.Okstad
C                                        Added lumped matrix option. The
C                                        diagonal element matrix is then
C                                        added DOF by DOF by calls to
C                                        SPRASM. Added call to ELMEQ2
C                                        and some consistency checking.
C                       03-08-01 / 2.1   K.M.Okstad
C                                        Added call to ELMEQ instead of
C                                        ELMEQ2 when MPAR(18) .EQ. 0.
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF,NSV
      INTEGER           IWORK(*)
      INTEGER           MADOF(*) ,MDOFNC(*),MEQN(*)  ,MMCEQ(*) ,
     $                  MMNPC(*) ,MPAR(50) ,MPMCEQ(*),MPMNPC(*),
     $                  MSICA(*) ,MSIFA(*) ,MSPAR(*) ,MTREES(*)
C
      DOUBLE PRECISION  EM(*)    ,SM(*)    ,SV(*)    ,TTCC(*)
C
C                                                ! local variables
C
      INTEGER           NCEQ,NMMCEQ
      INTEGER           NDOF,NEL,NEQ,NIWORK,NOFSUB,NSUPER,NZEROL
      INTEGER           SPRCON,SPRNOD
      INTEGER           ELNOD,LINDX,NODES,XELNOD,XLINDX,XLNZ
      INTEGER           XNODES,XSUPER
      INTEGER           I,IDOF,IPNOD,NELDOF,NENOD,NESLV,NEPRD,NNDOF
C
C ----------------------------------------------------------------------
C
      IF ( IEL .LT. 1 .OR. IEL .GT. MPAR(2) ) THEN
         IERR = -1
         CALL SPRER1 (22,'SPRADM',IEL,MPAR(2),0,LPU,IERR)
         RETURN
      ELSE IF ( ABS(NSV) .GT. 1 .AND. ABS(NSV) .NE. MPAR(9) ) THEN
         IERR = -2
         CALL SPRER1 (57,'SPRADM',IEL,NSV,MPAR(9),LPU,IERR)
         RETURN
      ELSE IF ( MSPAR(51) .EQ. 1 .AND. MSPAR(52) .EQ. 1 ) THEN
         IERR = -3
         CALL SPRER1 (23,'SPRADM',IEL,MSPAR(51),MSPAR(52),LPU,IERR)
         RETURN
      ENDIF

      IPNOD = MPMNPC(IEL)
      NENOD = MPMNPC(IEL+1) - IPNOD

      IF ( MPAR(18) .EQ. 0 ) THEN
C        Assumes that NNDOF = MADOF(I+I) - MADOF(I) for node I
         NNDOF = MPAR(3)
      ELSE
C        Assumes that NEDOF = NNDOF*NENOD, NNDOF is constant
         NNDOF = NEDOF/NENOD
      ENDIF
      CALL SPREEQ (MADOF,MMNPC(IPNOD),MPMCEQ,MEQN,NENOD,NNDOF,
     $             IWORK,NELDOF,NESLV,NEPRD)
      IF ( NELDOF .NE. NEDOF ) THEN
         IERR = -3
         CALL SPRER1 (22,'SPRADM',IEL,NELDOF,NEDOF,LPU,IERR)
         RETURN
      ENDIF

C                                                ! extract parameters
      NEL    = MPAR(2)
      NDOF   = MPAR(3)
      NCEQ   = MPAR(7)
      NEQ    = MPAR(11)
      NMMCEQ = MPAR(16)
      NIWORK = MSPAR(39)
      NOFSUB = MSPAR(15)
      NSUPER = MSPAR(11)
      NZEROL = MSPAR(16)
      SPRNOD = MSPAR(6)
      SPRCON = MSPAR(7)
C                                                ! extract pointers
      XELNOD = MSPAR(41)
      ELNOD  = MSPAR(42)
      XNODES = MSPAR(43)
      NODES  = MSPAR(44)
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)
C
C                                                ! assemble EM
C
      IF ( MSPAR(51) .EQ. 1 ) THEN
C                                                ! diagonal matrix
         IF ( MPAR(18) .EQ. 0 ) THEN
            IERR = -3
            CALL SPRER1 (21,'SPRADM',IEL,MSPAR(51),MPAR(18),LPU,IERR)
            RETURN
         ENDIF

         DO 100 I = 1, NEDOF
            IDOF = MADOF(MMNPC(IPNOD+(I-1)/NNDOF)) + MOD(I-1,NNDOF)
            CALL SPRASM (EM(I) , TTCC  , IWORK(I),MEQN , MMCEQ , MPMCEQ,
     $                   IDOF  , NCEQ  , NDOF  , 1     , NDOF  , NESLV ,
     $                   NEPRD , NEQ   , NMMCEQ, NSV   , MSICA (XELNOD),
     $                   MSICA (ELNOD) , MSICA (XNODES), MSICA (NODES) ,
     $                   MTREES(XSUPER), MTREES(XLINDX), MSIFA (LINDX) ,
     $                   MSIFA (XLNZ)  , MDOFNC(1)     , MDOFNC(NEQ+1) ,
     $                   NSUPER, SPRNOD, SPRCON, NOFSUB, NZEROL,
     $                   IWORK(1+NEDOF), NIWORK-NEDOF  , SM(NEQ+1), SV ,
     $                   LPU   , IERR    )
  100    CONTINUE

      ELSE

         CALL SPRASM (EM    , TTCC  , IWORK , MEQN  , MMCEQ , MPMCEQ,
     $                IEL   , NCEQ  , NDOF  , NEDOF , NEL   , NESLV ,
     $                NEPRD , NEQ   , NMMCEQ, NSV   , MSICA (XELNOD),
     $                MSICA (ELNOD) , MSICA (XNODES), MSICA (NODES) ,
     $                MTREES(XSUPER), MTREES(XLINDX), MSIFA (LINDX) ,
     $                MSIFA (XLNZ)  , MDOFNC(1)     , MDOFNC(NEQ+1) ,
     $                NSUPER, SPRNOD, SPRCON, NOFSUB, NZEROL,
     $                IWORK(1+NEDOF), NIWORK-NEDOF  , SM(NEQ+1), SV ,
     $                LPU   , IERR    )

      ENDIF
C
      RETURN
      END

      SUBROUTINE SPREEQ (MADOF,MNPC,MPMCEQ,MEQN,NENOD,NNDOF,
     +                   MEEN,NEDOF,NESLV,NEPRD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPREEQ                   GROUP 9 / PRIVATE
C
C     T A S K :  TO FIND  MEEN ('MATRIX OF' ELEMNET EQUATION NUMBERS),
C     NEDOF (NO. OF ELEMENT  D O F S), NESLV (NO. OF ELEMENT SLAVE
C     D O F S) AND NEPRD (NO. OF ELEMENT PRESCRIBED  D O F S) FOR AN
C     ELEMENT WHOSE  MNPC(NENOD)  IS INPUT
C
C     This is the  s p a r s e  version of ELMEQ/ELMEQ2.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-25 / 1.0
C                       02-07-11 / 2.0   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           NEDOF,NENOD,NNDOF,NEPRD,NESLV
      INTEGER           MADOF(*),MEEN(*),MEQN(*),MNPC(NENOD),MPMCEQ(*)
C
      INTEGER           I,ICEQ,IE,INOD,IS,NMST,NN
C ----------------------------------------------------------------------
      NEDOF = 0
      NESLV = 0
      NEPRD = 0
      DO 100 INOD=1,NENOD
         NN = MNPC(INOD)
         IS = MADOF(NN)
         IE = MIN(IS+NNDOF,MADOF(NN+1)) - 1
         DO 50 I=IS,IE
            NEDOF       = NEDOF+1
            MEEN(NEDOF) = MEQN(I)
            IF (MEQN(I).GE.0)  GO TO 50
C
            ICEQ =-MEQN(I)
            NMST = MPMCEQ(ICEQ+1) - MPMCEQ(ICEQ) - 1
C                                               ** PRESCRIBED  D O F
            IF (NMST.EQ.0) NEPRD = NEPRD+1
C                                               ** DEPENDENT  D O F
            IF (NMST.GT.0) NESLV = NESLV+1
   50    CONTINUE
  100 CONTINUE
C
      RETURN
      END


      SUBROUTINE   SPRASM   ( EM    , TTCC  ,
     $                        MEEN  , MEQN  , MMCEQ , MPMCEQ,
     $                        IEL   , NCEQ  , NDOF  , NEDOF , NEL   ,
     $                        NESLV , NEPRD , NEQ   , NMMCEQ, NSV   ,
     $                        XELNOD, ELNOD , XNODES, NODES , XSUPER,
     $                        XLINDX, LINDX , XLNZ  , NODMAP, SUPMAP,
     $                        NSUPER, SPRNOD, SPRCON, NOFSUB, NZEROL,
     $                        IWORK , NIWORK, SM    , SV    ,
     $                        LPU   , IERR                            )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           IEL   , NCEQ  , NDOF  , NEDOF , NEL   , NESLV ,
     $                  NEPRD , NEQ   , NMMCEQ, NSV   ,
     $                  NSUPER, SPRNOD, SPRCON, NOFSUB, NZEROL, NIWORK,
     $                  LPU   , IERR

      INTEGER           MEEN(NEDOF)           , MEQN(NDOF)            ,
     $                  MMCEQ(NMMCEQ)         , MPMCEQ(NCEQ+1)

      INTEGER           XELNOD(NEL+1)         , ELNOD(SPRCON)         ,
     $                  XNODES(SPRNOD+1)      , NODES(NEQ)            ,
     $                  XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  LINDX(NOFSUB)         , XLNZ(NSUPER+1)        ,
     $                  NODMAP(NEQ)           , SUPMAP(NEQ)           ,
     $                  IWORK(NIWORK)

      DOUBLE PRECISION  EM(*)                 , TTCC(*)               ,
     $                  SM(*)                 , SV(*)

C     ------------------------------------------------------------------
C
C     Purpose
C     -------
C
C     SPRASM
C
C --- To assemble an element matrix into the system matrix/vector.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : Mar. 03, 2003 (kmo)
C                 Removed arguments MPAR, MADOF, MPMNPC and MMNPC.
C                 Added new argument MEEN.
C
C
C     EM     : DOUBLE PRECISION(NEDOF,NEDOF)
C      Entry : Element stiffness matrix to be assembled into SM/SV.
C      Exit  : Unchanged.
C     TTCC   : DOUBLE PRECISION(NMMCEQ)
C      Entry : Table of constraint coefficients.
C      Exit  : Unchanged.
C     MEEN   : INTEGER(NEDOF)
C      Entry : Matrix of equation numbers for the element DOFs.
C      Exit  : Unchanged.
C     MEQN   : INTEGER(NDOF)
C      Entry : Matrix of equation numbers as output from SPRTRS.
C      Exit  : Unchanged.
C     MPMCEQ : INTEGER(NCEQ+1)
C      Entry : Matrix of pointers to MMCEQ.
C      Exit  : Unchanged.
C     MMCEQ  : INTEGER(NMMCEQ)
C      Entry : Matrix of constraint equation definition for all the
C              constraint equations.
C      Exit  : Unchanged.
C     IEL    : INTEGER
C      Entry : The current element.
C      Exit  : Unchanged.
C     NEDOF  : INTEGER
C      Entry : Number of element degrees of freedom.
C      Exit  : Unchanged.
C     LPU    : INTEGER
C      Entry : Output unit for error messages.
C      Exit  : Unchanged.
C     NSV    : INTEGER
C      Entry : Specifies by its absolute value the number of system
C              vectors to be assembled.
C      Exit  : Unchanged.
C     SM     : DOUBLE PRECISION(NZEROL)
C      Entry : Contains the current contents of the lower triangular
C              part of a symmetric system matrix stored on compressed
C              sparse form.
C      Exit  : For NSV .GE. 0, SPRASM modifies SM by adding to it
C              appropriate elements of EM. If NSV .LT. 0, SM is dummy.
C     SV     : DOUBLE PRECISION(NEQ*ABS(NSV))
C      Entry : Containts the current contents of a system vector.
C      Exit  : For NSV .NE. 0, SPRASM modifies SV. If NSV .EQ. 0, SV
C              is dummy.
C     XELNOD : INTEGER(NEL+1)
C      Entry : The SPR version of MPMNPC.
C      Exit  : Unchanged.
C     ELNOD  : INTEGER(SPRCON)
C      Entry : The SPR version of MMNPC.
C      Exit  : Unchanged.
C     XNODES : INTEGER(SPRNOD+1)
C      Entry : The SPR version of MADOF.
C      Exit  : Unchanged.
C     NODES  : INTEGER(NEQ)
C      Entry : For each variable in a SPR-node, this array holds the
C              SAM-library defined degree of freedom number.
C      Exit  : Unchanged.
C     XSUPER : INTEGER(NSUPER+1)
C      Entry : The supernode version of MADOF.
C      Exit  : Unchanged.
C     XLINDX : INTEGER(NSUPER+1)
C      Entry : Pointers to the index lists stored in LINDX.
C      Exit  : Unchanged.
C     LINDX  : INTEGER(NOFSUB)
C      Entry : The set of index lists for each supernode.
C      Exit  : Unchanged.
C     XLNZ   : INTEGER(NSUPER+1)
C      Entry : Pointer into SM for the the start of each supernode.
C      Exit  : Unchanged.
C     NODMAP : INTEGER(NEQ)
C      Entry : Maps each variable to its internal SPR-node.
C      Exit  : Unchanged.
C     SUPMAP : INTEGER(NEQ)
C      Entry : Maps each variable to its supernode.
C      Exit  : Unchanged.
C     NDOF   : INTEGER
C      Entry : Number of degrees of freedom.
C      Exit  : Unchanged.
C     NEL    : INTEGER
C      Entry : Number of elements.
C      Exit  : Unchanged.
C     NEQ    : INTEGER
C      Entry : Number of free degrees of freedom.
C      Exit  : Unchanged.
C     NSUPER : INTEGER
C      Entry : Number of supernodes.
C      Exit  : Unchanged.
C     SPRNOD : INTEGER
C      Entry : Number of internal SPR-nodes.
C      Exit  : Unchanged.
C     SPRCON : INTEGER
C      Entry : Length of the SPR connectivity array ELNOD.
C      Exit  : Unchanged.
C     NOFSUB : INTEGER
C      Entry : Number of indices stored to represent SM.
C      Exit  : Unchanged.
C     NZEROL : INTEGER
C      Entry : Length of SM.
C      Exit  : Unchanged.
C     NIWORK : INTEGER
C      Entry : Length of the internal work array IWORK.
C              MAX(MAXDOF,MSPAR(24))+NEQ+MSPAR(24)*(MSPAR(24)+1)/2+
C              MSPAR(25)*(MSPAR(25)+1)/2. MAXDOF is the order of the
C              largest original finite-element matrix.
C      Exit  : Unchanged.
C     IERR   : INTEGER
C      Entry : Need not be set.
C      Exit  : Error flag, if zero a normal return.
C
C
C     Working arrays
C     --------------
C
C     IWORK  : INTEGER(NIWORK)
C      Entry : Not defined.
C      Exit  : Need not be saved.
C
C     Procedures
C     ----------
C
C     SPRT48
C     SPRAD1
C     SPRAD2
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

      INTEGER           ELSIZE, I     , IDOF  , IPNODE, IWUSED, LOCIND,
     $                  NNODEL, NODE  , NUMNOD, PSM   , POINTR
      LOGICAL           LERROR
      EXTERNAL          SPRT48, SPRAD1, SPRAD2, SPRER1

C     -----------
C     INITIALIZE.
C     -----------
      IERR   = 0
      ELSIZE = 0

C     ==================================================================

      IF ( NSV .GE. 0 ) THEN

C        ----------------------------
C        SET SPACE IN IWORK FOR MEEN.
C        ----------------------------
         NUMNOD = XELNOD(IEL+1) - XELNOD(IEL)
         IWUSED = 1 + NUMNOD
         IF ( (IWUSED-1) .GT. NIWORK ) THEN
            IERR = -1
            CALL SPRER1 ( 34, 'SPRASM', IWUSED-1, NIWORK, 0, LPU, IERR )
            RETURN
         ENDIF

C        -----------------------------------------------------
C        BUILD THE GLOBAL INDEX LIST FOR THE EMBEDDED ELEMENT.
C        -----------------------------------------------------
         NNODEL = 0
         DO 100 I = XELNOD(IEL), XELNOD(IEL+1)-1
            NNODEL = NNODEL + 1
            NODE   = ELNOD(I)
            IPNODE = XNODES(NODE)
            IDOF   = NODES(IPNODE)
            IWORK(NNODEL) = MEQN(IDOF)
            IPNODE = XNODES(NODE+1) - IPNODE
            ELSIZE = ELSIZE + IPNODE
  100    CONTINUE

C        ----------------------------------------------
C        SORT THE INDICES IN MEEN INTO ASCENDING ORDER.
C        ----------------------------------------------
         IF ( NNODEL .GT. 1 ) THEN
            CALL SPRT48 ( NNODEL, IWORK, LERROR )
            IF ( LERROR ) THEN
               IERR = -2
               CALL SPRER1 ( 23, 'SPRASM', 0, 0, 0, LPU, IERR )
               RETURN
            ENDIF
         ENDIF

C        ------------------------------------------------
C        PARTITION THE WORK ARRAY FOR POINTR, PSM, LOCIND
C        ------------------------------------------------
         LOCIND = IWUSED
         PSM    = LOCIND + NEQ
         POINTR = PSM    + (ELSIZE*(ELSIZE+1)/2)
         IWUSED = POINTR + (NUMNOD*(NUMNOD+1)/2)
         IF ( (IWUSED-1) .GT. NIWORK ) THEN
            IERR = -1
            CALL SPRER1 ( 34, 'SPRASM', IWUSED-1, NIWORK, 0, LPU, IERR )
            RETURN
         ENDIF

C        -----------------------------------
C        CALL SPRAD1 TO COMPUTE PSM, LOCIND.
C        -----------------------------------
         CALL SPRAD1 ( XNODES, XSUPER, XLINDX, LINDX, XLNZ,
     $                 NODMAP, SUPMAP, NEQ, NSUPER, SPRNOD, NOFSUB,
     $                 NUMNOD, ELSIZE, IWORK(PSM), IWORK(LOCIND),
     $                 IWORK , IWORK(POINTR) )
      ELSE

C        -------------------------------
C        ONLY ASSEMBLE SYSTEM VECTOR SV.
C        -------------------------------
         LOCIND = 1
         PSM    = 1
         IF ( NIWORK .LT. 0 ) THEN
            IERR = -1
            CALL SPRER1 ( 34, 'SPRASM', 0, NIWORK, 0, LPU, IERR )
            RETURN
         ENDIF

      ENDIF

C     -------------------------------------------
C     ASSEMBLE EM INTO SM/SV BY A CALL TO SPRAD2.
C     -------------------------------------------
      CALL SPRAD2 ( EM, MEEN, IWORK(LOCIND), IWORK(PSM),
     $              TTCC, MEQN, MMCEQ, MPMCEQ,
     $              ELSIZE, NCEQ, NDOF, NESLV+NEPRD, NEDOF, NEPRD,
     $              NEQ, NMMCEQ, NSV, NZEROL, SM, SV )

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRASM
      END


      SUBROUTINE   SPRAD1   ( XNODES, XSUPER, XLINDX, LINDX ,
     $                        XLNZ  , NODMAP, SUPMAP, NEQ   , NSUPER,
     $                        SPRNOD, NOFSUB, NUMNOD, ELSIZE, PSM   ,
     $                        LOCIND, VARIND, POINTR                  )

C     ------------------------------------------------------------------

      IMPLICIT          NONE

      INTEGER           NEQ   , NSUPER, SPRNOD, NOFSUB, NUMNOD, ELSIZE
      INTEGER           XNODES(SPRNOD+1)      ,
     $                  XSUPER(NSUPER+1)      , XLINDX(NSUPER+1)      ,
     $                  LINDX(NOFSUB)         , XLNZ(NSUPER+1)
      INTEGER           NODMAP(NEQ)           , SUPMAP(NEQ)           ,
     $                  POINTR(NUMNOD*(NUMNOD+1)/2), LOCIND(NEQ)      ,
     $                  PSM(ELSIZE*(ELSIZE+1)/2), VARIND(NUMNOD)

C     ------------------------------------------------------------------
C
C     Purpose
C     -------
C
C     SPRAD1
C
C --- To assemble an element matrix into the system matrix.
C
C     Created   : Jan. 11, 1994 (acd)
C     Revisions : mnt. xx, 199x (acd)
C
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

      INTEGER           COLLEN, FSTJ  , I     , II    , IEQ   , IJUMP ,
     $                  INOD  , IPLOC , ISKIP1, ISKIP2, ISUP  , ISZE  ,
     $                  IVAR  , IPLNZ , IPLX  , IPJS  , J     , JJ    ,
     $                  JEQ   , JJUMP , JNOD  , JSLEN , JSUP  , JSZE  ,
     $                  JVAR  , K     , KK    , LPOINT, NODE  , NODSZE,
     $                  NXTVAR, OFFSET, POINT , XLPNT

      POINT(I,J,K) = ((J-1)*K) - ((J-1)*(J-2)/2) + I - J + 1

C     ==================================================================

C     ---------------------------------------
C     COMPUTE POINTERS INTO THE SYSTEM MATRIX
C     FOR THE FIRST VARIABLE IN EACH NODE IN
C     THE EMBEDDED MATRIX.
C     ---------------------------------------
      K = NUMNOD
      DO 400 J = 1, NUMNOD
C        -----------------------------------------
C        GET THE FIRST VARIABLE AND ITS SUPERNODE.
C        -----------------------------------------
         JVAR   = VARIND(J)
         JSUP   = SUPMAP(JVAR)
C        ----------------------------------------
C        GET FIRST VARIABLE IN THE SUPERNODE
C        AND THE POINTER TO THE START OF INDICES.
C        ----------------------------------------
         FSTJ = XSUPER(JSUP)
         IPLX = XLINDX(JSUP)
C        --------------------------------
C        GET THE SIZE OF THE SUPERNODE
C        ASSOCIATED WITH JVAR AND COMPUTE
C        THE LENGTH OF THE SUPERNODE.
C        --------------------------------
         ISKIP1 = XSUPER(JSUP+1) - JVAR
         ISKIP2 = XSUPER(JSUP+1) - FSTJ
         COLLEN = XLINDX(JSUP+1) - IPLX
C        ------------------------------------------
C        COMPUTE THE COLUMN IN JSUP ASSOCIATED WITH
C        JVAR AND COMPUTE THE OFFSET INTO SM.
C        ------------------------------------------
         JJUMP  = JVAR  - FSTJ
         OFFSET = (COLLEN*JJUMP) - (JJUMP*(JJUMP-1)/2)
C        -----------------------
C        POINTER FOR DIAG IN SM.
C        -----------------------
         LPOINT = XLNZ(JSUP) + OFFSET
         POINTR(POINT(J,J,K)) = LPOINT
         XLPNT  = LPOINT
         LPOINT = LPOINT + ISKIP1
C        ----------------------------
C        INCREMENT POINTER TO INDICES
C        TO FIRST OFF-SUPER INDEX.
C        ----------------------------
         IPLX = IPLX + ISKIP2
         DO 300 I = J+1, NUMNOD
C           -------------------------------
C           GET INFO FOR THE NEXT VARIABLE.
C           -------------------------------
            IVAR   = VARIND(I)
            ISUP   = SUPMAP(IVAR)
            IF ( ISUP .EQ. JSUP ) THEN
C              ----------------------------------
C              IVAR IS IN SAME SUPERNODE AS JVAR,
C              NO INDEX SEARCH IS NECESSARY.
C              ----------------------------------
               IJUMP  = IVAR  - JVAR
               POINTR(POINT(I,J,K)) = XLPNT + IJUMP
            ELSE
C              -----------------------------------
C              IVAR IS NOT IN SAME SUPERNODE,
C              SEARCH IN LINDX FOR MATCHING INDEX.
C              -----------------------------------
  200          NXTVAR = LINDX(IPLX)
               IF ( NXTVAR .LE. IVAR ) THEN
                  IF ( NXTVAR .EQ. IVAR ) POINTR(POINT(I,J,K)) = LPOINT
                  NODE   = NODMAP(NXTVAR)
                  NODSZE = XNODES(NODE+1) - XNODES(NODE)
                  IPLX = IPLX + NODSZE
                  LPOINT = LPOINT + NODSZE
                  IF ( NXTVAR .LT. IVAR ) GOTO 200
               ENDIF
            ENDIF
  300    CONTINUE
  400 CONTINUE

C     -----------------------------------------------
C     COMPUTE LOCAL INDICES INTO THE EMBEDDED MATRIX.
C     -----------------------------------------------
      IPLOC = 0
      DO 600 I = 1, NUMNOD
         IVAR   = VARIND(I)
         NODE   = NODMAP(IVAR)
         NODSZE = XNODES(NODE+1) - XNODES(NODE)
         DO 500 J = 1, NODSZE
            IPLOC = IPLOC + 1
            LOCIND(IVAR) = IPLOC
            IVAR = IVAR + 1
  500    CONTINUE
  600 CONTINUE

C     -----------------------------------------
C     COMPUTE THE REST OF THE POINTERS INTO SM.
C     -----------------------------------------
      KK = ELSIZE
      DO 1200 J = 1, NUMNOD
C        ----------------------------
C        GET VARIABLE AND BASIC INFO.
C        ----------------------------
         JVAR = VARIND(J)
         JSUP = SUPMAP(JVAR)
         FSTJ = XSUPER(JSUP)
         JJUMP = JVAR - FSTJ
         COLLEN = XLINDX(JSUP+1) - XLINDX(JSUP) - JJUMP
         JNOD = NODMAP(JVAR)
         JSZE = XNODES(JNOD+1) - XNODES(JNOD)
         IPLNZ = POINTR(POINT(J,J,K))
C        --------------------------------
C        POINTERS FOR THE DIAGONAL BLOCK.
C        --------------------------------
         JSLEN = COLLEN
         IPJS = IPLNZ
         DO 800 JJ = JVAR, JVAR+JSZE-1
            JEQ = LOCIND(JJ)
            DO 700 II = JJ, JVAR+JSZE-1
               IEQ = LOCIND(II)
               PSM(POINT(IEQ,JEQ,KK)) = IPJS
               IPJS = IPJS + 1
  700       CONTINUE
            IPLNZ = IPLNZ + JSLEN
            IPJS = IPLNZ
            JSLEN = JSLEN - 1
  800    CONTINUE
C        ------------------------------------
C        POINTERS FOR THE OFF-DIAGONAL BLOCK.
C        ------------------------------------
         COLLEN = COLLEN - 1
         DO 1100 I = J+1, NUMNOD
C           ----------------------------
C           GET VARIABLE AND BASIC INFO.
C           ----------------------------
            IVAR = VARIND(I)
            INOD = NODMAP(IVAR)
            ISZE = XNODES(INOD+1) - XNODES(INOD)
            IPLNZ = POINTR(POINT(I,J,K))
            JSLEN = COLLEN
            IPJS = IPLNZ
            DO 1000 JJ = JVAR, JVAR+JSZE-1
               JEQ = LOCIND(JJ)
               DO 900 II = IVAR, IVAR+ISZE-1
                  IEQ = LOCIND(II)
                  PSM(POINT(IEQ,JEQ,KK)) = IPJS
                  IPJS = IPJS + 1
  900          CONTINUE
               IPLNZ = IPLNZ + JSLEN
               IPJS = IPLNZ
               JSLEN = JSLEN - 1
 1000       CONTINUE
 1100    CONTINUE
 1200 CONTINUE

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRAD1
      END


      SUBROUTINE   SPRAD2   ( EM    , MEEN  , LOCIND, PSM   ,
     $                        TTCC  , MEQN  , MMCEQ , MPMCEQ,
     $                        ELSIZE, NCEQ  , NDOF  , NECEQ , NEDOF ,
     $                        NEPRD , NEQ   , NMMCEQ, NSV   , NZEROL,
     $                        SM    , SV      )

      IMPLICIT NONE

      INTEGER           ELSIZE, NCEQ  , NDOF  , NECEQ , NEDOF , NEPRD ,
     $                  NEQ   , NMMCEQ, NSV   , NZEROL
      INTEGER           MEEN(NEDOF)           , LOCIND(NEQ)           ,
     $                  PSM(ELSIZE*(ELSIZE+1)/2)                      ,
     $                  MEQN(NDOF)            , MMCEQ(NMMCEQ)         ,
     $                  MPMCEQ(NCEQ+1)

      DOUBLE PRECISION  EM(NEDOF,NEDOF)       , TTCC(NMMCEQ)          ,
     $                  SM(NZEROL)            , SV(NEQ)

C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRAD2                   GROUP 7 / PRIVATE
C
C     T A S K :  TO ADD/SUBTRACT THE ELEMENTS CORRESPONDING TO FREE AND
C     CONSTRAINED  D O F S  OF THE ELEMENT MATRIX  EM  TO/FROM THE
C     CURRENT CONTENT OF THE APPROPRIATE ELEMENTS OF
C       MATRIX     SM  IF  NSV.GE.0  (ADDITION WITHOUT/WITH WEIGHTING)
C       VECTOR(S)  SV  IF  NSV.NE.0  (SUBTRACTION WITH WEIGHTING)
C     SM IS A SYMMETRIC SYSTEM MATRIX STORED IN 'COMPRESSED SPARSE' FORM
C     AND SV IS A SYSTEM VECTOR (ABS(NSV)=1) OR A SET OF SYSTEM VECTORS
C     (ABS(NSV)=NPDOF) CONTAINING THE CONTRIBUTIONS TO THE RIGHT-HAND
C     SIDE(S) FROM PRESCRIBED AND (POSSIBLY) DEPENDENT D O F S.
C
C     ROUTINES CALLED/REFERENCED :  ABS    (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-01-14 / 1.0
C                       94-01-11 / 2.0   A.C.Damhaug
C                                  Changed ADDEM to SPRAD2 and adapted
C                                  to new definition of storage for SM.
C                       02-07-20 / 2.1   K.M.Okstad
C                                  Replaced ELMEQ by ELMEQ2 to allow
C                                  elements having fewer nodal DOFs than
C                                  given by the system array MADOF.
C                       02-09-06 / 2.2   K.M.Okstad
C                                  Removed the local arrays MMST[IJ] and
C                                  call to SLVEQ, such that there is no
C                                  limit on the number of master DOFs
C                                  for a constraint equation.
C                       03-03-03 / 3.0   K.M.Okstad
C                                  Moved call to ELMEQ2 and some error
C                                  checking to the calling routine.
C                                  Removed MPAR,MADOF,MPMNPC and MMNPC.
C                                  Added MEEN as an argument instead.
C                       05-10-10 / 3.1   K.M.Okstad
C                                  Check for MMCEQ.GT.0 and TTCC.NE.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************

      INTEGER           I     , ICEQ  , IEI   , IEQ   , II    , IP    ,
     $                  J     , JCEQ  , JEJ   , JEQ   , JJ    , JP    ,
     $                  K     , KEQ   , KP    , NMSTI , NMSTJ , PPSM
      DOUBLE PRECISION  C0,ZERO
      PARAMETER       ( ZERO = 0.0D0 )
      INTRINSIC         ABS

      PPSM(I,J,K) = ((J-1)*K) - ((J-1)*(J-2)/2) + I - J + 1

C ----------------------------------------------------------------------

      IF (NSV.LT.0)            GO TO 220

C ----------------------------------------------------------------------
C  ADD ELEMENTS CORRESPONDING TO FREE  D O F S  IN  EM  INTO  SM
C ----------------------------------------------------------------------
      DO 200 J=1,NEDOF
         JEQ = MEEN(J)
         IF (JEQ.LT.1)         GO TO 200
         DO 150 I=1,J
            IEQ = MEEN(I)
            IF (IEQ.LT.1)      GO TO 150

            IF     ( IEQ.LT.JEQ ) THEN
C              -----------------
C              SWAP THE INDICES.
C              -----------------
               IEI = LOCIND(JEQ)
               JEJ = LOCIND(IEQ)

            ELSEIF ( IEQ.GT.JEQ ) THEN
C              -----------------
C              KEEP THE INDICES.
C              -----------------
               IEI = LOCIND(IEQ)
               JEJ = LOCIND(JEQ)

            ELSE
C              ---------------
C              DIAGONAL ENTRY.
C              ---------------
               IEI = LOCIND(IEQ)
               JEJ = IEI

            ENDIF

C           ------------------------
C           FIND THE LOCATION IN SM.
C           ------------------------
            K = PSM(PPSM(IEI,JEJ,ELSIZE))

C           --------------------------
C           ASSEMBLE THE CONTRIBUTION.
C           --------------------------
            SM(K) = SM(K) + EM(I,J)

  150    CONTINUE
  200 CONTINUE

  220 IF (NECEQ.EQ.0)                  GO TO 1000
      IF (NSV.LT.(-1).AND.NEPRD.EQ.0)  GO TO 1000

C ----------------------------------------------------------------------
C  ADD (APPROPRIATELY WEIGHTED) ELEMENTS CORRESPONDING TO CONSTRAINED
C  (DEPENDENT AND PRESCRIBED)  D O F S  IN  EM  INTO  SM  AND/OR  SV
C ----------------------------------------------------------------------
      DO 600 J=1,NEDOF
         IF (MEEN(J).GE.0)     GO TO 600
         JCEQ =-MEEN(J)
         JP   = MPMCEQ(JCEQ)
         IF (NSV.EQ.0)         GO TO 320
         C0   = TTCC(JP)
         IF (C0.EQ.ZERO)       GO TO 310
         KP   = 0
         IF (ABS(NSV).EQ.1)    GO TO 250
         K    = 0
         IF (MPMCEQ(JCEQ+1)-MPMCEQ(JCEQ).GT.1) GO TO 310
         DO 240 I=1,JCEQ
            IF (MPMCEQ(I+1)-MPMCEQ(I).EQ.1) K = K + 1
  240    CONTINUE
         IF (K.EQ.0)           GO TO 310
         KP   = (K-1)*NEQ

C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SV  (R-H SIDE)
  250    DO 300 I=1,NEDOF
            IEQ = MEEN(I)
            IF (IEQ.EQ.0)      GO TO 300
            IF (IEQ.LT.0)      GO TO 260

            KEQ = KP+IEQ
            SV(KEQ) = SV(KEQ) - C0*EM(I,J)
            GO TO 300

  260       ICEQ =-IEQ
            NMSTI = MPMCEQ(ICEQ+1) - MPMCEQ(ICEQ) - 1
            IF (NMSTI.LE.0)    GO TO 300
            IP = MPMCEQ(ICEQ)
            DO 280 II=1,NMSTI
               IP  = IP+1
               IF (MMCEQ(IP).LE.0 .OR. TTCC(IP).EQ.ZERO) GO TO 280
               KEQ = KP+MEQN(MMCEQ(IP))
               SV(KEQ) = SV(KEQ) - C0*TTCC(IP)*EM(I,J)
  280       CONTINUE

  300    CONTINUE

  310    IF (NSV.LT.0)         GO TO 600
  320    NMSTJ = MPMCEQ(JCEQ+1) - MPMCEQ(JCEQ) - 1
         IF (NMSTJ.LE.0)       GO TO 600

C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SM
         DO 500 JJ=1,NMSTJ
            JP  = JP+1
            IF (MMCEQ(JP).LE.0 .OR. TTCC(JP).EQ.ZERO) GO TO 500
            JEQ = MEQN(MMCEQ(JP))
            DO 400 I=1,NEDOF
               IEQ = MEEN(I)
               IF (IEQ.EQ.0)   GO TO 400
               IF (IEQ.LT.0)   GO TO 360
               IF (IEQ.NE.JEQ) THEN

                  IF (IEQ.LT.JEQ) THEN
C                    -----------------
C                    SWAP THE INDICES.
C                    -----------------
                     IEI = LOCIND(JEQ)
                     JEJ = LOCIND(IEQ)

                  ELSE
C                    -----------------
C                    KEEP THE INDICES.
C                    -----------------
                     IEI = LOCIND(IEQ)
                     JEJ = LOCIND(JEQ)

                  ENDIF

C                 ------------------------
C                 FIND THE LOCATION IN SM.
C                 ------------------------
                  K = PSM(PPSM(IEI,JEJ,ELSIZE))

C                 --------------------------
C                 ASSEMBLE THE CONTRIBUTION.
C                 --------------------------
                  SM(K) = SM(K) + TTCC(JP)*EM(I,J)

               ELSE
C                 ---------------
C                 DIAGONAL ENTRY.
C                 ---------------
                  IEI = LOCIND(IEQ)
                  JEJ = IEI

C                 ------------------------
C                 FIND THE LOCATION IN SM.
C                 ------------------------
                  K = PSM(PPSM(IEI,JEJ,ELSIZE))

C                 --------------------------
C                 ASSEMBLE THE CONTRIBUTION.
C                 --------------------------
                  SM(K) = SM(K) + (TTCC(JP)+TTCC(JP))*EM(I,J)

               ENDIF

               GOTO 400

  360          ICEQ =-IEQ
               NMSTI = MPMCEQ(ICEQ+1) - MPMCEQ(ICEQ) - 1
               IF (NMSTI.LE.0) GO TO 400

               IP    = MPMCEQ(ICEQ)
               DO 380 II=1,NMSTI
                  IP  = IP+1
                  IF (MMCEQ(IP).LE.0 .OR. TTCC(IP).EQ.ZERO) GO TO 380
                  IEQ = MEQN(MMCEQ(IP))
                  IF (IEQ.LE.JEQ) THEN

C                    -----------------
C                    SWAP THE INDICES.
C                    -----------------
                     IEI = LOCIND(JEQ)
                     JEJ = LOCIND(IEQ)

C                    ------------------------
C                    FIND THE LOCATION IN SM.
C                    ------------------------
                     K = PSM(PPSM(IEI,JEJ,ELSIZE))

C                    --------------------------
C                    ASSEMBLE THE CONTRIBUTION.
C                    --------------------------
                     SM(K) = SM(K) + TTCC(IP)*TTCC(JP)*EM(I,J)
                  ENDIF
  380          CONTINUE
  400       CONTINUE
  500    CONTINUE
  600 CONTINUE
C ----------------------------------------------------------------------
 1000 RETURN
      END
