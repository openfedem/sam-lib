C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE ELWNUM (MPMNPC,MMNPC,MINEX,
     +                   MSC,MADOF,
     +                   MTMP,MNNN,
     +                   NANOD,NEL,NWTMP,IPSW,LPU,IFLAG,
     +                   NSTART,NEFF,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ELWNUM                  GROUP 9 / PUBLIC
C
C     TASK :  Renumber a set of NANOD active nodal points interconnected
C             by NEL elements such as to make the profile length of the
C             corresponding stiffness matrix as small as possible.
C             A modified version of an algorithm due E.L.Wilson (used
C             in EDSAP) and similar to reverse Cuthill-McKee is used.
C             As an option (IFLAG=1), retained nodes in substructure
C             analysis (or static condensation) may be accounted for.
C             These nodes, identified by having at least one 2-status
C             dof and no 1-status dofs recorded in MSC/MADOF, are in-
C             cluded in the renumbering, but they are "moved down" such
C             as to receive numbers after the internal (local) nodes,
C             in a sequence implied by the old (input) numbering.
C             The effectiveness of the algorithm is dependent on the
C             starting node.  This node may be supplied by the user
C             (in argument NSTART), or it may be determined auto-
C             matically (if NSTART=0 on input).
C             The effectiveness of the renumbering is returned in
C             (integer) argument NEFF which gives the length of the new
C             (renumbered) profile in per cent of the "old" profile
C             (corresponding to the input node numbering).
C
C
C     ROUTINES CALLED/REFERENCED :   NNNINT, NODELC, NUMERR  (SAM-9)
C                                    PFLNGT, ELWALG          (SAM-9)
C                                    PRITAB                  (SAM-8)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-25 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,IFLAG,IPSW,LPU,NANOD,NEFF,NEL,NSTART,NWTMP,
     +                  MADOF(*),MINEX(*),MMNPC(*),MNNN(*),MPMNPC(*),
     +                  MSC(*),MTMP(*)
C                                                ! local variables
C
      INTEGER           IEL,INOD,IPCE,IPED,IPMP,IPND,IPNN,L,LE,LS,
     +                  NBW,NHIGH,NPLNEW,NPLOLD,NW
C ----------------------------------------------------------------------
      IERR = 0
C
C --- check integrety of input information
C
      DO 20 IEL=1,NEL
         LS = MPMNPC(IEL)
         LE = MPMNPC(IEL+1) - 1
         DO 10 L=LS,LE
            IF (MMNPC(L) .GT. NANOD) THEN
               CALL NUMERR (1,IEL,MMNPC(L),NANOD,LPU,IERR)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      IF (IERR .LT. 0)  GO TO 100
C
      IF (NSTART .GT. 0) THEN
C                                                ! start node
         DO 30 INOD=1,NANOD
            IF (MINEX(INOD) .EQ. NSTART) THEN
               NSTART = INOD
               GO TO 40
            ENDIF
   30    CONTINUE
         NSTART = 0
      ENDIF
C
C --- set pointers in MTMP and establish MPMNCE / MMNCE
C
   40 IPED = 1
      IPND = IPED+NEL
      IPNN = IPND+NANOD
      IPMP = IPNN+NANOD
      IPCE = IPMP+NANOD+1
      NW   = NWTMP-IPCE+1
C
      IF (NW .LT. 0) THEN
         CALL NUMERR (2,NW,NW,NW,LPU,IERR)
         GO TO 100
      ENDIF
C
      CALL NODELC (MINEX,MPMNPC,MMNPC,MTMP(IPND),MTMP(IPMP),
     +             MTMP(IPCE),NANOD,NEL,LPU,NW,IERR)
      IF (IERR .LT. 0)  GO TO 100
C ------------------------------------------------ print
      IF (IPSW.GT.1 .AND. LPU.GT.0) THEN
         NW = MTMP(IPMP+NANOD) - 1
         WRITE (LPU,630) NW
      ENDIF
C ------------------------------------------------
C
C --- initiate MNNN and determine "old" profile and bandwidth
C
      CALL NNNINT (MSC,MADOF,MNNN,NANOD,IFLAG)
      CALL PFLNGT (MPMNPC,MMNPC,MNNN,MTMP(IPND),
     +             NANOD,NEL,NPLOLD,NBW)
      DO 50 INOD=1,NANOD
         IF (MNNN(INOD) .GT. 0)  MNNN(INOD) = 0
   50 CONTINUE
C ------------------------------------------------ print
      IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
         WRITE (LPU,610) NPLOLD,NBW
      ENDIF
C ------------------------------------------------
C
C --- renumber
C
      CALL ELWALG (MPMNPC,MMNPC,MTMP(IPMP),MTMP(IPCE),MNNN,
     +             MTMP(IPND),MTMP(IPED),MTMP(IPNN),
     +             NANOD,NEL,IPSW,LPU,NSTART,NHIGH,IERR)
      IF (IERR .LT. 0)  GO TO 100
      IF (IFLAG.NE.1 .AND. NHIGH.NE.NANOD) THEN
         CALL NUMERR (4,NHIGH,NANOD,NANOD,LPU,IERR)
         GO TO 100
      ENDIF
C
C --- determine "new" profile and bandwidth and compute "effectiveness"
C     of renumbering
C
      CALL PFLNGT (MPMNPC,MMNPC,MNNN,MTMP(IPND),
     +             NANOD,NEL,NPLNEW,NBW)
      NEFF   = (100*NPLNEW)/NPLOLD
      NSTART = MINEX(NSTART)
C ------------------------------------------------ print
      IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
         WRITE (LPU,620) NPLNEW,NBW,NSTART
      ENDIF
C ------------------------------------------------
C
      IF (NHIGH .LT. NANOD) THEN
C                                                ! number retained nodes
         L = NHIGH
         DO 60 INOD=1,NANOD
            IF (MNNN(INOD) .LT. 0) THEN
               L = L+1
               MNNN(INOD) = L
            ENDIF
   60    CONTINUE
      ENDIF
C ------------------------------------------------ print
      IF (IPSW.GT.2 .AND. LPU.GT.0) THEN
         WRITE (LPU,640)
         CALL PRITAB (MNNN,NANOD,LPU)
      ENDIF
C ------------------------------------------------
C
  100 RETURN
C ------------------------------------------------ formats
  610 FORMAT (///5X,'"OLD" PROFILE-LENGTH =',I8,
     +          /5X,'"OLD" BANDWIDTH      =',I8 )
  620 FORMAT (///5X,'"NEW" PROFILE-LENGTH =',I8,
     +          /5X,'"NEW" BANDWIDTH      =',I8,
     +          /5X,'START NODE NUMBER    =',I8 )
  630 FORMAT (///5X,'Size of table MMNCE  =',I8 )
  640 FORMAT ('1',4X,'Table of "new" node numbers:' / )
C
      END
      SUBROUTINE PFMNUM (MPMNPC,MMNPC,MINEX,
     +                   MSC,MADOF,
     +                   MTMP,MNNN,
     +                   NANOD,NEL,NWTMP,IPSW,LPU,IFLAG,
     +                   NSTART,NEFF,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PFMNUM                  GROUP 9 / PUBLIC
C
C     TASK :  Renumber a set of NANOD active nodal points interconnected
C             by NEL elements such as to make the profile length of the
C             corresponding stiffness matrix as small as possible.
C             A modified version of the PFM algorithm due to Hoit and
C             Wilson is used (modification due to Hoit and Garcelon).
C             As an option (IFLAG=1) retained nodes in substructure
C             analysis (or static condensation) may be accounted for.
C             These nodes, identified by having at least one 2-status
C             dof and no 1-status dofs recorded in MSC/MADOF, are not
C             included in the renumbering.  They are numbered after the
C             internal (local) nodes are numbered, and in the sequence
C             implied by the old (input) numbering.
C             The effectiveness of the algorithm is dependent on the
C             starting node.  This node may be supplied by the user
C             (in argument NSTART), or it may be determined auto-
C             matically (if NSTART=0 on input).
C             The effectiveness of the renumbering is returned in
C             (integer) argument NEFF which gives the length of the new
C             (renumbered) profile in per cent of the "old" profile
C             (corresponding to the input node numbering).
C
C
C     ROUTINES CALLED/REFERENCED :   NNNINT, NODELC, NUMERR  (SAM-9)
C                                    PFLNGT, PFMALG          (SAM-9)
C                                    PRITAB                  (SAM-8)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-23 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IFLAG,IPSW,LPU,NANOD,NEFF,NEL,NSTART,NWTMP,
     +                  MADOF(*),MINEX(*),MMNPC(*),MNNN(*),MPMNPC(*),
     +                  MSC(*),MTMP(*)
C                                                ! local variables
C
      INTEGER           IEL,INOD,IPCE,IPED,IPMP,IPND,IPNF,L,LE,LS,
     +                  NBW,NHIGH,NPLNEW,NPLOLD,NW
C ----------------------------------------------------------------------
      IERR = 0
C
C --- check integrety of input information
C
      DO 20 IEL=1,NEL
         LS = MPMNPC(IEL)
         LE = MPMNPC(IEL+1) - 1
         DO 10 L=LS,LE
            IF (MMNPC(L) .GT. NANOD) THEN
               CALL NUMERR (1,IEL,MMNPC(L),NANOD,LPU,IERR)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      IF (IERR .LT. 0)  GO TO 100
C
      IF (NSTART .GT. 0) THEN
C                                                ! start node
         DO 30 INOD=1,NANOD
            IF (MINEX(INOD) .EQ. NSTART) THEN
               NSTART = INOD
               GO TO 40
            ENDIF
   30    CONTINUE
         NSTART = 0
      ENDIF
C
C --- set pointers in MTMP and establish MPMNCE / MMNCE
C
   40 IPED = 1
      IPND = IPED+NEL
      IPNF = IPND+NANOD
      IPMP = IPNF+NANOD
      IPCE = IPMP+NANOD+1
      NW   = NWTMP-IPCE+1
C
      IF (NW .LT. 0) THEN
         CALL NUMERR (2,NW,NW,NW,LPU,IERR)
         GO TO 100
      ENDIF
C
      CALL NODELC (MINEX,MPMNPC,MMNPC,MTMP(IPND),MTMP(IPMP),
     +             MTMP(IPCE),NANOD,NEL,LPU,NW,IERR)
      IF (IERR .LT. 0)  GO TO 100
C ------------------------------------------------ print
      IF (IPSW.GT.1 .AND. LPU.GT.0) THEN
         NW = MTMP(IPMP+NANOD) - 1
         WRITE (LPU,630) NW
      ENDIF
C ------------------------------------------------
C
C --- initiate MNNN and determine "old" profile and bandwidth
C
      CALL NNNINT (MSC,MADOF,MNNN,NANOD,IFLAG)
      CALL PFLNGT (MPMNPC,MMNPC,MNNN,MTMP(IPND),
     +             NANOD,NEL,NPLOLD,NBW)
      DO 50 INOD=1,NANOD
         IF (MNNN(INOD) .GT. 0)  MNNN(INOD) = 0
   50 CONTINUE
C ------------------------------------------------ print
      IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
         WRITE (LPU,610) NPLOLD,NBW
      ENDIF
C ------------------------------------------------
C
C --- renumber
C
      CALL PFMALG (MPMNPC,MMNPC,MTMP(IPMP),MTMP(IPCE),MNNN,
     +             MTMP(IPND),MTMP(IPED),MTMP(IPNF),
     +             NANOD,NEL,IPSW,LPU,NSTART,NHIGH,IERR)
      IF (IERR .LT. 0)  GO TO 100
      IF (IFLAG.NE.1 .AND. NHIGH.NE.NANOD) THEN
         CALL NUMERR (4,NHIGH,NANOD,NANOD,LPU,IERR)
         GO TO 100
      ENDIF
C
C --- determine "new" profile and bandwidth and compute "effectiveness"
C     of renumbering
C
      CALL PFLNGT (MPMNPC,MMNPC,MNNN,MTMP(IPND),
     +             NANOD,NEL,NPLNEW,NBW)
      NEFF   = (100*NPLNEW)/NPLOLD
      NSTART = MINEX(NSTART)
C ------------------------------------------------ print
      IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
         WRITE (LPU,620) NPLNEW,NBW,NSTART
      ENDIF
C ------------------------------------------------
C
      IF (NHIGH .LT. NANOD) THEN
C                                                ! number retained nodes
         L = NHIGH
         DO 60 INOD=1,NANOD
            IF (MNNN(INOD) .LT. 0) THEN
               L = L+1
               MNNN(INOD) = L
            ENDIF
   60    CONTINUE
      ENDIF
C ------------------------------------------------ print
      IF (IPSW.GT.2 .AND. LPU.GT.0) THEN
         WRITE (LPU,640)
         CALL PRITAB (MNNN,NANOD,LPU)
      ENDIF
C ------------------------------------------------
C
  100 RETURN
C ------------------------------------------------ formats
  610 FORMAT (///5X,'"OLD" PROFILE-LENGTH =',I8,
     +          /5X,'"OLD" BANDWIDTH      =',I8 )
  620 FORMAT (///5X,'"NEW" PROFILE-LENGTH =',I8,
     +          /5X,'"NEW" BANDWIDTH      =',I8,
     +          /5X,'START NODE NUMBER    =',I8 )
  630 FORMAT (///5X,'Size of table MMNCE  =',I8 )
  640 FORMAT ('1',4X,'Table of "new" node numbers:' / )
C
      END
      SUBROUTINE SWSNUM (MPMNPC,MMNPC,MINEX,
     +                   MSC,MADOF,
     +                   MTMP,MNNN,
     +                   NANOD,NEL,NWTMP,IPSW,LPU,IFLAG,
     +                   NSTART,NEFF,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SWSNUM                  GROUP 9 / PUBLIC
C
C     TASK :  Renumber a set of NANOD active nodal points interconnected
C             by NEL elements such as to make the profile length of the
C             corresponding stiffness matrix as small as possible.
C             An algorithm due to S.W.Sloan is used, and the code is
C             a slightly rewritten version of Sloan's code (so as to fit
C             the SAM data structures).
C             As an option (IFLAG=1) retained nodes in substructure
C             analysis (or static condensation) may be accounted for.
C             These nodes, identified by having at least one 2-status
C             dof and no 1-status dofs recorded in MSC/MADOF, are not
C             included in the renumbering.  They are numbered after the
C             internal (local) nodes are numbered, and in the sequence
C             implied by the old (input) numbering.
C             The initial start node may be supplied by the user (in
C             argument NSTART), or it may be determined automatically
C             (if NSTART=0 on input).  The start node actually used for
C             the renumbering is obtained through an iterative process;
C             its "value" (in terms of user specified numbers) is re-
C             turned in argument NSTART.
C             The effectiveness of the renumbering is returned in
C             (integer) argument NEFF which gives the length of the new
C             (renumbered) profile in per cent of the "old" profile
C             (corresponding to the input node numbering).
C
C
C     ROUTINES CALLED/REFERENCED :   INNNI,  AGRAPH, NUMERR  (SAM-9)
C                                    PPNODE, RLVSTR, INSORT  (SAM-9)
C                                    SWSALG, PFRLL           (SAM-9)
C                                    PRITAB                  (SAM-8)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-11-18 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,IFLAG,IPSW,LPU,NANOD,NEFF,NEL,NSTART,NWTMP,
     +                  MADOF(*),MINEX(*),MMNPC(*),MNNN(*),MPMNPC(*),
     +                  MSC(*),MTMP(*)
C                                                ! local variables
C
      INTEGER           IEL,INOD,IP1,IP2,IP3,IP4,IP5,L,LE,LS,LSTNUM,
     +                  NBWNEW,NBWOLD,NEND,NHIGH,NMMAN,NNC,NPLNEW,
     +                  NPLOLD,NW
C ----------------------------------------------------------------------
      IERR = 0
C
C --- check integrety of input information
C
      DO 20 IEL=1,NEL
         LS = MPMNPC(IEL)
         LE = MPMNPC(IEL+1) - 1
         DO 10 L=LS,LE
            IF (MMNPC(L) .GT. NANOD) THEN
               CALL NUMERR (1,IEL,MMNPC(L),NANOD,LPU,IERR)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      IF (IERR .LT. 0) GO TO 100
C
C --- determine adjacency list (MMAN) and associated pointer array
C     (MPMAN) for the graph defined by MMNPC and MPMNPC
C     MPMAN is stored in MTMP, from index IP1
C     MMAN  is stored in MTMP, from index IP2
C
      IP1 = 1
      IP2 = IP1 + NANOD + 1
      NW  = NWTMP - IP2 + 1
      CALL AGRAPH (MPMNPC,MMNPC,MINEX,MTMP(IP1),MTMP(IP2),
     +             NANOD,NEL,LPU,NW,IERR)
      IF (IERR .LT. 0)  GO TO 100
      NMMAN = MTMP(NANOD+1) - 1
C ----------------------------------------------- print
      IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
         WRITE (LPU,600)
         WRITE (LPU,610) NW,NMMAN
      ENDIF
C -----------------------------------------------
C
C --- allocate remaining local storage
C
      IP3 = IP2 + NMMAN
      IP4 = IP3 + NANOD + 1
      IP5 = IP4 + NANOD
      NW  = IP5 + NANOD - 1
      IF (NW .GT. NWTMP) THEN
         NW = NW - NWTMP
         CALL NUMERR (3,NWTMP,NW,NW,LPU,IERR)
         GO TO 100
      ENDIF
C
C --- initiate MNNN and determine the number of nodes to be
C     renumbered (NHIGH)
C
      CALL INNNI (MSC,MADOF,MNNN,NANOD,IFLAG,NHIGH)
C
      LSTNUM = 0
      IF (NSTART .GT. 0) THEN
         DO 30 INOD=1,NANOD
            IF (MINEX(INOD) .EQ. NSTART) THEN
               NSTART = INOD
               GO TO 50
            ENDIF
   30    CONTINUE
         NSTART = 0
      ENDIF
C
C === while some nodes remain unnumbered do ============================
C
   50 IF (LSTNUM .LT. NHIGH) THEN
         IF (LSTNUM .GT. 0)  NSTART=0
C
C ------ find start and end node in current component of graph, and
C        compute distance of nodes from end node
C
         CALL PPNODE (MTMP(IP1),MTMP(IP2),MNNN,
     +                MTMP(IP3),MTMP(IP4),MTMP(IP5),
     +                NANOD,NSTART,NEND,NNC)
C
C ----------------------------------------------- print
         IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
            WRITE (LPU,620) MINEX(NSTART),MINEX(NEND)
         ENDIF
C -----------------------------------------------
C
C ------ number nodes in this component of the graph
C
         CALL SWSALG (MTMP(IP1),MTMP(IP2),MNNN,MTMP(IP4),
     +                MTMP(IP5),NNC,NSTART,LSTNUM)
         GO TO 50
      ENDIF
C ======================================================================
C
C --- compute new and old profile lengths and bandwidths, and effective-
C     ness of renumbering
C
      CALL PRFLL (MTMP(IP1),MTMP(IP2),MNNN,NANOD,
     +            NPLOLD,NPLNEW,NBWOLD,NBWNEW)
      NEFF   = (100*NPLNEW)/NPLOLD
      NSTART = MINEX(NSTART)
C ------------------------------------------------ print
      IF (IPSW.GT.0 .AND. LPU.GT.0) THEN
         WRITE (LPU,630) NPLOLD,NBWOLD,NPLNEW,NBWNEW
      ENDIF
C ------------------------------------------------
  100 RETURN
C ------------------------------------------------ formats
C
  600 FORMAT (///5X,'LOCAL PRINT FROM RENUMBERING ROUTINE SWSNUM:')
  610 FORMAT (///5X,'Initial size of adjacency array =',I8,
     +          /5X,'Final size of adjacency array   =',I8 )
  620 FORMAT (///5X,'Start node for current comp. of graph =',I8,
     +          /5X,'End node for current comp. of graph   =',I8 )
  630 FORMAT (///5X,'Profile length for "old" numbering =',I8,
     +          /5X,'Bandwidth for "old" numbering      =',I8,
     +          /5X,'Profile length for "new" numbering =',I8,
     +          /5X,'Bandwidth for "new" numbering      =',I8 )
C
      END
      SUBROUTINE RESHF1 (MNNN,X,AUX,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RESHF1                  GROUP 9 / PUBLIC
C
C     TASK :  Rearrange the N elements in X according to the
C             "renumbering" array MNNN by means of auxiliary array AUX
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S=Single precision      D=Double precision
C **********************************************************************
C     IMPLICIT          NONE
      INTEGER           N,     MNNN(N)
      DOUBLE PRECISION  AUX(N),X(N)
C                                                ! local variables
      INTEGER           I,J
C ----------------------------------------------------------------------
      DO 20 I=1,N
         J      = MNNN(I)
         AUX(J) = X(I)
   20 CONTINUE
C
      DO 40 I=1,N
         X(I) = AUX(I)
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE RESHF2 (MNNN,MMNPC,NMMNPC)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RESHF2                  GROUP 9 / PUBLIC
C
C     TASK :  Rearrange the content of the connectivity array MMNPC
C             according to "new" node numbering - recorded in MNNN.
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-05-13 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NMMNPC, MMNPC(*),MNNN(*)
C                                                ! local variables
      INTEGER           I,IOLD
C ----------------------------------------------------------------------
      DO 10 I=1,NMMNPC
         IOLD     = MMNPC(I)
         MMNPC(I) = MNNN(IOLD)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE RESHF3 (MNNN,MADOF,MSC,MTMP1,MTMP2,NANOD,NDOF)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RESHF3                  GROUP 9 / PUBLIC
C
C     TASK :  Rearrange the content of "control matrices"
C
C                 MADOF and MSC
C
C             due to "new" node numbering - recorded in MNNN.
C             The "reshuffling" is accomplished by means of two
C             temporary (scratch) arrays of at least NANOD and NDOF
C             locations, respectively.
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-10-17 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NANOD,NDOF,
     +                  MADOF(*),MNNN(*),MSC(*),MTMP1(*),MTMP2(*)
C
C                                                ! local variables
C
      INTEGER           IDOF,INEW,INOD,IOLD,IP,L,LE,LS,N1,N2
C ----------------------------------------------------------------------
C
C --- "New" (temporary) array MADOF, containing only the first NANOD
C     entries, stored in MTMP1
C
      DO 5 IOLD=1,NANOD
         INEW = MNNN(IOLD)
         MTMP1(INEW) = MADOF(IOLD+1) - MADOF(IOLD)
    5 CONTINUE
      IF (NANOD .GT. 1) THEN
         N1       = MTMP1(1)
         MTMP1(1) = 1
         DO 10 INEW=2,NANOD
            N2          = MTMP1(INEW)
            MTMP1(INEW) = MTMP1(INEW-1) + N1
            N1          = N2
   10    CONTINUE
      ENDIF
C                                                ! MSC - first in MTMP2
      DO 20 IOLD=1,NANOD
         INEW = MNNN(IOLD)
         LS   = MADOF(IOLD)
         LE   = MADOF(IOLD+1) - 1
         IP   = MTMP1(INEW)
         DO 15 L=LS,LE
            MTMP2(IP) = MSC(L)
            IP        = IP+1
   15    CONTINUE
   20 CONTINUE
      DO 30 IDOF=1,NDOF
         MSC(IDOF) = MTMP2(IDOF)
   30 CONTINUE
C                                                ! MADOF (from MTMP1)
      DO 40 INOD=1,NANOD
         MADOF(INOD) = MTMP1(INOD)
   40 CONTINUE
C
      RETURN
      END
      SUBROUTINE RESHF4 (MNNN,MINEX,MTMP,NANOD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RESHF4                  GROUP 9 / PUBLIC
C
C     TASK :  Rearrange the content of "control matrix" MINEX,
C             due to "new" node numbering - recorded in MNNN.
C             The "reshuffling" is accomplished by means of a
C             temporary (scratch) array of at least NANOD locations.
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-10-17 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NANOD,    MINEX(*),MNNN(*),MTMP(*)
C
C                                                ! local variables
      INTEGER           INEW,INOD,IOLD
C ----------------------------------------------------------------------
      DO 10 IOLD=1,NANOD
         INEW       = MNNN(IOLD)
         MTMP(INEW) = MINEX(IOLD)
   10 CONTINUE
      DO 20 INOD=1,NANOD
         MINEX(INOD) = MTMP(INOD)
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE AGRAPH (MPMNPC,MMNPC,MINEX,MPMAN,MMAN,
     +                   NANOD,NEL,LPU,NMX,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  AGRAPH                 GROUP 9 / PRIVATE
C
C     TASK :  To form the adjacency list for a graph corresponding to
C             a finite element mesh defined by connectivity information
C             stored in MPMNPC and MMNPC.
C             NMX contains on entry the maximum number of storage
C             locations available for the adjacency list MMAN, and on
C             exit it contains the initial size of MMAN.
C             The final size of MMAN is the exit value of MPMAN(NANOD+1)
C             minus 1.
C             The difference in initial and final size corresponds to
C             the necessary "elbow room" required by the assembling
C             algorithm.
C
C     This is a rewritten version of Sloan's subroutine GRAPH
C
C
C     ROUTINES CALLED/REFERENCED :  NUMERR     (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell   (after Scott W. Sloan)
C     DATE/VERSION  :   90-11-07 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LPU,NANOD,NEL,NMX
      INTEGER           MINEX(*),MMAN(*),MMNPC(*),MPMAN(*),MPMNPC(*)
C
C                                                ! local variables
C
      INTEGER           I,INODE,J,JE,JNODE,JS,K,L,LE,LS,M,ME,MS,NEN
C ----------------------------------------------------------------------
      DO 10 I=1,NANOD+1
         MPMAN(I) = 0
   10 CONTINUE
C
C --- estimate the degree of each node (always an overestimate)
C
      DO 20 I=1,NEL
         JS  = MPMNPC(I)
         JE  = MPMNPC(I+1) - 1
         NEN = JE-JS
         DO 15 J=JS,JE
            INODE = MMNPC(J)
            MPMAN(INODE) = MPMAN(INODE) + NEN
   15    CONTINUE
   20 CONTINUE
C
C --- reconstruct MPMAN to point to start of each set of neighbours
C
      J = 1
      DO 30 I=1,NANOD
         J = J + MPMAN(I)
         MPMAN(I) = J - MPMAN(I)
   30 CONTINUE
      MPMAN(NANOD+1) = J
C
C --- check storage, reset NMX, and initialize adjacency list
C
      IF (J-1 .GT. NMX) THEN
         J = J-NMX-1
         CALL NUMERR (3,NMX,J,J,LPU,IERR)
         GO TO 100
      ENDIF
      NMX = J-1
      DO 40 I=1,NMX
         MMAN(I) = 0
   40 CONTINUE
C
C --- form adjacency list (which may contain zeros)
C
      DO 70 I=1,NEL
         JS = MPMNPC(I)
         JE = MPMNPC(I+1) - 1
         DO 60 J=JS,JE-1
            INODE = MMNPC(J)
            LS    = MPMAN(INODE)
            LE    = MPMAN(INODE+1) - 1
            DO 55 K=J+1,JE
               JNODE = MMNPC(K)
               DO 50 L=LS,LE
                  IF (MMAN(L) .EQ. JNODE)  GO TO 55
                  IF (MMAN(L) .EQ. 0) THEN
                     MMAN(L) = JNODE
                     MS = MPMAN(JNODE)
                     ME = MPMAN(JNODE+1) - 1
                     DO 45 M=MS,ME
                        IF (MMAN(M) .EQ. 0) THEN
                           MMAN(M) = INODE
                           GO TO 55
                        ENDIF
   45                CONTINUE
                     CALL NUMERR (8,I,I,I,LPU,IERR)
                     GO TO 100
                  ENDIF
   50          CONTINUE
               CALL NUMERR (8,I,I,I,LPU,IERR)
               GO TO 100
   55       CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C --- strip any zeros from adjacency list
C
      K  = 0
      JS = 1
      DO 85 INODE=1,NANOD
         JE = MPMAN(INODE+1) - 1
         DO 75 J=JS,JE
            IF (MMAN(J) .EQ. 0) GO TO 80
            K = K+1
            MMAN(K) = MMAN(J)
   75    CONTINUE
   80    MPMAN(INODE+1) = K+1
         JS = JE+1
         IF (MPMAN(INODE+1) .EQ. MPMAN(INODE)) THEN
            CALL NUMERR (7,MINEX(INODE),I,I,LPU,IERR)
            GO TO 100
         ENDIF
   85 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE ELWALG (MPMNPC,MMNPC,MPMNCE,MMNCE,
     +                   MNNN,MGND,MGED,MNN,
     +                   NANOD,NEL,IPSW,LPU,NSTART,NHIGH,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ELWALG                 GROUP 9 / PRIVATE
C
C     TASK :  Renumber nodes by a modified simple ("brute force")
C             algorithm due to E.L.Wilson (EDSAP) - algorithm is
C             similar to reverse Cuthill-McKee.
C             Retained nodes (in substructure analysis), marked by a
C             -2 entry in the MNNN array, are not numbered.
C             "New" node number, INCUR, of "old" node INOD is returned
C             in MNNN(INOD).
C
C             Some key local variables:
C
C             NODE1 - number of start node (in terms of "old" numbers)
C             NELEL - number of elements eliminated
C             NNL   - number of nodes on list
C             IPASS - current "pass through" number
C
C
C     ROUTINES CALLED/REFERENCED :  NUMERR,  WENDEG     (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-25 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,IPSW,LPU,NANOD,NEL,NHIGH,NSTART,
     +                  MGED(NEL),MGND(NANOD),MMNCE(*),MMNPC(*),
     +                  MNNN(NANOD),MNN(NANOD),MPMNCE(*),MPMNPC(*)
C
C                                                ! local variables
C
      INTEGER           IEL,INCUR,INOD,IPASS,J,JE,JNOD,JS,
     +                  L,LE,LS,NDEEP,NELEL,NLARGE,NNL,NODE1
C
      PARAMETER         ( NDEEP=3 , NLARGE=100000000 )
C ----------------------------------------------------------------------
C                                                ! initiate
      IERR  = 0
      INCUR = 0
      IPASS = 0
      NELEL = 0
      NNL   = 0
      DO 20 IEL=1,NEL
         MGED(IEL) = 1
   20 CONTINUE
C ------------------------------------------------ print
      IF (IPSW.GT.1 .AND. LPU.GT.0) THEN
         WRITE (LPU,6000) NDEEP
      ENDIF
      IF (IPSW.GT.3 .AND. LPU.GT.0) THEN
         WRITE (LPU,6100)
      ENDIF
C ------------------------------------------------
C
C --- form weighted nodal and element "degrees" and determine start node
C
  100 NODE1 = 0
      IPASS = IPASS+1
      CALL WENDEG (MPMNPC,MMNPC,MGND,MGED,
     +             NANOD,NEL,NLARGE,NDEEP,NODE1)
C
      IF (NODE1 .EQ. 0)  THEN
         CALL NUMERR (6,NODE1,NODE1,NODE1,LPU,IERR)
         GO TO 1000
      ENDIF
C
      IF (IPASS .EQ. 1) THEN
         IF (NSTART .GT. 0) THEN
            NODE1 = NSTART
         ELSE
            NSTART = NODE1
         ENDIF
      ENDIF
C
C --- put start node on list
C
      NNL         = NNL+1
      MNN(NNL)    = NODE1
      MGND(NODE1) = 0
C
C --- add all nodes connected to current node to list
C
  200 INCUR = INCUR+1
C                                                ! disconnected ?
      IF (INCUR .GT. NNL)  GO TO 100
C
      INOD  = MNN(INCUR)
      LS    = MPMNCE(INOD)
      LE    = MPMNCE(INOD+1) - 1
C
C --- loop over elements connected to node
C
      DO 400 L=LS,LE
         IEL = MMNCE(L)
         IF (MGED(IEL) .GT. 0) THEN
            JS = MPMNPC(IEL)
            JE = MPMNPC(IEL+1) - 1
C
C --------- loop over nodes in element
C
            DO 300 J=JS,JE
               JNOD = MMNPC(J)
               IF (MGND(JNOD) .GT. 0) THEN
C                                                ! add node to list
C                                                  of numbered nodes
                  NNL        = NNL+1
                  MNN(NNL)   = JNOD
                  MGND(JNOD) = 0
               ENDIF
  300       CONTINUE
            MGED(IEL) = 0
            NELEL     = NELEL+1
C -------------------------------------------------------------- print
            IF (IPSW.GT.3 .AND. LPU.GT.0)  WRITE (LPU,6200) IEL
C -------------------------------------------------------------
         ENDIF
  400 CONTINUE
C
      IF (NELEL .LT. NEL)  GO TO 200
C
C --- reverse numbering (store temporarily in MGND)
C
      L     = NANOD+1
      INCUR = 0
  500 L = L-1
      INOD = MNN(L)
      IF (MNNN(INOD) .GE. 0) THEN
         INCUR = INCUR+1
         MGND(INCUR) = INOD
      ENDIF
      IF (L .GT. 1)  GO TO 500
C
      NHIGH = INCUR
C
C --- establish MNNN by converting the numbering
C
      DO 700 INCUR=1,NHIGH
         INOD = MGND(INCUR)
         MNNN(INOD) = INCUR
  700 CONTINUE
C
 1000 RETURN
C ----------------------------------------------- formats
 6000 FORMAT (5X,'Parameter NDEEP      =',I8 )
 6100 FORMAT ('1',
     +         5X,'Elements "eliminated" in the following order'/
     +         5X,'during node renumbering (by ELWNUM) :' / )
 6200 FORMAT (I12)
C
      END
      SUBROUTINE INNNI (MSC,MADOF,MNNN,NANOD,IFLAG,NHIGH)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  INNNI                  GROUP 9 / PRIVATE
C
C     TASK :  Initiate array MNNN ("Matrix of New Node Numbers") by
C             setting the entries of all NHIGH local (or internal) nodes
C             (having only local, status-1, and possibly specified dofs)
C             equal to 0 (null).  Retained (or external) nodes (having
C             at least one status-2 dof and no status-1 dofs) are given
C             new (final) numbers from NHIGH+1 to NANOD.
C             The above action is taken for IFLAG=1.
C             For all other values of IFLAG:
C             Retained nodes are disregarded, and all nodes are given
C             a null entry in MNNN and NHIGH is set equal to NANOD.
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   90-11-18 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IFLAG,NANOD,NHIGH,    MADOF(*),MNNN(*),MSC(*)

C                                                ! local variables
      INTEGER           IDOF,INOD,LDOF,LE,LS,NN
C ----------------------------------------------------------------------
C
      IF (IFLAG .EQ. 1) THEN
         NHIGH = 0
         DO 100 INOD=1,NANOD
            MNNN(INOD) = 0
            LS = MADOF(INOD)
            LE = MADOF(INOD+1) - 1
            DO 50 IDOF=LS,LE
               IF (MSC(IDOF) .EQ. 2) THEN
                  DO 25 LDOF=LS,LE
                     IF (MSC(LDOF) .EQ. 1)  GO TO 75
   25             CONTINUE
                  MNNN(INOD) =-1
                  GO TO 75
               ENDIF
   50       CONTINUE
   75       IF (MNNN(INOD) .EQ. 0)  NHIGH = NHIGH+1
  100    CONTINUE
         IF (NHIGH .LT. NANOD) THEN
            NN = NHIGH
            DO 200 INOD=1,NANOD
               IF (MNNN(INOD) .EQ. -1) THEN
                  NN = NN+1
                  MNNN(INOD) = NN
               ENDIF
  200       CONTINUE
         ENDIF
C
      ELSE
C
         NHIGH = NANOD
         DO 300 INOD=1,NANOD
            MNNN(INOD) = 0
  300    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE INSORT (LIST,KEY,NL,NK)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  INSORT                 GROUP 9 / PRIVATE
C
C     TASK :  To order a list of integers (LIST) in ascending order
C             of their keys (stored in KEY) using insertion sort
C
C     This is a slightly rewritten version of Sloan's subroutine ISORTI
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell  (after Scott W. Sloan)
C     DATE/VERSION  :   90-11-04 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NK,NL,    LIST(NL),KEY(NK)
C
C                                                ! local variables
      INTEGER           I,J,K,L
C ----------------------------------------------------------------------
      DO 20 I=2,NL
         J = LIST(I)
         K = KEY(J)
         DO 10 L=I-1,1,-1
            IF (K .GE. KEY(LIST(L))) THEN
               LIST(L+1) = J
               GO TO 20
            ENDIF
            LIST(L+1) = LIST(L)
   10    CONTINUE
         LIST(1) = J
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE NNNINT (MSC,MADOF,MNNN,NANOD,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NNNINT                 GROUP 9 / PRIVATE
C
C     TASK :  Initiate array MNNN ("Matrix of New Node Numbers") by
C             setting the entries of all retained nodes (having at
C             least one status-2 dof and no status-1 dofs) equal to -2.
C             All other (local) nodes (having only local, status-1, and
C             possibly specified dofs) are given numbers consecutively
C             from 1 as they are encountered when marching from 1 to
C             NANOD.
C             The above action is taken for IFLAG=1.
C             For all other values of IFLAG:
C             Retained nodes are disregarded, and all nodes are numbered
C             consecutively from 1 to NANOD.
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-16 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IFLAG,NANOD,    MADOF(*),MNNN(*),MSC(*)

C                                                ! local variables
      INTEGER           I,IDOF,INOD,LDOF,LE,LS
C ----------------------------------------------------------------------
C
      IF (IFLAG .EQ. 1) THEN
         I=1
         DO 100 INOD=1,NANOD
            MNNN(INOD) = I
            LS = MADOF(INOD)
            LE = MADOF(INOD+1) - 1
            DO 50 IDOF=LS,LE
               IF (MSC(IDOF) .EQ. 2) THEN
                  DO 25 LDOF=LS,LE
                     IF (MSC(LDOF) .EQ. 1)  GO TO 75
   25             CONTINUE
                  MNNN(INOD) =-2
                  GO TO 75
               ENDIF
   50       CONTINUE
   75       IF (MNNN(INOD) .GT. 0)  I = I+1
  100    CONTINUE
C
      ELSE
C
         DO 200 INOD=1,NANOD
            MNNN(INOD) = INOD
  200    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE NODELC (MINEX,MPMNPC,MMNPC,MTMP,MPMNCE,MMNCE,
     +                   NANOD,NEL,LPU,NW,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NODELC                 GROUP 9 / PRIVATE
C
C     TASK :  Establish the "matrix of node-connected elements",
C             MMNCE, and its associated pointer array, MPMNCE
C
C
C     ROUTINES CALLED/REFERENCED :  NUMERR         (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-21 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LPU,NANOD,NEL,NW,
     +                  MINEX(*),MMNCE(*),MMNPC(*),MPMNCE(*),MPMNPC(*),
     +                  MTMP(*)
C                                                ! local variables
C
      INTEGER           IEL,INOD,L,LE,LP,LS,NMMNCE
C ----------------------------------------------------------------------
C
      IERR = 0
C
      DO 10 INOD=1,NANOD
         MTMP(INOD) = 0
   10 CONTINUE
C
C --- count elements attached to nodes
C
      DO 50 IEL=1,NEL
         LS = MPMNPC(IEL)
         LE = MPMNPC(IEL+1) - 1
         DO 40 L=LS,LE
            INOD = MMNPC(L)
            MTMP(INOD) = MTMP(INOD) + 1
   40    CONTINUE
   50 CONTINUE
C
C --- form pointer array
C
      MPMNCE(1) = 1
      DO 100 INOD=2,NANOD+1
         IF (MTMP(INOD-1) .GT. 0) THEN
            MPMNCE(INOD) = MPMNCE(INOD-1) + MTMP(INOD-1)
         ELSE
            L = MINEX(INOD-1)
            CALL NUMERR (7,L,L,L,LPU,IERR)
         ENDIF
  100 CONTINUE
      IF (IERR .LT. 0)  GO TO 500
      NMMNCE = MPMNCE(NANOD+1) - 1
      IF (NMMNCE .GT. NW) THEN
         L  = NEL + 3*NANOD + NW + 1
         NW = NMMNCE-NW
         CALL NUMERR (3,L,NW,NW,LPU,IERR)
         GO TO 500
      ENDIF
C
C --- form "matrix of node-connected elements"
C
      DO 200 IEL=1,NEL
         LS = MPMNPC(IEL)
         LE = MPMNPC(IEL+1) - 1
         DO 150 L=LS,LE
            INOD = MMNPC(L)
            LP   = MPMNCE(INOD+1) - MTMP(INOD)
            MMNCE(LP)  = IEL
            MTMP(INOD) = MTMP(INOD) - 1
  150    CONTINUE
  200 CONTINUE
C
  500 RETURN
      END
      SUBROUTINE NUMERR (NUM,IVAL1,IVAL2,IVAL3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NUMERR                 GROUP 9 / PRIVATE
C
C     TASK :  Print error messages and set the error flag for the
C             renumbering routines
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-24 / 1.0
C                       99-04-26 / 1.1   K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IVAL1,IVAL2,IVAL3,LPU,NUM
C ----------------------------------------------------------------------
      IERR =-1
      IF (LPU .GT. 0) THEN
         WRITE (LPU,600)
         IF (NUM .EQ. 1) THEN
            WRITE (LPU,650)
            WRITE (LPU,601) IVAL1,IVAL2,IVAL3
         ELSEIF (NUM .EQ. 2) THEN
            WRITE (LPU,602)
         ELSEIF (NUM .EQ. 3) THEN
            WRITE (LPU,603) IVAL1,IVAL2
         ELSEIF (NUM .EQ. 4) THEN
            WRITE (LPU,650)
            WRITE (LPU,604) IVAL1,IVAL2
         ELSEIF (NUM .EQ. 5) THEN
            WRITE (LPU,650)
            WRITE (LPU,605)
         ELSEIF (NUM .EQ. 6) THEN
            WRITE (LPU,650)
            WRITE (LPU,606)
         ELSEIF (NUM .EQ. 7) THEN
            WRITE (LPU,607) IVAL1
         ELSEIF (NUM .EQ. 8) THEN
            WRITE (LPU,608)
            WRITE (LPU,650)
         ENDIF
      ENDIF
C
      RETURN
C -------------------------------------------------------------- formats
  600 FORMAT (///' *** ERROR RETURN DURING NODE RENUMBERING')
  601 FORMAT (5X,'Element',I8,'  connected to node',I8 /
     +        5X,'in problem with only',I8,'  nodes' )
  602 FORMAT (5X,'Temporary scratch storage far too small')
  603 FORMAT (5X,'Temporary scratch storage of',I8,'  locations' /
     +        5X,'exceeded by',I8,'  locations' )
  604 FORMAT (5X,I8,'  nodes out of',I8,'  are numbered' )
  605 FORMAT (5X,'No element connected to nodes on front')
  606 FORMAT (5X,'Cannot find starting node')
  607 FORMAT (5X,'No elements attached to node',I8 )
  608 FORMAT (5X,'Cannot assemble adjacency list')
C
  650 FORMAT (5X,'Inconsistent (corrupt?) data')
C
      END
      SUBROUTINE PFLNGT (MPMNPC,MMNPC,MNNN,MBW,NANOD,NEL,NPL,NBW)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PFLNGT                 GROUP 9 / PRIVATE
C
C     TASK :  Determine profile length, NPL, and max. bandwidth, NBW,
C             - in terms of nodes - of current node numbering contained
C             in MNNN
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-16 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NANOD,NBW,NEL,NPL,
     +                  MBW(*),MNNN(*),MMNPC(*),MPMNPC(*)
C
C                                                ! local variables
C
      INTEGER           IEL,INCUR,INOD,INOLD,KNCUR,KNOLD,L,LE,LL,LS
C ----------------------------------------------------------------------
C                                                ! initiate MBW
      DO 20 INOD=1,NANOD
         MBW(INOD) = INOD
   20 CONTINUE
C
C --- for each local (non-retained) node find its lowest numbered
C     adjacent node
C
      DO 100 IEL=1,NEL
         LS = MPMNPC(IEL)
         LE = MPMNPC(IEL+1) - 1
         DO 80 L=LS,LE
            INOLD = MMNPC(L)
            IF (MNNN(INOLD) .GT. 0) THEN
               INCUR = MNNN(INOLD)
               DO 60 LL=LS,LE
                  KNOLD = MMNPC(LL)
                  IF (MNNN(KNOLD) .GT. 0) THEN
                     KNCUR = MNNN(KNOLD)
                     IF (KNCUR .LT. MBW(INCUR)) THEN
                        MBW(INCUR) = KNCUR
                     ENDIF
                  ENDIF
   60          CONTINUE
            ENDIF
   80    CONTINUE
  100 CONTINUE
C
C --- determine profile and largest bandwidth
C
      NPL = 0
      NBW = 1
      DO 200 INOD=1,NANOD
         IF (MNNN(INOD) .GT. 0) THEN
            KNCUR      = MNNN(INOD)
            MBW(KNCUR) = KNCUR - MBW(KNCUR) + 1
            NPL        = NPL + MBW(KNCUR)
            IF (MBW(KNCUR) .GT. NBW)  NBW = MBW(KNCUR)
         ENDIF
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE PFMALG (MPMNPC,MMNPC,MPMNCE,MMNCE,
     +                   MNNN,MGND,MGED,MNF,
     +                   NANOD,NEL,IPSW,LPU,NSTART,NHIGH,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PFMALG                 GROUP 9 / PRIVATE
C
C     TASK :  Renumber nodes by a modified PFM-algorithm (due to Hoit,
C             Wilson and Garcelon).
C             Retained nodes (in substructure analysis), marked by a
C             -2 entry in the MNNN array, are not numbered.
C             "New" node number, INCUR, of "old" node INOD is returned
C             in MNNN(INOD).
C
C             Some key local variables:
C
C             NODE1 - number of start node (in terms of "old" numbers)
C             INCUR - current "new" node number
C             IPASS - current "pass-through" number
C             NNF   - current number of nodes in front
C
C
C     ROUTINES CALLED/REFERENCED :  NUMERR, PFRONT, WENDEG  (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-21 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IPSW,LPU,NANOD,NEL,NHIGH,NSTART,
     +                  MGED(NEL),MGND(NANOD),MMNCE(*),MMNPC(*),
     +                  MNNN(NANOD),MNF(NANOD),MPMNCE(*),MPMNPC(*)
C
C                                                ! local variables
C
      INTEGER           IEL,INCUR,INOD,INNF,IPASS,J,JE,JNOD,JS,JSUM,
     +                  KEL,L,LE,LS,NDEEP,NLARGE,NMIN,NNF,NODE1
C
      PARAMETER         ( NDEEP=3 , NLARGE=100000000 )
C ----------------------------------------------------------------------
C                                                ! initiate
      IERR  = 0
      INCUR = 0
      IPASS = 0
      DO 20 IEL=1,NEL
         MGED(IEL) = 1
   20 CONTINUE
C ------------------------------------------------ print
      IF (IPSW.GT.1 .AND. LPU.GT.0) THEN
         WRITE (LPU,6000) NDEEP
      ENDIF
      IF (IPSW.GT.3 .AND. LPU.GT.0) THEN
         WRITE (LPU,6100)
      ENDIF
C ------------------------------------------------
C
C --- form weighted nodal and element "degrees" and determine start node
C
  100 NODE1 = 0
      IPASS = IPASS+1
      CALL WENDEG (MPMNPC,MMNPC,MGND,MGED,
     +             NANOD,NEL,NLARGE,NDEEP,NODE1)
C                                                ! exit ?
      IF (NODE1 .EQ. 0)  GO TO 800
C                                                ! continue numbering
      IF (IPASS .EQ. 1) THEN
         IF (NSTART .GT. 0) THEN
            NODE1 = NSTART
         ELSE
            NSTART = NODE1
         ENDIF
      ENDIF
C
C --- put start node on front and mark it
C
      MNF(1) = NODE1
      NNF    = 1
      IF (MNNN(NODE1) .EQ. 0) THEN
         MNNN(NODE1) =-1
      ELSEIF (MNNN(NODE1) .EQ. -2) THEN
         MNNN(NODE1) =-3
      ENDIF
C
C --- find next element, KEL, connected to front to include (eliminate)
C --- choose element with "least connection"
C
  200 KEL  = 0
      NMIN = NLARGE
C
C --- loop over nodes on front
C
      DO 500 INNF=1,NNF
         INOD = MNF(INNF)
         LS   = MPMNCE(INOD)
         LE   = MPMNCE(INOD+1) - 1
C
C ------ loop over elements connected to node
C
         DO 400 L=LS,LE
            IEL = MMNCE(L)
            IF (MGED(IEL) .GT. 0) THEN
               JS   = MPMNPC(IEL)
               JE   = MPMNPC(IEL+1) - 1
               JSUM = 0
C
C ------------ loop over nodes connected to element
C
               DO 300 J=JS,JE
                  JNOD = MMNPC(J)
C
C --------------- decrease if only one element remains connected
C
                  IF (MGND(JNOD) .EQ. MGED(IEL))  JSUM = JSUM-3
C
C --------------- decrease if node is on front
C
                  IF (MNNN(JNOD) .EQ. -1)  JSUM = JSUM-1
C
C --------------- increase if node is not on front
C
                  IF (MNNN(JNOD) .EQ. 0)   JSUM = JSUM+3
  300          CONTINUE
               IF (JSUM .LT. NMIN) THEN
                  NMIN = JSUM
                  KEL  = IEL
               ENDIF
            ENDIF
  400    CONTINUE
  500 CONTINUE
C
      IF (KEL .EQ. 0)  THEN
         CALL NUMERR (5,KEL,KEL,KEL,LPU,IERR)
         GO TO 1000
      ENDIF
C
C --- move (eliminate) element - subtract its degree from the degrees
C     of the nodes connected to it - and number "completed" nodes
C
      LS = MPMNPC(KEL)
      LE = MPMNPC(KEL+1) - 1
      DO 600 L=LS,LE
         INOD = MMNPC(L)
         MGND(INOD) = MGND(INOD) - MGED(KEL)
C                                                ! add to front if new
C                                                  node
         CALL PFRONT (MNF,MNNN,INOD,NNF,1)
         IF (MGND(INOD) .EQ. 0) THEN
C                                                ! number node (if not
C                                                  retained)
            IF (MNNN(INOD) .EQ. -1) THEN
               INCUR      = INCUR+1
               MNNN(INOD) = INCUR
            ENDIF
C                                                ! remove node from
C                                                  front
            CALL PFRONT (MNF,MNNN,INOD,NNF,-1)
         ENDIF
  600 CONTINUE
      MGED(KEL) = 0
C ------------------------------------------------ print
C
      IF (IPSW.GT.3 .AND. LPU.GT.0) WRITE (LPU,6200) KEL
C ------------------------------------------------------
C
C --- if no more nodes on front, go back to 100 and check that all
C     elements have been eliminated (they may not if structure is
C     "disconnected") and exit from there (or continue)
C     else go to 200 and continue
C
      IF (NNF .EQ. 0)  GO TO 100
      GO TO 200
C
  800 NHIGH = INCUR
C
C --- reverse numbering
C
      DO 850 INOD=1,NANOD
         INCUR = MNNN(INOD)
         IF (INCUR .GT. 0)  MNNN(INOD) = NHIGH-INCUR+1
  850 CONTINUE
C
C
 1000 RETURN
C ----------------------------------------------- formats
 6000 FORMAT (5X,'Parameter NDEEP      =',I8 )
 6100 FORMAT ('1',
     +         5X,'Elements "eliminated" in the following order'/
     +         5X,'during node renumbering (by PFMNUM) :' / )
 6200 FORMAT (I12)
C
      END
      SUBROUTINE PFRONT (MNF,MNNN,INOD,NNF,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PFRONT                 GROUP 9 / PRIVATE
C
C     TASK : Update nodes on front
C            IFLAG > 0 : Add node to front
C            IFLAG < 0 : Remove node from front
C
C            Convention:
C            MNNN(INOD) > 0 : node is numbered
C                       = 0 : node is not numbered and not on front
C                       =-1 : node is on front
C                       =-2 : node is flagged as retained but not
C                             on front
C                       =-3 : node is flagged as retained and on front
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-16 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IFLAG,INOD,NNF,   MNF(*),MNNN(*)
C
C                                                ! local variables
      INTEGER           I,J,KODE
C ----------------------------------------------------------------------
      IF (IFLAG .GT. 0) THEN
C                                                ! add to front
         KODE = MNNN(INOD)
         IF (KODE .EQ. 0) THEN
C                                                ! node not on front and
C                                                  not numbered
            NNF        = NNF+1
            MNF(NNF)   = INOD
            MNNN(INOD) =-1
         ELSEIF (KODE .EQ. -2) THEN
C                                                ! retained node
            NNF        = NNF+1
            MNF(NNF)   = INOD
            MNNN(INOD) =-3
         ENDIF
C
      ELSEIF (IFLAG .LT. 0) THEN
C                                                ! remove from front
         I = 0
   10    I = I+1
         IF (INOD .EQ. MNF(I)) THEN
            NNF = NNF-1
            IF (NNF .GE. I) THEN
C                                                ! pack
               DO 20 J=I,NNF
                  MNF(J) = MNF(J+1)
   20          CONTINUE
            ENDIF
            GO TO 100
         ENDIF
         IF (I .LT. NNF)  GO TO 10
      ENDIF
C
  100 RETURN
      END
      SUBROUTINE PPNODE (MPMAN,MMAN,MASK,MPLS,MLS,MHLN,
     +                   NANOD,NSTART,NEND,NNC)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PPNODE                 GROUP 9 / PRIVATE
C
C     TASK :  To find nodes which define a pseudo-diameter of a graph,
C             and store distance from end node (NEND) in MASK.
C             Initial start node may be input (in NSTART) or else is
C             chosen as node of minimum degree.
C
C     This is a rewritten version of Sloan's subroutine DIAMTR
C
C
C     ROUTINES CALLED/REFERENCED :  RLVSTR and INSORT     (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell   (after Scott W. Sloan)
C     DATE/VERSION  :   90-11-05 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NANOD,NEND,NNC,NSTART
      INTEGER           MASK(*),MHLN(*),MLS(*),MMAN(*),MPMAN(*),MPLS(*)
C
C                                                ! local variables
C
      INTEGER           I,IDEG,IE,IH,INODE,IS,L,LD,LW,
     +                  N1,NDEPTH,NWIDTH,MINDEG
C ----------------------------------------------------------------------
      IF (NSTART .EQ. 0) THEN
C
C ------ choose first guess for starting node by min degree - ignore
C        nodes that are invisible (MASK .ne. 0)
C
         MINDEG = NANOD
         DO 10 INODE=1,NANOD
            IF (MASK(INODE) .EQ. 0) THEN
               IDEG = MPMAN(INODE+1) - MPMAN(INODE)
               IF (IDEG .LT. MINDEG) THEN
                  NSTART = INODE
                  MINDEG = IDEG
               ENDIF
            ENDIF
   10    CONTINUE
      ENDIF
C
C --- generate level structure for first start node
C
      N1  = NANOD+1
      CALL RLVSTR (MPMAN,MMAN,MASK,MPLS,MLS,NSTART,N1,NDEPTH,NWIDTH)
      NNC = MPLS(NDEPTH+1) - 1
C ----------------------------------------------------------------------
C     iterate to find start (NSTART) and end (NEND) node
C ----------------------------------------------------------------------
C
C --- store list of nodes that are at max distance from starting node,
C     and store their degrees in MPLS
C
   20 IH = 0
      IS = MPLS(NDEPTH)
      IE = MPLS(NDEPTH+1) - 1
      DO 30 I=IS,IE
         INODE = MLS(I)
         IH    = IH+1
         MHLN(IH) = INODE
         MPLS(INODE) = MPMAN(INODE+1) - MPMAN(INODE)
   30 CONTINUE
C
C --- sort list of nodes in ascending sequence of their degree
C
      IF (IH .GT. 1)  CALL INSORT (MHLN,MPLS,IH,NANOD)
C
C --- remove nodes with duplicate degree
C
      IE = IH
      IH = 1
      IDEG = MPLS(MHLN(1))
      DO 40 I=2,IE
         INODE = MHLN(I)
         IF (MPLS(INODE) .NE. IDEG) THEN
            IDEG     = MPLS(INODE)
            IH       = IH+1
            MHLN(IH) = INODE
         ENDIF
   40 CONTINUE
C
C --- loop over nodes in shrunken level
C
      LW = NNC+1
      DO 50 I=1,IH
         INODE = MHLN(I)
         CALL RLVSTR (MPMAN,MMAN,MASK,MPLS,MLS,INODE,LW,LD,NWIDTH)
         IF (NWIDTH .LT. LW) THEN
C
C --------- level structure was not aborted during assembly
C
            IF (LD .GT. NDEPTH) THEN
C
C ------------ level structure of greater depth found - store new
C              starting node, new max depth, and begin new iteration
C
               NSTART = INODE
               NDEPTH = LD
               GO TO 20
            ENDIF
C
C --------- level structure width for this end node is smallest so far
C           store end node and new min width
C
            NEND = INODE
            LW   = NWIDTH
         ENDIF
   50 CONTINUE
C
C --- generate level structure rooted at end node if necessary
C
      IF (INODE .NE. NEND) THEN
         CALL RLVSTR (MPMAN,MMAN,MASK,MPLS,MLS,NEND,N1,LD,NWIDTH)
      ENDIF
C
C --- store distance of each node from end node
C
      DO 70 L=1,LD
         IS = MPLS(L)
         IE = MPLS(L+1) - 1
         DO 60 I=IS,IE
            MASK(MLS(I)) = L-1
   60    CONTINUE
   70 CONTINUE
C
      RETURN
      END
      SUBROUTINE PRFLL (MPMAN,MMAN,MNNN,NANOD,
     +                  LOLD,LNEW,NBOLD,NBNEW)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRFLL                  GROUP 9 / PRIVATE
C
C     TASK :  To determine profile length and bandwidth corresponding
C             to both old and new numbering.
C
C     This is a slightly rewritten version of Sloan's subroutine PROFIL
C
C
C     ROUTINES CALLED/REFERENCED :  MIN    (Fortran library)
C
C     PROGRAMMED BY :   Kolbein Bell  (after Scott W. Sloan)
C     DATE/VERSION  :   90-11-08 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           LNEW,LOLD,NANOD,NBNEW,NBOLD,
     +                  MMAN(*),MNNN(*),MPMAN(*)
C
C                                                ! local variables
C
      INTEGER           INODE,J,JE,JS,KN,KO,MOLD,MNEW
C ----------------------------------------------------------------------
      LOLD  = 0
      LNEW  = 0
      NBOLD = 0
      NBNEW = 0
C
C --- do for each node in graph
C
      DO 20 INODE=1,NANOD
         JS   = MPMAN(INODE)
         JE   = MPMAN(INODE+1) - 1
         MOLD = MMAN(JS)
         MNEW = MNNN(MMAN(JS))
C
C ------ find lowest numbered neighbour of node INODE
C
         DO 10 J=JS+1,JE
            MOLD = MIN(MOLD,MMAN(J))
            MNEW = MIN(MNEW,MNNN(MMAN(J)))
   10    CONTINUE
C
C ------ update profiles and if appropriate bandwidths
C
         KO = INODE - MOLD
         KN = MNNN(INODE) - MNEW
         IF (KO .GT. NBOLD) NBOLD = KO
         IF (KN .GT. NBNEW) NBNEW = KN
         IF (KO .GT. 0)     LOLD  = LOLD+KO
         IF (KN .GT. 0)     LNEW  = LNEW+KN
   20 CONTINUE
C
C --- add diagonal terms
C
      LOLD  = LOLD+NANOD
      LNEW  = LNEW+NANOD
      NBOLD = NBOLD+1
      NBNEW = NBNEW+1
C
      RETURN
      END
      SUBROUTINE RLVSTR (MPMAN,MMAN,MASK,MPLS,MLS,NROOT,MAXWDT,
     +                   NDEPTH,NWIDTH)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  RLVSTR                 GROUP 9 / PRIVATE
C
C     TASK :  To generate a rooted level structure for the "active"
C             component of the graph defined by MPMAN and MMAN.
C             The nodes of each level I (= 1,..,NDEPTH) of the graph
C             (starting at node NROOT) is stored consecutively in MLS
C             from MPLS(I) onwards.
C
C     This is a rewritten version of Sloan's subroutine ROOTLS
C
C     Some local variables:
C     NNC - final value is the number of nodes in current component
C           og graph
C     LS  - points to start of current level
C     LE  - points to end of current level
C     LW  - width of current level
C
C
C     ROUTINES CALLED/REFERENCED :  MAX     (Fortran library)
C
C     PROGRAMMED BY :   Kolbein Bell   (after Scott W. Sloan)
C     DATE/VERSION  :   90-11-04 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           MAXWDT,NDEPTH,NROOT,NWIDTH
      INTEGER           MASK(*),MLS(*),MMAN(*),MPMAN(*),MPLS(*)
C
C                                                ! local variables
C
      INTEGER           INODE,J,JE,JNODE,JS,L,LE,LS,LW,NNC
C ----------------------------------------------------------------------
      MASK(NROOT) = 1
      MLS(1)      = NROOT
      NNC    = 1
      NWIDTH = 1
      NDEPTH = 0
      LW     = 1
      LE     = 0
C
   10 IF (LW .GT. 0) THEN
         LS     = LE+1
         LE     = NNC
         NDEPTH = NDEPTH+1
         MPLS(NDEPTH) = LS
C
C ------ generate next level by finding all visible neighbours of nodes
C        in current level
C
         DO 30 L=LS,LE
            INODE = MLS(L)
            JS    = MPMAN(INODE)
            JE    = MPMAN(INODE+1) - 1
            DO 20 J=JS,JE
               JNODE = MMAN(J)
               IF (MASK(JNODE) .EQ. 0) THEN
                  NNC         = NNC+1
                  MLS(NNC)    = JNODE
                  MASK(JNODE) = 1
               ENDIF
   20       CONTINUE
   30    CONTINUE
C
C ------ compute width of level just assembled and the width of the
C        level structure so far
C
         LW     = NNC-LE
         NWIDTH = MAX(LW,NWIDTH)
C
C ------ abort assembly if level structure is too wide
C
         IF (NWIDTH .GE. MAXWDT)  GO TO 40
         GO TO 10
      ENDIF
      MPLS(NDEPTH+1) = LE+1
C
C --- reset MASK=0 for nodes in the level structure
C
   40 DO 50 L=1,NNC
         MASK(MLS(L)) = 0
   50 CONTINUE
C
      RETURN
      END
      SUBROUTINE SWSALG (MPMAN,MMAN,MNNN,MNQ,MNP,
     +                   NNC,NSTART,LSTNUM)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SWSALG                 GROUP 9 / PRIVATE
C
C     TASK :  To (re)number the nodes in the cuurrent component of the
C             graph by an algorithm due to S.W.Sloan.
C             The graph is defined by arrays MPMAN and MMAN.
C             The NNC nodes in the current component of the graph are
C             contained in array MNQ (this array is also used to store
C             queue of active and preactive nodes).
C             MNNN contains on input the new numbers of nodes already
C             numbered and the distance of each node in this component
C             from the end node.
C             NSTART is the ("old") node at which numbering starts.
C             LSTNUM is the count of nodes which have alredy been
C             numbered (is updated by the routine).
C
C             During the numbering process, MNNN serves as a list
C             giving the status of the nodes:
C             MNNN(I) gt  0 indicates node I is postactive
C             MNNN(I)  =  0 indicates node I is active
C             MNNN(I)  = -1 indicates node I is preactive
C             MNNN(I)  = -2 indicates node I is inactive
C
C     This is a rewritten version of Sloan's subroutine NUMBER
C
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell   (after Scott W. Sloan)
C     DATE/VERSION  :   90-11-18 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           LSTNUM,NNC,NSTART
      INTEGER           MMAN(*),MNNN(*),MNP(*),MNQ(*),MPMAN(*)
C
C                                                ! local variables
C
      INTEGER           I,IE,INDX,INODE,IPRI,IS,IW1,IW2,J,JE,JS,
     +                  MXPRI,NBHR,NBHR1,NEXT,NQ
C
      PARAMETER         (IW1=1 , IW2=2)
C ----------------------------------------------------------------------
C
C --- initialise priorities and status for each node in this component:
C     = iw1*(distance from end node) - iw2*(initial current degree)
C
      DO 10 I=1,NNC
         INODE       = MNQ(I)
         MNP(INODE)  = IW1*MNNN(INODE)
     +               - IW2*(MPMAN(INODE+1) - MPMAN(INODE) + 1)
         MNNN(INODE) = -2
   10 CONTINUE
C
C --- insert starting node in queue and assign it a preactive status
C     NQ is the size of queue
C
      NQ = 1
      MNQ(NQ)      = NSTART
      MNNN(NSTART) = -1
C
C === loop while queue is not empty ====================================
C
   20 IF (NQ .GT. 0) THEN
C
C ------ scan queue for node with max priority (MXPRI)
C
         INDX  = 1
         MXPRI = MNP(MNQ(1))
         IF (NQ .GT. 1) THEN
            DO 30 I=2,NQ
               IPRI = MNP(MNQ(I))
               IF (IPRI .GT. MXPRI) THEN
                  INDX  = I
                  MXPRI = IPRI
               ENDIF
   30       CONTINUE
         ENDIF
C
C ------ NEXT is the node to be numbered next - delete it from queue
C
         NEXT      = MNQ(INDX)
         MNQ(INDX) = MNQ(NQ)
         NQ        = NQ-1
         IS        = MPMAN(NEXT)
         IE        = MPMAN(NEXT+1) - 1
C
         IF (MNNN(NEXT) .EQ. -1) THEN
C
C --------- node NEXT is preactive - examine its neighbours;
C           decrease current degree of neighbour (NBHR) by -1;
C           add neighbour to queue if it is inactive and assign it a
C           preactive status
C
            DO 40 I=IS,IE
               NBHR      = MMAN(I)
               MNP(NBHR) = MNP(NBHR) + IW2
               IF (MNNN(NBHR) .EQ. -2) THEN
                  NQ         = NQ+1
                  MNQ(NQ)    = NBHR
                  MNNN(NBHR) = -1
               ENDIF
   40       CONTINUE
         ENDIF
C
C ------ store new number for node NEXT (gives it postactive status)
C
         LSTNUM     = LSTNUM+1
         MNNN(NEXT) = LSTNUM
C
C ------ search for preactive neighbours of node NEXT
C
         DO 60 I=IS,IE
            NBHR = MMAN(I)
            IF (MNNN(NBHR) .EQ. -1) THEN
C
C ------------ decrease current degree of preactive neighbour by -1 and
C              assign neighbour an active status
C
               MNP(NBHR)  = MNP(NBHR) + IW2
               MNNN(NBHR) = 0
C
C ------------ loop over nodes adjacent to preactive neighbour and
C              decrease decrease their current degrees by -1
C
               JS = MPMAN(NBHR)
               JE = MPMAN(NBHR+1) - 1
               DO 50 J=JS,JE
                  NBHR1      = MMAN(J)
                  MNP(NBHR1) = MNP(NBHR1) + IW2
                  IF (MNNN(NBHR1) .EQ. -2) THEN
C
C ------------------ insert inactive node in queue with preactive status
C
                     NQ          = NQ+1
                     MNQ(NQ)     = NBHR1
                     MNNN(NBHR1) = -1
                  ENDIF
   50          CONTINUE
            ENDIF
   60    CONTINUE
C
         GO TO 20
C ======================================================================
      ENDIF
C
      RETURN
      END
      SUBROUTINE WENDEG (MPMNPC,MMNPC,MGND,MGED,
     +                   NANOD,NEL,NLARGE,NDEEP,NSTART)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  WENDEG                 GROUP 9 / PRIVATE
C
C     TASK :  Form weighted ("global") nodal and element degrees -
C             compounded to a level of NDEEP - to obtain indication of
C             interconnectivity.
C             "Least connected" node is chosen as starting node.
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-11-18 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           NANOD,NDEEP,NEL,NLARGE,NSTART,
     +                  MGED(*),MGND(*),MMNPC(*),MPMNPC(*)
C
C                                                ! local variables
C
      INTEGER           IEL,INOD,ISUM,L,LE,LEVEL,LS,NMIN
C ----------------------------------------------------------------------
      DO 10 INOD=1,NANOD
         MGND(INOD) = 1
   10 CONTINUE
C
C --- "compounded" loop to evaluate weigthed nodal and element degrees:
C      element degree = sum of nodal degrees of nodes in element
C      nodal degree   = sum of element degrees of all elements con-
C                       nected to node
C
      DO 200 LEVEL=1,NDEEP
         DO 40 IEL=1,NEL
            IF (MGED(IEL) .GT. 0) THEN
               LS = MPMNPC(IEL)
               LE = MPMNPC(IEL+1) - 1
               ISUM = 0
               DO 20 L=LS,LE
                  INOD = MMNPC(L)
                  ISUM = ISUM + MGND(INOD)
   20          CONTINUE
C                                                ! element degree
               MGED(IEL) = ISUM
            ENDIF
   40    CONTINUE
C                                                ! init. nodal degrees
         DO 60 INOD=1,NANOD
            MGND(INOD) = 0
   60    CONTINUE
C                                                ! comp. nodal degrees
         DO 100 IEL=1,NEL
            IF (MGED(IEL) .GT. 0) THEN
               LS = MPMNPC(IEL)
               LE = MPMNPC(IEL+1) - 1
               DO 80 L=LS,LE
                  INOD = MMNPC(L)
                  MGND(INOD) = MGND(INOD) + MGED(IEL)
   80          CONTINUE
            ENDIF
  100    CONTINUE
  200 CONTINUE
C
C --- find node with minimum degree
C
      NMIN = NLARGE
      DO 300 INOD=1,NANOD
         IF (MGND(INOD) .LE. NMIN) THEN
            IF (MGND(INOD) .GT. 0) THEN
               NMIN   = MGND(INOD)
               NSTART = INOD
            ENDIF
         ENDIF
  300 CONTINUE
C
      RETURN
      END

