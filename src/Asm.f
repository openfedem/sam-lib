C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE PRASS (MADOF,MPMNPC,MMNPC,
     +                  MSC,MPMCEQ,MMCEQ,LPU,
     +                  MPAR,MEEN,MEQN,MSKY,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRASS                   GROUP 7 / PUBLIC
C
C     T A S K :  TO PREPARE (CHECK AND DETERMINE) CONTROL INFORMATION
C     FOR THE ASSEMBLY PROCESS (PRE-ASSEMBLY)
C
C     ROUTINES CALLED/REFERENCED :  SYSEQ, SKYLIN,
C                                   BANDW AND ASMERR     (SAM - 7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-24 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*),MSC(*),MSKY(*)
C
      INTEGER           NDOF,NEL,NANOD
C
      EXTERNAL          ASMERR,BANDW,SKYLIN,SYSEQ
C ----------------------------------------------------------------------
      IERR  = 0
      NANOD = MPAR(1)
      NEL   = MPAR(2)
      NDOF  = MADOF(NANOD+1) - 1
      IF (MPAR(3).NE.NDOF)   CALL ASMERR (7,NDOF,NDOF,LPU,IERR)
      IF (IERR.LT.0)         GO TO 100
      MPAR(15) = MPMNPC(NEL+1) - 1
C
      CALL SYSEQ (MSC,MPMCEQ,MMCEQ,LPU,MPAR,MEQN,IERR)
      IF (IERR.LT.0)         GO TO 100
C
      CALL SKYLIN (MADOF,MPMNPC,MMNPC,MPMCEQ,MMCEQ,MEQN,MPAR,MEEN,MSKY)
C
      CALL BANDW (MADOF,MEQN,MSKY,MPAR)
C
  100 RETURN
      END
      SUBROUTINE PRASSX (MADOF,MINEX,MPMNPC,MMNPC,
     +                   MSC,MPMCEQ,MMCEQ,MNNN,LPU,
     +                   MPAR,MWORK,MEQN,MSKY,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRASSX                  GROUP 7 / PUBLIC
C
C     TASK :  To prepare (check and determine) control information for
C             the assembly process (pre-assembly) assuming a "skyline"
C             storage format - control arrays MEQN and MSKY are de-
C             termined, as well as a number of parameters in MPAR.
C             This is a special version of PRASS in which MEQN, and
C             consequently MSKY, are determined with due consideration
C             to "new" node numbers (recorded in MNNN).
C             The relationship between internal node numbers and
C             external node identifiers is recorded in MINEX.
C
C
C     ROUTINES CALLED/REFERENCED :  SYSEQX, SKYLIN,
C                                   BANDW AND ASMERR     (SAM - 7)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-03-14 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU
      INTEGER           MADOF(*),MEQN(*),MINEX(*),MMCEQ(*),MMNPC(*),
     &                  MNNN(*),MPAR(*),MPMCEQ(*),MPMNPC(*),MSC(*),
     &                  MSKY(*),MWORK(*)
C
      INTEGER           NDOF,NEL,NANOD
C ----------------------------------------------------------------------
      IERR  = 0
      NANOD = MPAR(1)
      NEL   = MPAR(2)
      NDOF  = MADOF(NANOD+1) - 1
      IF (MPAR(3).NE.NDOF) THEN
         CALL ASMERR (27,NDOF,NDOF,LPU,IERR)
         GO TO 100
      ENDIF
      MPAR(15) = MPMNPC(NEL+1) - 1
C
      CALL SYSEQX (MADOF,MINEX,MSC,MPMCEQ,MMCEQ,
     &             MNNN,LPU,MPAR,MEQN,MWORK,IERR)
      IF (IERR.LT.0)  GO TO 100
C
      CALL SKYLIN (MADOF,MPMNPC,MMNPC,MPMCEQ,MMCEQ,MEQN,MPAR,MWORK,MSKY)
C
      CALL BANDW (MADOF,MEQN,MSKY,MPAR)
C
  100 RETURN
      END
      SUBROUTINE ADDEM (EM,TTCC,MPAR,MADOF,MEQN,
     +                  MPMNPC,MMNPC,MPMCEQ,MMCEQ,MSKY,
     +                  IEL,NEDOF,LPU,NSV,
     +                  SM,SV,MEEN,IERR )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDEM                   GROUP 7 / PUBLIC
C
C     T A S K :  TO ADD/SUBTRACT THE ELEMENTS CORRESPONDING TO FREE
C     AND CONSTRAINED  D O F S  OF THE ELEMENT MATRIX  EM ( OF ELEMENT
C     NO. IEL) TO/FROM THE CURRENT CONTENT OF THE APPROPRIATE ELEMENTS
C     OF
C       MATRIX     SM  IF  NSV.GE.0  (ADDITION WITHOUT/WITH WEIGHTING)
C       VECTOR(S)  SV  IF  NSV.NE.0  (SUBTRACTION WITH WEIGHTING)
C     SM IS A SYMMETRIC SYSTEM MATRIX STORED IN 'SKYLINE' FORM, AND
C     SV IS A SYSTEM VECTOR (ABS(NSV)=1) OR A SET OF SYSTEM VECTORS
C     (ABS(NSV)=NPDOF) CONTAINING THE CONTRIBUTIONS TO THE RIGHT-HAND
C     SIDE(S) FROM PRESCRIBED AND (POSSIBLY) DEPENDENT D O F S.
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ,SLVEQ,NOPRED, ASMERR (SAM-7)
C                                   ABS              (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-01-14 / 1.0
C                       89-01-13 / 1.1   K.Bell
C                       00-11-13 / 1.2   B.Haugen
C                       04-01-15 / 1.3   K.M.Okstad
C                       05-10-10 / 1.4   K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF,NSV
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*),MSKY(*)
      DOUBLE PRECISION  TTCC(*),EM(NEDOF,NEDOF),SM(*),SV(*)
C
      INTEGER           I,ICEQ,IEQ,II,IP,J,JCEQ,JEQ,JJ,JP,JV,K,KEQ,KP,
     +                  NECEQ,NELDOF,NENOD,NEPRD,NESLV,NMSTI,NMSTJ,NV,
     +                  NOPRED
      DOUBLE PRECISION  C0,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR,ELMEQ,NOPRED
C ----------------------------------------------------------------------
      IERR = 0
      IF (IEL.GT.0 .AND. IEL.LE.MPAR(2))  GO TO 10
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
C
   10 NV = ABS(NSV)
      IF (NV.LT.2)             GO TO 20
      IF (NV.EQ.MPAR(9))       GO TO 20
      CALL ASMERR (17,IEL,NSV,LPU,IERR)
      GO TO 1000
C
   20 NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
      IP    = MPMNPC(IEL)
      CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +            MEEN,NELDOF,NESLV,NEPRD )
      IF (NELDOF.EQ.NEDOF)     GO TO 40
      CALL ASMERR (11,IEL,IEL,LPU,IERR)
      GO TO 1000
   40 NECEQ = NESLV+NEPRD
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
            K   = MSKY(IEQ)
            IF (IEQ.LT.JEQ)    K = MSKY(JEQ)-JEQ+IEQ
            IF (IEQ.GT.JEQ)    K = K-IEQ+JEQ
            SM(K) = SM(K) + EM(I,J)
  150    CONTINUE
  200 CONTINUE
C
  220 IF (NECEQ.EQ.0)                  GO TO 1000
      IF (NSV.LT.-1 .AND. NEPRD.EQ.0)  GO TO 1000
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
         IF (NV.EQ.1)          GO TO 250
         JV   = NOPRED(MPMCEQ,JCEQ)
         IF (JV.EQ.0)          GO TO 310
         KP   = (JV-1)*MPAR(11)
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
            nmsti = mpmceq(jceq+1) - mpmceq(jceq) - 1
            IF (NMSTI.EQ.0)    GO TO 300
            IP = MPMCEQ(ICEQ)
            DO 280 II=1,NMSTI
               IP  = IP+1
               if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) go to 280
               KEQ = KP+meqn(mmceq(ip))
               SV(KEQ) = SV(KEQ) - C0*TTCC(IP)*EM(I,J)
  280       CONTINUE
  300    CONTINUE
C
  310    IF (NSV.LT.0)         GO TO 600
  320    nmstj = mpmceq(jceq+1) - mpmceq(jceq) - 1
         IF (NMSTJ.EQ.0)       GO TO 600
C                                               ** ADD CONTRIBUTIONS
C                                                  TO  SM
         DO 500 JJ=1,NMSTJ
            JP  = JP+1
            if (mmceq(jp).le.0 .or. ttcc(jp).eq.zero) go to 500
            jeq = meqn(mmceq(jp))
            DO 400 I=1,NEDOF
               IEQ = MEEN(I)
               IF (IEQ.EQ.0)   GO TO 400
               IF (IEQ.LT.0)   GO TO 360
               IF (IEQ.LT.JEQ) GO TO 330
               IF (IEQ.GT.JEQ) GO TO 340
               K     = MSKY(IEQ)
               SM(K) = SM(K) + (TTCC(JP)+TTCC(JP))*EM(I,J)
               GO TO 400
  330          K     = MSKY(JEQ)-JEQ+IEQ
               GO TO 350
  340          K     = MSKY(IEQ)-IEQ+JEQ
  350          SM(K) = SM(K) + TTCC(JP)*EM(I,J)
               GO TO 400
C
  360          ICEQ  =-IEQ
               nmsti = mpmceq(iceq+1) - mpmceq(iceq) - 1
               IF (NMSTI.EQ.0) GO TO 400
               IP    = MPMCEQ(ICEQ)
               DO 380 II=1,NMSTI
                  IP  = IP+1
                  if (mmceq(ip).le.0 .or. ttcc(ip).eq.zero) go to 380
                  ieq = meqn(mmceq(ip))
                  IF (IEQ.GT.JEQ)  GO TO 380
                  K   = MSKY(JEQ)-JEQ+IEQ
                  SM(K) = SM(K) + TTCC(IP)*TTCC(JP)*EM(I,J)
  380          CONTINUE
  400       CONTINUE
  500    CONTINUE
  600 CONTINUE
C ----------------------------------------------------------------------
 1000 RETURN
      END
      SUBROUTINE ADDEV (EV,TTCC,MPAR,MADOF,MEQN,
     +                  MPMNPC,MMNPC,MPMCEQ,MMCEQ,
     +                  IEL,NEDOF,LPU,SV,MEEN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDEV                   GROUP 7 / PUBLIC
C
C     T A S K :  TO ADD THE ELEMENTS CORRESPONDING TO FREE AND DE-
C     PENDENT  D O F S  OF THE ELEMENT VECTOR  EV  (ELEMENT NO. IEL)
C     TO THE CURRENT CONTENT OF THE APPROPRIATE ELEMENTS OF THE SYSTEM
C     VECTOR  SV
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ, SLVEQ AND ASMERR  (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-04-01 / 1.0
C                       89-01-13 / 1.1   K.Bell
C                       00-11-13 / 1.2   B.Haugen
C                       04-01-15 / 1.3   K.M.Okstad
C                       05-10-10 / 1.4   K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IEL,IERR,LPU,NEDOF
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*)
      DOUBLE PRECISION  TTCC(*),EV(NEDOF),SV(*)
C
      INTEGER           I,ICEQ,IEQ,II,IP,NELDOF,NENOD,NEPRD,NESLV,NMST
C
      EXTERNAL          ASMERR,ELMEQ
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (IEL.GT.0 .AND. IEL.LE.MPAR(2))  GO TO 10
      CALL ASMERR (12,IEL,IEL,LPU,IERR)
      GO TO 100
C
   10 NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
      IP    = MPMNPC(IEL)
      CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +            MEEN,NELDOF,NESLV,NEPRD)
      IF (NELDOF.EQ.NEDOF)  GO TO 20
      CALL ASMERR (12,IEL,IEL,LPU,IERR)
      WRITE(LPU,"(5X,'NELDOF, NEDOF =',2I6/)") NELDOF, NEDOF
      GO TO 100
C
   20 DO 50 I=1,NEDOF
         IEQ = MEEN(I)
         IF (IEQ.EQ.0)      GO TO 50
         IF (IEQ.LT.0)      GO TO 30
         SV(IEQ) = SV(IEQ) + EV(I)
         GO TO 50
C
   30    ICEQ =-IEQ
         nmst = mpmceq(iceq+1) - mpmceq(iceq) - 1
         IF (NMST.EQ.0)     GO TO 50
         IP   = MPMCEQ(ICEQ)
         DO 40 II=1,NMST
            IP  = IP+1
            if (mmceq(ip).le.0) go to 40
            ieq = meqn(mmceq(ip))
            SV(IEQ) = SV(IEQ) + TTCC(IP)*EV(I)
   40    CONTINUE
   50 CONTINUE
C ----------------------------------------------------------------------
  100 RETURN
      END
      SUBROUTINE ADDNV ( VNOD, TTCC,   MPAR, MADOF,
     &                   MEQN, MPMCEQ, MMCEQ,
     &                   INOD, LPU,    SV,   IERR )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ADDNV                   GROUP 7 / PUBLIC
C
C     TASK :  To add the elements corresponding to free and dependent
C             d o f s  of the vector VNOD, which contains one element
C             per  d o f  of node INOD, to the current content of the
C             appropriate elements of the system vector SV.
C             The elements of VNOD are assumed to be expressed in the
C             same coordinates as the corresponding elements of the
C             system vector SV.
C
C     ROUTINES CALLED/REFERENCED :    SLVEQ AND ASMERR  (SAM-7)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-03-26 / 1.0
C                       04-01-15 / 1.1   K.M.Okstad
C                       05-10-10 / 1.2   K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           INOD,IERR,LPU
      INTEGER           MADOF(*),MMCEQ(*),MEQN(*),MPAR(*),MPMCEQ(*)
      DOUBLE PRECISION  TTCC(*),VNOD(*),SV(*)
C
      INTEGER           ICEQ,IDOF,IEQ,II,IP,LDOF,NMST,NNDOF
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (INOD.LT.1 .OR. INOD.GT.MPAR(1)) THEN
         CALL ASMERR (18,INOD,INOD,LPU,IERR)
         GO TO 100
      ENDIF
C
      NNDOF = MADOF(INOD+1) - MADOF(INOD)
      IDOF  = MADOF(INOD) - 1
C
      DO 50 LDOF=1,NNDOF
         IDOF = IDOF+1
         IEQ  = MEQN(IDOF)
         IF (IEQ.GT.0) THEN
C                                                ! free  d o f
            SV(IEQ) = SV(IEQ) + VNOD(LDOF)
         ELSEIF (IEQ.LT.0) THEN
C                                                ! slave  d o f
            ICEQ =-IEQ
            nmst = mpmceq(iceq+1) - mpmceq(iceq) - 1
            IF (NMST.GT.0) THEN
               IP = MPMCEQ(ICEQ)
               DO 25 II=1,NMST
                  IP  = IP+1
                  if (mmceq(ip).le.0) go to 25
                  ieq = meqn(mmceq(ip))
                  SV(IEQ) = SV(IEQ) + TTCC(IP)*VNOD(LDOF)
   25          CONTINUE
            ENDIF
         ENDIF
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE EXPAND (SVEQ,TTCC,MPMCEQ,MMCEQ,MEQN,FF,FS,NDOF,NEQ,SV)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EXPAND                  GROUP 7 / PUBLIC
C
C     T A S K :  TO EXPAND AND REARRANGE A SYSTEM VECTOR EXPRESSED IN
C     EQUATION ORDER, SVEQ, INTO A VECTOR EXPRESSED IN NODAL POINT
C     ORDER, MULTIPLY THIS BY A SCALAR AND ADD THE RESULT TO THE CON-
C     TENT OF AN INPUT VECTOR, SV.
C     ONE SCALAR (FF) IS USED FOR FREE  D O F S, ANOTHER (FS) FOR
C     SPECIFIED  D O F S.
C
C     ROUTINES CALLED/REFERENCED :   SLVEQ   (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-04-01 / 1.0
C                       89-01-13 / 1.1   K.Bell
C                       00-11-16 / 1.2   B.Haugen
C                       03-07-30 / 1.3   K.M.Okstad
C                       05-10-10 / 1.4   K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           NDOF,NEQ, MMCEQ(*),MEQN(*),MPMCEQ(*)
      DOUBLE PRECISION  FF,FS, TTCC(*),SV(*),SVEQ(*)
C
      INTEGER           I,ICEQ,IDOF,IEQ,IP,NMST
C
C ----------------------------------------------------------------------
      DO 100 IDOF=1,NDOF
         IEQ = MEQN(IDOF)
         IF (IEQ.EQ.0 .OR. IEQ.GT.NEQ) GO TO 100
         IF (IEQ.LT.0)    GO TO  20
         SV(IDOF) = SV(IDOF) + FF*SVEQ(IEQ)
         GO TO 100
C
   20    ICEQ =-IEQ
         IP   = MPMCEQ(ICEQ)
         SV(IDOF) = SV(IDOF) + FS*TTCC(IP)
         nmst = mpmceq(iceq+1) - mpmceq(iceq) - 1
         IF (NMST.EQ.0)   GO TO 100
         DO 50 I=1,NMST
            IP  = IP+1
            if (mmceq(ip).le.0) go to 50
            ieq = meqn(mmceq(ip))
            IF (IEQ.LT.1 .OR. IEQ.GT.NEQ) GO TO 50
            SV(IDOF) = SV(IDOF) + FF*TTCC(IP)*SVEQ(IEQ)
   50    CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE EXTEV (SV,MADOF,MPMNPC,MMNPC,
     +                  IEL,EV,NEDOF )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EXTEV                   GROUP 7 / PUBLIC
C
C     T A S K :  TO EXTRACT A VECTOR (EV), FOR ELEMENT IEL, FROM A
C     SYSTEM VECTOR (SV) IN WHICH THE ELEMENTS ARE ARRANGED IN NODAL
C     POINT ORDER  -  THE NUMBER OF ELEMENT  D O F S  (NEDOF) IS ALSO
C     DETERMINED
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IEL,NEDOF,  MADOF(*),MMNPC(*),MPMNPC(*)
      DOUBLE PRECISION  EV(*),SV(*)
C
      INTEGER           I,IE,INOD,IS,J,JE,JP,JS
C ----------------------------------------------------------------------
      JP    = 0
      IS    = MPMNPC(IEL)
      IE    = MPMNPC(IEL+1) - 1
      DO 100 I=IS,IE
         INOD = MMNPC(I)
         JS   = MADOF(INOD)
         JE   = MADOF(INOD+1) - 1
         DO 50 J=JS,JE
            JP     = JP+1
            EV(JP) = SV(J)
   50    CONTINUE
  100 CONTINUE
      NEDOF = JP
C
      RETURN
      END
      SUBROUTINE EXSEM (SM,MSKY,NEQ,NEQ2,LPU,SEM,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EXSEM                   GROUP 7 / PUBLIC
C
C     T A S K :  TO EXTRACT THE FULL, BUT SYMMETRIC SUPERELEMENT MATRIX
C     SEM  FROM THE PARTIALLY FACTORIZED TOTAL SUBSTRUCTURE MATRIX  SM.
C     THE SUBSTRUCTURE MATRIX (SM) IS A SYMMETRIC MATRIX STORED IN
C     'SKYLINE' FORM, AND THE PARTIAL FACTORIZATION HAS BEEN PERFORMED
C     BY  S A M  ROUTINE -SKYSOL-  (WITH  NEQ1 = abs(NEQ)-NEQ2)
C
C     ROUTINES CALLED/REFERENCED :  ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-29 / 1.0
C                       11-04-11 / 1.1   K.M.Okstad (added NEQ < 0)
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU,NEQ,NEQ2,   MSKY(*)
      DOUBLE PRECISION  SEM(NEQ2,NEQ2),SM(*)
C
      INTEGER           I,IP,I1,J,JJ
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR = 0
      IF (NEQ2.LT.1 .OR. NEQ2.GT.abs(NEQ)) GO TO 90
C
      if (NEQ.LT.0) THEN
         JJ = 1
      else
         JJ = NEQ-NEQ2+1
      end if
      IP = MSKY(JJ)
      SEM(1,1) = SM(IP)
      IF (NEQ2.EQ.1)       GO TO 100
      DO 50 J=2,NEQ2
         JJ = JJ+1
         I1 = MSKY(JJ-1) + 1
         IP = MSKY(JJ) - J
         DO 30 I=1,J
            IP = IP+1
            SEM(I,J) = SM(IP)
            IF (IP.LT.I1)  SEM(I,J) = ZERO
            SEM(J,I) = SEM(I,J)
   30    CONTINUE
   50 CONTINUE
      GO TO 100
C
   90 CALL ASMERR (13,NEQ2,NEQ2,LPU,IERR)
C
  100 RETURN
      END
      SUBROUTINE EXSEV (SV,NEQ,NEQ2,LPU,SEV,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EXSEV                   GROUP 7 / PUBLIC
C
C     T A S K :  TO EXTRACT THE SUPERELEMENT 'LOAD' VECTOR  SEV  FROM
C     THE PARTIALLY REDUCED SUBSTRUCTURE (LOAD) VECTOR  SV.
C     SV  IS PARTIALLY REDUCED BY  S A M  ROUTINE  -SKYSOL-
C
C     ROUTINES CALLED/REFERENCED :  ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU,NEQ,NEQ2
      DOUBLE PRECISION  SEV(NEQ2),SV(NEQ)
C
      INTEGER           I,II
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR = 0
      IF (NEQ2.LT.1 .OR. NEQ2.GT.NEQ)  GO TO 90
C
      II = NEQ-NEQ2
      DO 50 I=1,NEQ2
         II = II+1
         SEV(I) = SV(II)
   50 CONTINUE
      GO TO 100
C
   90 CALL ASMERR (13,NEQ2,NEQ2,LPU,IERR)
C
  100 RETURN
      END
      SUBROUTINE INSEV (SEV,NEQ,NEQ2,LPU,SV,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  INSEV                   GROUP 7 / PUBLIC
C
C     T A S K :  TO INSERT THE 'EXTERNAL' SUBSTRUCTURE DISPLACEMENT
C     VECTOR  SEV (= THE SUPERELEMENT DISPLACEMENT VECTOR) INTO THE
C     PARTIALLY REDUCED SUBSTRUCTURE  LOAD/DISPLACEMENT  VECTOR  SV
C
C     ROUTINES CALLED/REFERENCED :  ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU,NEQ,NEQ2
      DOUBLE PRECISION  SEV(NEQ2),SV(NEQ)
C
      INTEGER           I,II
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR = 0
      IF (NEQ2.LT.1 .OR. NEQ2.GT.NEQ)  GO TO 90
C
      II = NEQ-NEQ2
      DO 50 I=1,NEQ2
         II = II+1
         SV(II) = SEV(I)
   50 CONTINUE
      GO TO 100
C
   90 CALL ASMERR (13,NEQ2,NEQ2,LPU,IERR)
C
  100 RETURN
      END
      SUBROUTINE GETS22 (SM,MSKY,NEQ,NDOF2,LPU,IFLAG,S22,MSKY22,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  GETS22                  GROUP 7 / PUBLIC
C
C     T A S K :  TO EXTRACT THE UPPER TRIANGULAR PART (S22) OF THE
C     SYMMETRIC SUBMATRIX CORRESPONDING TO 2-STATUS DOFS, FROM A SYM-
C     METRIC SYSTEM MATRIX (SM) STORED IN 'SKYLINE' FORMAT DEFINED BY
C     MSKY.
C     S22 IS ALSO STORED IN A ONE-DIM. ARRAY, COLUMN BY COLUMN,
C       - IN ACTUAL SKYLINE FORMAT   IF  IFLAG.EQ.0
C       - OR AS A FULL MATRIX (WITH HORIZONTAL SKYLINE)  IF  IFLAG.NE.0
C     THE SKYLINE DEFINITION TABLE IS PREPARED (IRRESPECTIVE OF IFLAG)
C     AND RETURNED IN ARGUMENT MSKY22
C
C     ROUTINES CALLED/REFERENCED :  ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   85-12-12 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,IFLAG,LPU,NEQ,NDOF2,    MSKY(*),MSKY22(*)
      DOUBLE PRECISION  SM(*),S22(*)
C
      INTEGER           I,IE,IP,IS,J,JJ,LP
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR = 0
      IF (NDOF2.LT.1 .OR. NDOF2.GT.NEQ)  GO TO 90
C
      LP        = 1
      JJ        = NEQ-NDOF2+1
      IP        = MSKY(JJ)
      S22(LP)   = SM(IP)
      MSKY22(1) = 1
      IF (NDOF2.EQ.1)    GO TO 100
C
      DO 50 J=2,NDOF2
         JJ = JJ + 1
         IS = J - MSKY(JJ) + MSKY(JJ-1) + 1
         IF (IS.GT.J)    GO TO 95
         IF (IS.LT.1)    IS = 1
         IF (IFLAG.EQ.0) GO TO 20
         IF (IS.EQ.1)    GO TO 20
C
         IE = IS-1
         DO 10 I=1,IE
            LP = LP+1
            S22(LP) = ZERO
   10 CONTINUE
C
   20    IP = MSKY(JJ) - J + IS
         DO 30 I=IS,J
            LP = LP+1
            S22(LP) = SM(IP)
            IP = IP+1
   30    CONTINUE
         MSKY22(J) = LP
   50 CONTINUE
      GO TO 100
C
   90 CALL ASMERR (14,NEQ,NDOF2,LPU,IERR)
      GO TO 100
   95 CALL ASMERR (15,NEQ,NEQ,LPU,IERR)
C
  100 RETURN
      END
      SUBROUTINE GETS12 (SM,MSKY,NDOF1,NDOF2,NEQ,LPU,S12,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY :  GETS12                          GROUP 7 / PUBLIC
C
C     T A S K :  TO EXTRACT THE RECTANGULAR SUBMATRIX (S12) CORRE-
C     SPONDING TO 1-STATUS  D O F S  ROW-WISE AND 2-STATUS  D O F S
C     COLUMN-WISE, FROM A SYMMETRIC SYSTEM MATRIX (SM) STORED IN 'SKY-
C     LINE' FORMAT DEFINED BY MSKY.
C     S12 IS RETURNED AS A FULL TWO-DIM. (NDOF1 BY NDOF2) ARRAY.
C
C     ROUTINES CALLED/REFERENCED :  ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   85-12-12 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU,NEQ,NDOF1,NDOF2,    MSKY(*)
      DOUBLE PRECISION  SM(*),S12(NDOF1,NDOF2)
C
      INTEGER           I,IP,IS,J,JJ
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR = 0
      IF (NDOF1.LT.1 .OR. NDOF2.LT.1)  GO TO 90
      IF (NDOF1+NDOF2.NE.NEQ)          GO TO 90
C
      DO 20 J=1,NDOF2
         DO 10 I=1,NDOF1
            S12(I,J) = ZERO
   10    CONTINUE
   20 CONTINUE
C
      JJ = NDOF1
      DO 50 J=1,NDOF2
         JJ = JJ + 1
         IS = JJ - MSKY(JJ) + MSKY(JJ-1) + 1
         IF (IS.LT.1 .OR. IS.GT.JJ)  GO TO 95
         IF (IS.GT.NDOF1)            GO TO 50
         IP = MSKY(JJ-1) + 1
         DO 40 I=IS,NDOF1
            S12(I,J) = SM(IP)
            IP = IP+1
   40    CONTINUE
   50 CONTINUE
      GO TO 100
C
   90 CALL ASMERR (16,NDOF1,NDOF2,LPU,IERR)
      GO TO 100
   95 CALL ASMERR (15,NEQ,NEQ,LPU,IERR)
C
  100 RETURN
      END
      SUBROUTINE ASMDOF (MADOF,MINEX,NANOD,IDOF,LDOF,INOD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ASMDOF                 GROUP 7 / PRIVATE
C
C     T A S K :  To determine the node identifier and local (node) dof
C                number corresponding to (global) dof number IDOF
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-03-13 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IDOF,INOD,LDOF,NANOD,
     &                  MADOF(*),MINEX(*)
C                                                ! Local variables
      INTEGER           I
C ----------------------------------------------------------------------
C
      INOD = 0
      LDOF = 0
      DO 10 I=1,NANOD
         IF (IDOF.LT.MADOF(I+1)) THEN
            INOD = MINEX(I)
            LDOF = IDOF - MADOF(I) + 1
            GO TO 100
         ENDIF
   10 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE ASMEQN (MPAR,MADOF,MINEX,MEQN,IEQ,LDOF,INOD,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ASMEQN                 GROUP 7 / PUBLIC
C
C     T A S K :  To determine the node identifier and local (node) dof
C                number corresponding to a particular equation number
C                IEQ
C
C
C     ROUTINES CALLED/REFERENCED :  INDXTV     (SAM-8)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   97-02-15 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IEQ,IERR,INOD,LDOF,
     &                  MADOF(*),MEQN(*),MINEX(*),MPAR(50)
C
C                                                ! Local variables
      INTEGER           I,IDOF,N,     INDXTV
      EXTERNAL          INDXTV
C ----------------------------------------------------------------------
C
      INOD = 0
      N    = MPAR(3)
      IDOF = INDXTV(MEQN,IEQ,N)
      IF (IDOF.EQ.0) THEN
         IERR = -1
         GO TO 100
      ENDIF
C
      N = MPAR(1) + 1
      DO 20 I=2,N
         IF (IDOF.LT.MADOF(I)) THEN
            INOD = MINEX(I-1)
            LDOF = IDOF - MADOF(I-1) + 1
            GO TO 100
         ENDIF
   20 CONTINUE
      IF (IDOF.EQ.0)  IERR = -1
C
  100 RETURN
      END
      SUBROUTINE SYSEQ (MSC,MPMCEQ,MMCEQ,LPU,MPAR,MEQN,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SYSEQ                   GROUP 7 / PRIVATE
C
C     T A S K :  TO CHECK CONTROL MATRICES  MSC  AND  MMCEQ, DETERMINE
C     PARAMETERS  NDOF1, NDOF2, NSPDOF, NSDOF, NPDOF, NDDOF, NEQ  AND
C     NMMCEQ  AND STORE THEM IN  MPAR  AND FINALLY, DETERMINE CONTROL
C     MATRIX  MEQN
C
C     ROUTINES CALLED/REFERENCED :  ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-26 / 1.0
C                       87-04-04 / 1.1
C                       87-08-25 / 1.2
C                       89-01-13 / 1.3   K.Bell
C                       00-11-13 / 1.4   B.Haugen
C                       04-01-15 / 1.5   K.M.Okstad
C                       05-10-10 / 1.6   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU
      INTEGER           MMCEQ(*),MEQN(*),MPAR(*),MPMCEQ(*),MSC(*)
C
      INTEGER           I,IDOF,IP,J,JDOF,NCEQ,NDDOF,NDOF,NDOF1,NDOF2,
     +                  NM,NPDOF,NSPDOF
C
      EXTERNAL          ASMERR
C ----------------------------------------------------------------------
      IERR   = 0
      NDOF   = MPAR(3)
      NCEQ   = MPAR(7)
C
      NDOF1  = 0
      NDOF2  = 0
      NSPDOF = 0
      NPDOF  = 0
      NDDOF  = 0
C ----------------------------------------------------------------------
C  CHECK  MSC  AND DETERMINE  NDOF1, NDOF2 AND NSPDOF
C ----------------------------------------------------------------------
      DO 100 I=1,NDOF
         IF (MSC(I).LT.0)  CALL ASMERR (1,I,I,LPU,IERR)
         IF (MSC(I).EQ.0)  THEN
            NSPDOF  = NSPDOF+1
            MEQN(I) = 0
         ENDIF
         IF (MSC(I).EQ.1)  NDOF1  = NDOF1 +1
         IF (MSC(I).EQ.2)  NDOF2  = NDOF2 +1
         IF (MSC(I).GT.2)  CALL ASMERR (1,I,I,LPU,IERR)
  100 CONTINUE
      IF (IERR.LT.0)       GO TO 1000
      MPAR( 4) = NDOF1
      MPAR( 5) = NDOF2
      MPAR( 6) = NSPDOF
      MPAR( 8) = NSPDOF-NCEQ
      MPAR(11) = NDOF1 +NDOF2
C
      MPAR(16) = 0
      IF (NCEQ.EQ.0)       GO TO 300
C ----------------------------------------------------------------------
C  CHECK  MMCEQ, DETERMINE  NPDOF AND NDDOF, AND SET NEGATIVE ELEMENTS
C  IN MEQN
C ----------------------------------------------------------------------
      DO 200 I=1,NCEQ
         NM   = MPMCEQ(I+1) - MPMCEQ(I) - 1
         IP   = MPMCEQ(I)
         IDOF = MMCEQ(IP)
C
         IF (IDOF.GT.0 .AND. IDOF.LE.NDOF) GO TO 120
         CALL ASMERR (3,IDOF,IDOF,LPU,IERR)
         GO TO 200
C
  120    IF (MSC(IDOF).EQ.0) GO TO 130
         CALL ASMERR (4,IDOF,IDOF,LPU,IERR)
         GO TO 200
C
  130    MEQN(IDOF) = -I
         IF (NM.EQ.0)      GO TO 160
         IF (NM.GT.0)      GO TO 140
         CALL ASMERR (5,I,I,LPU,IERR)
         GO TO 200
C                                               ** SLAVE  D O F
  140    NDDOF = NDDOF+1
         DO 150 J=1,NM
            IP   = IP+1
            JDOF = MMCEQ(IP)
            IF (JDOF.LE.0) GO TO 150
            IF (MSC(JDOF).LT.1)  CALL ASMERR (6,IDOF,IDOF,LPU,IERR)
  150    CONTINUE
         GO TO 200
C                                               ** PRESCRIBED  D O F
  160    NPDOF = NPDOF+1
  200 CONTINUE
      IF (IERR.LT.0)       GO TO 1000
      MPAR(16) = MPMCEQ(NCEQ+1) - 1
C
  300 MPAR( 9) = NPDOF
      MPAR(10) = NDDOF
C ----------------------------------------------------------------------
C  COMPLETE  MEQN
C ----------------------------------------------------------------------
      I = 1
      J = NDOF1+1
      DO 400 IDOF=1,NDOF
         IF (MSC(IDOF).EQ.1)  GO TO 320
         IF (MSC(IDOF).EQ.2)  GO TO 340
         GO TO 400
  320    MEQN(IDOF) = I
         I = I+1
         GO TO 400
  340    MEQN(IDOF) = J
         J = J+1
  400 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE SYSEQX (MADOF,  MINEX, MSC,
     &                   MPMCEQ, MMCEQ, MNNN,  LPU,
     &                   MPAR,   MEQN,  MWORK, IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SYSEQX                 GROUP 7 / PRIVATE
C
C     TASK :  To check control matrices MSC and and MMCEQ, to determine
C             parameters NDOF1, NDOF2, NSPDOF, NSDOF, NPDOF, NDDOF, NEQ
C             and NMMCEQ and store them in MPAR and, finally, determine
C             control matrix MEQN with due consideration to "new" node
C             numbers (recorded in MNNN).
C             This is a modified version of SYSEQ (SAM-7).
C
C     ROUTINES CALLED/REFERENCED :  ASMDOF and ASMERR    (SAM-7)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-03-13 / 1.0
C                       04-01-15 / 1.1   K.M.Okstad
C                       05-10-10 / 1.2   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,LPU
      INTEGER           MADOF(*),MEQN(*),MINEX(*),MMCEQ(*),MNNN(*),
     &                  MPAR(*),MPMCEQ(*),MSC(*),MWORK(*)
C
      INTEGER           I,IDOF,IE,IEQ1,IEQ2,INOD,IP,IS,
     &                  J,JDOF,KNOD,LDOF,NANOD,NCEQ,
     &                  NDDOF,NDOF,NDOF1,NDOF2,NM,NPDOF,NSPDOF
C
C ----------------------------------------------------------------------
      IERR   = 0
      NANOD  = MPAR(1)
      NDOF   = MPAR(3)
      NCEQ   = MPAR(7)
C
      NDOF1  = 0
      NDOF2  = 0
      NSPDOF = 0
      NPDOF  = 0
      NDDOF  = 0
C
C --- Check  MSC  and determine  NDOF1, NDOF2 and NSPDOF
C
      DO 100 I=1,NDOF
         MEQN(I) = 0
         IF (MSC(I).EQ.0) THEN
            NSPDOF  = NSPDOF+1
         ELSEIF (MSC(I).EQ.1) THEN
            NDOF1  = NDOF1 +1
         ELSEIF (MSC(I).EQ.2) THEN
            NDOF2  = NDOF2 +1
         ELSE
            CALL ASMDOF (MADOF,MINEX,NANOD,I,LDOF,INOD)
            CALL ASMERR (21,LDOF,INOD,LPU,IERR)
         ENDIF
  100 CONTINUE
      IF (IERR.LT.0)  GO TO 1000
      MPAR( 4) = NDOF1
      MPAR( 5) = NDOF2
      MPAR( 6) = NSPDOF
      MPAR( 8) = NSPDOF-NCEQ
      MPAR(11) = NDOF1 +NDOF2
C
      MPAR(16) = 0
C
      IF (NCEQ.GT.0) THEN
C
C ------ Check MMCEQ, determine NPDOF and NDDOF, and set negative
C        elements in MEQN
C
         DO 200 I=1,NCEQ
            NM   = MPMCEQ(I+1) - MPMCEQ(I) - 1
            IP   = MPMCEQ(I)
            IDOF = MMCEQ(IP)
C
            IF (IDOF.LT.1 .OR. IDOF.GT.NDOF) THEN
               CALL ASMERR (23,IDOF,IDOF,LPU,IERR)
               GO TO 200
            ENDIF
C
            IF (MSC(IDOF).NE.0) THEN
               CALL ASMDOF (MADOF,MINEX,NANOD,IDOF,LDOF,INOD)
               CALL ASMERR (24,LDOF,INOD,LPU,IERR)
               GO TO 200
            ENDIF
C
            MEQN(IDOF) = -I
            IF (NM.GT.0) THEN
C                                                ! slave  d o f
               NDDOF = NDDOF+1
               DO 150 J=1,NM
                  IP   = IP+1
                  JDOF = MMCEQ(IP)
                  IF (JDOF.LE.0) GO TO 150
                  IF (MSC(JDOF).LT.1) THEN
                     CALL ASMDOF (MADOF,MINEX,NANOD,IDOF,LDOF,INOD)
                     CALL ASMERR (26,LDOF,INOD,LPU,IERR)
                  ENDIF
  150          CONTINUE
            ELSEIF (NM.EQ.0) THEN
C                                                ! prescribed  d o f
               NPDOF = NPDOF+1
            ELSE
               CALL ASMERR (25,I,I,LPU,IERR)
            ENDIF
  200    CONTINUE
         IF (IERR.LT.0)       GO TO 1000
         MPAR(16) = MPMCEQ(NCEQ+1) - 1
      ENDIF
C
      MPAR( 9) = NPDOF
      MPAR(10) = NDDOF
C
C --- Establish a temporary table (in MWORK) indicating the "old" node
C     number corresponding to the new number (which is the index in
C     MWORK)
C
      DO 300 KNOD=1,NANOD
         INOD = MNNN(KNOD)
         IF (INOD.GT.NANOD) THEN
            IE = MINEX(KNOD)
            CALL ASMERR (28,IE,INOD,LPU,IERR)
            GO TO 1000
         ENDIF
         MWORK(INOD) = KNOD
  300 CONTINUE
C
C --- Complete  MEQN
C
      IEQ1 = 1
      IEQ2 = NDOF1 + 1
      DO 400 KNOD=1,NANOD
         INOD = MWORK(KNOD)
         IS   = MADOF(INOD)
         IE   = MADOF(INOD+1) - 1
         DO 350 IDOF=IS,IE
            IF (MSC(IDOF).EQ.1) THEN
               IF (MEQN(IDOF).EQ.0) THEN
                  MEQN(IDOF) = IEQ1
                  IEQ1       = IEQ1 + 1
               ELSE
                  CALL ASMERR (29,I,I,LPU,IERR)
                  GO TO 1000
               ENDIF
            ELSEIF (MSC(IDOF).EQ.2) THEN
               IF (MEQN(IDOF).EQ.0) THEN
                  MEQN(IDOF) = IEQ2
                  IEQ2       = IEQ2 + 1
               ELSE
                  CALL ASMERR (29,I,I,LPU,IERR)
               ENDIF
            ENDIF
  350    CONTINUE
  400 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE SKYLIN (MADOF,MPMNPC,MMNPC,MPMCEQ,MMCEQ,MEQN,
     +                   MPAR,MEEN,MSKY )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SKYLIN                  GROUP 7 / PRIVATE
C
C     T A S K :  DETERMINE THE 'SKYLINE' OF A SYMMETRIC SYSTEM MATRIX
C     (STORED ACTIVE COLUMN BY ACTIVE COLUMN), I.E., DETERMINE THE
C     ELEMENTS OF CONTROL MATRIX  MSKY  -  THE NUMBER OF ELEMENTS UNDER
C     THE SKYLINE, NESKY, IS SET IN  MPAR(12)
C
C     ROUTINES CALLED/REFERENCED :  ELMEQ AND SLVEQ   (SAM-7)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-25 / 1.0
C                       89-01-13 / 1.1   K.Bell
C                       00-11-13 / 1.2   B.Haugen
C                       04-01-15 / 1.3   K.M.Okstad
C                       05-10-10 / 1.4   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           MADOF(*),MMCEQ(*),MEEN(*),MEQN(*),MMNPC(*),
     +                  MPAR(*),MPMCEQ(*),MPMNPC(*),MSKY(*)
C
      INTEGER           I,ICEQ,IEL,IEQ,IP,J,LEQ,NEDOF,NEL,NMST,NENOD,
     +                  NEPRD,NEQ,NESLV
C
      EXTERNAL          ELMEQ
C ----------------------------------------------------------------------
      NEL = MPAR( 2)
      NEQ = MPAR(11)
C                                               ** INITIATE  MSKY
      DO 100 I=1,NEQ
         MSKY(I) = I
  100 CONTINUE
C ----------------------------------------------------------------------
C  FIND EQ. NO. OF FIRST NON-ZERO ELEMENT OF EACH COLUMN (STORE IN MSKY)
C ----------------------------------------------------------------------
      DO 500 IEL=1,NEL
         NENOD = MPMNPC(IEL+1) - MPMNPC(IEL)
         IP    = MPMNPC(IEL)
C                                               ** DETERMINE  MEEN
         CALL ELMEQ (MADOF,MMNPC(IP),MPMCEQ,MEQN,NENOD,
     +               MEEN,NEDOF,NESLV,NEPRD)
C                                               ** LOWEST EQ.NO.  LEQ
         LEQ = NEQ
         DO 200 I=1,NEDOF
            IF (MEEN(I).LT.1)        GO TO 200
            IF (MEEN(I).LT.LEQ)      LEQ = MEEN(I)
  200    CONTINUE
C                                               ** CHECK  MSKY
         DO 250 I=1,NEDOF
            IEQ = MEEN(I)
            IF (IEQ.LT.1)            GO TO 250
            IF (MSKY(IEQ).GT.LEQ)    MSKY(IEQ) = LEQ
  250    CONTINUE
C
         IF (NESLV.EQ.0)             GO TO 500
C                                               ** MODIFY FOR LINEAR
C                                                  CONSTRAINT EQS.
         DO 300 I=1,NEDOF
            IF (MEEN(I).GE.0)        GO TO 300
            ICEQ = -MEEN(I)
            nmst = mpmceq(iceq+1) - mpmceq(iceq) - 1
            IF (NMST.EQ.0)           GO TO 300
C                                               ** CHECK/ADJUST  LEQ
            ip = mpmceq(iceq)
            DO 280 J=1,NMST
               ip  = ip+1
               if (mmceq(ip) .le. 0) go to 280
               ieq = meqn(mmceq(ip))
               IF (IEQ.LT.LEQ)       LEQ = IEQ
  280       CONTINUE
  300    CONTINUE
C                                               ** CHECK  MSKY  AGAIN
         DO 400 I=1,NEDOF
            IF (MEEN(I).EQ.0)        GO TO 400
            IF (MEEN(I).GT.0)        GO TO 350
            ICEQ = -MEEN(I)
            nmst = mpmceq(iceq+1) - mpmceq(iceq) - 1
            IF (NMST.EQ.0)           GO TO 400
            ip = mpmceq(iceq)
            DO 320 J=1,NMST
               ip  = ip+1
               if (mmceq(ip) .le. 0) go to 320
               ieq = meqn(mmceq(ip))
               IF (MSKY(IEQ).GT.LEQ) MSKY(IEQ) = LEQ
  320       CONTINUE
            GO TO 400
  350       IEQ = MEEN(I)
            IF (MSKY(IEQ).GT.LEQ)    MSKY(IEQ) = LEQ
  400    CONTINUE
  500 CONTINUE
C
      IF (NEQ.EQ.1)                  GO TO 700
C ----------------------------------------------------------------------
C  CONVERT CURRENT CONTENT OF  MSKY  TO FINAL CONTENT
C ----------------------------------------------------------------------
      DO 600 I=2,NEQ
         J       = I - MSKY(I) + 1
         MSKY(I) = MSKY(I-1) + J
  600 CONTINUE
C
  700 MPAR(12) = MSKY(NEQ)
C ----------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE ELMEQ (MADOF,MNPC,MPMCEQ,MEQN,NENOD,
     +                  MEEN,NEDOF,NESLV,NEPRD )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ELMEQ                   GROUP 7 / PRIVATE
C
C     T A S K :  TO FIND  MEEN ('MATRIX OF' ELEMNET EQUATION NUMBERS),
C     NEDOF (NO. OF ELEMENT  D O F S), NESLV (NO. OF ELEMENT SLAVE
C     D O F S) AND NEPRD (NO. OF ELEMENT PRESCRIBED  D O F S) FOR AN
C     ELEMENT WHOSE  MNPC(NENOD)  IS INPUT
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-25 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           NEDOF,NENOD,NEPRD,NESLV
      INTEGER           MADOF(*),MEEN(*),MEQN(*),MNPC(NENOD),MPMCEQ(*)
C
      INTEGER           I,ICEQ,IDOF,IE,INOD,IS,NMST,NN
C ----------------------------------------------------------------------
      IDOF  = 0
      NESLV = 0
      NEPRD = 0
      DO 100 INOD=1,NENOD
         NN = MNPC(INOD)
         IS = MADOF(NN)
         IE = MADOF(NN+1) - 1
         DO 50 I=IS,IE
            IDOF       = IDOF+1
            MEEN(IDOF) = MEQN(I)
            IF (MEQN(I).GE.0)  GO TO 50
C
            ICEQ =-MEQN(I)
            NMST = MPMCEQ(ICEQ+1) - MPMCEQ(ICEQ) - 1
C                                               ** PRESCRIBED  D O F
            IF (NMST.EQ.0)     NEPRD = NEPRD+1
C                                               ** DEPENDENT  D O F
            IF (NMST.GT.0)     NESLV = NESLV+1
   50    CONTINUE
  100 CONTINUE
C
      NEDOF = IDOF
C
      RETURN
      END
      SUBROUTINE SLVEQ (MPMCEQ,MMCEQ,MEQN,ICEQ,MMST,NMST)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SLVEQ                   GROUP 7 / PRIVATE
C
C     T A S K :  TO FIND AND RETURN THE NO. OF MASTERS (NMST) AND THEIR
C     EQUATION NOS. (LISTED IN MMST) FOR CONSTRAINT EQUATION NO. ICEQ
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-25 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           ICEQ,NMST,  MMCEQ(*),MEQN(*),MMST(*),MPMCEQ(*)
C
      INTEGER           I,IDOF,IP
C ----------------------------------------------------------------------
      NMST = MPMCEQ(ICEQ+1) - MPMCEQ(ICEQ) - 1
      IF (NMST.LT.1)  GO TO 100
      IP   = MPMCEQ(ICEQ)
      DO 50 I=1,NMST
         IP      = IP+1
         IDOF    = MMCEQ(IP)
         MMST(I) = MEQN(IDOF)
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE BANDW (MADOF,MEQN,MSKY,MPAR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  BANDW                   GROUP 7 / PRIVATE
C
C     T A S K :  TO DETERMINE THE LARGEST SEMI-BANDWIDTH  NLBW  AND THE
C     CORRESPONDING NODE NUMBER NODLBW, OF THE FIRST NEQ1 EQUATIONS,
C     AND STORE THEM IN  MPAR
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-03-25 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           MADOF(*),MEQN(*),MPAR(*),MSKY(*)
C
      INTEGER           I,IDF,IDOF,IEQ,INOD,NLBW,NDOF,NEQ1,NANOD
C ----------------------------------------------------------------------
      NANOD = MPAR(1)
      NDOF  = MPAR(3)
      NEQ1  = MPAR(4)
      IF (NEQ1.GT.0)             GO TO 10
      NLBW  = 0
      INOD  = 0
      GO TO 100
   10 NLBW  = 1
      IEQ   = 1
      IF (NEQ1.EQ.1)             GO TO 30
C
      DO 20 I=2,NEQ1
         IDF = MSKY(I)-MSKY(I-1)
         IF (IDF.LE.NLBW)        GO TO 20
         NLBW = IDF
         IEQ = I
   20 CONTINUE
C
   30 DO 40 IDOF=1,NDOF
         IF (MEQN(IDOF).EQ.IEQ)  GO TO 60
   40 CONTINUE
      IDOF = NDOF
C
   60 DO 80 I=1,NANOD
         IF (IDOF.GE.MADOF(I))   GO TO 80
         INOD = I-1
         GO TO 100
   80 CONTINUE
      INOD = NANOD
C
  100 MPAR(13) = NLBW
      MPAR(14) = INOD
C
      RETURN
      END
      INTEGER FUNCTION NOPRED(MPMCEQ,I)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NOPRED                  GROUP 7 / PRIVATE
C
C     T A S K :  TO DETERMINE IF THE I'TH CONSTRAINT EQUATION DEFINES A
C     PRESCRIBED  D O F,  AND IF SO WHICH NUMBER OF PRESCRIBED  D O F
C     IT IS  (COUNTING FROM THE START OF MPMCEQ/MMCEQ)
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-01-14 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           I,MPMCEQ(*)
      INTEGER           K,N
C ----------------------------------------------------------------------
      N = 0
      IF (I.LT.1)                           GO TO 20
      IF ((MPMCEQ(I+1)-MPMCEQ(I)).GT.1)     GO TO 20
C
      DO 10 K=1,I
         IF ((MPMCEQ(K+1)-MPMCEQ(K)).EQ.1)  N = N+1
   10 CONTINUE
C
   20 NOPRED = N
C
      RETURN
      END
      SUBROUTINE ASMERR (NUM,IVAL1,IVAL2,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ASMERR                  GROUP 7 / PRIVATE
C
C     T A S K :  TO PRINT ERROR MESSAGES AND DECREMENT THE ERROR FLAG
C     FOR OTHER ROUTINES IN THE  A S M  PACKAGE
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   85-12-12 / 1.0
C                       94-03-14 / 1.1   K.Bell (major rewrite)
C
C **********************************************************************
C
      IMPLICIT          NONE
C
      INTEGER           IERR,IVAL1,IVAL2,LPU,NUM
C ----------------------------------------------------------------------
      IF (LPU.LE.0)  GO TO 100
C
      WRITE (LPU,600)
C
      IF (NUM.LE.10) THEN
         IF (NUM.EQ. 1)  WRITE (LPU,601)  IVAL1
         IF (NUM.EQ. 2)  WRITE (LPU,602)  IVAL1,IVAL2
         IF (NUM.EQ. 3)  WRITE (LPU,603)  IVAL1
         IF (NUM.EQ. 4)  WRITE (LPU,604)  IVAL1
         IF (NUM.EQ. 5)  WRITE (LPU,605)
         IF (NUM.EQ. 6)  WRITE (LPU,606)  IVAL1
         IF (NUM.EQ. 7)  WRITE (LPU,607)  IVAL1
         WRITE (LPU,650)
      ELSEIF (NUM.LE.20) THEN
         IF (NUM.EQ.11) THEN
            WRITE (LPU,611)  IVAL1
            WRITE (LPU,651)
         ENDIF
         IF (NUM.EQ.12) THEN
            WRITE (LPU,611)  IVAL1
            WRITE (LPU,652)
         ENDIF
         IF (NUM.EQ.13)  WRITE (LPU,613)  IVAL1
         IF (NUM.EQ.14)  WRITE (LPU,614)  IVAL1,IVAL2
         IF (NUM.EQ.15)  WRITE (LPU,615)
         IF (NUM.EQ.16)  WRITE (LPU,616)  IVAL1,IVAL2
         IF (NUM.EQ.17)  WRITE (LPU,617)  IVAL2,IVAL1
         IF (NUM.EQ.18)  THEN
            WRITE (LPU,618)  IVAL1
            WRITE (LPU,653)
         ENDIF
      ELSEIF (NUM.LE.30) THEN
         IF (NUM.EQ.21)  WRITE (LPU,621)  IVAL1,IVAL2
         IF (NUM.EQ.22)  WRITE (LPU,622)  IVAL1,IVAL2
         IF (NUM.EQ.23)  WRITE (LPU,623)  IVAL1
         IF (NUM.EQ.24)  WRITE (LPU,624)  IVAL1,IVAL2
         IF (NUM.EQ.25)  WRITE (LPU,625)
         IF (NUM.EQ.26)  WRITE (LPU,626)  IVAL1,IVAL2
         IF (NUM.EQ.27)  WRITE (LPU,627)  IVAL1
         IF (NUM.EQ.28)  WRITE (LPU,628)  IVAL1,IVAL2
         IF (NUM.EQ.29)  WRITE (LPU,629)
         WRITE (LPU,660)
      ENDIF
C
  100 IERR = IERR-1
#ifdef FT_DEBUG
      IERR = IERR/0 ! Force a core dump for easy debugging
#endif
C
      RETURN
C ----------------------------------------------------------------------
  600 FORMAT(///' *** ERROR return from a SAM library routine')
C
  601 FORMAT(/5X,'Illegal status code encountered for  d o f  no.',I8)
  602 FORMAT(/5X,'D o f  no.',I8,' is coupled to too many (=',I4,
     +           ' master  d o f s')
  603 FORMAT(/5X,'Illegal  d o f  (no.=',I8,' ) is constrained')
  604 FORMAT(/5X,'D o f  no.',I8,' is constrained but not specified',
     +           ' (in MSC)')
  605 FORMAT(/5X,'Control matrix  MPMCEQ  in error')
  606 FORMAT(/5X,'D o f  no.',I8,' is coupled to a specified  d o f')
  607 FORMAT(/5X,'Illegal or inconsistent no. of  d o f s (=',I8,' )')
C
  611 FORMAT(/5X,'Inconsistent control information encountered'
     +       /5X,'during assembly of element no.',I8 )
  613 FORMAT(/5X,'Incorrect superelement dimension  (=',I8,' )'/
     +       /5X,'Detected by subroutine  EXSEM, EXSEV or INSEV')
  614 FORMAT(/5X,'Incorrect dimension (NEQ =',I8,'  or  NDOF2 =',I8,')'/
     +       /5X,'Detected by subroutine GETS22' )
  615 FORMAT(/5X,'Incorrect profile definition  (MSKY in error)'/
     +       /5X,'Detected by subroutine GETS22 OR GETS12' )
  616 FORMAT(/5X,'Incorrect dimension parameter(s)'/
     +        5X,'NEQ  and/or  NDOF1 =',I8,'  and/or  NDOF2 =',I8,/
     +       /5X,'Detected by subroutine GETS12' )
  617 FORMAT(/5X,'Incorrect value of argument NSV (=',I8,')'/
     +        5X,'when calling  ADDEM  for element ',I8 )
  618 FORMAT(/5X,'Illegal node number (=',I8,' ) or inconsistent'
     +       /5X,'control information encountered during assembly')
C
  621 FORMAT(/5X,'Illegal status code encontered for  d o f  no.',I3,
     +           ' of node',I8 )
  622 FORMAT(/5X,'D o f  no.',I3,' of node',I8,' is coupled to too ',
     +           'many master  d o f s')
  623 FORMAT(/5X,'Nonexisting  d o f (no.',I8,') is constrained')
  624 FORMAT(/5X,'D o f  no.',I3,' of node',I8,' is constrained,'/
     +        5X,'but not specified (in MSC)')
  625 FORMAT(/5X,'Control matrix MPMCEQ appears to be in error')
  626 FORMAT(/5X,'D o f  no.',I3,' of node',I8,' is coupled to a ',
     +           'specified  d o f' )
  627 FORMAT(/5X,'Illegal or inconsistent number of  d o f s (=',I8,')')
  628 FORMAT(/5X,'Active node',I8,' has been assigned to an internal'/
     +        5X,'number (=',I8,' ) which is too large')
  629 FORMAT(/5X,'Control matrix MNNN (new node numbers) appears'/
     +        5X,'to be incorrect/corrupt' )
C
  650 FORMAT(/5X,'Detected by subroutine PRASS, directly or indirectly')
  651 FORMAT(/5X,'Detected by subroutine ADDEM')
  652 FORMAT(/5X,'Detected by subroutine ADDEV')
  653 FORMAT(/5X,'Detected by subroutine ADDNV')
C
  660 FORMAT(/5X,'Detected by subroutine PRASSX, directly or ',
     +           'indirectly')
      END
