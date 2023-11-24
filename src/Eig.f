C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE LANCZ1 (A,B,TOL,MSKY,MIP,
     +                   MOP,EVL,EVERR,V,R,RWA,LWA,
     +                   SHIFT,N,NEVAL,NEVEC,MAXNIV,LPU,IPSW,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LANCZ1                  GROUP 5 / PUBLIC
C
C     T A S K :  TO DETERMINE THE NEVAL EIGENVALUES LAMBDAi OF THE
C                GENERALIZED, SYMMETRIC EIGENPROBLEM
C                  (A - LAMBDAi*B)*Qi = 0
C                THAT ARE ALGEBRAICALLY CLOSEST TO A SPECIFIED SHIFT
C                OF ORIGON.
C                THE ROUTINE ALSO RETURNS NIV (=MOP(3)) ORTHOGONAL
C                LANCZOS VECTORS (IN V) OR, OPTIONALLY, APPROXIMATIONS
C                TO NEVEC EIGENVECTORS (CORRESPONDING TO THE NEVEC FIRST
C                EIGENVALUES).  THE EIGENVECTORS, WHICH ARE B-ORTHOGONAL
C                ARE RETURNED (AS THE FIRST VECTORS) IN V.
C     A IS A SYMMETRIC SKYLINE MATRIX, AND B IS A SYMMETRIC SKYLINE
C     MATRIX (SAME SKYLINE AS A) OR A DIAGONAL MATRIX.  FOR BOTH A AND B
C     THE ORIGINAL MATRICES OR THEIR FACTORS ARE ACCEPTED AS INPUT ( IN
C     ARRAYS A AND B).
C     A TRUNCATED LANCZOS ITERATION, ASSUMING ALL MATRICES IN PRIMARY
C     STORAGE, IS USED TO TRANSFORM THE PROBLEM TO TRIDIAGONAL, SPECIAL
C     FORM.  THE SMALL, TRIDIAGONAL EIGENPROBLEM IS SOLVED BY IMPLICIT
C     QL TRANSFORMATION.
C
C     ROUTINES CALLED/REFERENCED :  IMP, RMINT, IMINT, PRACC    (SAM-0)
C                                   SCALE, VASCA, ICOPY, RCOPY  (SAM-0)
C                                   BIGMUL                      (SAM-3)
C                                   BIGMUL, NORVEC              (SAM-3)
C                                   QLS3D, QLS3DV, EVARR        (SAM-5)
C                                   MDFCTA, BLBLBT, SPSMUL      (SAM-5)
C                                   VRORT, ABRLER, TOLCHK       (SAM-5)
C                                   BCKSUB, LCZPR1, LCZPR2      (SAM-5)
C                                   LCZPR3, LCZER1              (SAM-5)
C                                   RANVEC, REARRV              (SAM-8)
C                                   ABS, SQRT         (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-04-01 / 1.0
C                       01-05-24 / 1.1  K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IPSW,LPU,MAXNIV,N,NEVAL,NEVEC,
     +                  LWA(*),MIP(5),MOP(8),MSKY(*)
      DOUBLE PRECISION  SHIFT, A(*),B(*),EVERR(MAXNIV),EVL(MAXNIV),
     +                  R(*),RWA(*),TOL(3),V(N,MAXNIV)
C
      INTEGER           I,IP,IPA,IPB,IPD,IPE,IPS,J,K,KEX,KODEX,KOUNT,
     +                  KSA,KSB,KSR,LLP,LPC,LPP,MAXIT,NACEV,NEVEXT,
     +                  NGDEV,NIV,NN,NREP,NSDIG,NV
      INTEGER           IMP
      DOUBLE PRECISION  ALFA,AUX,BETA,EPS1,EPS2,EPS3,ONE,RAN,SMALL,
     +                  TEN,TINY,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 , TEN = 10.0D0 )
C
      EXTERNAL          ABRLER,BCKSUB,BIGMUL,BLBLBT,EVARR,ICOPY,IMINT,
     +                  IMP,LCZER1,LCZPR1,LCZPR2,LCZPR3,MDFCTA,NORVEC,
     +                  PRACC,QLS3D,QLS3DV,RANVEC,RCOPY,REARRV,RMINT,
     +                  SCALE,SPSMUL,TOLCHK,VASCA,VRORT
C
C ----------------------------------------------------------------------
C   CHECK INPUT PARAMETERS
C ----------------------------------------------------------------------
      IF (IPSW.GE.1.AND.NEVAL.GT.0) THEN
         CALL LCZPR1 (MSKY,MIP,SHIFT,N,NEVAL,NEVEC,MAXNIV,LPU)
      ENDIF
C
      IERR  = 0
      LLP   = LPU
      IF (IPSW.LT.0)  LLP =-1
      KSA   = MIP(1)
      KSB   = MIP(2)
      KSR   = ABS(MIP(3))
      MAXIT = MIP(4)
      KEX   = MIP(5)
      IF (KEX.LT.1.OR.KEX.GT.5)  KEX = 2
C
      IF (N.LT.2)           CALL LCZER1 (11,N,N,N,LLP,IERR)
      IF (ABS(KSA).NE.1)    CALL LCZER1 (12,KSA,KSB,KSB,LLP,IERR)
      IF (NEVAL.EQ.0) THEN
         IF (KSA.EQ.(-1))   CALL LCZER1 (13,KSA,NEVAL,NEVAL,LLP,IERR)
      ELSE
         IF (IERR.EQ.0)     GO TO 100
      ENDIF
      IF (NEVAL.LT.0)       CALL LCZER1 (14,NEVAL,NEVEC,MAXNIV,LLP,IERR)
      IF (NEVAL.GT.MAXNIV)  CALL LCZER1 (14,NEVAL,NEVEC,MAXNIV,LLP,IERR)
      IF (MAXNIV.LT.1)      CALL LCZER1 (14,NEVAL,NEVEC,MAXNIV,LLP,IERR)
      IF (MAXNIV.GT.N)      CALL LCZER1 (14,NEVAL,NEVEC,MAXNIV,LLP,IERR)
      IF (NEVEC.GT.NEVAL)   CALL LCZER1 (14,NEVAL,NEVEC,MAXNIV,LLP,IERR)
      IF (ABS(KSB).LT.1)    CALL LCZER1 (12,KSA,KSB,KSB,LLP,IERR)
      IF (ABS(KSB).GT.2)    CALL LCZER1 (12,KSA,KSB,KSB,LLP,IERR)
C
      IF (IERR.LT.0)        GO TO 1000
C ----------------------------------------------------------------------
C   SET/INITIATE PARAMETERS, VARIABLES AND POINTERS
C ----------------------------------------------------------------------
  100 NEVEXT = 0
      NACEV  = 0
      NGDEV  = 0
      NIV    = 0
      NSDIG  = IMP(3)
      I      = NSDIG - NSDIG/5
      SMALL  = TEN**(-I)
      TINY   = SMALL*SMALL
      RAN    = 0.1234567891234567D0
C
      IF (ABS(SHIFT).LT.TINY)  SHIFT = ZERO
      IF (KSA.LT.0)            SHIFT = ZERO
C                                               ** DEFAULT VALUES ?
      I    = 2*NSDIG/3
      EPS1 = TOL(1)
      IF (EPS1.LE.ZERO) THEN
         EPS1   = TEN**(-I)
         TOL(1) = EPS1
      ENDIF
      EPS2 = TOL(2)
      IF (EPS2.LE.ZERO) THEN
         EPS2   = TEN**(-I)
         TOL(2) = EPS2
      ENDIF
      EPS3 = TOL(3)
      IF (EPS3.LE.ZERO) THEN
         EPS3   = SMALL
         TOL(3) = EPS3
      ENDIF
C
C                                               ** POINTERS IN RWA AND LWA
      IPA = 1
      LPP = 1
      IF (NEVAL.GT.0) THEN
         IPB = IPA+MAXNIV
         IPD = IPB+MAXNIV
         IPS = IPD+MAXNIV
         IPE = IPS+MAXNIV
         LPC = LPP+MAXNIV
C                                               ** INITIALIZE (TO ZERO)
         NN  = MAXNIV*(MAXNIV+4)
         CALL RMINT (RWA,NN,1,ZERO)
         CALL RMINT (EVL,MAXNIV,1,ZERO)
         NN  = 2*MAXNIV
         CALL IMINT (LWA,NN,1,0)
         CALL IMINT (MOP,8,1,0)
      ELSE
         IPB = 1
         IPD = 1
         IPS = 1
         IPE = 1
         LPC = 1
      ENDIF
C ----------------------------------------------------------------------
C     IF NECESSARY, MODIFY AND FACTORIZE MATRIX A, FACTORIZE MATRIX  B
C     AND PREPARE START VECTOR
C ----------------------------------------------------------------------
      IF (KSA.EQ.1) THEN
         CALL MDFCTA (A,B,MSKY,EPS2,SHIFT,N,KSB,LLP,NN,IERR)
         IF (IERR.LT.0) THEN
            CALL LCZER1 (21,I,I,I,LLP,IERR)
            MOP(1) = -NN
            GO TO 1000
         ENDIF
         KSA    =-1
         MOP(4) = NN
         IF (NEVAL.EQ.0)  GO TO 1000
      ENDIF
C                                               ** START VECTOR
      IF (KSR.EQ.1)  CALL RANVEC (R,RAN,N)
      IF (KSR.EQ.3)  CALL RMINT (R,N,1,ONE)
      IF (KSR.EQ.2)  THEN
         K = ABS(KSB)
         DO 120 I=1,N
            J = I
            IF (K.EQ.1) J = MSKY(I)
            R(I) = B(J)
  120    CONTINUE
      ENDIF
C                                               ** MATRIX B
      IF (KSB.GT.0) THEN
         CALL BLBLBT (B,MSKY,EPS2,N,KSB,LLP,NN,IERR)
         IF (IERR.LT.0) THEN
            CALL LCZER1 (22,I,I,I,LLP,IERR)
            GO TO 1000
         ENDIF
         KSB =-KSB
      ELSE
         NN  = 0
         K   = ABS(KSB)
         DO 140 I=1,N
            J = I
            IF (K.EQ.1) J = MSKY(I)
            IF (B(J).EQ.ZERO)  NN = NN+1
  140    CONTINUE
      ENDIF
      MOP(6) = NN
      IF (MAXNIV.GT.(N-NN)) THEN
         NN = N-NN
         CALL LCZER1 (15,MAXNIV,NN,NN,LLP,IERR)
         GO TO 1000
      ENDIF
C ======================================================================
C
C   GENERATE LANCZOS VECTORS AND REDUCE THE INVERSE EIGENPROBLEM TO
C   TRIDIAGONAL FORM, OF DIMENSION K BY K, AND SOLVE THIS WHEN
C   K = NEVAL, NEVAL+INT(3*NEVAL/2) AND THEN FOR EACH KEX'TH STEP,
C   UNTIL NEVAL EIGENVALUES ARE ACCEPTED (TWICE!)
C   K = 1,2,.....,M     WHERE  M  IS LESS THAN OR EQUAL TO MAXNIV
C
C ======================================================================
      K     = 0
      KODEX = NEVAL
      NREP  = 0
C                                               ** PREPARE START VECTOR
      IF (MIP(3).GT.0) THEN
         CALL SPSMUL (A,B,R(1),R(N+1),R(N+1),MSKY,
     +                EPS2,N,KSA,KSB,3,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL LCZER1 (23,K,K,K,LLP,IERR)
            GO TO 1000
         ENDIF
      ENDIF
C                                               ** NORMALIZE START
C                                                  VECTOR
      CALL NORVEC (R,N,1,IERR)
      IF (IERR.LT.0) THEN
         CALL RANVEC (R,RAN,N)
         CALL NORVEC (R,N,1,IERR)
      ENDIF
C
      BETA = ONE
C ================================================ START LANCZOS LOOP
C
  200 K     = K+1
      KOUNT = 0
      KODEX = KODEX-1
C                                               ** FORM LANCZOS VECTOR
      IP      = IPB-1+K
      RWA(IP) = BETA
      AUX     = ONE/BETA
      CALL SCALE (R,V(1,K),N,1,AUX)
      IF (K.GT.1) THEN
C                                               ** REORTHOGONALIZE
         NN = K-1
         CALL VRORT (V,V(1,K),EPS3,N,NN,MAXIT,IPSW,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL LCZER1 (24,K,K,K,LLP,IERR)
            GO TO 1000
         ENDIF
      ENDIF
C                                               ** PREPARE NEW R
C
      CALL SPSMUL (A,B,V(1,K),R(1),R(N+1),
     +             MSKY,EPS2,N,KSA,KSB,4,LLP,IERR)
      IF (IERR.LT.0) THEN
         CALL LCZER1 (23,K,K,K,LLP,IERR)
         GO TO 1000
      ENDIF
      IF (K.GT.1) THEN
         AUX =-BETA
         CALL VASCA (R,V(1,K-1),1,1,AUX,N)
      ENDIF
C                                               ** ELEMENT ALPHA-K OF
C                                                  TRIDIAG. FORM
  250 ALFA    = PRACC(V(1,K),R,1,1,N)
      IP      = IPA+K-1
      RWA(IP) = ALFA
C
      NIV = K
C
      IF (KODEX.EQ.0.OR.K.EQ.MAXNIV) THEN
C ------------------------------------------------ SOLVE REDUCED, TRI-
C                                                  DIAG. EIGENPROBLEM
C
         IF (NEVEXT.GT.0) CALL ICOPY (K,LWA(LPC),1,LWA(LPP),1)
         CALL RCOPY (RWA(IPA),RWA(IPD),K,1,K)
         CALL RCOPY (RWA(IPB),RWA(IPS),K,1,K)
C
         CALL QLS3D (RWA(IPD),RWA(IPS),K,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL LCZER1 (26,K,K,K,LLP,IERR)
            GO TO 1000
         ENDIF
         NEVEXT = NEVEXT+1
         CALL REARRV (RWA(IPD),K,4)
C                                               ** RELATIVE ERROR
         CALL ABRLER (EVL,RWA(IPD),EVERR,K)
         CALL RCOPY (RWA(IPD),EVL,K,1,K)
C                                               ** CHECK TOLERANCE
C
         CALL TOLCHK (EVERR,LWA(LPP),LWA(LPC),EPS1,K,NGDEV,NACEV)
C
         IF (IPSW.GE.2) THEN
C                                               ** PRINT
            NN = 1
            IF (IPSW.GE.3) NN=2
            CALL LCZPR3 (RWA(IPA),RWA(IPB),EVL,EVERR,
     +                   LWA(LPC),EPS1,NIV,NGDEV,LPU,NN)
         ENDIF
C                                               ** NEXT EIGENVALUE
C                                                  EXTRACTION
         KODEX = KEX
         IF (NEVEXT.EQ.1) THEN
            NN = NEVAL/2
            IF (NN.GT.0)  KODEX = NN
         ENDIF
C                                               ** LOOP EXIT
C
         IF (NGDEV.GE.NEVAL.OR.K.EQ.MAXNIV)  GO TO 300
      ENDIF
C
      AUX =-ALFA
      CALL VASCA (R,V(1,K),1,1,AUX,N)
C                                               ** LENGTH OF R
      AUX = PRACC(R,R,1,1,N)
      AUX = SQRT(AUX)
C                                               ** BREAKDOWN ?
      IF (K.GT.1) THEN
         IF (AUX.LT.BETA*EPS1 .OR. AUX.LT.ABS(ALFA)*EPS3) THEN
            NREP = NREP+1
            IF (NREP.GT.3 .OR. KOUNT.GT.0) THEN
               CALL LCZER1 (30,K,KOUNT,NREP,LLP,IERR)
               GO TO 1000
            ENDIF
            KOUNT = 1
            GO TO 250
         ENDIF
      ENDIF
C
      BETA = AUX
      GO TO 200
C ================================================ END LANCZOS LOOP
  300 MOP(1) = NGDEV
      MOP(2) = NACEV
      MOP(3) = NIV
      MOP(7) = NEVEXT
      IF (NGDEV.EQ.0) THEN
         CALL LCZER1 (32,MAXNIV,K,K,LLP,IERR)
         GO TO 1000
      ENDIF
      IF (NGDEV.LT.NEVAL) THEN
         CALL LCZER1 (1,NGDEV,MAXNIV,NEVAL,LLP,IERR)
      ENDIF
C
      IF (NEVEC.LT.1) THEN
C                                                  RESTORE EIGENVALUES
C                                                  LAMBDAi
         DO 350 I=1,NIV
            EVL(I) = SHIFT + ONE/EVL(I)
  350    CONTINUE
C
      ELSE
C ----------------------------------------------------------------------
C   DETERMINE (B-ORTHOGONAL) APPROXIMATIONS TO THE NV FIRST EIGENVECTORS
C ----------------------------------------------------------------------
         NV = NEVEC
         IF (NV.GT.NGDEV)  NV = NGDEV
         MOP(5) = NV
C                                               ** RE-SOLVE SMALL
C                                                  EIGENPROBLEM
C
         CALL QLS3DV (RWA(IPA),RWA(IPB),RWA(IPE),NIV,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL LCZER1 (26,NIV,NIV,NIV,LLP,IERR)
            GO TO 1000
         ENDIF
C                                               ** REARRANGE ORDERING
C
         CALL EVARR (RWA(IPA),RWA(IPE),NIV,NIV,NV,4)
C
C                                               ** RESTORE EIGENVECTORS
C                                                  OF SPECIAL PROBLEM
         CALL BIGMUL (V,RWA(IPE),V,R,N,NIV,NV,0)
C                                               ** RESTORE DESIRED
C                                                  EIGENVECTORS
         IF (MOP(6).EQ.0) THEN
            DO 400 J=1,NV
               CALL BCKSUB (B,V(1,J),MSKY,N,KSB,IERR)
               IF (IERR.LT.0) THEN
                  CALL LCZER1 (27,J,J,J,LLP,IERR)
                  GO TO 1000
               ENDIF
  400       CONTINUE
C
         ELSE
C
            DO 450 J=1,NV
               CALL SPSMUL (A,B,V(1,J),R,R,
     +                      MSKY,EPS2,N,KSA,KSB,1,LLP,IERR)
               IF (IERR.LT.0) THEN
                  CALL LCZER1 (27,J,J,J,LLP,IERR)
                  GO TO 1000
               ENDIF
               AUX = ONE/EVL(J)
               CALL SCALE (V(1,J),V(1,J),N,1,AUX)
  450       CONTINUE
         ENDIF
C                                               ** RESTORE EIGENVALUES
C                                                  LAMBDAi
         DO 475 I=1,NIV
            EVL(I) = SHIFT + ONE/EVL(I)
  475    CONTINUE
C
      ENDIF
C
      IF (IPSW.GE.1) THEN
C                                               ** PRINT
C
         CALL LCZPR2 (EVL,EVERR,TOL,LWA(LPC),MOP,NIV,LPU)
      ENDIF
C ------------------------------------------------ EXIT
C
 1000 RETURN
      END
      SUBROUTINE LCZPR1 (MSKY,MIP,SHIFT,N,NEVAL,NEVEC,MAXNIV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LCZPR1                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT KEY INPUT INFORMATION TO THE EIGENSOLUTION
C                IN SUBROUTINE LANCZ1
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-10-06 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           LPU,MAXNIV,N,NEVAL,NEVEC,    MSKY(N),MIP(5)
      DOUBLE PRECISION  SHIFT
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         WRITE (LPU,600)
         I = MSKY(N)/N
         WRITE (LPU,610) N,I
         I = ABS(MIP(2))
         IF (I.EQ.1) WRITE (LPU,620)
         IF (I.EQ.2) WRITE (LPU,630)
         WRITE (LPU,640) MIP(3),NEVAL,NEVEC,MAXNIV,MIP(4),SHIFT
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT('1'/5X,43('*')//
     +           5X,'PRINT FROM  S A M  LIBRARY ROUTINE LANCZ1'/
     +           5X,'(GENERAL EIGENPROBLEM BY TRUNCATED LANCZOS)'//
     +           5X,43('*')//)
  610 FORMAT(5X,'MATRIX DIMENSION  :',I5/
     +       5X,'AVERAGE BANDWIDTH :',I5)
  620 FORMAT(5X,'FORM OF B-MATRIX  :  SKYLINE')
  630 FORMAT(5X,'FORM OF B-MATRIX  :  DIAGONAL')
  640 FORMAT(/5X,'START VECTOR CODE        :',I5/
     +        5X,'REQUESTED EIGENVALUES    :',I5/
     +        5X,'REQUESTED EIGENVECTORS   :',I5/
     +        5X,'MAX. NUMBER OF STEPS     :',I5/
     +        5X,'REORTHOGONALIZATION CODE :',I5//
     +        5X,'SHIFT VALUE              :',1PE11.3)
C
      END
      SUBROUTINE LCZPR2 (EVL,EVERR,TOL,MANCUR,MOP,NIV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LCZPR2                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT KEY INFORMATION ON COMPLETION OF EIGENSOLUTION
C                IN LANCZ1, INCLUDING EIGENVALUES AND THEIR RELATIVE
C                ERRORS
C
C     ROUTINES CALLED/REFERENCED :  LCZPR3     (SAM-5)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-04-04 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           LPU,NIV,       MANCUR(NIV),MOP(8)
      DOUBLE PRECISION  EVL(NIV),EVERR(NIV),TOL(3)
C
      EXTERNAL          LCZPR3
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         WRITE (LPU,600) MOP(1),MOP(2),MOP(3),MOP(4),MOP(5),
     +                   MOP(6),MOP(7)
         WRITE (LPU,610) TOL(1),TOL(2),TOL(3)
         CALL LCZPR3 (EVL,EVL,EVL,EVERR,MANCUR,TOL(1),NIV,MOP(2),LPU,0)
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT(////5X,14('*')/5X,'LEAVING LANCZ1'/5X,14('*')//
     +           5X,'NO. OF "GOOD" EIGENVALUES             :',I5/
     +           5X,'NO. OF ACCEPTABLE EIGENVALUES         :',I5/
     +           5X,'NO. OF LANCZOS STEPS (VECTORS)        :',I5/
     +           5X,'NO. OF EIGENVALUES SMALLER THAN SHIFT :',I5/
     +           5X,'NO. OF EIGENVECTORS RETURNED          :',I5/
     +           5X,'NO. OF ZERO ELEMENTS IN DIAG.(B)      :',I5/
     +           5X,'NO. OF TRIDIAG. EGV.PROBLEMS SOLVED   :',I5)
  610 FORMAT(/   5X,'RELATIVE EIGENVALUE TOLERANCE    :',1PE11.3 /
     +           5X,'RELATIVE SINGULARITY TOLERANCE   :',1PE11.3 /
     +           5X,'RELATIVE ORTHOGONALITY TOLERANCE :',1PE11.3 )
C
      END
      SUBROUTINE LCZPR3 (ALPHA,BETA,EVL,EVERR,MANCUR,
     +                   EPS,NIV,NGDEV,LPU,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LCZPR3                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT COMPUTED EIGENVALUES AND THEIR CORRESPONDING
C                RELATIVE ERRORS AFTER A GIVEN LANCZOS STEP (NO. NIV).
C                THE STEP AT WHICH THE EIGENVALUES WERE ACCEPTED IS
C                ALSO PRINTED
C                OPTIONALLY (IFLAG.GT.1), THE TRIDIAGONAL FORM, STORED
C                IN ALPHA AND BETA, IS ALSO PRINTED
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-04-01 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,LPU,NGDEV,NIV,       MANCUR(NIV)
      DOUBLE PRECISION  EPS,  ALPHA(NIV),BETA(NIV),EVERR(NIV),EVL(NIV)
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         IF (IFLAG.GT.0) THEN
            WRITE (LPU,600) NIV
            IF (IFLAG.GT.1) THEN
               WRITE (LPU,610)
               DO 25 I=1,NIV
                  WRITE (LPU,620) I,BETA(I),ALPHA(I)
   25          CONTINUE
            ENDIF
            WRITE (LPU,630) NGDEV,EPS
         ENDIF
         WRITE (LPU,640)
         DO 50 I=1,NIV
            WRITE (LPU,650) I,EVL(I),EVERR(I),MANCUR(I)
   50    CONTINUE
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT(////5X,21('*')/5X,'LANCZOS STEP NO.',I5 / 5X,21('*'))
  610 FORMAT(//  5X,'TRIDIAGONAL FORM :' ///
     +          12X,'SUBDIAGONAL',6X,'DIAGONAL'/)
  620 FORMAT(I7,1P2D16.5)
  630 FORMAT(//  5X,'NUMBER OF "GOOD" EIGENVALUES :',I5,'   (TOL =',
     +           1PE11.3,' )')
  640 FORMAT(// 14X,'COMPUTED',13X,'RELATIVE',6X,'ACCEPTED'/
     +          14X,'EIGENVALUES',10X,'ERRORS',8X,'AT STEP'/)
  650 FORMAT(I7,1PD23.14,D15.4,I9)
C
      END
      SUBROUTINE LCZER1 (N,I1,I2,I3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LCZER1                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT ERROR/WARNING MESSAGES AND SET THE ERROR FLAG
C                FOR EIGENPROBLEM ROUTINE LANCZ1
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-04-04 / 1.0
C
C **********************************************************************
C
      INTEGER           IERR,I1,I2,I3,LPU,N
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         IF (N.LT.10) THEN
C                                               ** WARNINGS
            WRITE (LPU,600)
            IF (N.EQ.1)  WRITE (LPU,601) I1,I2,I3
         ELSE
C                                               ** ERRORS
            WRITE (LPU,690)
            IF (N.GT.10.AND.N.LT.20) THEN
               WRITE (LPU,610)
               IF (N.EQ.11)  WRITE (LPU,611) I1
               IF (N.EQ.12)  WRITE (LPU,612) I1,I2
               IF (N.EQ.13)  WRITE (LPU,613) I1,I2
               IF (N.EQ.14)  WRITE (LPU,614) I1,I2,I3
               IF (N.EQ.15)  WRITE (LPU,615) I1,I2
            ENDIF
            IF (N.EQ.21)  WRITE (LPU,621)
            IF (N.EQ.22)  WRITE (LPU,622)
            IF (N.EQ.23)  WRITE (LPU,623) I1
            IF (N.EQ.24)  WRITE (LPU,624) I1
            IF (N.EQ.26)  WRITE (LPU,626) I1
            IF (N.EQ.27)  WRITE (LPU,627) I1
            IF (N.EQ.30)  WRITE (LPU,630) I1,I2,I3
            IF (N.EQ.32)  WRITE (LPU,632) I1
         ENDIF
      ENDIF
C
      IF (N.EQ.1)   IERR = 1
      IF (N.GT.10.AND.N.LT.20)  IERR =-1
      IF (N.EQ.21)  IERR =-3
      IF (N.EQ.22)  IERR =-2
      IF (N.EQ.23)  IERR =-4
      IF (N.EQ.24)  IERR =-5
      IF (N.EQ.26)  IERR =-6
      IF (N.EQ.27)  IERR =-7
      IF (N.EQ.30)  IERR =-8
      IF (N.EQ.32)  IERR =-9
C
      RETURN
C ------------------------------------------------ FORMATS
  600 FORMAT(///' *** WARNING FROM  S A M  LIBRARY ROUTINE LANCZ1')
  601 FORMAT(5X,'ONLY',I5,'  EIGENVALUES ACCEPTED'/
     +       5X,'IN',I5,'  LANCZOS STEPS'/
     +       5X,'(THE REQUESTED NUMBER WAS',I5,'  )')
  610 FORMAT(5X,'ILLEGAL OR INCONSISTENT INPUT PARAMETERS')
  611 FORMAT(5X,'DIMENSION TOO SMALL ( N =',I4,' )')
  612 FORMAT(5X,'STORAGE CODES :  KSA =',I5,5X,'KSB =',I5)
  613 FORMAT(5X,'KSA =',I5,'  AND  NEVAL =',I5,'  ARE INCONSISTENT')
  614 FORMAT(5X,'NEVAL =',I5,5X,'NEVEC =',I5,5X,'MAXNIV =',I5)
  615 FORMAT(5X,'PARAMETER MAXNIV (=',I5,' ) EXCEEDS THE NUMBER OF'/
     +       5X,'FINITE EIGENVALUES (=',I5,' )')
C
  621 FORMAT(5X,'ERROR DURING FACTORIZATION OF MATRIX A')
  622 FORMAT(5X,'MATRIX B CANNOT BE FACTORIZED')
  623 FORMAT(5X,'ERROR DURING PREMULTIPLICATION (BY MATRIX H)'/
     +       5X,'LANCZOS STEP NO.',I5)
  624 FORMAT(5X,'ERROR DURING REORTHOGONALIZATION OF LANCZOS VECTOR NO.'
     +         ,I5)
  626 FORMAT(5X,'ERROR DURING SOLUTION OF SMALL, TRIDIAG. EIGENPROBLEM'
     +      /5X,'(OF DIMENSION',I4,' )')
  627 FORMAT(5X,'ERROR DURING RECOVERY OF EIGENVECTOR NO.',I5)
  630 FORMAT(5X,'BREAKDOWN OF ALGORITHM DURING LANCZOS STEP NO.',I5,/
     +       5X,'( KOUNT =',I3,5X,'NREP =',I3,' )')
  632 FORMAT(5X,'NO EIGENVALUES ACCEPTED (AFTER',I5,'  LANCZOS STEPS')
C
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE LANCZ1')
C
      END
      SUBROUTINE SSSIT (A,B,TOL,MSKY,MIP,
     +                  MOP,EVL,EVERR,EVEC,RWA,
     +                  SHIFT,N,NEVAL,NIV,LPU,IPSW,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SSSIT                   GROUP 5 / PUBLIC
C
C     T A S K :  TO DETERMINE THE NEVAL EIGENVALUES LAMBDAi AND CORRE-
C                SPONDING EIGENVECTORS Qi OF THE GENERAL, SYMMETRIC
C                EIGENPROBLEM
C                        (A - LAMBDAi*B)*Qi = 0
C                THAT ARE ALGEBRAICALLY CLOSEST TO A SPECIFIED SHIFT
C                OF ORIGON.
C     A IS A SYMMETRIC SKYLINE MATRIX, AND B IS A SYMMETRIC SKYLINE
C     MATRIX (SAME SKYLINE AS A) OR A DIAGONAL MATRIX.  BOTH THE
C     ORIGINAL MATRIX  A  AND ITS FACTORS  LA/DA  (A = LA*DA*LAT) ARE
C     ACCEPTED AS INPUT.
C     A STANDARD SUBSPACE ITERATION METHOD (WITH NIV ITERATION VECTORS)
C     ASSUMING ALL MATRICES IN PRIMARY STORAGE IS USED.
C     NOTE: IF THE PROBLEM IS SMALL (N.LE.6), IT IS SOLVED DIRECTLY
C           (BY GENJAC), THAT IS, WITHOUT ITERATION.
C
C     ROUTINES CALLED/REFERENCED :  IMP, IMINT, RMINT          (SAM-0)
C                                   RCOPY, SCALE, SCADD        (SAM-0)
C                                   PRESKY, SKYTRA, SMATB      (SAM-3)
C                                   PSTMUL                     (SAM-3)
C                                   SKYSOL                     (SAM-4)
C                                   MDFCTA, CHBDG, SSSSVC      (SAM-5)
C                                   GENJAC, EVARR, SSSER1      (SAM-5)
C                                   SSSPR1, SSSPR2, SSSPR3     (SAM-5)
C                                   SKYJAC                     (SAM-5)
C                                   RANVEC                     (SAM-8)
C                                   ABS, SQRT        (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-10-10 / 1.0
C                       88-01-30 / 1.1    K.Bell
C                       89-01-31 / 1.2    K.Bell
C                       92-03-04 / 1.3    K.Bell
C                       92-12-05 / 1.4    K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IPSW,LPU,N,NEVAL,NIV, MIP(5),MOP(4),MSKY(*)
      DOUBLE PRECISION  SHIFT,A(*),B(*),EVEC(N,NIV),EVERR(NIV),EVL(NIV),
     +                  RWA(*),TOL(3)
C
C                                               ** LOCAL VARIABLES
C
      INTEGER           I,INRAN,IP,IPAB,IPBB,IPSC,IPW,IPX,IT,J,KSA,KSB,
     +                  KSKY,KSR,LLP,MAXIT,NACEV,NFEV,NSDIG,NZDB,
     +                  IMP
      DOUBLE PRECISION  EPS1,EPS2,EPS3,ERR,ONE,RAN,S,SHIFTL,TEN,ZERO
C
      DOUBLE PRECISION  SIGMA
C
      PARAMETER         (ZERO = 0.0D0 , ONE = 1.0D0 , TEN = 10.0D0)
C
      EXTERNAL          CHBDG,EVARR,GENJAC,IMINT,IMP,MDFCTA,PRESKY,
     +                  PSTMUL,RANVEC,RCOPY,RMINT,SCADD,SCALE,SKYSOL,
     +                  SKYTRA,SMATB,SSSER1,SSSPR1,SSSPR2,SSSPR3,SSSSVC,
     +                  SKYJAC
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (IPSW.GT.0.AND.NEVAL.GT.0) THEN
         CALL SSSPR1 (MSKY,MIP,SHIFT,N,NEVAL,NIV,LPU)
      ENDIF
C                                               ** INITIATE PARAMETERS
      LLP = LPU
      IF (IPSW.LT.0)   LLP =-1
      KSA   = MIP(1)
      KSB   = MIP(2)
      KSR   = MIP(3)
      MAXIT = MIP(4)
      INRAN = 0
      IF (MAXIT.LT.1)  MAXIT = 15
      IF (KSR.LT.0) THEN
         INRAN = 1
         KSR   =-KSR
      ENDIF
      CALL IMINT (MOP,4,1,0)
C                                               ** CHECK INPUT
C
      IF (N.LT.1)           CALL SSSER1 (11,N,N,N,LLP,IERR)
      IF (ABS(KSA).NE.1)    CALL SSSER1 (12,KSA,KSB,KSB,LLP,IERR)
      IF (KSB.LT.1)         CALL SSSER1 (12,KSA,KSB,KSB,LLP,IERR)
      IF (KSB.GT.2)         CALL SSSER1 (12,KSA,KSB,KSB,LLP,IERR)
C
      CALL CHBDG (B,MSKY,N,KSB,NZDB,I)
      IF (I.LT.0) THEN
         IF (MIP(5) .EQ. 2) THEN
            CALL SSSER1 (22,I,I,I,LLP,IERR)
            GO TO 1000
         ENDIF
      ENDIF
      IF (KSB .EQ. 2) THEN
         NFEV = N-NZDB
      ELSE
         NFEV = N
      ENDIF
      MOP(3) = NZDB
C
      IF (NEVAL.EQ.0) THEN
         IF (KSA.EQ.(-1)) THEN
            CALL SSSER1 (13,KSA,NEVAL,NEVAL,LLP,IERR)
            GO TO 1000
         ENDIF
         GO TO 100
      ENDIF
C
      IF (NEVAL.LT.0)    CALL SSSER1 (14,NEVAL,NIV,NIV,LLP,IERR)
      IF (NEVAL.GT.NIV)  CALL SSSER1 (14,NEVAL,NIV,NIV,LLP,IERR)
      IF (NIV.GT.NFEV)   CALL SSSER1 (15,NIV,NFEV,NFEV,LLP,IERR)
C
C ----------------------------------------------------------------------
C  SMALL PROBLEM  (N.LE.6)  - SEQUENCE INTRODUCED IN VERSION 1.2
C ----------------------------------------------------------------------
      IF (N .LE. 6) THEN
         IF (KSA .EQ. 1) THEN
            NSDIG = IMP(3)
            EPS3  = TOL(3)
            IF (EPS3.LE.ZERO) THEN
               I      = NSDIG - NSDIG/5
               EPS3   = TEN**(-I)
               TOL(3) = EPS3
            ENDIF
            CALL SKYJAC (A,B,MSKY,RWA(1),RWA(37),
     +                   RWA(73),RWA(79),EPS3,N,KSB,
     +                   MIP(5),2,LPU,NFEV,IERR)
            IF (IERR.LT.0) THEN
               CALL SSSER1 (30,I,I,I,LLP,IERR)
               GO TO 1000
            ENDIF
            CALL RCOPY (RWA(73),EVL,NIV,1,1)
            CALL RCOPY (RWA(79),EVEC,N,NIV,1)
            CALL RMINT (EVERR,NIV,1,ZERO)
            IF (SHIFT .NE. ZERO) THEN
               DO 20 I=1,NFEV
                  IF (RWA(72+I) .LT. SHIFT)  MOP(2) = MOP(2)+1
   20          CONTINUE
               DO 30 I=1,NIV
                  EVL(I) = EVL(I) - SHIFT
   30          CONTINUE
               CALL EVARR (EVL,EVEC,N,NIV,NIV,2)
               DO 40 I=1,NIV
                  EVL(I) = EVL(I) + SHIFT
   40          CONTINUE
            ENDIF
            MOP(1) = NIV
            GO TO 800
         ENDIF
      ENDIF
C -------------------------------------------------------------------
C
      IF (NEVAL.GT.1) THEN
         IF (NEVAL.EQ.NIV)  CALL SSSER1 (14,NEVAL,NIV,NIV,LLP,IERR)
      ENDIF
C
      IF (IERR.LT.0)        GO TO 1000
C
C                                               ** SPECIAL CASE:  N=1
      IF (N.EQ.1) THEN
         EVL(1)    = A(1)/B(1)
         EVEC(1,1) = ONE/SQRT(B(1))
         MOP(1)    = 1
         IF (EVL(1).LT.SHIFT)  MOP(3) = 1
         GO TO 1000
      ENDIF
C                                               ** SET PARAMETERS
  100 NSDIG = IMP(3)
      RAN   = 0.1234567891234567D0
      KSKY  = KSB
      IF (KSB.EQ.2)  KSKY  = 0
      IF (KSA.LT.0)  SHIFT = ZERO
      IF (SHIFT.GT.ZERO) THEN
         SHIFTL = SHIFT+SHIFT
      ELSE
         SHIFTL = ZERO
      ENDIF
C                                               ** DEFAULT VALUES ?
      I    = 2*NSDIG/3
      EPS1 = TOL(1)
      IF (EPS1.LE.ZERO) THEN
         EPS1   = TEN**(-I)
         TOL(1) = EPS1
      ENDIF
      EPS2 = TOL(2)
      IF (EPS2.LE.ZERO) THEN
         EPS2   = TEN**(-I)
         TOL(2) = EPS2
      ENDIF
      EPS3 = TOL(3)
      IF (EPS3.LE.ZERO) THEN
         I      = NSDIG - NSDIG/5
         EPS3   = TEN**(-I)
         TOL(3) = EPS3
      ENDIF
C
C                                               ** POINTERS IN RWA
      IPAB = 1
      IF (NEVAL.GT.0) THEN
         IPSC = IPAB + NIV*NIV
         IPBB = IPSC + N
         IPX  = IPBB + NIV*NIV
         IPW  = IPSC
C                                               ** INITIALIZE (TO ZERO)
         CALL RMINT (EVL,NIV,1,ZERO)
      ELSE
         IPSC = 1
         IPBB = 1
         IPX  = 1
         IPW  = 1
      ENDIF
C ----------------------------------------------------------------------
C  IF NECESSARY, PREPARE START VECTORS AND MODIFY/FACTORIZE MATRIX  A
C ----------------------------------------------------------------------
      IF (KSR.EQ.1) THEN
         DO 210 J=1,NIV
            CALL RANVEC (EVEC(1,J),RAN,N)
  210    CONTINUE
      ELSEIF (KSR.EQ.2) THEN
         IF (KSA.EQ.1) THEN
            CALL SSSSVC (A,B,EVEC,RWA(IPSC),MSKY,RAN,N,NIV,KSB)
         ELSE
            DO 220 J=1,NIV
               CALL RANVEC (EVEC(1,J),RAN,N)
  220       CONTINUE
         ENDIF
      ENDIF
C
      IF (KSA.EQ.1) THEN
         CALL MDFCTA (A,B,MSKY,EPS2,SHIFT,N,KSB,LLP,I,IERR)
         IF (IERR.LT.0) THEN
            CALL SSSER1 (21,J,J,J,LLP,IERR)
            GO TO 1000
         ENDIF
         KSA    =-1
         MOP(2) = I
         IF (NEVAL.EQ.0)  GO TO 1000
      ENDIF
C
C ------------------------------------------------ START ITERATION LOOP
      DO 500 IT=1,MAXIT
C                                               ** PREMULTIPLY BY B
         CALL PRESKY (B,EVEC,EVEC,RWA(IPSC),
     +                MSKY,N,NIV,KSKY,1,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL SSSER1 (23,IT,IT,IT,LLP,IERR)
            GO TO 1000
         ENDIF
C                                               ** COPY INTO  W
         CALL RCOPY (EVEC,RWA(IPW),N,NIV,1)
C                                               ** SOLVE FOR RITZ VECTS.
         CALL SKYSOL (A,EVEC,MSKY,EPS2,N,N,
     +                NIV,LLP,4,I,IERR)
         IF (IERR.LT.0) THEN
            CALL SSSER1 (23,IT,IT,IT,LLP,IERR)
            GO TO 1000
         ENDIF
C                                               ** SUBSPACE MATRICES
         CALL SMATB (EVEC,RWA(IPW),RWA(IPAB),
     +               N,NIV,0)
         CALL SKYTRA (B,EVEC,RWA(IPBB),RWA(IPSC),
     +                MSKY,N,NIV,KSKY,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL SSSER1 (23,IT,IT,IT,LLP,IERR)
            GO TO 1000
         ENDIF
C                                               ** SOLVE REDUCED
C                                                  EIGENPROBLEM
         IF (MIP(5).EQ.2) THEN
C                                               ** ORIGINAL PROBLEM
C
            CALL GENJAC (RWA(IPAB),RWA(IPBB),EVERR,
     +                   RWA(IPX),EPS3,NIV,NIV,2,LLP,IERR)
            IF (IERR.LT.0) THEN
               CALL SSSER1 (24,IT,IT,IT,LLP,IERR)
               GO TO 1000
            ENDIF
         ELSE
C                                               ** INVERSE PROBLEM
            IF (SHIFTL.GT.ZERO) THEN
               CALL SCADD (RWA(IPAB),RWA(IPBB),
     +                     RWA(IPAB),NIV,NIV,SHIFTL)
            ENDIF
            CALL GENJAC (RWA(IPBB),RWA(IPAB),EVERR,
     +                   RWA(IPX),EPS3,NIV,NIV,0,LLP,IERR)
            IF (IERR.LT.0) THEN
               CALL SSSER1 (24,IT,IT,IT,LLP,IERR)
               GO TO 1000
            ENDIF
C
            SIGMA = EPS1
            DO 250 I=1,NIV
               IF (EVERR(I) .LT. SIGMA)  SIGMA = EVERR(I)
  250       CONTINUE
            IF (SIGMA .LT. EPS1) THEN
               SIGMA = ABS(SIGMA) + EPS1
               DO 260 I=1,NIV
                  EVERR(I) = EVERR(I) + SIGMA
  260          CONTINUE
            ELSE
               SIGMA = ZERO
            ENDIF
C
C                                               ** SCALE EIGENVECTORS
            IP = IPX
            DO 300 I=1,NIV
               S  = ONE/SQRT(EVERR(I))
               CALL SCALE (RWA(IP),RWA(IP),NIV,1,S)
               IP = IP+NIV
  300       CONTINUE
C
            DO 320 I=1,NIV
               EVERR(I) = ONE/(EVERR(I) - SIGMA) - SHIFTL
  320       CONTINUE
            CALL EVARR (EVERR,RWA(IPX),NIV,NIV,NIV,2)
         ENDIF
C                                               ** EIGENVAL. ESTIMATES
C                                                  AND THEIR REL. ERRORS
         DO 350 I=1,NIV
            IF (EVERR(I).EQ.ZERO) THEN
               IF (EVL(I).EQ.ZERO) THEN
                  ERR = ZERO
               ELSE
                  ERR = ONE
               ENDIF
            ELSE
               ERR = ABS(EVERR(I)-EVL(I))/ABS(EVERR(I))
            ENDIF
            EVL(I)   = EVERR(I)
            EVERR(I) = ERR
  350    CONTINUE
C                                               ** IMPROVED EIGENVECTORS
         CALL PSTMUL (RWA(IPX),EVEC,
     +                RWA(IPSC),N,NIV,1)
C                                               ** CHECK CONVERGENCE
         NACEV = 0
         DO 400 I=1,NIV
            IF (EVERR(I).GT.EPS1)  GO TO 410
            NACEV = NACEV+1
  400    CONTINUE
  410    MOP(1) = NACEV
         MOP(4) = IT
C                                               ** EXIT IF ENOUGH EIGEN-
C                                                  VALUES ARE ACCEPTED
         IF (NACEV.GE.NEVAL)  GO TO 600
C
         IF (IPSW.GT.1) THEN
C                                               ** PRINT
            CALL SSSPR3 (EVL,EVERR,EPS1,
     +                   NIV,IT,NACEV,LLP,1)
         ENDIF
C
         IF (INRAN.EQ.1.AND.NIV.GT.2) THEN
C                                               ** INCL. A RANDOM VECTOR
            CALL RANVEC (EVEC(1,NIV),RAN,N)
         ENDIF
  500 CONTINUE
C ------------------------------------------------ END OF ITERATION LOOP
      IF (NACEV.EQ.0) THEN
         CALL SSSER1 (25,MAXIT,MAXIT,MAXIT,LLP,IERR)
      ELSE
         CALL SSSER1 (1,NACEV,MAXIT,NEVAL,LLP,IERR)
      ENDIF
C
  600 IF (SHIFT.NE.ZERO) THEN
C                                               ** RESTORE EIGENVALUES
C                                                  LAMBDAi
         DO 700 I=1,NIV
            EVL(I) = EVL(I) + SHIFT
  700    CONTINUE
      ENDIF
C
  800 IF (IPSW.GT.0) THEN
C                                               ** PRINT
         CALL SSSPR2 (EVL,EVERR,TOL,MOP,NIV,LPU)
      ENDIF
C
 1000 RETURN
      END
      SUBROUTINE SSSPR1 (MSKY,MIP,SHIFT,N,NEVAL,NIV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SSSPR1                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT KEY INPUT INFORMATION TO THE EIGENSOLUTION
C                IN SUBROUTINE SSSIT
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-08-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           LPU,N,NEVAL,NIV,    MSKY(N),MIP(5)
      DOUBLE PRECISION  SHIFT
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         WRITE (LPU,600)
         I = MSKY(N)/N
         WRITE (LPU,610) N,I
         I = MIP(2)
         IF (I.EQ.1) WRITE (LPU,620)
         IF (I.EQ.2) WRITE (LPU,630)
         WRITE (LPU,640) MIP(3),NEVAL,NIV,MIP(4),MIP(5),SHIFT
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT('1'/5X,43('*')//
     +           5X,'PRINT FROM  S A M  LIBRARY ROUTINE SSSIT'/
     +          13X,'GENERAL EIGENPROBLEM BY'/
     +          11X,'STANDARD SUBSPACE ITERATION'//
     +           5X,43('*')//)
  610 FORMAT(5X,'MATRIX DIMENSION  :',I5/
     +       5X,'AVERAGE BANDWIDTH :',I5)
  620 FORMAT(5X,'FORM OF B-MATRIX  :  SKYLINE')
  630 FORMAT(5X,'FORM OF B-MATRIX  :  DIAGONAL')
  640 FORMAT(/5X,'START VECTOR CODE         :',I5/
     +        5X,'REQUESTED EIGENVALUES     :',I5/
     +        5X,'NO. OF ITERATION VECTORS  :',I5/
     +        5X,'MAX. NUMBER OF ITERATIONS :',I5/
     +        5X,'POS.DEF. INDICATOR        :',I5//
     +        5X,'SHIFT VALUE               :',1PE11.3)
C
      END
      SUBROUTINE SSSPR2 (EVL,EVERR,TOL,MOP,NIV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SSSPR2                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT KEY INFORMATION ON COMPLETION OF EIGENSOLUTION
C                IN SSSIT, INCLUDING EIGENVALUES AND THEIR RELATIVE
C                ERRORS
C
C     ROUTINES CALLED/REFERENCED :  SSSPR3     (SAM-5)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-10-10 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           LPU,NIV,            MOP(4)
      DOUBLE PRECISION  EVL(NIV),EVERR(NIV),TOL(3)
C
      EXTERNAL          SSSPR3
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         WRITE (LPU,600) MOP(1),MOP(2),MOP(3),MOP(4)
         WRITE (LPU,610) TOL(1),TOL(2),TOL(3)
         CALL SSSPR3 (EVL,EVERR,TOL(1),NIV,MOP(4),MOP(1),LPU,0)
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT(////5X,13('*')/5X,'LEAVING SSSIT'/5X,13('*')//
     +           5X,'NO. OF EIGENVALUES ACCEPTED           :',I5/
     +           5X,'NO. OF EIGENVALUES SMALLER THAN SHIFT :',I5/
     +           5X,'NO. OF ZERO ELEMENTS IN DIAG.(B)      :',I5/
     +           5X,'NO. OF ITERATIONS                     :',I5)
  610 FORMAT(/   5X,'RELATIVE EIGENVALUE TOLERANCE    :',1PE11.3 /
     +           5X,'RELATIVE SINGULARITY TOLERANCE   :',1PE11.3 /
     +           5X,'RELATIVE EIGENVALUE TOLERANCE'/
     +           5X,'FOR THE REDUCED EIGENPROBLEM     :',1PE11.3 )
C
      END
      SUBROUTINE SSSPR3 (EVL,EVERR,EPS,NIV,IT,NACEV,LPU,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SSSPR3                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT COMPUTED EIGENVALUES AND THEIR CORRESPONDING
C                RELATIVE ERRORS AFTER
C                 - ITERATION NO. IT     (IF IFLAG.GT.0)
C                 - THE LAST ITERATION   (IF IFLAG.EQ.0)
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-10-10 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,IT,LPU,NACEV,NIV
      DOUBLE PRECISION  EPS,  EVERR(NIV),EVL(NIV)
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         IF (IFLAG.GT.0) THEN
            WRITE (LPU,600) IT
            WRITE (LPU,610) NACEV,EPS
         ENDIF
         WRITE (LPU,620)
         DO 10 I=1,NIV
            WRITE (LPU,630) I,EVL(I),EVERR(I)
   10    CONTINUE
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT(////5X,18('*')/5X,'ITERATION NO.',I5 / 5X,18('*'))
  610 FORMAT(//  5X,'NUMBER OF ACCEPTED EIGENVALUES :',I5,
     +       '   (TOL =',1PE11.3,' )')
  620 FORMAT(// 14X,'COMPUTED',13X,'RELATIVE'/
     +          14X,'EIGENVALUES',10X,'ERRORS'/)
  630 FORMAT(I7,1PD23.14,D15.4)
C
      END
      SUBROUTINE SSSER1 (N,I1,I2,I3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SSSER1                 GROUP 5 / PRIVATE
C
C     T A S K :  TO PRINT ERROR/WARNING MESSAGES AND SET THE ERROR FLAG
C                FOR EIGENPROBLEM ROUTINE SSSIT
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-10-10 / 1.0
C
C **********************************************************************
C
      INTEGER           IERR,I1,I2,I3,LPU,N
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         IF (N.LT.10) THEN
C                                               ** WARNINGS
            WRITE (LPU,600)
            IF (N.EQ.1)  WRITE (LPU,601) I1,I2,I3
         ELSE
C                                               ** ERRORS
            WRITE (LPU,690)
            IF (N.GT.10.AND.N.LT.20) THEN
               WRITE (LPU,610)
               IF (N.EQ.11)  WRITE (LPU,611) I1
               IF (N.EQ.12)  WRITE (LPU,612) I1
               IF (N.EQ.13)  WRITE (LPU,613) I1,I2
               IF (N.EQ.14)  WRITE (LPU,614) I1,I2
               IF (N.EQ.15)  WRITE (LPU,615) I1,I2
            ENDIF
            IF (N.EQ.21)  WRITE (LPU,621)
            IF (N.EQ.22)  WRITE (LPU,622)
            IF (N.EQ.23)  WRITE (LPU,623) I1
            IF (N.EQ.24)  WRITE (LPU,624) I1
            IF (N.EQ.25)  WRITE (LPU,625) I1
            IF (N.EQ.30)  WRITE (LPU,630)
         ENDIF
      ENDIF
C
      IF (N.EQ.1)   IERR = 1
      IF (N.GT.10.AND.N.LT.20)  IERR =-1
      IF (N.EQ.21)  IERR =-2
      IF (N.EQ.22)  IERR =-3
      IF (N.EQ.23)  IERR =-4
      IF (N.EQ.24)  IERR =-5
      IF (N.EQ.25)  IERR =-6
      IF (N.EQ.30)  IERR =-5
C
      RETURN
C ----------------------------------------------------------------------
  600 FORMAT(///' *** WARNING FROM  S A M  LIBRARY ROUTINE SSSIT')
  601 FORMAT(5X,'ONLY',I5,'  EIGENVALUES ACCEPTED'/
     +       5X,'IN',I5,'  ITERATIONS'/
     +       5X,'(THE REQUESTED NUMBER WAS',I5,'  )')
  610 FORMAT(5X,'ILLEGAL OR INCONSISTENT INPUT PARAMETERS')
  611 FORMAT(5X,'DIMENSION TOO SMALL ( N =',I4,' )')
  612 FORMAT(5X,'STORAGE CODE :  KSA =',I5)
  613 FORMAT(5X,'KSA =',I5,'  AND  NEVAL =',I5,'  ARE INCONSISTENT')
  614 FORMAT(5X,'NEVAL =',I5,5X,'NIV =',I5)
  615 FORMAT(5X,'PARAMETER  NIV (=',I5,' ) EXCEEDS THE NUMBER OF'/
     +       5X,'FINITE EIGENVALUES (=',I5,' )')
C
  621 FORMAT(5X,'ERROR DURING FACTORIZATION OF MATRIX A')
  622 FORMAT(5X,'MATRIX B IS NEGATIVE DEFINITE')
  623 FORMAT(5X,'ERROR DURING MATRIX MULTIPLICATION'/
     +       5X,'OR SUBSTITUTION IN ITERATION NO.',I5)
  624 FORMAT(5X,'ERROR DURING SOLUTION OF REDUCED EIGENPROBLEM'
     +      /5X,'IN ITERATION NO.',I5)
  625 FORMAT(5X,'NO EIGENVALUES ACCEPTED  (AFTER',I4,'ITERATIONS')
  630 FORMAT(5X,'ERROR DURING DIRECT (JACOBI) SOLUTION')
C
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE SSSIT')
C
      END
      SUBROUTINE LANCZ2 (A,B,TOL,MSKY,MIP,
     &                   MOP,EVL,EVERR,V,R,RWA,LWA,
     &                   SHIFT,N,NEVAL,NEVEC,MAXLAN,LPU,IPSW,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LANCZ2                  GROUP 5 / PUBLIC
C
C     TASK :  To determine the NEVAL eigenvalues Lambda_i of the
C             generalized , symmetric eigenproblem
C                  (A - Lambda_i*B)*q_i = 0
C             that are closest to a specified SHIFT of origin.
C             The routine also returns NLV (=MOP(3)) ortogonal Lanczos
C             vectors (in V) or, optionally, approximations to NEVEC
C             eigenvectors (corresponding to the NEVEC first eigen-
C             values).  The eigenvectors, which are B-orthonormal, are
C             returned as the first NEVEC vectors in V.
C     A is a symmetric skyline matrix, and B is a symmetric skyline
C     matrix (same skyline as A) or a diagonal matrix.  For both A and B
C     the original matrices or their appropriate factors are accepted as
C     input (in arrays A and B).
C     A truncated Lanczos method, assuming all matrices in primary
C     storage, is used to transform the problem to tridiagonal, special
C     form.  The reduced, tridiagonal eigenproblem is solved by implicit
C     QL transformation.
C     Various reorthogonalization schemes may be specified.
C
C
C     ROUTINES CALLED/REFERENCED :  IMP, RMINT, IMINT, PRACC    (SAM-0)
C                                   SCALE, VASCA, ICOPY, RCOPY  (SAM-0)
C                                   BIGMUL, NORVEC              (SAM-3)
C                                   FACSKY, FSBSKY, BCKSKY      (SAM-3)
C                                   SKYSOL                      (SAM-4)
C                                   QLS3D,  QLS3DV, EVARR       (SAM-5)
C                                   MDSHFT, SPSMM,  REORL2      (SAM-5)
C                                   MGSOR1, ABRLER, TOLCHK      (SAM-5)
C                                   LANPR1, LANPR2              (SAM-5)
C                                   LANPR3, LANERR              (SAM-5)
C                                   RANVEC, REARRV              (SAM-8)
C                                   ABS, SQRT         (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-06-15 / 1.0
C                       96-12-29 / 1.1  K.Bell
C                       01-05-24 / 1.2  K.M.Okstad
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IPSW,LPU,MAXLAN,N,NEVAL,NEVEC,
     &                  LWA(*),MIP(6),MOP(10),MSKY(*)
      DOUBLE PRECISION  SHIFT, A(*),B(*),EVERR(MAXLAN),EVL(MAXLAN),
     &                  R(*),RWA(*),TOL(2),V(N,MAXLAN)
C
      INTEGER           I,IFLAG,IP,IPA,IPB,IPD,IPE,IPK,IPKM1,IPKP1,IPS,
     &                  J,JORT,JP,K,KEIG,KEVEX,KEX,KM1,KORT,KPD,KRST,
     &                  KSA,KSB,KSR,LLP,LPC,LPP,NACEV,NEVEXT,
     &                  NGDEV,NLV,NN,NSDIG,NVC
      INTEGER           IMP
      DOUBLE PRECISION  ALFA,AUX,BETA,EPS1,EPS2,ONE,RAN,RLPR,
     &                  S,SIGMA,SMALL,TEN,TINY,TRSH,VPMAX,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 , TEN = 10.0D0 )
C
      EXTERNAL          IMP,PRACC
C ----------------------------------------------------------------------
C
      IERR = 0
C
      IF (IPSW.GT.0 .AND. NEVAL.GT.0) THEN
         CALL LANPR1 (MSKY,MIP,SHIFT,N,NEVAL,NEVEC,MAXLAN,LPU)
      ENDIF
C
C --- Check input parameters -------------------------------------------
C
      LLP   = LPU
      IF (IPSW .LT. 0)  LLP =-1
      KSA  = MIP(1)
      KSB  = MIP(2)
      KSR  = ABS(MIP(3))
      KORT = MIP(4)
      KEX  = MIP(5)
      KPD  = MIP(6)
C
      IF (KSA.LT.1 .OR. KSA.GT.3)  CALL LANERR (11,KSA,KSB,I,LLP,IERR)
      I = ABS(KSB)
      IF (I.LT.1   .OR.   I.GT.2)  CALL LANERR (11,KSA,KSB,I,LLP,IERR)
      IF (ABS(KORT) .GT. 2)        CALL LANERR (12,KORT,I,I,LLP,IERR)
      IF (KEX.LT.1 .OR. KEX.GT.5)  CALL LANERR (13,KEX,I,I,LLP,IERR)
      IF (KPD.LT.1 .OR. KPD.GT.2)  CALL LANERR (14,KPD,I,I,LLP,IERR)
      IF (KPD .EQ. 1) THEN
         IF (KSA.EQ.2 .OR. KSB.LT.0) THEN
            CALL LANERR (15,KPD,KSA,KSB,LLP,IERR)
         ENDIF
      ENDIF
      IF (KSB.LT.0 .AND. KSR.EQ.2) THEN
         CALL LANERR (16,KSB,MIP(3),I,LLP,IERR)
      ENDIF
      IF (NEVAL.EQ.0 .AND. KSA.GT.1) THEN
         CALL LANERR (17,KSA,NEVAL,NEVAL,LLP,IERR)
      ENDIF
      IF (NEVAL .GT. MAXLAN) THEN
         CALL LANERR (18,NEVAL,NEVEC,MAXLAN,LLP,IERR)
      ENDIF
      IF (MAXLAN.LT.1 .AND. NEVAL.GT.0) THEN
         CALL LANERR (18,NEVAL,NEVEC,MAXLAN,LLP,IERR)
      ENDIF
      IF (MAXLAN.GT.N)      CALL LANERR (18,NEVAL,NEVEC,MAXLAN,LLP,IERR)
      IF (NEVEC.GT.NEVAL)   CALL LANERR (18,NEVAL,NEVEC,MAXLAN,LLP,IERR)
C
      IF (IERR.LT.0)        GO TO 1000
      IF (N .EQ. 1) THEN
C
C ------ Special case:  N = 1 ------------------------------------------
C
         IF (B(1) .EQ. ZERO) THEN
            CALL LANERR (23,I,I,I,LLP,IERR)
            GO TO 1000
         ENDIF
         EVL(1)    = A(1) / B(1)
         IF (NEVEC .GT. 0) THEN
            IF (B(1) .GT. ZERO) THEN
               V(1,1) = ONE / SQRT(B(1))
            ELSE
               V(1,1) = ONE
            ENDIF
         ENDIF
         MOP(1)    = 1
         IF (EVL(1) .LT. SHIFT)  MOP(2) = 1
         GO TO 1000
      ENDIF
C
C --- Set/initiate parameters, variables and pointers ------------------
C
      NEVEXT = 0
      NACEV  = 0
      NGDEV  = 0
      NLV    = 0
      NSDIG  = IMP(3)
      RLPR   = TEN**(-NSDIG)
      TRSH   = SQRT(RLPR)
      I      = NSDIG - NSDIG/5
      SMALL  = TEN**(-I)
      TINY   = SMALL*SMALL
      RAN    = 0.1234567891234567D0
C
      IF (ABS(SHIFT) .LT. TINY)  SHIFT = ZERO
      IF (SHIFT .NE. ZERO) THEN
         IF (KPD .EQ. 1) THEN
            CALL LANERR (2,I,I,I,LLP,IERR)
         ENDIF
         IF (KSA.GT.1 .OR. KSB.LT.0) THEN
            CALL LANERR (19,KSA,KSB,I,LLP,IERR)
            GO TO 1000
         ENDIF
      ENDIF
C                                                !  default values ?
      I    = 2*NSDIG/3
      AUX  = TEN**(-I)
      EPS1 = TOL(1)
      IF (EPS1 .LE. ZERO) THEN
         EPS1   = AUX
         TOL(1) = EPS1
      ENDIF
      EPS2 = TOL(2)
      IF (EPS2 .LE. ZERO) THEN
         EPS2   = AUX
         TOL(2) = EPS2
      ENDIF
C
C                                                ! pointers in RWA and LWA
      IPA = 1
      LPP = 1
      IF (NEVAL .GT. 0) THEN
         IPB   = IPA + MAXLAN
         IPD   = IPB + MAXLAN
         IPS   = IPD + MAXLAN
         IPE   = IPS + MAXLAN
         IPK   = IPE + MAXLAN*MAXLAN
         IPKM1 = IPK + MAXLAN
         IPKP1 = IPKM1 + MAXLAN
         LPC   = LPP + MAXLAN
C                                                ! initialize (to zero)
         NN  = MAXLAN*(MAXLAN+7)
         CALL RMINT (RWA,NN,1,ZERO)
         CALL RMINT (EVL,MAXLAN,1,ZERO)
         NN  = 2*MAXLAN
         CALL IMINT (LWA,NN,1,0)
         CALL IMINT (MOP,10,1,0)
      ELSE
         IPB   = 1
         IPD   = 1
         IPS   = 1
         IPE   = 1
         IPK   = 1
         IPKM1 = 1
         IPKP1 = 1
         LPC   = 1
      ENDIF
C
C --- If necessary, modify and factorize matrix A,
C     if KPD=2 factorize matrix B,
C     and prepare start vector
C
      IF (KSA .EQ. 1) THEN
C
         CALL MDSHFT (A,B,MSKY,SHIFT,N,KSB)
C
         IF (KPD.EQ.1 .AND. NEVAL.GT.0) THEN
            CALL FACSKY (A,MSKY,EPS2,N,LLP,I,IERR)
            IF (IERR .LT. 0) THEN
               CALL LANERR (21,I,I,I,LLP,IERR)
               MOP(1) = -I
               GO TO 1000
            ENDIF
            KSA = 3
         ELSE
            CALL SKYSOL (A,B,MSKY,EPS2,N,N,N,LLP,1,NN,IERR)
            IF (IERR .LT. 0) THEN
               CALL LANERR (21,I,I,I,LLP,IERR)
               MOP(1) = -NN
               GO TO 1000
            ENDIF
            KSA = 2
            MOP(4) = NN
C                                                ! exit if NEVAL=0
            IF (NEVAL .EQ. 0)  GO TO 1000
         ENDIF
      ENDIF
C                                                ! start vector
      IF (KSR .EQ. 1)  CALL RANVEC (R,RAN,N)
      IF (KSR .EQ. 3)  CALL RMINT (R,N,1,ONE)
      IF (KSR .EQ. 2)  THEN
         DO 20 I=1,N
            J = I
            IF (KSB .EQ. 1)  J = MSKY(I)
            R(I) = B(J)
   20    CONTINUE
      ENDIF
C                                                ! matrix B
      NN = 0
      IF (KPD.EQ.2 .AND. KSB.GT.0) THEN
         IF (KSB .EQ. 1) THEN
            CALL FACSKY (B,MSKY,EPS2,N,LLP,I,IERR)
            IF (IERR .LT. 0) THEN
               CALL LANERR (22,KSB,I,I,LLP,IERR)
               MOP(1) = -I
               GO TO 1000
            ENDIF
            KSB = -1
         ELSE
C                                                ! B is diagonal
            DO 40 I=1,N
               IF (B(I) .GT. ZERO) THEN
                  B(I) = SQRT(B(I))
               ELSEIF (B(I) .EQ. ZERO) THEN
                  NN = NN+1
               ELSE
                  CALL LANERR (22,KSB,I,I,LLP,IERR)
                  MOP(1) = -I
                  GO TO 1000
               ENDIF
   40       CONTINUE
            KSB = -2
         ENDIF
      ELSEIF (ABS(KSB) .EQ. 2) THEN
         DO 60 I=1,N
            IF (B(I) .EQ. ZERO)  NN = NN+1
   60    CONTINUE
      ENDIF
      MOP(6) = NN
C ----------------------------------------------------------------------
C
C   Generate Lanczos vectors and reduce the inverse eigenproblem to
C   tridiagonal form, of dimension K by K, and if KEIG=0 solve this when
C      K = NEVAL, NEVAL+INT(3*NEVAL/2) and then for each KEX'th step
C   until NEVAL eigenvalues are accepted (twice!)
C   If  KEIG=1  tridiagonalization is completed and the tridiagonal
C   problem is solved only once (i.e. complete solution)
C
C   K = 1,2,......,M    where  M  is less than or equal to MAXLAN
C
C ----------------------------------------------------------------------
      K    = 0
C
      KEIG = 0
      IF (MAXLAN.EQ.N .AND. N.LT.25)  KEIG=1
C
      KEVEX = NEVAL
C
      IF (MIP(3) .GT. 0) THEN
C                                               ! premultiply start
C                                                 vector by matrix H
         IF (KPD .EQ. 1) THEN
            IFLAG = 5
         ELSE
            IFLAG = 3
         ENDIF
         CALL SPSMM (A,B,R(1),R(N+1),R(N+1),
     &               MSKY,N,KSA,KSB,IFLAG,LLP,IERR)
         IF (IERR .LT. 0) THEN
            CALL LANERR (24,K,K,K,LLP,IERR)
            GO TO 1000
         ENDIF
      ENDIF
C                                                ! length of R
      AUX = PRACC(R,R,1,1,N)
      IF (AUX .LT. SMALL) THEN
         CALL RANVEC (R,RAN,N)
         AUX = PRACC(R,R,1,1,N)
      ENDIF
      BETA = SQRT(AUX)
      IF (KPD .EQ. 1) THEN
         IFLAG = 6
      ELSE
         IFLAG = 4
      ENDIF
      RWA(IPB) = ONE
C
      KRST = 0
      JORT = 0
      IF (KORT .EQ. -2)  JORT = -1
C
C ================================================ Start Lanczos loop
C
  100 K     = K+1
      IF (KEIG .EQ. 0)  KEVEX = KEVEX-1
C                                                ! form Lanczos vector
      AUX = ONE/BETA
      CALL SCALE (R,V(1,K),N,1,AUX)
      IF (KRST.EQ.0 .AND. K.GT.1) THEN
C
C --- Reorthogonalization - if necessary and/or if specified -----------
C
         KM1 = K-1
         IF (KORT .GT. 0) THEN
            MOP(8) = MOP(8)+1
            DO 150 I=1,KORT
               CALL MGSOR1 (V,V(1,K),N,KM1)
               IF (IPSW .GT. 1)  WRITE (LPU,610) K
               MOP(9) = MOP(9)+KM1
  150       CONTINUE
            IF (IPSW .GT. 1)  WRITE (LPU,600)
         ELSEIF (KORT.EQ.(-1) .AND. K.GT.2) THEN
            IF (JORT .EQ. 1) THEN
               CALL MGSOR1 (V,V(1,K),N,KM1)
               IF (IPSW .GT. 1) THEN
                  WRITE (LPU,610) K
                  WRITE (LPU,600)
               ENDIF
               MOP(8) = MOP(8)+1
               MOP(9) = MOP(9)+KM1
               CALL RMINT (RWA(IPK),K,1,RLPR)
               CALL RMINT (RWA(IPKM1),K,1,RLPR)
               JORT = 0
            ELSE
               CALL REORL2 (RWA(IPA),RWA(IPB),RWA(IPK),
     &                      RWA(IPKM1),RWA(IPKP1),BETA,
     &                      RLPR,TRSH,KM1,VPMAX)
               IF (IPSW .GT. 1)  WRITE (LPU,620) K,VPMAX
               IF (VPMAX .GT. TRSH) THEN
                  CALL MGSOR1 (V,V(1,K),N,KM1)
                  IF (IPSW .GT. 1)  WRITE (LPU,610) K
                  MOP(8) = MOP(8)+1
                  MOP(9) = MOP(9)+KM1
                  JORT   = 1
               ENDIF
               IF (IPSW .GT. 1)  WRITE (LPU,600)
            ENDIF
         ELSEIF (KORT .EQ. -2) THEN
            JORT = JORT+1
            IF (JORT .GT. 0) THEN
               CALL MGSOR1(V,V(1,K),N,KM1)
               IF (IPSW .GT. 1) THEN
                  WRITE (LPU,610) K
                  WRITE (LPU,600)
               ENDIF
               MOP(8) = MOP(8)+1
               MOP(9) = MOP(9)+KM1
            ENDIF
            IF (JORT .EQ. 2)  JORT = -2
         ENDIF
      ENDIF
C                                                ! prepare new R
      CALL SPSMM (A,B,V(1,K),R(1),R(N+1),
     &            MSKY,N,KSA,KSB,IFLAG,LLP,IERR)
      IF (IERR.LT.0) THEN
         CALL LANERR (24,K,K,K,LLP,IERR)
         GO TO 1000
      ENDIF
      IF (KRST.EQ.0 .AND. K.GT.1) THEN
         AUX =-BETA
         CALL VASCA (R,V(1,K-1),1,1,AUX,N)
      ENDIF
      KRST = 0
C                                                ! element alpha_k of
C                                                  tridiag. form
      ALFA    = PRACC(V(1,K),R,1,1,N)
      IP      = IPA+K-1
      RWA(IP) = ALFA
C
      NLV = K
C
      IF (KEVEX.EQ.0 .OR. K.EQ.MAXLAN) THEN
C
C --- Solve reduced, tridiagonal eigenproblem --------------------------
C
         IF (NEVEXT .GT. 0) CALL ICOPY (K,LWA(LPC),1,LWA(LPP),1)
         CALL RCOPY (RWA(IPA),RWA(IPD),K,1,K)
         CALL RCOPY (RWA(IPB),RWA(IPS),K,1,K)
C
         CALL QLS3D (RWA(IPD),RWA(IPS),K,LLP,IERR)
         IF (IERR.LT.0) THEN
            CALL LANERR (25,K,K,K,LLP,IERR)
            GO TO 1000
         ENDIF
         NEVEXT = NEVEXT+1
         CALL REARRV (RWA(IPD),K,4)
C                                                ! relative error
         CALL ABRLER (EVL,RWA(IPD),EVERR,K)
         CALL RCOPY (RWA(IPD),EVL,K,1,K)
C                                                ! check tolerance
C
         CALL TOLCHK (EVERR,LWA(LPP),LWA(LPC),EPS1,K,NGDEV,NACEV)
C
         IF (IPSW .GT. 2) THEN
C                                                ! Print
            NN = 1
            IF (IPSW .GT. 3) NN=2
            CALL LANPR3 (RWA(IPA),RWA(IPB),EVL,EVERR,
     &                   LWA(LPC),EPS1,NLV,NGDEV,LPU,NN)
         ENDIF
C                                                ! next eigenvalue
C                                                  extraction
         KEVEX = KEX
         IF (NEVEXT .EQ. 1) THEN
            NN = NEVAL/2
            IF (NN .GT. 0)  KEVEX = NN
         ENDIF
C                                                ! loop exit
C
         IF (NGDEV.GE.NEVAL .OR. K.EQ.MAXLAN)  GO TO 300
      ENDIF
C
C --- Next (unscaled) Lanczos vector -----------------------------------
C
      AUX =-ALFA
      CALL VASCA (R,V(1,K),1,1,AUX,N)
C                                                ! length of R
      AUX  = PRACC(R,R,1,1,N)
      BETA = SQRT(AUX)
      AUX  = ABS(RWA(IPA))
      IF (BETA .LT. RLPR*AUX) THEN
C
C ------ Breakdown; generate new (random) start vector and make it
C        ortonormal to previous Lanczos vectors, and start again
C
         CALL RANVEC (R,RAN,N)
         AUX = PRACC(R,R,1,1,N)
         AUX = ONE/SQRT(AUX)
         CALL SCALE (R,V(1,K),N,1,AUX)
         CALL MGSOR1 (V,R,N,K)
         BETA = ONE
         IP =IPB+K
         RWA(IP) = ZERO
C
         MOP(10) = MOP(10) + 1
         IF (MOP(10) .GT. 10) THEN
            CALL LANERR (31,K,K,K,LLP,IERR)
            GO TO 100
         ENDIF
         KRST    = 1
         GO TO 100
      ENDIF
C
      RWA(IPB+K) = BETA
      GO TO 100
C
C ================================================ End of Lanczos loop
C
  300 SMALL = RLPR*ABS(EVL(1))
C
      IF (NGDEV .EQ. 0) THEN
         CALL LANERR (32,MAXLAN,K,K,LLP,IERR)
         GO TO 1000
      ENDIF
      IF (NGDEV .LT. NEVAL) THEN
         CALL LANERR (1,NGDEV,MAXLAN,NEVAL,LLP,IERR)
      ENDIF
C                                                ! restore eigenvalues
C                                                  Lambda_i
      DO 350 I=1,NLV
         IF (ABS(EVL(I)) .GT. SMALL) THEN
            EVL(I) = SHIFT + ONE/EVL(I)
         ELSE
            NN = I-1
            IF (NN .LT. NGDEV) NGDEV=NN
            IF (NN .LT. NEVAL) THEN
               CALL LANERR (3,NN,I,I,LLP,IERR)
               GO TO 400
            ENDIF
         ENDIF
  350 CONTINUE
C
  400 MOP(1) = NGDEV
      MOP(2) = NACEV
      MOP(3) = NLV
      MOP(7) = NEVEXT
C
      IF (NEVEC .GT. 0) THEN
C
C --- Determine (B-orthonormal) approx. to the NVC first eigenvectors
C
         NVC = NEVEC
         IF (NVC .GT. NGDEV)  NVC = NGDEV
         MOP(5) = NVC
C                                                ! re-solve small
C                                                  eigenproblem
C
         CALL QLS3DV (RWA(IPA),RWA(IPB),RWA(IPE),NLV,LLP,IERR)
         IF (IERR .LT. 0) THEN
            CALL LANERR (25,NLV,NLV,NLV,LLP,IERR)
            GO TO 1000
         ENDIF
C                                                ! rearrange ordering
C
         CALL EVARR (RWA(IPA),RWA(IPE),NLV,NLV,NVC,4)
C
C ------ Restore eigenvectors of special problem -----------------------
C
         IF (KPD .EQ. 1) THEN
C                                                ! scale eigenvectors
            SIGMA = TRSH
            NN    = IPA+NVC-1
            DO 410 IP=IPA,NN
               IF (RWA(IP) .LT. SIGMA)  SIGMA = RWA(IP)
  410       CONTINUE
            IF (SIGMA .LT. TRSH) THEN
               SIGMA = ABS(SIGMA) + TRSH
               DO 420 IP=IPA,NN
                  RWA(IP) = RWA(IP) + SIGMA
  420          CONTINUE
            ELSE
               SIGMA = ZERO
            ENDIF
            IP = IPA
            JP = IPE
            DO 430 I=1,NVC
               S = ONE / SQRT(RWA(IP))
               CALL SCALE (RWA(JP),RWA(JP),NLV,1,S)
               RWA(IP) = RWA(IP) - SIGMA
               IP = IP+1
               JP = JP+NLV
  430       CONTINUE
         ENDIF
C
         CALL BIGMUL (V,RWA(IPE),V,R,N,NLV,NVC,0)
C
C ------ Restore desired eigenvectors ----------------------------------
C
         IF (KPD .EQ. 1) THEN
            CALL BCKSKY (A,V,MSKY,N,NVC)
         ELSE
            IF (KSB .EQ. -1) THEN
C                                                ! skyline B
               CALL BCKSKY (B,V,MSKY,N,NVC)
            ELSE
C                                                ! diagonal B
               IF (MOP(6) .EQ. 0) THEN
                  DO 450 J=1,NVC
                     DO 440 I=1,N
                        V(I,J) = V(I,J) / B(I)
  440                CONTINUE
  450             CONTINUE
               ELSE
                  DO 460 J=1,NVC
                     CALL SPSMM (A,B,V(1,J),R,R,MSKY,
     &                           N,KSA,KSB,1,LLP,IERR)
                     IF (IERR .LT. 0) THEN
                        CALL LANERR (26,J,J,J,LLP,IERR)
                        GO TO 1000
                     ENDIF
                     AUX = ONE / RWA(IPA-1+J)
                     CALL SCALE (V(1,J),V(1,J),N,1,AUX)
  460             CONTINUE
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
      IF (IPSW .GT. 0) THEN
C                                                ! Print
C
         CALL LANPR2 (EVL,EVERR,TOL,LWA(LPC),MOP,NLV,LPU)
      ENDIF
C ------------------------------------------------ Exit
C
 1000 RETURN
C
  600 FORMAT (' ')
  610 FORMAT (5X,'Lanczos step',I5,' :  reorthogonalize')
  620 FORMAT (5X,'Lanczos step',I5,' :  max. inner product =',1PE12.4)
C
      END
      SUBROUTINE LANPR1 (MSKY,MIP,SHIFT,N,NEVAL,NEVEC,MAXLAN,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LANPR1                  GROUP 5 / PRIVATE
C
C     TASK :  To print key input information to the eigensolution in
C             in subroutine LANCZ2
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-07-03 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           LPU,MAXLAN,N,NEVAL,NEVEC,    MSKY(N),MIP(6)
      DOUBLE PRECISION  SHIFT
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         WRITE (LPU,600)
         I = MSKY(N)/N
         WRITE (LPU,610) N,I
         I = ABS(MIP(2))
         IF (I .EQ. 1)  WRITE (LPU,620)
         IF (I .EQ. 2)  WRITE (LPU,630)
         WRITE (LPU,640) MIP(3),MIP(6),MIP(5),MIP(4),
     &                   NEVAL,NEVEC,MAXLAN,SHIFT
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT('1'/5X,43('*')//
     &           5X,'Print from  S A M  library routine LANCZ2'/
     &           5X,'(General eigenproblem by truncated Lanczos)'//
     &           5X,43('*')//)
  610 FORMAT(5X,'Matrix dimension  :',I5/
     &       5X,'Average bandwidth :',I5)
  620 FORMAT(5X,'Form of B-matrix  :  skyline')
  630 FORMAT(5X,'Form of B-matrix  :  diagonal')
  640 FORMAT(/5X,'Start vector code        :',I5/
     &        5X,'Factorization indicator  :',I5/
     &        5X,'Parameter KEX            :',I5/
     &        5X,'Reorthogonalization code :',I5/
     &        5X,'Requested eigenvalues    :',I5/
     &        5X,'Requested eigenvectors   :',I5/
     &        5X,'Max. number of steps     :',I5/
     &        5X,'Shift value              :',1PE11.3)
C
      END
      SUBROUTINE LANPR2 (EVL,EVERR,TOL,MANCUR,MOP,NLV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LANPR2                 GROUP 5 / PRIVATE
C
C     TASK :  To print key information on completion of eigensolution
C             in LANCZ2, including eigenvalues and their relative errors
C
C
C     ROUTINES CALLED/REFERENCED :  LANPR3     (SAM-5)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-07-03 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           LPU,NLV,       MANCUR(NLV),MOP(10)
      DOUBLE PRECISION  EVL(NLV),EVERR(NLV),TOL(2)
C
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         WRITE (LPU,600) MOP(1),MOP(2),MOP(3),MOP(4),MOP(5),
     &                   MOP(6),MOP(7),MOP(8),MOP(9),MOP(10)
         WRITE (LPU,610) TOL(1),TOL(2)
         CALL LANPR3 (EVL,EVL,EVL,EVERR,MANCUR,TOL(1),NLV,MOP(2),LPU,0)
      ENDIF
C
      RETURN
C ----------------------------------------------------------------------
C
  600 FORMAT(////5X,14('*')/5X,'Leaving LANCZ2'/5X,14('*')//
     &           5X,'No. of "good" eigenvalues             :',I5/
     &           5X,'No. of acceptable eigenvalues         :',I5/
     &           5X,'No. of Lanczos steps (vectors)        :',I5/
     &           5X,'No. of eigenvalues smaller than shift :',I5/
     &           5X,'No. of eigenvectors returned          :',I5/
     &           5X,'No. of zero elements in diag.(B)      :',I5/
     &           5X,'No. of tridiag. eigv.problems solved  :',I5/
     &           5X,'No. of Lanczos vectors for which'/
     &           5X,'reorthogonalization was carried out   :',I5/
     &           5X,'No. of vector dot products performed'/
     &           5X,'during reorthogonalization            :',I5/
     &           5X,'No. of restarts                       :',I5)
  610 FORMAT(/   5X,'Relative eigenvalue tolerance    :',1PE11.3 /
     &           5X,'Diagonal decay tolerance         :',1PE11.3 )
C
      END
      SUBROUTINE LANPR3 (ALPHA,BETA,EVL,EVERR,MANCUR,
     +                   EPS,NLV,NGDEV,LPU,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LANPR3                 GROUP 5 / PRIVATE
C
C     TASK :  To print computed eigenvalues and their corresponding
C             relative errors after a given Lanczos step (no. NLV).
C             The step at which the eigenvalues were accepted is also
C             printed.
C             Optionally (IFLAG > 1) the tridiagonal form, stored in
C             ALPHA and BETA, is also printed
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-07-031 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IFLAG,LPU,NGDEV,NLV,       MANCUR(NLV)
      DOUBLE PRECISION  EPS,  ALPHA(NLV),BETA(NLV),EVERR(NLV),EVL(NLV)
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         IF (IFLAG .GT. 0) THEN
            WRITE (LPU,600)  NLV
            IF (IFLAG .GT. 1) THEN
               WRITE (LPU,610)
               DO 25 I=1,NLV
                  WRITE (LPU,620) I,BETA(I),ALPHA(I)
   25          CONTINUE
            ENDIF
            WRITE (LPU,630) NGDEV,EPS
         ENDIF
         WRITE (LPU,640)
         DO 50 I=1,NLV
            WRITE (LPU,650) I,EVL(I),EVERR(I),MANCUR(I)
   50    CONTINUE
      ENDIF
C
      RETURN
C ----------------------------------------------------------------------
C
  600 FORMAT(////5X,21('*')/5X,'Lanczos step no.',I5 / 5X,21('*'))
  610 FORMAT(//  5X,'Tridiagonal form :' ///
     &          12X,'Subdiagonal',6X,'diagonal'/)
  620 FORMAT(I7,1P2D16.5)
  630 FORMAT(//  5X,'Number of "good" eigenvalues :',I5,'   (tol =',
     &           1PE11.3,' )')
  640 FORMAT(// 15X,'COMPUTED',13X,'RELATIVE',6X,'ACCEPTED'/
     &          15X,'EIGENVALUES',10X,'ERRORS',8X,'AT STEP'/)
  650 FORMAT(I7,1PD24.15,D15.4,I9)
C
      END
      SUBROUTINE LANERR (N,I1,I2,I3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  LANERR                  GROUP 5 / PRIVATE
C
C     TASK :  To print error/warning messages and set the error flag
C             for eigenproblem routine LANCZ2
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-07-01 / 1.0
C                       96-12-29 / 1.1  K.Bell
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,I1,I2,I3,LPU,N
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         IF (N .LT. 10) THEN
C                                                ! warnings
            WRITE (LPU,600)
            IF (N .EQ. 1) THEN
               WRITE (LPU,601) I1,I2,I3
            ELSEIF (N .EQ. 2) THEN
               WRITE (LPU,602)
            ELSEIF (N .EQ. 3) THEN
               WRITE (LPU,603) I1
            ENDIF
         ELSE
C                                                ! errors
            WRITE (LPU,690)
            IF (N.GT.9 .AND. N.LT.21) THEN
               WRITE (LPU,610)
               IF (N .EQ. 11)  WRITE (LPU,611) I1,I2
               IF (N .EQ. 12)  WRITE (LPU,612) I1
               IF (N .EQ. 13)  WRITE (LPU,613) I1
               IF (N .EQ. 14)  WRITE (LPU,614) I1
               IF (N .EQ. 15)  WRITE (LPU,615) I1,I2,I3
               IF (N .EQ. 16)  WRITE (LPU,616) I1,I2
               IF (N .EQ. 17)  WRITE (LPU,617) I1,I2
               IF (N .EQ. 18)  WRITE (LPU,618) I1,I2,I3
               IF (N .EQ. 19)  WRITE (LPU,619) I1,I2
            ELSEIF (N .EQ. 21) THEN
               WRITE (LPU,621)
            ELSEIF (N .EQ. 22) THEN
               WRITE (LPU,622) I1
            ELSEIF (N .EQ. 23) THEN
               WRITE (LPU,623)
            ELSEIF (N .EQ. 24) THEN
               WRITE (LPU,624) I1
            ELSEIF (N .EQ. 25) THEN
               WRITE (LPU,625) I1
            ELSEIF (N .EQ. 26) THEN
               WRITE (LPU,626) I1
            ELSEIF (N .EQ. 31) THEN
               WRITE (LPU,631) I1
            ELSEIF (N .EQ. 32) THEN
               WRITE (LPU,632) I1
            ENDIF
         ENDIF
      ENDIF
C
      IF (N .LT. 10)             IERR = 1
      IF (N.GT.9 .AND. N.LT.21)  IERR =-1
      IF (N.EQ.21)               IERR =-2
      IF (N.EQ.22)               IERR =-3
      IF (N.EQ.23)               IERR =-4
      IF (N.EQ.24)               IERR =-5
      IF (N.EQ.25)               IERR =-6
      IF (N.EQ.26)               IERR =-7
      IF (N.EQ.31)               IERR =-8
      IF (N.EQ.32)               IERR =-9
C
      RETURN
C ----------------------------------------------------------------------
  600 FORMAT(///' *** WARNING from  S A M  library routine LANCZ2')
  601 FORMAT(5X,'Only',I5,'  eigenvalues accepted'/
     &       5X,'in',I5,'  Lanczos steps'/
     &       5X,'(the required number was',I5,' )')
  602 FORMAT(5X,'Nonzero shift may render A non-positive definite'/
     &       5X,'and cause an error during Cholesky factorization' )
  603 FORMAT(5X,'Problem has only',I5,'  finite eigenvalues')
  610 FORMAT(5X,'Illegal or inconsistent input parameters')
  611 FORMAT(5X,'Storage codes :  KSA =',I5,5X,'KSB =',I5)
  612 FORMAT(5X,'Reorthogonalization code =',I5)
  613 FORMAT(5X,'Parameter KEX (=',I5,' ) out of range')
  614 FORMAT(5X,'Parameter KPD (=',I5,' ) out of range')
  615 FORMAT(5X,'Check parameters KPD, KSA and KSB =',3I6)
  616 FORMAT(5X,'Check parameters KSB and KSR =',2I6)
  617 FORMAT(5X,'Check parameters KSA and NEVAL =',2I6)
  618 FORMAT(5X,'Check parameters NEVAL, NEVEC and MAXLAN =',3I6)
  619 FORMAT(5X,'KSA =',I4,' and/or  KSB =',I4,'  are inconsistent'/
     &       5X,'with non-zero shift' )
C
  621 FORMAT(5X,'Error during factorization of matrix A')
  622 FORMAT(5X,'Matrix B (KSB =',I4,') cannot be factorized')
  623 FORMAT(5X,'Problem has no finite eigenvalues (N=1)')
  624 FORMAT(5X,'Error during premultiplication (by matrix H)'/
     &       5X,'Lanczos step no.',I5)
  625 FORMAT(5X,'Error during solution of small, tridiag. eigenproblem'
     &      /5X,'(of dimension',I4,' )')
  626 FORMAT(5X,'Error during recovery of eigenvector no.',I5)
  631 FORMAT(5X,'Breakdown! - 11th restart required in step no.',I5)
  632 FORMAT(5X,'No eigenvalues accepted (after',I5,'  Lanczos steps')
C
  690 FORMAT(///' *** ERROR return from  S A M  library routine LANCZ2')
C
      END
      SUBROUTINE MGSOR1 (G,V,N,M)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MGSOR1                 GROUP 5 / PRIVATE
C
C     T A S K :  To make vector V orthonormal to all M, mutually
C                orthonormal, column-vectors of G, by a modified Gram-
C                Schmidt procedure.
C
C
C     ROUTINES CALLED/REFERENCED :    PRACC, SCALE   (SAM-0)
C                                     SQRT           (Fortran library)
C
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-06-15 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           M,N
      DOUBLE PRECISION  G(N,M),V(N)
C
      INTEGER           I,J
      DOUBLE PRECISION  ONE,X,    PRACC
C
      PARAMETER         ( ONE = 1.0D0 )
C ----------------------------------------------------------------------
C
      DO 20 J=1,M
         X = PRACC(V,G(1,J),1,1,N)
         DO 10 I=1,N
            V(I) = V(I) - X*G(I,J)
   10    CONTINUE
   20 CONTINUE
      X = PRACC(V,V,1,1,N)
      X = ONE / SQRT(X)
      CALL SCALE (V,V,N,1,X)
C
      RETURN
      END
      SUBROUTINE REORL2 (ALFA,BETA,VVK,VVKM1,VVKP1,
     &                   BTAKP1,RLPR,TRSH,K,VPMAX)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  REORL2                 GROUP 5 / PRIVATE
C
C     T A S K :  To estimate the largest inner-product (VPMAX) of the
C                current (K+1)th Lanczos vector and all previous Lanczos
C                vectors, using Simon's recursion formula (see Sehmi's
C                book, p 64)
C
C
C     ROUTINES CALLED/REFERENCED :  ABS     (FORTRAN library)
C
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-06-15 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           K
      DOUBLE PRECISION  BTAKP1,RLPR,TRSH,VPMAX
      DOUBLE PRECISION  ALFA(*),BETA(*),VVK(*),VVKM1(*),VVKP1(*)
C
      INTEGER           I
      DOUBLE PRECISION  ONE,X
C
      PARAMETER         ( ONE = 1.0D0 )
C ----------------------------------------------------------------------
C
      VVK(K)     = ONE
      VVKM1(K-1) = ONE
C
      IF (K .EQ. 2) THEN
C                                               ! initialize VVK & VVKM1
         VVK(1)   = RLPR
         VVK(2)   = ONE
         VVKM1(1) = ONE
      ENDIF
      VPMAX = RLPR
C
      X = ((ALFA(1)-ALFA(K))*VVK(1) + BETA(2)*VVK(2)
     &     - BETA(K)*VVKM1(1)) / BTAKP1
      IF (ABS(X) .GT. VPMAX)  VPMAX = ABS(X)
      VVKP1(1) = X
C
      IF (K .GT. 2) THEN
         DO 20 I=2,K-1
            X = (  BETA(I)*VVK(I-1)
     &           + (ALFA(I)-ALFA(K))*VVK(I)
     &           + BETA(I+1)*VVK(I+1)
     &           - BETA(K)*VVKM1(I) ) / BTAKP1
            IF (ABS(X) .GT. VPMAX)  VPMAX = ABS(X)
            VVKP1(I) = X
   20    CONTINUE
      ENDIF
      VVKP1(K) = RLPR
      IF (VPMAX .LE. TRSH) THEN
C                                               ! prepare for next step
         DO 40 I=1,K
            VVKM1(I) = VVK(I)
            VVK(I)   = VVKP1(I)
   40    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE SPSMM (A,B,X,Y,W,MSKY,N,KSA,KSB,IFLAG,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPSMM                  GROUP 5 / PRIVATE
C
C     T A S K :  To determine
C
C                   X := (A-1)*UbT*X         IFLAG = 1
C                   Y  = (A-1)*UbT*X         IFLAG = 2
C                   X := Ub*(A-1)*UbT*X      IFLAG = 3
C                   Y  = Ub*(A-1)*UbT*X      IFLAG = 4
C                   X := (Ua-T)*B*(Ua-1)*X   IFLAG = 5
C                   Y  = (Ua-T)*B*(Ua-1)*X   IFLAG = 6
C     where
C           X and Y are vectors,
C           A is a symmetric skyline matrix, represented by its factors
C           La*Da*LaT, stored in A (KSA = 2),
C           Ub is the (Cholesky) factor of B,stored in B (KSB<0), and
C           Ua is the Cholesky factor of A, stored in A (KSA=3)
C     If  KSB=1 or -1, B / Ub is a symmetric / upper triangular skyline
C     matrix, and
C     if  KSB=2 or -2, B / Ub is a diagonal matrix
C
C
C     ROUTINES CALLED/REFERENCED :  RCOPY                    (SAM-0)
C                                   TSKYMB, PRESKY           (SAM-3)
C                                   BCKSKY, FSBSKY, SKYSOL   (SAM-4)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-06-12 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IFLAG,LPU,KSA,KSB,N,     MSKY(N)
      DOUBLE PRECISION  A(*),B(*),Y(N),X(N),W(N)
C
      INTEGER           I,LSB
      DOUBLE PRECISION  DUM
C ----------------------------------------------------------------------
      IERR = 0
      DUM = 0.0
C
      LSB = 1
      IF (ABS(KSB) .EQ. 2)  LSB = -1
C
      IF (IFLAG .EQ. 1) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GO TO 90
         CALL TSKYMB (B,X,X,W,MSKY,N,1,LSB,2,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL SKYSOL (A,X,MSKY,DUM,N,N,1,LPU,4,I,IERR)
         IF (IERR .LT. 0)   GO TO 95
      ELSEIF (IFLAG .EQ. 2) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GO TO 90
         CALL TSKYMB (B,X,W,Y,MSKY,N,1,LSB,22,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL SKYSOL (A,Y,MSKY,DUM,N,N,1,LPU,4,I,IERR)
         IF (IERR .LT. 0)   GO TO 95
      ELSEIF (IFLAG .EQ. 3) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GO TO 90
         CALL TSKYMB (B,X,Y,W,MSKY,N,1,LSB,22,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL SKYSOL (A,W,MSKY,DUM,N,N,1,LPU,4,I,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL TSKYMB (B,W,Y,X,MSKY,N,1,LSB,11,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
      ELSEIF (IFLAG .EQ. 4) THEN
         IF (KSA.NE.2 .OR. KSB.GT.0)  GO TO 90
         CALL TSKYMB (B,X,Y,W,MSKY,N,1,LSB,22,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL SKYSOL (A,W,MSKY,DUM,N,N,1,LPU,4,I,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL TSKYMB (B,W,X,Y,MSKY,N,1,LSB,11,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
      ELSEIF (IFLAG .EQ. 5) THEN
         IF (KSA.NE.3 .OR. KSB.LT.0)  GO TO 90
         CALL BCKSKY (A,X,MSKY,N,1)
         CALL PRESKY (B,X,X,W,MSKY,N,1,LSB,1,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL FSBSKY (A,X,MSKY,N,1)
      ELSEIF (IFLAG .EQ. 6) THEN
         IF (KSA.NE.3 .OR. KSB.LT.0)  GO TO 90
         CALL RCOPY (X,W,N,1,N)
         CALL BCKSKY (A,W,MSKY,N,1)
         CALL PRESKY (B,W,Y,Y,MSKY,N,1,LSB,88,LPU,IERR)
         IF (IERR .LT. 0)   GO TO 95
         CALL FSBSKY (A,Y,MSKY,N,1)
      ENDIF
      GO TO 100
C                                                ! error return
   90 IERR = -1
      IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
         WRITE (LPU,691) IFLAG,KSA,KSB
      ENDIF
      GO TO 100
   95 IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
      ENDIF
C
  100 RETURN
C
  690 FORMAT(///' *** ERROR return from  S A M library routine  SPSMM')
  691 FORMAT(5X,'Illegal parameter(s) - IFLAG, KSA, KSB = ',3I5)
C
      END
      SUBROUTINE MSSIT (A,B,TOL,MSKY,MIP,
     +                  MOP,EVL,EVERR,EVEC,RWA,
     +                  SHIFT,N,NEVAL,NIV,LPU,IPSW,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MSSIT                   GROUP 5 / PUBLIC
C
C     TASK :  To determine NEVAL eigenvalues LAMBDAi and corresponding
C             eigenvectors Qi of the general, symmetric eigenproblem
C
C                     (A -LAMBDAi*B)Qi = 0
C
C             that are closest to a user-specified "shift" of origin.
C             A is a symmetric "skyline" matrix, and B is a symmetric
C             skyline matrix (same skyline as A) or a diagonal matrix.
C             Both the original matrix A and its factors La/Da (where
C             A = La Da LaT) are accepted as input.
C             A fairly standard "subspace iteration method" (with NIV
C             iteration vectors), assuming all matrices in primary
C             storage, is used - matrix operations on converged eigen-
C             vectors may (as an option) be avoided.
C             NOTE: If the problem is small (N < 7), it is solved
C                   directly (by GENHQL), that is, without iteration.
C
C     SOME LOCAL VARIABLES:
C
C     INRAN  - "flag" signalling, if = 1, inclusion of a "new", random
C              vector in each iteration (starting with the 3rd)
C     IPAB   |
C     IPBB   |
C     IPEV   |
C     IPSC   > "pointers" in the local "work array" RWA
C     IPW    |
C     IPWJ   |
C     IPX    |
C     JS     - first vector to be "operated on"
C     KSA    - storage code for matrix A:
C              = 1 : A stores the matrix itself
C              =-1 : A stores the factors, LaT and Da
C     KSB    - storage code for matrix B:
C              = 1 : B stores a "skyline" matrix
C              = 2 : B stores a diagonal matrix
C     KSR    - start vector code:
C              = 1 : pseudo-random start vectors are used
C              = 2 : special start vectors, based on diag(A) and
C                    diag(B), are used
C              = all other values :  start vectors are input (in EVEC)
C     KSKY   - storage code for matrix B:
C              = 0 : B stores a diagonal matrix
C              = 1 : B stores a "skyline" matrix
C     LLP    - "local" print unit number (may be negative)
C     MAXIT  - maximum number of iteration cycles (default = 20)
C     NACEV  - number of accepted eigenpairs
C     NFEV   - number of finite eigenvalues
C     NJ     - number of iteration-vectors to be "operated on"
C     NZDB   - number of zero-elements on diagonal of B
C     SIGMA1 |
C     SIGMA2 > local "shift" values
C
C
C     ROUTINES CALLED/REFERENCED :  IMP,   IMINT, RMINT        (SAM-0)
C                                   RCOPY, SCALE, SCADD        (SAM-0)
C                                   PRESKY                     (SAM-3)
C                                   SKYSOL                     (SAM-4)
C                                   MDFCTA, CHBDG, SSSSVC      (SAM-5)
C                                   GENHQL, EVARR, MSSER1      (SAM-5)
C                                   MSSPR1, MSSPR2, MSSPR3     (SAM-5)
C                                   MPSTML, MSMATB, EVARR      (SAM-5)
C                                   SKYCNV, RANVEC             (SAM-8)
C                                   ABS, SQRT        (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-26 / 1.0
C                       97-10-17 / 1.1   K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,IPSW,LPU,N,NEVAL,NIV,
     &                  MIP(6),MOP(4),MSKY(*)
      DOUBLE PRECISION  SHIFT,A(*),B(*),EVEC(N,NIV),EVERR(NIV),EVL(NIV),
     &                  RWA(*),TOL(2)
C
C                                                ! local variables
C
      INTEGER           I,INRAN,IP,IPAB,IPBB,IPEV,IPSC,IPW,IPWJ,IPX,IT,
     &                  JP,JS,KSA,KSB,KSKY,KSR,LLP,MAXIT,NACEV,NFEV,NJ,
     &                  NSDIG,NZDB,
     &                  IMP
      DOUBLE PRECISION  EPS1,EPS2,ERR,ONE,RAN,S,SIGMA1,SIGMA2,TEN,ZERO
C
      PARAMETER         (ZERO = 0.0D0 , ONE = 1.0D0 , TEN = 10.0D0)
C
      EXTERNAL          IMP
C ----------------------------------------------------------------------
C
      IERR = 0
C
C --- Preliminaries ----------------------------------------------------
C
      IF (IPSW.GT.0 .AND. NEVAL.GT.0) THEN
C                                                ! print input param.
C
         CALL MSSPR1 (MSKY,MIP,SHIFT,N,NEVAL,NIV,LPU)
      ENDIF
C                                                ! initiate parameters
      LLP = LPU
      IF (IPSW .LT. 0)   LLP =-1
      KSA   = MIP(1)
      KSB   = MIP(2)
      KSR   = MIP(3)
      MAXIT = MIP(4)
      INRAN = 0
      IF (KSR .LT. 0) THEN
         INRAN = 1
         KSR   =-KSR
      ENDIF
      NSDIG = IMP(3)
      RAN   = 0.1234567891234567D0
      KSKY  = KSB
      IF (KSB .EQ. 2)  KSKY  = 0
      IF (KSA .LT. 0)  SHIFT = ZERO
      CALL IMINT (MOP,4,1,0)
C                                                ! default values ?
      IF (MAXIT .LT. 1)  MAXIT = 20
      I    = 2*NSDIG/3
      EPS1 = TOL(1)
      IF (EPS1 .LE. ZERO) THEN
         EPS1   = TEN**(-I)
         TOL(1) = EPS1
      ENDIF
      EPS2 = TOL(2)
      IF (EPS2 .LE. ZERO) THEN
         EPS2   = TEN**(-I)
         TOL(2) = EPS2
      ENDIF
C                                                ! check input
C
      IF (N .LT. 1)           CALL MSSER1 (11,N,N,N,LLP,IERR)
      IF (ABS(KSA) .NE. 1)    CALL MSSER1 (12,KSA,KSB,KSB,LLP,IERR)
      IF (IERR .LT. 0)  GO TO 1000
      IF (NEVAL .LT. 0) THEN
         CALL MSSER1 (14,NEVAL,NIV,NIV,LLP,IERR)
         GO TO 1000
      ELSEIF (NEVAL .EQ. 0) THEN
         IF (KSA .EQ. (-1)) THEN
            CALL MSSER1 (13,KSA,NEVAL,NEVAL,LLP,IERR)
            GO TO 1000
         ENDIF
         IPW  = 1
         IPX  = 1
         IPEV = 1
         IPAB = 1
         IPBB = 1
         IPSC = 1
         GO TO 200
      ENDIF
      IF (KSB.LT.1 .OR. KSB.GT.2) THEN
         CALL MSSER1 (12,KSA,KSB,KSB,LLP,IERR)
         GO TO 1000
      ENDIF
C
      CALL CHBDG (B,MSKY,N,KSB,NZDB,I)
      IF (I .LT. 0) THEN
         IF (MIP(5) .EQ. 2) THEN
            CALL MSSER1 (22,I,I,I,LLP,IERR)
            GO TO 1000
         ENDIF
      ENDIF
      IF (KSB .EQ. 2) THEN
         NFEV = N-NZDB
      ELSE
         NFEV = N
      ENDIF
      IF (NFEV .EQ. 0) THEN
         CALL MSSER1 (26,I,I,I,LLP,IERR)
         GO TO 1000
      ENDIF
      MOP(3) = NZDB
C
      IF (N .EQ. 1) THEN
C
C ------ Special case:  N = 1 ------------------------------------------
C
         EVL(1)    = A(1) / B(1)
         EVEC(1,1) = ONE / SQRT(B(1))
         MOP(1)    = 1
         IF (EVL(1) .LT. SHIFT)  MOP(2) = 1
         GO TO 1000

      ELSEIF (N .LT. 7) THEN
C
C ------ Small problem -  N < 7 ----------------------------------------
C
         IF (NIV .GT. N) NIV = N
         IF (KSA .EQ. 1) THEN
C                                                ! convert matrices
            CALL SKYCNV (A,MSKY,RWA(7),N,KSA)
            CALL SKYCNV (B,MSKY,RWA(43),N,KSB)
            IF (MIP(5) .EQ. 2) THEN
               IF (NZDB .EQ. 0) THEN
C                                                ! solve orign. problem
C
                  CALL GENHQL (RWA(7),RWA(43),RWA(1),RWA(79),
     &                         RWA(115),EPS2,N,N,2,LLP,IERR)
                  IF (IERR .LT. 0) THEN
                     CALL MSSER1 (30,I,I,I,LLP,IERR)
                     GO TO 1000
                  ENDIF
               ELSE
                  CALL MSSER1 (27,I,I,I,LLP,IERR)
                  GO TO 1000
               ENDIF
            ELSE
C                                                ! solve inverse problem
C
               CALL GENHQL (RWA(43),RWA(7),RWA(1),RWA(79),
     &                      RWA(115),EPS2,N,N,4,LLP,IERR)
               IF (IERR .LT. 0) THEN
                  CALL MSSER1 (30,I,I,I,LLP,IERR)
                  GO TO 1000
               ENDIF
               SIGMA1 = EPS1
               DO 10 I=1,N
                  IF (RWA(I) .LT. SIGMA1)  SIGMA1 = RWA(I)
   10          CONTINUE
               IF (SIGMA1 .LT. EPS1) THEN
                  SIGMA1 = ABS(SIGMA1) + EPS1
                  DO 20 I=1,N
                     RWA(I) = RWA(I) + SIGMA1
   20             CONTINUE
               ELSE
                  SIGMA1 = ZERO
               ENDIF
C                                                ! scale eigenvectors
               IP = 79
               DO 30 I=1,N
                  S = ONE / SQRT(RWA(I))
                  CALL SCALE (RWA(IP),RWA(IP),N,1,S)
                  IP = IP+N
   30          CONTINUE
               DO 40 I=1,N
                  RWA(I) = ONE / (RWA(I) - SIGMA1)
   40          CONTINUE
            ENDIF
            CALL RCOPY (RWA(1),EVL,NIV,1,1)
            CALL RCOPY (RWA(79),EVEC,N,NIV,1)
            CALL RMINT (EVERR,NIV,1,ZERO)
            IF (SHIFT .NE. ZERO) THEN
               DO 50 I=1,NIV
                  IF (RWA(I) .LT. SHIFT)  MOP(2) = MOP(2)+1
   50          CONTINUE
               DO 60 I=1,NIV
                  EVL(I) = EVL(I) - SHIFT
   60          CONTINUE
               CALL EVARR (EVL,EVEC,N,NIV,NIV,2)
               DO 70 I=1,NIV
                  EVL(I) = EVL(I) + SHIFT
   70          CONTINUE
            ENDIF
            MOP(1) = NIV
            GO TO 800
         ENDIF
      ENDIF
C -------------------------------------------------------------------
C                                                ! more checking
      IF (NEVAL.GT.1) THEN
         IF (NEVAL.EQ.NIV)  CALL MSSER1 (14,NEVAL,NIV,NIV,LLP,IERR)
      ENDIF
      IF (NEVAL .GT. NIV)   CALL MSSER1 (14,NEVAL,NIV,NIV,LLP,IERR)
      IF (NIV .GT. NFEV)    CALL MSSER1 (15,NIV,NFEV,NFEV,LLP,IERR)
      IF (IERR.LT.0)        GO TO 1000

C                                                ! pointers in RWA
      IPW  = 1
      IPX  = IPW + N*NIV
      IPEV = IPX + NIV*NIV
      IPAB = IPEV + NIV
      IPBB = IPAB + NIV*NIV
      IPSC = IPBB + NIV*NIV
C                                                ! initialize (to zero)
      CALL RMINT (EVL,NIV,1,ZERO)
C
      IF (KSR .EQ. 1) THEN
C
C ------ Prepare random start vectors
C
         DO 100 I=1,NIV
            CALL RANVEC (EVEC(1,I),RAN,N)
  100    CONTINUE
C
      ELSEIF (KSR .EQ. 2) THEN
C
         IF (KSA.EQ.1) THEN
C
C ------ Prepare start vectors based on diag(A) and diag(B)
C
            CALL SSSSVC (A,B,EVEC,RWA(IPSC),MSKY,RAN,N,NIV,KSB)
         ELSE
C                                                ! random vectors
            DO 110 I=1,NIV
               CALL RANVEC (EVEC(1,I),RAN,N)
  110       CONTINUE
         ENDIF
      ENDIF
C
  200 CONTINUE
C
      IF (KSA .EQ. 1) THEN
C
C ------ Modify matrix A (with respect to "shift") and factorize it
C
         CALL MDFCTA (A,B,MSKY,EPS2,SHIFT,N,KSB,LLP,I,IERR)
         IF (IERR .LT. 0) THEN
            CALL MSSER1 (21,I,I,I,LLP,IERR)
            GO TO 1000
         ENDIF
         KSA    =-1
         MOP(2) = I
C                                                ! exit if NEVAL=0
         IF (NEVAL .EQ. 0)  GO TO 1000
      ENDIF
C
C --- S u b s p a c e   i t e r a t i o n ------------------------------
C
      NACEV = 0
      JS    = 1
C                                                ! compute W
      CALL PRESKY (B,EVEC,RWA(IPW),RWA(IPSC),
     &             MSKY,N,NIV,KSKY,1,LLP,IERR)
      IF (IERR .LT. 0) THEN
         IT = 0
         CALL MSSER1 (23,IT,IT,IT,LLP,IERR)
         GO TO 1000
      ENDIF
C
      DO 500 IT=1,MAXIT
C
         IPWJ = IPW + N*(JS-1)
         NJ   = NIV + 1 - JS
C                                                ! copy W into V (EVEC)
         CALL RCOPY (RWA(IPWJ),EVEC(1,JS),N,NJ,1)
C                                                ! solve for Ritz vects.
         CALL SKYSOL (A,EVEC(1,JS),MSKY,EPS2,N,N,
     &                NJ,LLP,4,I,IERR)
         IF (IERR .LT. 0) THEN
            CALL MSSER1 (23,IT,IT,IT,LLP,IERR)
            GO TO 1000
         ENDIF
C                                                ! subspace matrices
         CALL MSMATB (EVEC,RWA(IPW),RWA(IPAB),
     &                EVL,N,NIV,JS,0)
C                                                ! compute W
C
         CALL PRESKY (B,EVEC,RWA(IPW),RWA(IPSC),
     &                MSKY,N,NIV,KSKY,1,LLP,IERR)
         IF (IERR .LT. 0) THEN
            CALL MSSER1 (23,IT,IT,IT,LLP,IERR)
            GO TO 1000
         ENDIF
         CALL MSMATB (EVEC,RWA(IPW),RWA(IPBB),
     &                EVL,N,NIV,JS,1)
C
C ------ Solve reduced eigenproblem
C
         IF (MIP(5) .EQ. 2) THEN
C                                                ! original problem
C
            CALL GENHQL (RWA(IPAB),RWA(IPBB),RWA(IPEV),RWA(IPX),
     &                   RWA(IPSC),EPS2,NIV,NIV,2,LLP,IERR)
            IF (IERR.LT.0) THEN
               CALL MSSER1 (24,IT,IT,IT,LLP,IERR)
               GO TO 1000
            ENDIF
         ELSE
C                                                ! inverse problem
            IF (JS .GT. 1) THEN
               IP = IPAB
               JP = IPBB
               DO 240 I=1,JS-1
                  IF (RWA(IP) .LT. ZERO) THEN
C                                                ! change sign of
C                                                  diagonal elements
                     RWA(IP) = -RWA(IP)
                     RWA(JP) = -RWA(JP)
                  ENDIF
                  IP = IP+NIV+1
                  JP = JP+NIV+1
  240          CONTINUE
            ENDIF
C
            IF (SHIFT.NE.ZERO .AND. MOP(2).GT.0) THEN
C
C                                                  shift
               SIGMA2 = SHIFT
               CALL SCADD (RWA(IPAB),RWA(IPBB),
     &                     RWA(IPAB),NIV,NIV,SIGMA2)
            ELSE
               SIGMA2 = ZERO
            ENDIF
            CALL GENHQL (RWA(IPBB),RWA(IPAB),RWA(IPEV),RWA(IPX),
     &                   RWA(IPSC),EPS2,NIV,NIV,2,LLP,IERR)
            IF (IERR.LT.0) THEN
               CALL MSSER1 (28,IT,IT,IT,LLP,IERR)
               GO TO 1000
            ENDIF
C
            SIGMA1 = EPS1
            IP     = IPEV
            DO 250 I=1,NIV
               IF (RWA(IP) .LT. SIGMA1)  SIGMA1 = RWA(IP)
               IP = IP+1
  250       CONTINUE
            IF (SIGMA1 .LT. EPS1) THEN
               SIGMA1 = ABS(SIGMA1) + EPS1
               IP     = IPEV
               DO 260 I=1,NIV
                  RWA(IP) = RWA(IP) + SIGMA1
                  IP = IP+1
  260          CONTINUE
            ELSE
               SIGMA1 = ZERO
            ENDIF
C                                                ! scale eigenvectors
            JP = IPX
            IP = IPEV
            DO 300 I=1,NIV
               S  = ONE/SQRT(RWA(IP))
               CALL SCALE (RWA(JP),RWA(JP),NIV,1,S)
               IP = IP+1
               JP = JP+NIV
  300       CONTINUE
C
            IP = IPEV
            DO 320 I=1,NIV
               RWA(IP) = ONE/(RWA(IP) - SIGMA1) - SIGMA2
               IP = IP+1
  320       CONTINUE
            CALL EVARR (RWA(IPEV),RWA(IPX),NIV,NIV,NIV,2)
         ENDIF
C
C ------ Eigenvalue estimates and their relative errors
C
         IP = IPEV-1+JS
         DO 350 I=JS,NIV
            IF (RWA(IP) .EQ. ZERO) THEN
               IF (EVL(I) .EQ. ZERO) THEN
                  ERR = ZERO
               ELSE
                  ERR = ONE
               ENDIF
            ELSE
               ERR = ABS(RWA(IP)-EVL(I)) / ABS(RWA(IP))
            ENDIF
            EVL(I)   = RWA(IP)
            EVERR(I) = ERR
            IP = IP+1
  350    CONTINUE
C
C ------ Check convergence
C
         IF (MIP(6) .NE. 1)  NACEV = 0
         DO 400 I=JS,NIV
            IF (EVERR(I) .GT. EPS1)  GO TO 410
            NACEV = NACEV+1
  400    CONTINUE
  410    MOP(1) = NACEV
         MOP(4) = IT
C
C ------ Improved eigenvectors
C
         IF (MIP(6) .EQ. 1) THEN
            IF ((NACEV-JS) .GE. 0) THEN
               CALL MPSTML (RWA(IPX),EVEC,RWA(IPSC),N,NIV,JS,NACEV)
               IF (NACEV .GE. NEVAL) THEN
C                                                ! exit iteration loop
                  GO TO 600
               ELSE
                  NJ   = NACEV-JS+1
                  IPWJ = IPW + N*(JS-1)
                  CALL RCOPY (EVEC(1,JS),RWA(IPWJ),N,NJ,1)
                  JS   = NACEV+1
               ENDIF
            ENDIF
         ELSEIF (NACEV.GE.NEVAL .OR. IT.EQ.MAXIT) THEN
            IF (NACEV .GT. 0) THEN
               CALL MPSTML (RWA(IPX),EVEC,RWA(IPSC),N,NIV,1,NACEV)
               IF (NACEV .GE. NEVAL) THEN
C                                                ! exit iteration loop
                  GO TO 600
               ELSE
                  CONTINUE
               ENDIF
            ENDIF
         ENDIF
         IF (IT .LT. MAXIT) THEN
C
C ------ Compute new Ritz-vectors
C
            CALL MPSTML (RWA(IPX),RWA(IPW),RWA(IPSC),N,NIV,JS,NIV)
         ENDIF
C
         IF (IPSW .GT. 1) THEN
C                                                ! print
C
            CALL MSSPR3 (RWA(IPEV),EVERR,EPS1,NIV,IT,NACEV,LLP,1)
         ENDIF
C
         IF (INRAN.EQ.1 .AND. IT.GT.2) THEN
            IF (NIV .GT. 5) THEN
C
C --------- Include a random vector
C
               CALL RANVEC (EVEC(1,NIV),RAN,N)
            ENDIF
         ENDIF
  500 CONTINUE
C ------------------------------------------------ end of iteration loop
      IF (NACEV .EQ. 0) THEN
         CALL MSSER1 (25,MAXIT,MAXIT,MAXIT,LLP,IERR)
      ELSE
         CALL MSSER1 (1,NACEV,MAXIT,NEVAL,LLP,IERR)
      ENDIF
C
  600 IF (SHIFT .NE. ZERO) THEN
C
C ------ Restore eigenvalues (lambda-i)
C
         DO 700 I=1,NIV
            EVL(I) = EVL(I) + SHIFT
  700    CONTINUE
      ENDIF
C
  800 IF (IPSW .GT. 0) THEN
C                                                ! print
         CALL MSSPR2 (EVL,EVERR,TOL,MOP,NIV,LPU)
      ENDIF
C
 1000 RETURN
      END
      SUBROUTINE MPSTML (B,C,WA,M,N,JS,JE)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MPSTML                  GROUP 5 / PRIVATE
C
C     TASK : To determine columns JS to JE of matrix C postmultiplied
C            by a square matrix B; the result is stored in C (other
C            columns of C are not affected).
C     This subroutine is a special version of SAM/MAT routine PSTMUL.
C
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           JE,JS,M,N
      DOUBLE PRECISION  B(N,N),C(M,N),WA(N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      DO 100 I=1,M
         DO 25 J=JS,JE
            WA(J) = PRACC(C(I,1),B(1,J),M,1,N)
   25    CONTINUE
         DO 50 J=JS,JE
            C(I,J) = WA(J)
   50    CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE MSMATB (A,B,C,V,M,N,JS,KODE)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MSMATB                  GROUP 5 / PRIVATE
C
C     TASK : To multiply two rectangular matrices, A-transpose and B,
C            assuming the product C to be symmetric.
C            If  JS > 1 and KODE=0, the first (JS-1)*(JS-1) diagonal
C                                   submatrix is set equal to a diagonal
C                                   matrix such that
C                                      C(I,I) = V(I), I=1,..,JS-1
C            If  JS > 1 and KODE=1, the first (JS-1)*(JS-1) diagonal
C                                   submatrix is set equal to I
C     This subroutine is special version of SAM/MAT routine SMATB
C
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           JS,KODE,M,N
      DOUBLE PRECISION  A(M,N),B(M,N),C(N,N),V(N)
C
      INTEGER           I,J
      DOUBLE PRECISION  ONE,ZERO,     PRACC
C
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0 )
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1 .OR. N.LT.1)  GO TO 100
C
      DO 20 J=JS,N
         DO 10 I=1,J
            C(I,J) = PRACC(A(1,I),B(1,J),1,1,M)
   10    CONTINUE
   20 CONTINUE
      IF (JS .GT. 1) THEN
         DO 40 J=1,JS-1
            DO 30 I=1,J
               C(I,J) = ZERO
   30       CONTINUE
            IF (KODE .EQ. 0) THEN
               C(J,J) = V(J)
            ELSE
               C(J,J) = ONE
            ENDIF
   40    CONTINUE
       ENDIF
C                                                ! symmetry
       DO 60 J=1,N
         DO 50 I=J,N
            C(I,J) = C(J,I)
   50    CONTINUE
   60 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE MSSPR1 (MSKY,MIP,SHIFT,N,NEVAL,NIV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MSSPR1                  GROUP 5 / PRIVATE
C
C     TASK :  To print key input information to the eigensolution in
C             subroutine MSSIT.
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           LPU,N,NEVAL,NIV,    MSKY(N),MIP(6)
      DOUBLE PRECISION  SHIFT
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         WRITE (LPU,600)
         I = MSKY(N)/N
         WRITE (LPU,610) N,I
         IF (MIP(2) .EQ. 1) THEN
            WRITE (LPU,620)
         ELSEIF (MIP(2) .EQ. 2) THEN
            WRITE (LPU,630)
         ENDIF
         WRITE (LPU,640) MIP(3),NEVAL,NIV,MIP(4),MIP(5),SHIFT
         IF (MIP(6) .EQ. 1) THEN
            WRITE (LPU,650)
         ENDIF
      ENDIF
C
      RETURN
C ------------------------------------------------ formats
C
  600 FORMAT('1'/5X,43('*')//
     &           5X,'PRINT from  S A M  library routine MSSIT'/
     &          13X,'general eigenproblem by'/
     &          15X,'SUBSPACE ITERATION'//
     &           5X,43('*')//)
  610 FORMAT(5X,'Matrix dimension  :',I5/
     &       5X,'Average bandwidth :',I5)
  620 FORMAT(5X,'Form of B-matrix  :  skyline')
  630 FORMAT(5X,'Form of B-matrix  :  diagonal')
  640 FORMAT(/5X,'Start vector code         :',I5/
     &        5X,'Requested eignevalues     :',I5/
     &        5X,'No. of iteration vectors  :',I5/
     &        5X,'Max. number of iterations :',I5/
     &        5X,'Pos.def. indicator        :',I5//
     &        5X,'Shift value               :',1PE11.3)
  650 FORMAT(/5X,'Converged eigenvectors are excluded from iteration')
C
      END
      SUBROUTINE MSSPR2 (EVL,EVERR,TOL,MOP,NIV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MSSPR2                  GROUP 5 / PRIVATE
C
C     TASK :  To print key information on completion of eigensolution in
C             MSSIT, including eigenvalues and their relative errors.
C
C
C     ROUTINES CALLED/REFERENCED :  MSSPR3         (SAM-5)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           LPU,NIV,            MOP(4)
      DOUBLE PRECISION  EVL(NIV),EVERR(NIV),TOL(2)
C
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         WRITE (LPU,600) MOP(1),MOP(2),MOP(3),MOP(4)
         WRITE (LPU,610) TOL(1),TOL(2)
         CALL MSSPR3 (EVL,EVERR,TOL(1),NIV,MOP(4),MOP(1),LPU,0)
      ENDIF
C
      RETURN
C ------------------------------------------------ formats
C
  600 FORMAT(////5X,13('*')/5X,'Leaving MSSIT'/5X,13('*')//
     &           5X,'No. of eigenvalues accepted           :',I5/
     &           5X,'No. of eigenvalues smaller than SHIFT :',I5/
     &           5X,'No. of zero elements in diag.(B)      :',I5/
     &           5X,'No. of iterations                     :',I5)
  610 FORMAT(/   5X,'Limit for relative eigenvalue error   :',1PE11.3 /
     &           5X,'Limit value for diagonal decay test' /
     &           5X,'during factorization                  :',1PE11.3 )
C
      END
      SUBROUTINE MSSPR3 (EVL,EVERR,EPS,NIV,IT,NACEV,LPU,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MSSPR3                  GROUP 5 / PRIVATE
C
C     TASK :  To print computed eigenvalues and their corresponding
C             relative errors after
C               - iteration no. IT    (if IFLAG > 0)
C               - the last iteration  (if IFLAG = 0)
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IFLAG,IT,LPU,NACEV,NIV
      DOUBLE PRECISION  EPS,  EVERR(NIV),EVL(NIV)
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         IF (IFLAG .GT. 0) THEN
            WRITE (LPU,600) IT
            WRITE (LPU,610) NACEV,EPS
         ENDIF
         WRITE (LPU,620)
         DO 10 I=1,NIV
            WRITE (LPU,630) I,EVL(I),EVERR(I)
   10    CONTINUE
      ENDIF
C
      RETURN
C ------------------------------------------------ formats
C
  600 FORMAT(////5X,18('*')/5X,'Iteration no.',I5 / 5X,18('*'))
  610 FORMAT(//  5X,'Number of accepted eigenvalues :',I5,
     +       '   (max. rel. error =',1PE11.3,' )')
  620 FORMAT(// 15X,'Computed',12X,'Relative'/
     +          14X,'eigenvalues',11X,'error'/)
  630 FORMAT(I7,1PD23.14,D15.4)
C
      END
      SUBROUTINE MSSER1 (N,I1,I2,I3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MSSER1                 GROUP 5 / PRIVATE
C
C     TASK :  To print error/warning messages and set the error flag
C             for eigenvalue routine MSSIT
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-30 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,I1,I2,I3,LPU,N
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         IF (N .LT. 10) THEN
C                                                ! warnings
            WRITE (LPU,600)
            IF (N .EQ. 1)  WRITE (LPU,601) I1,I2,I3
         ELSE
C                                                ! errors
            WRITE (LPU,690)
            IF (N.GT.10 .AND. N.LT.20) THEN
               WRITE (LPU,610)
               IF (N .EQ. 11)  WRITE (LPU,611) I1
               IF (N .EQ. 12)  WRITE (LPU,612) I1
               IF (N .EQ. 13)  WRITE (LPU,613) I1,I2
               IF (N .EQ. 14)  WRITE (LPU,614) I1,I2
               IF (N .EQ. 15)  WRITE (LPU,615) I1,I2
            ENDIF
            IF (N .EQ. 21)  WRITE (LPU,621)
            IF (N .EQ. 22)  WRITE (LPU,622)
            IF (N .EQ. 23)  WRITE (LPU,623) I1
            IF (N .EQ. 24)  WRITE (LPU,624) I1
            IF (N .EQ. 25)  WRITE (LPU,625) I1
            IF (N .EQ. 26)  WRITE (LPU,626)
            IF (N .EQ. 27)  WRITE (LPU,627)
            IF (N .EQ. 28)  WRITE (LPU,628) I1
            IF (N. EQ. 30)  WRITE (LPU,630)
         ENDIF
      ENDIF
C
      IF (N .EQ. 1)   IERR = 1
      IF (N.GT.10 .AND. N.LT.20)  IERR =-1
      IF (N .EQ. 21)  IERR =-2
      IF (N .EQ. 22)  IERR =-3
      IF (N .EQ. 23)  IERR =-4
      IF (N .EQ. 24)  IERR =-5
      IF (N .EQ. 25)  IERR =-6
      IF (N .EQ. 26)  IERR =-7
      IF (N .EQ. 27)  IERR =-3
      IF (N .EQ. 28)  IERR =-5
      IF (N .EQ. 30)  IERR =-5
C
      RETURN
C ----------------------------------------------------------------------
  600 FORMAT(///' *** WARNING from  S A M  library routine MSSIT')
  601 FORMAT(5X,'Only',I5,'  eigenvalues accepted'/
     &       5X,'in',I5,'  iterations'/
     &       5X,'(the requested number was',I5,'  )')
  610 FORMAT(5X,'Illegal or inconsistent input parameters')
  611 FORMAT(5X,'Dimension too small ( N =',I5,' )')
  612 FORMAT(5X,'Storage code :  KSA =',I5)
  613 FORMAT(5X,'KSA =',I5,'  and  NEVAL =',I5,'  are inconsistent')
  614 FORMAT(5X,'NEVAL =',I5,5X,'NIV =',I5)
  615 FORMAT(5X,'Parameter  NIV (=',I5,' ) exceeds the number of'/
     &       5X,'finite eigenvalues (=',I4,' )')
C
  621 FORMAT(5X,'Error during factorization of matrix A')
  622 FORMAT(5X,'Matrix B is indefinite')
  623 FORMAT(5X,'Error during matrix multiplication'/
     &       5X,'or substitution in iteration no.',I4)
  624 FORMAT(5X,'Error during solution of reduced eigenproblem'
     &      /5X,'in iteration no.',I4)
  625 FORMAT(5X,'No eigenvalues accepted (after',I4,' iterations')
  626 FORMAT(5X,'Problem has no finite eigenvalues')
  627 FORMAT(5X,'Matrix B is not positive definite')
  628 FORMAT(5X,'Error during solution of reduced (inverse) eigenproblem
     &',   / 5X,'in iteration no.',I4)
  630 FORMAT(5X,'Error during direct (GENHQL) solution')
C
  690 FORMAT(///' *** ERROR return from  S A M  library routine MSSIT')
C
      END
      SUBROUTINE GENHQL (A,B,EVAL,EVEC,WA,EPS,N,NVEC,LORD,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  GENHQL                  GROUP 5 / PUBLIC
C
C     TASK:  To determine all eigenvalues and all or none eigenvectors
C            of the generalized symmetric eigenproblem
C               (A - EVAL*B) EVEC = 0
C            by
C               1) reduction to special form (Cholesky),
C               2) reduction to tridiagonal form (Householder),
C               3) eigenvalue/eigenvector "extraction" by implicit QL,
C               4) transformation of eigenvectors back to original
C                  problem, if relevant.
C            A and B are full symmetric matrices, and B is also positive
C            definite.
C            Eigenvalues/eigenvectors are ordered according to argument
C            LORD, and the eigenvectors are, if computed, B-orthonormal.
C
C     ROUTINES CALLED/REFERENCED :  PRACC                        (SAM-0)
C                                   FACHOL, FSCHOL, BSCHOL       (SAM-4)
C                                   HOUS3D, QLS3D, QLS3DV, EVARR (SAM-5)
C                                   REARRV                       (SAM-8)
C                                   SQRT               (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   93-04-13 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LORD,LPU,N,NVEC
      DOUBLE PRECISION  EPS, A(N,N),B(N,N),EVAL(N),EVEC(N,N),WA(*)
C
C                                                ! local variables
C
      INTEGER           I,IP1,IP2,J,K,KP1,L,NM1,NV
      DOUBLE PRECISION  ONE,S,ZERO,    PRACC
C
      PARAMETER         (ZERO = 0.0D0, ONE = 1.0D0)
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (N .GT. 0) THEN
         IF (N .EQ. 1) THEN
            IF (B(1,1) .LT. ZERO) THEN
               IERR = -2
               GO TO 90
            ELSEIF (B(1,1) .EQ. ZERO) THEN
               IERR = -3
               GO TO 90
            ENDIF
            EVAL(1) = A(1,1) / B(1,1)
            IF (NVEC .EQ. 1) THEN
               EVEC(1,1) = ONE / SQRT(B(1,1))
            ENDIF
         ELSE
C                                                ! preliminaries
            NV = 0
            IF (NVEC .EQ. N)  NV = N
            IP1 = 1
            IP2 = IP1+N
C
C --------- reduction to special form (Cholesky)
C
            CALL FACHOL (B,EPS,N,LPU,IERR)
            IF (IERR .LT. 0)  GO TO 90
            CALL FSCHOL (B,A,N,N,1)
            CALL FSCHOL (B,A,N,N,-1)
C
C --------- reduction to tridiagonal form (Householder)
C
            CALL HOUS3D (A,EVAL,WA(IP1),WA(IP2),N)
C
C --------- compute and arrange in specified order all eigenvalues and
C           if requested all eigenvectors
C
            IF (NV .EQ. N) THEN
               CALL QLS3DV (EVAL,WA(IP1),EVEC,N,LPU,IERR)
               IF (IERR .LT. 0) THEN
                  IERR = -4
                  GO TO 90
               ENDIF
               CALL EVARR (EVAL,EVEC,N,NV,NV,LORD)
            ELSE
               CALL QLS3D (EVAL,WA(IP1),N,LPU,IERR)
               IF (IERR .LT. 0) THEN
                  IERR = -4
                  GO TO 90
               ENDIF
               CALL REARRV (EVAL,N,LORD)
            ENDIF
C
            IF (NV .EQ. N) THEN
C
C ------------ transform eigenvectors
C
               IF (N .GT. 2) THEN
C
C --------------- to eigenvectors of special problem
C
                  NM1 = N-1
                  DO 30 J=1,NV
                     DO 20 I=2,NM1
                        K   = N-I
                        KP1 = K+1
                        S   = PRACC(A(KP1,K),EVEC(KP1,J),1,1,I)
                        DO 10 L=KP1,N
                           EVEC(L,J) = EVEC(L,J) - S*A(L,K)
   10                   CONTINUE
   20                CONTINUE
   30             CONTINUE
               ENDIF
C
C ------------ to eigenvectors of original (generalized) problem
C
               CALL BSCHOL (B,EVEC,N,N,1)
            ENDIF
         ENDIF
      ELSE
         IERR = -1
         GO TO 90
      ENDIF
      GO TO 100
C ------------------------------------------------ error exit
   90 IF (LPU .GT. 0) THEN
         WRITE (LPU,690)
         IF (IERR .EQ. -1)  WRITE (LPU,691)
         IF (IERR .EQ. -2)  WRITE (LPU,692)
         IF (IERR .EQ. -3)  WRITE (LPU,693)
         IF (IERR .EQ. -4)  WRITE (LPU,694)
      ENDIF
C
  100 RETURN
C ----------------------------------------------------------------------
  690 FORMAT(///' *** ERROR RETURN from  S A M  LIBRARY routine GENHQL')
  691 FORMAT(5X,'Illegal matrix dimension:  N =',I6)
  692 FORMAT(5X,'Matrix B is not pos.definite')
  693 FORMAT(5X,'Matrix B is singular or near-singular')
  694 FORMAT(5X,'An insufficient number of eigenvalues have converged')
C
      END
      SUBROUTINE GENJAC (A,B,EVAL,EVEC,TOL,N,NVEC,LORD,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  GENJAC                  GROUP 5 / PUBLIC
C
C     T A S K :  TO COMPUTE ALL EIGENVALUES (EVAL) AND, OPTIONALLY (IF
C                NVEC=N), ALL CORRESPONDING EIGENVECTORS (STORED IN
C                EVEC) OF THE GENERAL, SYMMETRIC EIGENPROBLEM
C                     (A - EVAL*B)*EVEC = 0
C                BY A GENERALIZED JACOBI PROCEDURE (DUE TO FALK AND
C                LANGEMEYER).
C     NOTE:  THE PROCEDURE MAY BREAK DOWN IF MATRIX  B  IS NOT POSITIVE
C            DEFINITE.
C
C     ROUTINES CALLED/REFERENCED :  IMP,RMINT,SCALE AND PRACC    (SAM-0)
C                                   EVARR                        (SAM-5)
C                                   REARRV                       (SAM-8)
C                                   ABS,MAX,SIGN,SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-08-14 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,LORD,LPU,N,NVEC
      DOUBLE PRECISION  TOL,A(N,N),B(N,N),EVAL(N),EVEC(N,*)
C
      INTEGER           I,IM1,IP1,ISWEEP,J,JM1,JP1,K,LR,NM1,NNH,NP1,
     +                  NSDIG,NSWEEP
      INTEGER           IMP
      DOUBLE PRECISION  AI,ALIM,AMX,AUX,BI,BLIM,BMX,DD,DELTA,ONE,PRECS,
     +                  P1,P2,Q,RAD,SMALL,TEN,TIJ,TJI,TWO,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, TEN=10.0D0 )
C
      EXTERNAL          EVARR,IMP,PRACC,REARRV,RMINT,SCALE
C ----------------------------------------------------------------------
C   PRELIMINARIES  -  CHECK INPUT AND SET PARAMETERS
C ----------------------------------------------------------------------
      BLIM   = ZERO
      IF (N.LT.1.OR.N.GT.1000)  GO TO 910
C
      IERR   = 0
      NSDIG  = IMP(3)
      NSWEEP = NSDIG+1
      I      = NSDIG - NSDIG/5
      SMALL  = TEN**(-I)
      PRECS  = TEN**(-NSDIG)
C
      IF (TOL.LT.PRECS)  TOL = SMALL
C
      NM1 = N-1
      NP1 = N+1
      NNH = N*NM1/2
C
      IF (NVEC.EQ.N) THEN
C                                               ** INITIALIZE EVEC
         CALL RMINT (EVEC,N,N,ZERO)
         DO 20 I=1,N
            EVEC(I,I) = ONE
   20    CONTINUE
      ENDIF
C ----------------------------------------------------------------------
C   TRANSFORMATION OF A AND B (TO DIAGONAL FORM) AND OF EVEC (IF NVEC=N)
C   THE TRANSFORMATIONS (ROTATIONS) ARE CARRIED OUT IN SWEEPS
C   (NOTE THAT EVAL IS USED AS A SCRATCH ARRAY IN THIS SEQUENCE)
C ----------------------------------------------------------------------
C
      DO 400 ISWEEP=1,NSWEEP
C
         AUX   = PRACC(B,B,NP1,NP1,N)
         DELTA = PRECS*SQRT(AUX/N)
C                                               ** INITIATE ALIM AND
C                                                  NORMALIZATION
         ALIM  = ZERO
         DO 100 I=1,N
            IF (B(I,I).LE.DELTA)  GO TO 920
            A(I,I)  = A(I,I)/B(I,I)
            ALIM    = ALIM + ABS(A(I,I))
            EVAL(I) = ONE/SQRT(B(I,I))
            B(I,I)  = ONE
  100    CONTINUE
C
         IF (NVEC.EQ.N) THEN
C                                               ** MODIFY EVEC
            DO 120 J=1,N
               CALL SCALE (EVEC(1,J),EVEC(1,J),N,1,EVAL(J))
  120       CONTINUE
         ENDIF
C                                               ** EXIT IF N=1
         IF (N.EQ.1)  GO TO 1000
C                                               ** SCALE A AND B
         AMX = ZERO
         BMX = ZERO
         DO 150 I=1,NM1
            IP1 = I+1
            DO 140 J=IP1,N
               AUX    = EVAL(I)*EVAL(J)
               A(I,J) = A(I,J)*AUX
               B(I,J) = B(I,J)*AUX
               AMX    = MAX(ABS(A(I,J)),AMX)
               BMX    = MAX(ABS(B(I,J)),BMX)
  140       CONTINUE
  150    CONTINUE
C                                               ** SET THRESHOLDS
C                                                  ALIM AND BLIM
         ALIM = ALIM/N
         IF (ALIM.LT.TOL)  ALIM = ONE
         AMX  = AMX/ALIM
         BLIM = MAX(AMX,BMX)
C                                               ** EXIT LOOP IF
C                                                  CONVERGED
         IF (BLIM.LT.TOL)  GO TO 500
C
         BLIM = BLIM/TEN
         ALIM = ALIM*BLIM
C                                               ** PREPARE SWEEP
         LR = 0
         I  = 0
C ------------------------------------------------ PERFORM SWEEP
C                                                  I = ROW NUMBER
C                                                  J = COLUMN NUMBER
  200    I = I+1
         IF (I.GT.NM1)  I=1
         J = I
  250    J = J+1
         IF (J.GT.N)    GO TO 200
C
         IF (ABS(A(I,J)).GT.ALIM.OR.ABS(B(I,J)).GT.BLIM) THEN
            LR =-1
C                                               ** PREPARE ROTATION
            P1  = A(I,I)*B(I,J) - B(I,I)*A(I,J)
            P2  = A(J,J)*B(I,J) - B(J,J)*A(I,J)
            Q   = A(I,I)*B(J,J) - B(I,I)*A(J,J)
            Q   = Q/TWO
            RAD = Q*Q + P1*P2
            IF (RAD.LT.ZERO)  GO TO 920
            AUX = SQRT(RAD)
            DD  = Q + SIGN(AUX,Q)
            IF (ABS(DD).LT.PRECS) THEN
               TIJ = ZERO
               TJI =-B(I,J)
            ELSE
               TIJ = P2/DD
               TJI =-P1/DD
            ENDIF
C                                               ** ROTATE A AND B
            IP1 = I+1
            IM1 = I-1
            JP1 = J+1
            JM1 = J-1
            IF (I.GT.1) THEN
               DO 300 K=1,IM1
                  AI = A(K,I)
                  BI = B(K,I)
                  A(K,I) = AI + TJI*A(K,J)
                  B(K,I) = BI + TJI*B(K,J)
                  A(K,J) = A(K,J) + TIJ*AI
                  B(K,J) = B(K,J) + TIJ*BI
  300          CONTINUE
            ENDIF
            IF (J.LT.N) THEN
               DO 320 K=JP1,N
                  AI = A(I,K)
                  BI = B(I,K)
                  A(I,K) = AI + TJI*A(J,K)
                  B(I,K) = BI + TJI*B(J,K)
                  A(J,K) = A(J,K) + TIJ*AI
                  B(J,K) = B(J,K) + TIJ*BI
  320          CONTINUE
            ENDIF
            IF (JM1.GT.I) THEN
               DO 340 K=IP1,JM1
                  AI = A(I,K)
                  BI = B(I,K)
                  A(I,K) = AI + TJI*A(K,J)
                  B(I,K) = BI + TJI*B(K,J)
                  A(K,J) = A(K,J) + TIJ*AI
                  B(K,J) = B(K,J) + TIJ*BI
  340          CONTINUE
            ENDIF
            AI = A(I,I)
            BI = B(I,I)
            A(I,I) = AI + TWO*TJI*A(I,J) + TJI*TJI*A(J,J)
            B(I,I) = BI + TWO*TJI*B(I,J) + TJI*TJI*B(J,J)
            A(J,J) = A(J,J) + TWO*TIJ*A(I,J) + TIJ*TIJ*AI
            B(J,J) = B(J,J) + TWO*TIJ*B(I,J) + TIJ*TIJ*BI
            A(I,J) = ZERO
            B(I,J) = ZERO
C
            IF (NVEC.EQ.N) THEN
C                                               ** ROTATE EVEC
               DO 350 K=1,N
                  AUX = EVEC(K,I)
                  EVEC(K,I) = AUX + TJI*EVEC(K,J)
                  EVEC(K,J) = EVEC(K,J) + TIJ*AUX
  350          CONTINUE
            ENDIF
         ENDIF
C
         LR = LR+1
         IF (LR.LT.NNH)  GO TO 250
C
  400 CONTINUE
C
      GO TO 930
C ----------------------------------------------------------------------
C   TRANSFER EIGENVALUES AND UPDATE LOWER TRIANGLE OF  A  AND  B
C ----------------------------------------------------------------------
  500 DO 600 I=1,N
         EVAL(I) = A(I,I)
         DO 550 J=I,N
            A(J,I) = A(I,J)
            B(J,I) = B(I,J)
  550    CONTINUE
  600 CONTINUE
C                                               ** REARRANGE ELEMENTS OF
C                                                  EVAL AND EVEC
      IF (NVEC.EQ.N) THEN
         CALL EVARR (EVAL,EVEC,N,N,N,LORD)
      ELSE
         CALL REARRV (EVAL,N,LORD)
      ENDIF
C
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      GO TO 950
  930 IERR =-3
C
  950 IF (LPU.GT.0) THEN
         WRITE (LPU,6000)
         IF (IERR.EQ.-1)  WRITE (LPU,6100) N
         IF (IERR.EQ.-2)  WRITE (LPU,6200)
         IF (IERR.EQ.-3)  WRITE (LPU,6300) NSWEEP,TEN*BLIM
      ENDIF
C ------------------------------------------------
 1000 RETURN
C ------------------------------------------------ FORMATS
C
 6000 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE GENJAC')
 6100 FORMAT( 5X,'MATRIX DIMENSION ( N =',I5,' ) OUT OF RANGE')
 6200 FORMAT( 5X,'BREAKDOWN OF ALGORITHM  -  IS  B  POS.DEFINITE ?')
 6300 FORMAT( 5X,'SPECIFIED TOLERANCE NOT ACHIEVED IN',I3,' SWEEPS' /
     +        5X,'ACHIEVED TOLERANCE =',1PE11.3)
C
      END
      SUBROUTINE SKYJAC (A,B,MSKY,
     +                   AA,BB,EVAL,EVEC,
     +                   TOL,N,KSB,KODE,LORD,LPU,NEV,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SKYJAC                  GROUP 5 / PUBLIC
C
C     T A S K :  To compute all (NEV) finite eigenvalues (EVAL) and
C                corresponding eigenvectors (EVEC) of the general,
C                symmetric eigenproblem
C                         (A - lamda*B)q = 0
C                by a generalized Jacobi procedure.
C     Matrix A is stored in skyline format (upper half) and so is B
C     (if KSB=1) unless it is a diagonal matrix (KSB=2) stored as a
C     vector.  These input matrices are converted into full, N by N
C     matices, AA and BB, before solution.
C     The procedure requires matrix B to be positive definite.  If it
C     is, parameter KODE should be set equal to 2 (and the original
C     problem is solved directly).  If B is not positive definite, or
C     if A is more likely to be so, KODE should be set equal to 1,
C     causing the inverse problem to be solved and then "inverted".
C
C
C     ROUTINES CALLED/REFERENCED :  SCALE           (SAM-0)
C                                   GENJAC, EVARR   (SAM-5)
C                                   SKYCNV          (SAM-8)
C                                   SQRT, ABS       (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   89-01-25 / 1.0
C                       92-03-04 / 1.1    K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,KODE,KSB,LORD,LPU,N,NEV,    MSKY(*)
      DOUBLE PRECISION  TOL, A(*),AA(N,N),B(*),BB(N,N),EVAL(N),EVEC(N,N)
C
C                                                ! Local variables
      INTEGER           I,J
      DOUBLE PRECISION  EPS,ONE,ZERO
C
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, EPS = 1.0D-10 )
C ----------------------------------------------------------------------
      IERR = 0
C                                                ! Convert A and B
      CALL SKYCNV (A,MSKY,AA,N,1)
      CALL SKYCNV (B,MSKY,BB,N,KSB)
C                                                ! Determine NEV - no.
C                                                  of finite eigenvalues
      NEV = 0
      DO 10 I=1,N
         IF (BB(I,I) .LT. ZERO) THEN
            IF (KODE .EQ. 2) THEN
               IERR =-4
               IF (LPU .GT. 0) THEN
                  WRITE (LPU,690)
                  WRITE (LPU,691)
               ENDIF
               GO TO 100
            ENDIF
         ENDIF
         IF (BB(I,I) .NE. ZERO)  NEV = NEV+1
   10 CONTINUE
      IF (KSB .EQ. 1)  NEV = N
C
      IF (NEV .EQ. 0) THEN
         IERR =-5
         IF (LPU .GT. 0) THEN
            WRITE (LPU,690)
            WRITE (LPU,692)
         ENDIF
         GO TO 100
      ENDIF
C
      IF (N.EQ.1) THEN
C                                                ! Special case:  N=1
         EVAL(1)   = A(1)/B(1)
         EVEC(1,1) = ONE/SQRT(B(1))
         GO TO 100
      ENDIF
C
      IF (KODE .EQ. 2) THEN
C                                                ! Solve orign. problem
C
         CALL GENJAC (AA,BB,EVAL,EVEC,TOL,N,N,LORD,LPU,IERR)
         IF (IERR .LT. 0) THEN
            IF (LPU .GT. 0) THEN
               WRITE (LPU,690)
               WRITE (LPU,693)
            ENDIF
            GO TO 100
         ENDIF
C
      ELSE
C                                                ! Solve inverse problem
C
         CALL GENJAC (BB,AA,EVAL,EVEC,TOL,N,N,4,LPU,IERR)
         IF (IERR .LT. 0) THEN
            IF (LPU .GT. 0) THEN
               WRITE (LPU,690)
               WRITE (LPU,694)
            ENDIF
            GO TO 100
         ENDIF
C
         J = 0
         DO 15 I=1,N
            IF (ABS(EVAL(I)) .GT. EPS)  J=J+1
   15    CONTINUE
         IF (J .LT. NEV)  NEV=J
         IF (NEV .EQ. 0) THEN
            IERR =-5
            IF (LPU .GT. 0) THEN
               WRITE (LPU,690)
               WRITE (LPU,692)
            ENDIF
            GO TO 100
         ENDIF
C
C                                                  Invert eigenvalues
         DO 30 I=1,NEV
            EVAL(I) = ONE/EVAL(I)
   30    CONTINUE
C                                                  Rearrange order
         CALL EVARR (EVAL,EVEC,N,NEV,NEV,LORD)
C
      ENDIF
C
  100 RETURN
C
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE SKYJAC')
  691 FORMAT(5X,'MATRIX B IS NEGATIVE DEFINITE')
  692 FORMAT(5X,'PROBLEM HAS NO FINITE EIGENVALUES')
  693 FORMAT(5X,'DURING SOLUTION OF ORIGINAL PROBLEM')
  694 FORMAT(5X,'DURING SOLUTION OF INVERSE PROBLEM')
C
      END
      SUBROUTINE SPSJAC (A,EVAL,EVEC,RWA,EPS,N,NVEC,LORD,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPSJAC                  GROUP 5 / PUBLIC
C
C     TASK:  To determine all eigenvalues and all or none eigenvectors
C            of the special, symmetric eigenproblem
C
C                  A*EVEC = EVAL*EVEC
C
C            a special, cyclic Jacobi procedure.  A is a full, real
C            symmetric matrix.
C            Eigenvalues/eigenvectors are ordered according to argument
C            LORD, and the eigenvectors are, if computed, orthonormal.
C
C     This subroutine is based on FEMPACK routine JACOBI (by Hansteen
C     and Bell).
C
C     ROUTINES CALLED/REFERENCED :   EVARR       (SAM-5)
C                                    REARRV      (SAM-8)
C
C
C
C     PROGRAMMED BY :  Kolbein Bell
C     DATE/VERSION  :  93-05-18 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           IERR,LORD,LPU,N,NVEC
      DOUBLE PRECISION  EPS,A(N,N),EVAL(N),EVEC(N,N),RWA(N)
C
C                                                ! local variables
C
      INTEGER           I,IEND,J,K,M,MAX,MROT,NM1,NN,NROT,NV
      DOUBLE PRECISION  AI,AIJ,AJ,AKK,AUX,C,ERR,HALF,ONE
      DOUBLE PRECISION  S,TLD,TWO,U,V,Y,Z,ZERO
C
      PARAMETER         (ZERO=0.0D0, HALF=0.5D0,ONE=1.0D0, TWO=2.0D0)
C ----------------------------------------------------------------------
C
      IERR = 0
C
      IF (N .EQ. 1) THEN
         EVAL(1) = A(1,1)
         IF (NVEC .EQ. 1)  EVEC(1,1) = ONE
      ELSEIF (N .GT. 1) THEN
         MAX  = 20
         NM1  = N-1
C                                                ! copy diagonal of A
         DO 10 I=1,N
            RWA(I) = A(I,I)
   10    CONTINUE
C                                                ! initialize
         NN   = N*(N-1)/2
         MROT = MAX*NN
         NROT = 0
         IEND = 0
         NV   = 0
         IF (NVEC .EQ. N)  NV = N
C
         IF (NV .GT. 0) THEN
            DO 30 J=1,N
               DO 20 I=1,N
                  EVEC(I,J) = ZERO
   20          CONTINUE
            EVEC(J,J) = ONE
   30       CONTINUE
         ENDIF
C
C------- Determine numerically largest on- and off-diagonal elements and
C        error
C
  100    AIJ = ZERO
         AKK = ABS(A(1,1))
         DO 120 J=2,N
            DO 110 I=1,J-1
               IF (ABS(A(I,J)) .GT. AIJ)  AIJ = ABS(A(I,J))
  110       CONTINUE
            IF (ABS(A(J,J)) .GT. AKK)  AKK = ABS(A(J,J))
  120    CONTINUE
C                                                ! relative error
         ERR = NM1*AIJ / AKK
C                                                ! exit loop ?
         IF (ERR .LT. EPS)  GO TO 700
C
         IF (IEND .EQ. 1) THEN
            IERR = 1
            IF (LPU .GT. 0) THEN
               WRITE (LPU,690) ERR,EPS
            ENDIF
            EPS = ERR
C                                                ! exit with warning
            GO TO 700
         ENDIF
C
C ------ Prepare first sweep
C
         TLD = HALF*AIJ
         AIJ = ZERO
C
         M = 0
         I = 0
C                                                ! start new row
  200    I = I+1
         IF (I .GT. NM1)  I = 1
         J = I
C                                                ! start new column
  300    J = J+1
         IF (J .GT. N)  GO TO 200
C                                                ! check off-diag. elem.
         IF (ABS(A(I,J)) .GT. TLD) THEN
C
C --------- A large off-diagonal element is encountered;
C           perform rotation in A
C
            Y = A(I,I) - A(J,J)
            IF (Y .LT. ZERO) THEN
               Y = -Y
               Z = -TWO*A(I,J)
            ELSE
               Z = TWO*A(I,J)
            ENDIF
            U =  Z*Z + Y*Y
            V = SQRT(U)
C                                                ! cosine and sine
            C = SQRT(HALF*(ONE + Y/V))
            S = HALF*Z / (C*V)
C
            DO 410 K=1,I-1
               AI = A(K,I)
               AJ = A(K,J)
               A(K,I) = C*AI + S*AJ
               A(K,J) = S*AI - C*AJ
  410       CONTINUE
            DO 420 K=J+1,N
               AI = A(I,K)
               AJ = A(J,K)
               A(I,K) = C*AI + S*AJ
               A(J,K) = S*AI - C*AJ
  420       CONTINUE
            DO 430 K=I+1,J-1
               AI = A(I,K)
               AJ = A(K,J)
               A(I,K) = C*AI + S*AJ
               A(K,J) = S*AI - C*AJ
  430       CONTINUE
C
            AJ = A(J,J)
            A(J,J) = C*C*AJ + S*S*A(I,I) - TWO*C*S*A(I,J)
            A(I,I) = S*S*AJ + C*C*A(I,I) + TWO*C*S*A(I,J)
            A(I,J) = ZERO
C                                                ! eigenvectors ?
            IF (NV .EQ. N) THEN
C
C ------------ Rotate EVEC
C
               DO 510 K=1,N
                  AUX = EVEC(K,I)
                  EVEC(K,I) = AUX*C + EVEC(K,J)*S
                  EVEC(K,J) = AUX*S - EVEC(K,J)*C
  510          CONTINUE
            ENDIF
C                                                ! prepare new search
            M    = 0
            AIJ  = ZERO
            NROT = NROT+1
C
            IF (NROT .GT. MROT)  THEN
C                                                ! all rotations used
               IEND = 1
               GO TO 100
            ENDIF
         ELSE
C
C --------- Off-diagonal element is small  (less than TLD)
C
            M = M+1
            AUX = ABS(A(I,J))
            IF (AUX .GT. AIJ)  AIJ = AUX
         ENDIF
C                                                ! continue
         IF (M .LT. NN)  GO TO 300
C
C ------ Check convergence
C
         IF (AIJ .GT. ZERO) THEN
            AKK = ZERO
            DO 610 K=1,N
               AUX = ABS(A(K,K))
               IF (AUX .GT. AKK)  AKK = AUX
  610       CONTINUE
            IF (NM1*AIJ .GT. AKK*EPS) THEN
C                                                ! prepare new sweep
               TLD = HALF*AIJ
               AIJ = ZERO
               M   = 0
               GO TO 300
            ENDIF
         ENDIF
C
C ------ Transfer eigenvalues and restore A
C
  700    DO 720 J=1,N
            EVAL(J) = A(J,J)
            A(J,J)  = RWA(J)
            DO 710 I=J+1,N
               A(J,I) = A(I,J)
  710       CONTINUE
  720    CONTINUE
C
         IF (LORD.GT.0 .AND. LORD.LT.5) THEN
C
C --------- Rearrange eigenvalues and (if NV=N) eigenvectors
C
            IF (NV .EQ. N) THEN
               CALL EVARR (EVAL,EVEC,N,N,N,LORD)
            ELSE
               CALL REARRV (EVAL,N,LORD)
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
C
  690 FORMAT(///' *** WARNING from  S A M library routine SPSJAC' /
     &       5X,'Relative eigenvalue error greater than specified'/
     &       5X,'Actual error =',1PE10.2,5X,'Specified error =',
     &       1PE10.2)
C
      END
      SUBROUTINE HQLS (A,EVAL,EVEC,WA,N,NVEC,LORD,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  HQLS                    GROUP 5 / PUBLIC
C
C     T A S K :  TO COMPUTE ALL EIGENVALUES (EVAL) AND ALL OR NONE
C                EIGENVECTORS (STORED IN EVEC) OF A SYMMETRIC, REAL
C                MATRIX (A) OF DIMENSION N
C     THE SOLUTION IS OBTAINED THROUGH TWO STEPS :
C         1)  HOUSEHOLDER REDUCTION TO TRIDIAGONAL FORM
C         2)  IMPLICIT QL ITERATION
C
C     ROUTINES CALLED/REFERENCED :  PRACC                     (SAM-0)
C                                   HOUS3D,QLS3D,QLS3DV,EVARR (SAM-5)
C                                   REARRV                    (SAM-8)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-06-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,LORD,LPU,N,NVEC
      DOUBLE PRECISION  A(N,N),EVAL(N),EVEC(N,N),WA(*)
C
      INTEGER           I,IP1,IP2,J,K,KP1,L,NM1,NV
      DOUBLE PRECISION  S,    PRACC
C
      EXTERNAL          EVARR,HOUS3D,PRACC,QLS3D,QLS3DV,REARRV
C ----------------------------------------------------------------------
      IF (N.LT.1)       GO TO 100
      NV = 0
      IF (NVEC.EQ.N)    NV = N
C                                               ** POINTERS IN WA
      IP1 = 1
      IP2 = IP1+N
C
C --- REDUCTION TO TRIDIAGONAL FORM
C
      CALL HOUS3D (A,EVAL,WA(IP1),WA(IP2),N)
C
C --- COMPUTE AND ARRANGE IN SPECIFIED ORDER ALL EIGENVALUES AND IF
C     REQUESTED ALL EIGENVECTORS
C
      IF (NV.EQ.N) THEN
         CALL QLS3DV (EVAL,WA(IP1),EVEC,N,LPU,IERR)
         IF (IERR.LT.0)  NV = -IERR
         CALL EVARR (EVAL,EVEC,N,NV,NV,LORD)
      ELSE
         K = N
         CALL QLS3D (EVAL,WA(IP1),N,LPU,IERR)
         IF (IERR.LT.0)  K = -IERR
         CALL REARRV (EVAL,K,LORD)
      ENDIF
C
C --- TRANSFORM TO EIGENVECTORS OF  A, IF RELEVANT
C
      IF (NV.GT.0) THEN
         IF (N.GT.2) THEN
            NM1 = N-1
            DO 30 J=1,NV
               DO 20 I=2,NM1
                  K   = N-I
                  KP1 = K+1
                  S   = PRACC(A(KP1,K),EVEC(KP1,J),1,1,I)
                  DO 10 L=KP1,N
                     EVEC(L,J) = EVEC(L,J) - S*A(L,K)
   10             CONTINUE
   20          CONTINUE
   30       CONTINUE
         ENDIF
      ENDIF
C
C --- RESTORE MATRIX A
C
      DO 50 J=1,N
         DO 40 I=J,N
            A(I,J) = A(J,I)
   40    CONTINUE
   50 CONTINUE
C
      IF (IERR.LT.0) THEN
         IF (LPU.GT.0)  WRITE (LPU,600)
      ENDIF
C
  100 RETURN
C
  600 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE HQLS')
      END
      SUBROUTINE QLS3D (D,SD,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  QLS3D                   GROUP 5 / PUBLIC
C
C     T A S K :  TO DETERMINE ALL EIGENVALUES OF A SYMMETRIC, TRI-
C                DIAGONAL MATRIX STORED IN  D (DIAGONAL) AND  SD
C                (SUBDIAGONAL)
C     THE EIGENVALUES ARE DETERMINED BY THE IMPLICIT QL METHOD (AND RE-
C     TURNED IN D)
C     THIS SUBROUTINE IS BASED ON THE  E I S P A C K  ROUTINE  IMTQL1
C
C     ROUTINES CALLED/REFERENCED :  IMP         (SAM-0)
C                                   ABS, SQRT   (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL  (AFTER EISPACK ROUTINE IMTQL1)
C     DATE/VERSION  :   86-04-15 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,LPU,N
      DOUBLE PRECISION  D(N),SD(N)
C
      INTEGER           I,II,ITR,J,L,LMJ,NM1,NSDIG
      INTEGER           IMP
      DOUBLE PRECISION  B,C,F,G,ONE,P,PRECS,R,S,TEN,TWO,ZERO
C
      PARAMETER         (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, TEN=10.0D0)
C
      EXTERNAL          IMP
C
C ----------------------------------------------------------------------
C   PRELIMINARIES
C ----------------------------------------------------------------------
      NSDIG = IMP(3)
      PRECS = TEN**(-NSDIG)
      IERR  = 0
C
      IF (N.EQ.1)  GO TO 1000
C
      DO 40 I=2,N
         SD(I-1) = SD(I)
   40 CONTINUE
      SD(N) = ZERO
      NM1   = N-1
C ----------------------------------------------------------------------
C     MAIN LOOP
C ----------------------------------------------------------------------
      DO 300 J =1,N
         ITR = 0
C                                               ** LOOK FOR SMALL SUB-
C                                                  DIAGONAL ELEMENT
  100    DO 120 L=J,NM1
            S = ABS(D(L)) + ABS(D(L+1))
            IF (ABS(SD(L)).LE.PRECS*S)   GO TO 140
  120    CONTINUE
         L = N
C
  140    P = D(J)
         IF (L.EQ.J)     GO TO 300
         IF (ITR.EQ.30)  GO TO 900
         ITR = ITR+1
C                                               ** FORM SHIFT
         G = (D(J+1)-P)/(TWO*SD(J))
         R = SQRT(G*G + ONE)
         IF (G.GT.ZERO) THEN
            G = D(L) - P + SD(J)/(G+R)
         ELSE
            G = D(L) - P + SD(J)/(G-R)
         ENDIF
         S   = ONE
         C   = ONE
         P   = ZERO
         LMJ = L-J
C
C    --- FOR I=L-1 STEP -1 UNTIL J DO
         DO 200 II=1,LMJ
            I = L-II
            F = S*SD(I)
            B = C*SD(I)
            IF (ABS(F).GE.ABS(G)) THEN
               C = G/F
               R = SQRT(C*C + ONE)
               S = ONE/R
               C = C*S
               SD(I+1) = F*R
            ELSE
               S = F/G
               R = SQRT(S*S + ONE)
               C = ONE/R
               S = S*C
               SD(I+1) = G*R
            ENDIF
            G = D(I+1) - P
            R = (D(I)-G)*S + TWO*C*B
            P = S*R
            D(I+1) = G+P
            G = C*R - B
  200    CONTINUE
C    --- ENDDO
         D(J)  = D(J) - P
         SD(J) = G
         SD(L) = ZERO
         GO TO 100
C
  300 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  900 IERR =-J
      IF (LPU.GT.0)  WRITE (LPU,6000) J,N
C
 1000 RETURN
C ------------------------------------------------ FORMAT
 6000 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE QLS3D'
     +         /5X,  'ONLY',I5,'  OF',I5,'  EIGENVALUES HAVE CONVERGED')
C
      END
      SUBROUTINE QLS3DV (D,SD,EVEC,N,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  QLS3DV                  GROUP 5 / PUBLIC
C
C     T A S K :  TO DETERMINE ALL EIGENVALUES AND CORRESPONDING EIGEN-
C                VECTORS OF A SYMMETRIC, TRIDIAGONAL MATRIX STORED IN
C                D (DIAGONAL) AND SD (SUBDIAGONAL)
C     THE EIGENVALUES ARE DETERMINED BY THE IMPLICIT QL METHOD (AND RE-
C     TURNED IN D), AND THE ACCUMULATED QL TRANSFORMATIONS ARE USED TO
C     COMPUTE THE EIGENVECTORS.
C     THIS SUBROUTINE IS BASED ON THE  E I S P A C K  ROUTINE  IMTQL2
C
C     ROUTINES CALLED/REFERENCED :  IMP         (SAM-0)
C                                   ABS, SQRT   (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL  (AFTER EISPACK ROUTINE IMTQL2)
C     DATE/VERSION  :   86-04-09 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,LPU,N
      DOUBLE PRECISION  D(N),EVEC(N,N),SD(N)
C
      INTEGER           I,II,ITR,J,K,L,LMJ,NM1,NSDIG
      INTEGER           IMP
      DOUBLE PRECISION  B,C,F,G,ONE,P,PRECS,R,S,TEN,TWO,ZERO
C
      PARAMETER         (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, TEN=10.0D0)
C
      EXTERNAL          IMP
C
C ----------------------------------------------------------------------
C   PRELIMINARIES
C ----------------------------------------------------------------------
      NSDIG = IMP(3)
      PRECS = TEN**(-NSDIG)
      IERR  = 0
C
      DO 20 J=1,N
         DO 10 I=1,N
            EVEC(I,J) = ZERO
   10    CONTINUE
         EVEC(J,J) = ONE
   20 CONTINUE
C
      IF (N.EQ.1)  GO TO 1000
C
      DO 40 I=2,N
         SD(I-1) = SD(I)
   40 CONTINUE
      SD(N) = ZERO
      NM1   = N-1
C ----------------------------------------------------------------------
C     MAIN LOOP
C ----------------------------------------------------------------------
      DO 300 J =1,N
         ITR = 0
C                                               ** LOOK FOR SMALL SUB-
C                                                  DIAGONAL ELEMENT
  100    DO 120 L=J,NM1
            S = ABS(D(L)) + ABS(D(L+1))
            IF (ABS(SD(L)).LE.PRECS*S)   GO TO 140
  120    CONTINUE
         L = N
C
  140    P = D(J)
         IF (L.EQ.J)     GO TO 300
         IF (ITR.EQ.30)  GO TO 900
         ITR = ITR+1
C                                               ** FORM SHIFT
         G = (D(J+1)-P)/(TWO*SD(J))
         R = SQRT(G*G + ONE)
         IF (G.GT.ZERO) THEN
            G = D(L) - P + SD(J)/(G+R)
         ELSE
            G = D(L) - P + SD(J)/(G-R)
         ENDIF
         S   = ONE
         C   = ONE
         P   = ZERO
         LMJ = L-J
C
C    --- FOR I=L-1 STEP -1 UNTIL J DO
         DO 200 II=1,LMJ
            I = L-II
            F = S*SD(I)
            B = C*SD(I)
            IF (ABS(F).GE.ABS(G)) THEN
               C = G/F
               R = SQRT(C*C + ONE)
               S = ONE/R
               C = C*S
               SD(I+1) = F*R
            ELSE
               S = F/G
               R = SQRT(S*S + ONE)
               C = ONE/R
               S = S*C
               SD(I+1) = G*R
            ENDIF
            G = D(I+1) - P
            R = (D(I)-G)*S + TWO*C*B
            P = S*R
            D(I+1) = G+P
            G = C*R - B
C                                               ** FORM VECTOR
            DO 175 K=1,N
               F = EVEC(K,I+1)
               EVEC(K,I+1) = S*EVEC(K,I) + C*F
               EVEC(K,I)   = C*EVEC(K,I) - S*F
  175       CONTINUE
  200    CONTINUE
C    --- ENDDO
         D(J)  = D(J) - P
         SD(J) = G
         SD(L) = ZERO
         GO TO 100
C
  300 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  900 IERR =-J
      IF (LPU.GT.0)  WRITE (LPU,6000) J,N
C
 1000 RETURN
C ------------------------------------------------ FORMAT
 6000 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE QLS3DV'
     +         /5X,  'ONLY',I5,'  OF',I5,'  EIGENVALUES HAVE CONVERGED')
C
      END
      SUBROUTINE HQRIS (A,EVAL,EVEC,RLWA,ILWA,N,NVEC,LORD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  HQRIS                   GROUP 5 / PUBLIC
C
C     T A S K :  TO COMPUTE ALL EIGENVALUES (EVAL) AND AS MANY EIGEN-
C                VECTORS (STORED IN EVEC) AS DESIRED OF A SYMMETRIC,
C                REAL MATRIX (A) OF DIMENSION N
C     THE SOLUTION IS OBTAINED THROUGH THREE STEPS :
C         1)  HOUSEHOLDER REDUCTION TO TRIDIAGONAL FORM
C         2)  KAHAN-VARAH  QR ITERATION
C         3)  INVERSE ITERATION
C
C     ROUTINES CALLED/REFERENCED :  PRACC                     (SAM-0)
C                                   HOUS3D,QRS3D AND IIS3D    (SAM-5)
C                                   REARRV                    (SAM-8)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-06-15 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           LORD,N,NVEC,     ILWA(*)
      DOUBLE PRECISION  A(N,N),EVAL(N),EVEC(N,*),RLWA(*)
C
      INTEGER           I,IP1,IP2,IP3,IP4,IP5,J,K,KP1,L,NM1
      DOUBLE PRECISION  S,    PRACC
C
      EXTERNAL          HOUS3D,IIS3D,PRACC,QRS3D,REARRV
C ----------------------------------------------------------------------
      IF (N.LT.1)       GO TO 100
C                                               ** POINTERS IN RLWA
      IP1 = 1
      IP2 = IP1+N
      IP3 = IP2+N
      IP4 = IP3+N
      IP5 = IP4+N
C
C --- REDUCTION TO TRIDIAGONAL FORM
C
      CALL HOUS3D (A,RLWA(IP1),RLWA(IP2),RLWA(IP3),N)
C
C --- COMPUTE AND ARRANGE IN SPECIFIED ORDER ALL EIGENVALUES
C
      CALL QRS3D  (RLWA(IP1),RLWA(IP2),RLWA(IP3),EVAL,N)
      CALL REARRV (EVAL,N,LORD)
C
C --- COMPUTE EIGENVECTORS, IF REQUESTED
C
      IF (NVEC.GT.0) THEN
         CALL IIS3D (RLWA(IP1),RLWA(IP2),EVAL,EVEC,
     +               RLWA(IP3),RLWA(IP4),RLWA(IP5),ILWA,N,NVEC)
         IF (N.GT.2) THEN
C                                               ** TRANSFORM TO
C                                                  EIGENVECTORS OF  A
            NM1 = N-1
            DO 30 J=1,NVEC
               DO 20 I=2,NM1
                  K   = N-I
                  KP1 = K+1
                  S   = PRACC(A(KP1,K),EVEC(KP1,J),1,1,I)
                  DO 10 L=KP1,N
                     EVEC(L,J) = EVEC(L,J) - S*A(L,K)
   10             CONTINUE
   20          CONTINUE
   30       CONTINUE
         ENDIF
      ENDIF
C
C --- RESTORE MATRIX A
C
      DO 50 J=1,N
         DO 40 I=J,N
            A(I,J) = A(J,I)
   40    CONTINUE
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE QRS3D (D,SD,V,EVAL,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE : QRS3D                    GROUP 5 / PUBLIC
C
C     T A S K :  TO COMPUTE ALL EIGENVALUES, EVAL, OF A SYMMETRIC, TRI-
C                DIAGONAL MATRIX OF DIMENSION N, WITH  D  CONTAINING
C                THE DIAGONAL AND THE LAST  N-1  ELEMENTS OF  SD  CON-
C                TAINING THE SUPER-/SUB-DIAGONAL
C     A QR-METHOD DUE TO KAHAN AND VARAH IS USED
C
C     ROUTINES CALLED/REFERENCED :  IMP              (SAM-0)
C                                   MINDX            (SAM-8)
C                                   ABS, SIGN, SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-04-01 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           N
      DOUBLE PRECISION  D(N),EVAL(N),SD(N),V(*)
C
      INTEGER           I,K,KP1,L,LM1,MAXEXP,NSDIG
      INTEGER           IMP,MINDX
      DOUBLE PRECISION  BMAX,C,DEL1,DEL2,E1,E2,F,G,ONE,PRECS,S,SCALE
      DOUBLE PRECISION  SHIFT,SQOV,T,TEN,TINY,TMAX,TOL,TWO,X,Y,ZERO
C
      PARAMETER         (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, TEN=10.0D0)
C
      EXTERNAL          IMP,MINDX
C
C ----------------------------------------------------------------------
C   PRELIMINARIES  -  SCALING AND TOLERANCE LIMITS
C ----------------------------------------------------------------------
      NSDIG  = IMP(3)
      MAXEXP = IMP(5)
      PRECS  = TEN**(-NSDIG)
      SQOV   = TEN**(MAXEXP/3)
      TINY   = TEN**(-MAXEXP+3)
C
      I      = MINDX(SD,1,N,4)
      BMAX   = ABS(SD(I))
      I      = MINDX(D,1,N,4)
      TMAX   = ABS(D(I))
      IF (BMAX.GT.TMAX)           TMAX = BMAX
C                                               ** NULL MATRIX ?
      IF (TMAX.LT.TINY)           GO TO 400
C
      SCALE  = ONE/TMAX
C                                               ** DIAGONAL MATRIX ?
      IF (BMAX*SCALE.LT.PRECS)    GO TO 400
C                                               ** SCALING FACTOR
      SCALE  = TEN
   20 SCALE  = SCALE/TEN
      IF (SCALE*TMAX.GT.SQOV)     GO TO 20
      DO 40 I=1,MAXEXP
         IF (SCALE*TMAX.GT.SQOV)  GO TO 60
         SCALE = SCALE*TEN
   40 CONTINUE
   60 SCALE = SCALE/TEN
C                                               ** SCALE TRIDIAG. FORM
      DO 80 I=1,N
         EVAL(I) = SCALE*D(I)
         V(I)    = (SCALE*SD(I))**2
   80 CONTINUE
      V(N+1) = ZERO
C                                               ** TOLERANCE LIMITS
      X    = N
      TOL  = PRECS/(TEN*X)
      DEL1 = TMAX*SCALE*TOL
      DEL2 = DEL1*DEL1
C ----------------------------------------------------------------------
C   COMPUTE THE EIGENVALUES
C ----------------------------------------------------------------------
      K = N
  100 L = K
      IF (L.LT.1) GO TO 300
         LM1 = L-1
C                                               ** NEGLIGIBLE OFF-
C                                                  DIAGONAL TERMS ?
         DO 120 I=1,L
            KP1 = K
            K   = K-1
            IF (V(KP1).LT.DEL2)   GO TO 140
  120    CONTINUE
  140    IF (KP1.LT.L)            GO TO 160
         V(L) = ZERO
      GO TO 100
C                                               ** COMPUTE SHIFT
  160    T = EVAL(L)-EVAL(LM1)
         X = V(L)
         Y = T/TWO
         IF (ABS(T).GT.DEL1) THEN
            S = (X/Y)/(ONE + SQRT(ONE + X/Y**2))
         ELSE
            S = SQRT(X)
         ENDIF
         E1 = EVAL(L)   + S
         E2 = EVAL(LM1) - S
C
         IF (KP1.LT.LM1)          GO TO 180
C
         EVAL(L)   = E1
         EVAL(LM1) = E2
         V(LM1)    = ZERO
      GO TO 100
C
  180    SHIFT = E1
         IF (ABS(T).GT.DEL1)      GO TO 200
         IF (ABS(E2).GT.ABS(E1))  GO TO 200
         SHIFT = E2
C                                               ** QR-TRANSFORMATION
  200    C = ONE
         S = ZERO
         G = EVAL(KP1) - SHIFT
  220    IF (ABS(G).LT.DEL1)      G = G + SIGN(C*DEL1,G)
         K    = KP1
         KP1  = K+1
         F    = (G**2)/C
         X    = V(KP1)
         T    = X+F
         V(K) = S*T
C
         IF (K.EQ.L)  GO TO 240
C
         C    = F/T
         S    = X/T
         X    = G
         G    = C*(EVAL(KP1)-SHIFT) - S*X
         EVAL(K) = X - G + EVAL(KP1)
         GO TO 220
C
  240    EVAL(K) = G + SHIFT
      GO TO 100
C                                               ** SCALE EIGENVALUES
  300 SCALE = ONE/SCALE
      DO 320 I=1,N
         EVAL(I) = SCALE*EVAL(I)
  320 CONTINUE
      GO TO 500
C ----------------------------------------------------------------------
C   DIAGONAL MATRIX
C ----------------------------------------------------------------------
  400 DO 420 I=1,N
         EVAL(I) = D(I)
  420 CONTINUE
C
C
  500 RETURN
      END
      SUBROUTINE IIS3D (D,SD,EVAL,EVEC,P,Q,U,LARR,N,NVEC)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  IIS3D                   GROUP 5 / PUBLIC
C
C     T A S K :  TO COMPUTE, BY INVERSE ITERATION, THE  NVEC  EIGEN-
C                VECTORS CORRESPONDING TO THE FIRST NVEC EIGENVALUES
C                EVAL(1),EVAL(2),.. OF A SYMMETRIC TRIDIAGONAL FORM
C                (D AND SD).
C     TWO STEPS OF INVERSE ITERATION WITH SHIFT ARE USED FOR EACH VECTOR
C
C     ROUTINES CALLED/REFERENCED :  IMP                (SAM-0)
C                                   RANVEC AND MINDX   (SAM-8)
C                                   DCOPY AND DSCAL    (BLAS)
C                                   ABS AND SQRT       (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL   (AFTER C.A.FELIPPA)
C     DATE/VERSION  :   86-06-15 / 1.0
C                       19-05-22 / 1.1  K.M.Okstad  (using BLAS)
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           N,NVEC,    LARR(N)
      DOUBLE PRECISION  D(N),EVAL(NVEC),EVEC(N,NVEC),P(N),Q(N),SD(N)
      DOUBLE PRECISION  U(N)
C
      INTEGER           I,ITR,J,L,NM1,NSDIG
      INTEGER           IMP,MINDX
C
      DOUBLE PRECISION  BMAX,C,DMAX,EV,PRECS,S,SCALE,SEP,SUM,TOL,X,Y,Z
      DOUBLE PRECISION  RAN,ONE,TEN,ZERO
C
      PARAMETER         ( ZERO = 0.0D0, ONE = 1.0D0, TEN = 10.0D0 )
C
      EXTERNAL          IMP,MINDX,RANVEC
C
C ----------------------------------------------------------------------
C   PRELIMINARIES  -  SCALING AND TOLERANCE LIMITS
C ----------------------------------------------------------------------
      IF (N.EQ.1)       GO TO 700
C
      RAN   = 0.1234567891234567D0
      NSDIG = IMP(3)
      PRECS = TEN**(-NSDIG)
      I     = MINDX(D,1,N,4)
      DMAX  = ABS(D(I))
      I     = MINDX(SD,1,N,4)
      BMAX  = ABS(SD(I))
C
      IF (DMAX.EQ.ZERO) GO TO 700
      SCALE = ONE/DMAX
C                                               ** DIAGONAL MATRIX ?
      IF (BMAX*SCALE.LT.PRECS)  GO TO 700
      X     = N
      TOL   = PRECS/(TEN*X)
      SEP   = TEN*(DMAX+BMAX)*PRECS
C                                               ** SCALE TRIDIAG. FORM
      CALL DSCAL (N,SCALE,D,1)
      CALL DSCAL (N,SCALE,SD,1)
C
C ----------------------------------------------------------------------
C   COMPUTE NVEC EIGENVECTORS OF TRIDIAGONAL FORM
C ----------------------------------------------------------------------
      EV = 0.0D0
      NM1 = N-1
      DO 600 J=1,NVEC
C                                               ** INITIATE (REDUCED)
C                                                  RIGHT-HAND SIDE
         IF (J.EQ.1 .OR. ABS(EVAL(J)-EV).GT.SEP) THEN
            CALL DCOPY (N,ONE,0,U,1)
         ELSE
            CALL RANVEC (U,RAN,N)
         ENDIF
C
         EV = SCALE*EVAL(J)
C                                               ** GAUSSIAN ELIMINATION
C                                                  WITH PART. PIVOTING
         X = D(1)-EV
         Y = SD(2)
         DO 200 I=1,NM1
            C = D(I+1)-EV
            S = SD(I+1)
            IF (ABS(X).GE.ABS(S))  GO TO 160
C                                               ** INTERCHANGE ROWS
            IF (ABS(S).LT.TOL)     S = TOL
            LARR(I) = 1
            P(I)    = S
            Q(I)    = C
            Z       =-X/S
            X       = Y + Z*C
            IF (I.LT.NM1)          Y = Z*SD(I+2)
            GO TO 180
C
  160       IF (ABS(X).LT.TOL)     X = TOL
            LARR(I) = 0
            P(I)    = X
            Q(I)    = Y
            Z       =-S/X
            X       = C + Z*Y
            IF (I.LT.NM1)          Y = SD(I+2)
C
  180       EVEC(I,J) = Z
  200    CONTINUE
         IF (ABS(X).LT.TOL)        X = TOL
C
         ITR = 0
  300    ITR = ITR+1
         IF (ITR.EQ.1)  GO TO 400
C                                               ** REDUCE (IMPROVE)
C                                                  RIGHT-HAND SIDE
         DO 350 I=1,NM1
            Z = EVEC(I,J)
            IF (LARR(I).EQ.1) THEN
               Y      = U(I)
               U(I)   = U(I+1)
               U(I+1) = Y + Z*U(I)
            ELSE
               U(I+1) = U(I+1) + Z*U(I)
            ENDIF
  350    CONTINUE
C                                               ** BACKSUBSTITUTION
  400    U(N) = U(N)/X
         DO 450 L=1,NM1
            I = N-L
            Y = U(I) - Q(I)*U(I+1)
            IF (LARR(I).EQ.1) THEN
               IF (I.LT.NM1)  Y = Y - SD(I+2)*U(I+2)
            ENDIF
            U(I) = Y/P(I)
  450    CONTINUE
C                                               ** NORMALIZE SOLUTION
         I   = MINDX(U,1,N,4)
         CALL DSCAL (N,ONE/U(I),U,1)
         SUM = ZERO
         DO 480 I=1,N
            SUM  = SUM + U(I)*U(I)
  480    CONTINUE
         CALL DSCAL (N,ONE/SQRT(SUM),U,1)
C
         IF (ITR.LT.2)  GO TO 300
C
         CALL DCOPY (N,U(1),1,EVEC(1,J),1)
         EV = EVAL(J)
  600 CONTINUE
      GO TO 1000
C ----------------------------------------------------------------------
C   DIAGONAL MATRIX
C ----------------------------------------------------------------------
  700 DO 800 J=1,NVEC
         CALL DCOPY (N,ZERO,0,EVEC(1,J),1)
         EVEC(J,J) = ONE
  800 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE HOUS3D (A,D,SD,V,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  HOUS3D                  GROUP 5/ PUBLIC
C
C     T A S K :  TO REDUCE A SYMMETRIC MATRIX  A  TO TRIDIAGONAL FORM
C                BY THE METHOD OF HOUSEHOLDER
C     THE DIAGONAL OF THE TRIDIAG. FORM IS RETURNED IN  D  AND THE
C     SUPER-/SUB-DIAGONAL IN THE LAST N-1 ELEMENTS OF  SD (SD(1)=0.).
C     THE UPPER TRIANGULAR PART OF  A  RETURNS UNALTERED (INCL. THE
C     DIAGONAL), WHEREAS THE SUB-DIAGONAL PART CONTAINS, IN EACH COLUMN,
C     THE NON-ZERO PART OF THE VECTORS (U-K) DEFINING THE REFLECTION
C     MATRICES  (P-K = I - (U-K)*(U-K)T)
C
C     ROUTINES CALLED/REFERENCED :  PRACC            (SAM-0)
C                                   ABS, SIGN, SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-03-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           N
      DOUBLE PRECISION  A(N,N),D(N),SD(N),V(N)
C
      INTEGER           I,J,K,KP1,KP2,NM1,NM2,NN
      DOUBLE PRECISION  AK1,ONE,TWO,X1,X2,X3,ZERO,       PRACC
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 , TWO = 2.0D0 )
C
      EXTERNAL          PRACC
C
C ----------------------------------------------------------------------
      SD(1) = ZERO
C                                               ** SAVE DIAGONAL OF  A
      DO 20 I=1,N
         D(I) = A(I,I)
   20 CONTINUE
C
      IF (N.EQ.1)  GO TO 1000
      IF (N.EQ.2)  GO TO  300
C ----------------------------------------------------------------------
C   REDUCTION OF  A  THROUGH N-2 ORTHOGONAL SIMILARITY TRANSFORMATIONS
C ----------------------------------------------------------------------
      NM1 = N-1
      NM2 = N-2
      DO 200 K=1,NM2
         KP1 = K+1
         KP2 = K+2
         AK1 = A(KP1,K)
         NN  = N-KP1
         X1  = PRACC(A(KP2,K),A(KP2,K),1,1,NN)
         IF (X1.GT.ZERO) THEN
            X2 = SQRT(X1 + AK1*AK1)
            SD(KP1) =-SIGN(X2,AK1)
C                                               ** VECTOR  U-K (STORED
C                                                  IN COLUMN K OF  A)
            X1 = SQRT(ONE + ABS(AK1)/X2)
            A(KP1,K) = X1
            X3 = SIGN(ONE/(X2*X1),AK1)
            DO 110 I=KP2,N
               A(I,K) = X3*A(I,K)
  110       CONTINUE
C                                               ** VECTOR  P-K
C                                                  (STORED IN  V)
            DO 120 I=KP1,NM1
               NN = I-K
               V(I) = PRACC(A(I,KP1),A(KP1,K),N,1,NN)
               NN = N-I
               J  = I+1
               V(I) = V(I) + PRACC(A(J,I),A(J,K),1,1,NN)
  120       CONTINUE
            NN = N-K
            V(N) = PRACC(A(N,KP1),A(KP1,K),N,1,NN)
C                                               ** VECTOR Q-K
C                                                  (STORED IN  V)
            X3 = PRACC(A(KP1,K),V(KP1),1,1,NN)/TWO
            DO 140 I=KP1,N
               V(I) = X3*A(I,K) - V(I)
  140       CONTINUE
C                                               ** THE MODIFIED MATRIX
C                                                  A-K (LOWER TRIANGLE)
            DO 160 J=KP1,N
               DO 150 I=J,N
                  A(I,J) = A(I,J) + V(I)*A(J,K) + V(J)*A(I,K)
  150          CONTINUE
  160       CONTINUE
C                                               ** ELEMENTS OF  U-K
C                                                  ALREADY ZERO
         ELSE
            A(KP1,K) = SQRT(TWO)
            SD(KP1)  =-AK1
            DO 180 I=KP2,N
               A(I,KP1) =-A(I,KP1)
  180       CONTINUE
         ENDIF
  200 CONTINUE
C ------------------------------------------------ DIAGONALS
  300 DO 400 I=1,N
         X1     = D(I)
         D(I)   = A(I,I)
         A(I,I) = X1
  400 CONTINUE
C                                               ** LAST ELEMENT OF SD
      SD(N) = A(N,N-1)
C
C
 1000 RETURN
      END
      SUBROUTINE ABRLER (V1,V2,E,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ABRLER                  GROUP 5 / PRIVATE
C
C     T A S K :  TO DETERMINE THE ABSOLUTE RELATIVE DIFFERENCE (ERROR)
C                BETWEEN CORRESPONDING ELEMENTS IN VECTORS V1 AND V2
C
C     ROUTINES CALLED/REFERENCED :  ABS  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-09-27 / 1.0
C                       99-10-14 / 1.1     K.Bell
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           N
      DOUBLE PRECISION  V1(N),V2(N),E(N)
C
      INTEGER           I
      DOUBLE PRECISION  ONE,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C ----------------------------------------------------------------------
      DO 10 I=1,N
         IF (V2(I).EQ.ZERO) THEN
            IF (V1(I).NE.ZERO) THEN
               E(I) = ONE
            ELSE
               E(I) = ZERO
            ENDIF
         ELSE
            E(I) = ABS((ABS(V1(I)) - ABS(V2(I)))/V2(I))
         ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BCKSUB (U,X,MSKY,N,KSU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  BCKSUB                  GROUP 5 / PRIVATE
C
C     T A S K :  TO PERFORME THE BACKSUBSTITUTION OPERATION
C                   X := (U-1)*X
C                WHERE U IS
C                - A DIAGONAL MATRIX    (KSU =-2 OR 2)   OR
C                - AN UPPER TRIANGULAR SKYLINE MATRIX (ALL OTHER VALUES
C                  OF KSU)
C
C     ROUTINES CALLED/REFERENCED :  ABS   (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-09-27 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KSU,N,     MSKY(N)
      DOUBLE PRECISION  U(*),X(N)
C
      INTEGER           I,IE,II,IP,IS,K,LC
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
C ----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1)       GO TO 100
      LC   = ABS(KSU)
C
      I = N+1
   10 I = I-1
C     DO WHILE I.GE.1
         IP = MSKY(I)
         IF (LC.EQ.2)        IP = I
         IF (U(IP).LE.ZERO)  GO TO 90
         X(I) = X(I)/U(IP)
C                                               ** LOOP EXIT
         IF (I.EQ.1)         GO TO 100
C
         IF (LC.EQ.2)        GO TO 10
C
         K  = I-MSKY(I)+MSKY(I-1)+1
         IS = MSKY(I-1)+1
         IE = MSKY(I)-1
         IF (IE.GE.IS) THEN
            DO 20 II=IS,IE
               X(K) = X(K)-U(II)*X(I)
               K    = K+1
   20       CONTINUE
         ENDIF
         GO TO 10
C     ENDDO
C
   90 IERR = -1
C
  100 RETURN
      END
      SUBROUTINE BLBLBT (B,MSKY,EPS,N,KSB,LPU,NZD,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  BLBLBT                  GROUP 5 / PRIVATE
C
C     T A S K :  TO FACTORIZE A MATRIX (B), WHICH MAY BE A SYMMETRIC
C                SKYLINE MATRIX (KSB=1) OR A DIAGONAL MATRIX (KSB=2),
C                INTO ITS CHOLESKY FACTOR (B = LB*LBT)
C                THE NUMBER OF ZERO DIAGONAL ELEMENTS IN B IS RETURNED
C                IN ARGUMENT NZD
C
C     ROUTINES CALLED/REFERENCED :  SKYSOL    (SAM-4)
C                                   SQRT      (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-07-13 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KSB,LPU,N,NZD,     MSKY(N)
      DOUBLE PRECISION  B(*),EPS
C
      INTEGER           I,IP,J,JE,JP,JS,NN
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          SKYSOL
C
C ----------------------------------------------------------------------
      IERR = 0
      NZD  = 0
      IF (KSB.EQ.1) THEN
C                                               ** SKYLINE MATRIX
C
         CALL SKYSOL (B,B,MSKY,EPS,N,N,N,LPU,1,NN,IERR)
         IF (IERR.LT.0)  GO TO 100
         IF (NN.GT.0)    GO TO  90
         B(1) = SQRT(B(1))
         IF (N.GT.1) THEN
            DO 50 I=2,N
               IP    = MSKY(I)
               B(IP) = SQRT(B(IP))
               NN    = IP-MSKY(I-1)-1
               IF (NN.GT.0) THEN
                  JS = I-NN
                  JE = I-1
                  IP = IP-NN
                  DO 25 J=JS,JE
                     JP    = MSKY(J)
                     B(IP) = B(IP)*B(JP)
                     IP    = IP+1
   25             CONTINUE
               ENDIF
   50       CONTINUE
         ENDIF
C
      ELSEIF (KSB.EQ.2) THEN
C                                               ** DIAGONAL MATRIX
         DO 75 I=1,N
            IF (B(I).LT.ZERO)  GO TO 90
            IF (B(I).EQ.ZERO)  NZD = NZD+1
            IF (B(I).GT.ZERO)  B(I) = SQRT(B(I))
   75    CONTINUE
C
      ELSE
C                                               ** UNKNOWN MATRIX
         IERR =-6
      ENDIF
      GO TO 100
C
C
   90 IERR =-5
C
C
  100 RETURN
C
      END
      SUBROUTINE SPSMUL (A,U,X,Y,W,MSKY,EPS,N,KSA,KSU,IFLAG,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPSMUL                  GROUP 5 / PRIVATE
C
C     T A S K :  TO DETERMINE
C
C                   X := (A-1)*UT*X       IFLAG = 1
C                   Y  = (A-1)*UT*X       IFLAG = 2
C                   X := U*(A-1)*UT*X     IFLAG = 3
C                   Y  = U*(A-1)*UT*X     IFLAG = 4
C
C     WHERE X AND Y ARE VECTORS AND  A  IS A SYMMETRIC SKYLINE MATRIX -
C     DEPENDING ON THE STORAGE CODE KSA, A (KSA=1) OR ITS FACTORS LT
C     AND D (KSA=-1) ARE STORED IN ARRAY A (A = L*D*LT) - U IS AN UPPER
C     TRIANGULAR SKYLINE MATRIX (ABS(KSU)=1) OR A DIAGONAL MATRIX
C     (ABS(KSB)=2)
C
C     ROUTINES CALLED/REFERENCED :  TSKYMB  (SAM-3)
C                                   SKYSOL  (SAM-4)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-09-27 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,LPU,KSA,KSU,N,     MSKY(N)
      DOUBLE PRECISION  EPS,A(*),U(*),Y(N),X(N),W(N)
C
      INTEGER           LCA,LCU,NN
C
      EXTERNAL          SKYSOL,TSKYMB
C ----------------------------------------------------------------------
      IF (IFLAG.LT.1.OR.IFLAG.GT.4)  GO TO 100
C
      LCA = 3
      IF (KSA.LT.0)       LCA = 4
      LCU = 1
      IF (ABS(KSU).EQ.2)  LCU =-1
C
      IF (IFLAG.EQ.1) THEN
         CALL TSKYMB (U,X,X,Y,MSKY,N,1,LCU,2,LPU,IERR)
         IF (IERR.LT.0)   GO TO 100
         CALL SKYSOL (A,X,MSKY,EPS,N,N,1,LPU,LCA,NN,IERR)
         IF (IERR.LT.0)   GO TO 100
      ELSEIF (IFLAG.EQ.4) THEN
         CALL TSKYMB (U,X,Y,W,MSKY,N,1,LCU,22,LPU,IERR)
         IF (IERR.LT.0)   GO TO 100
         CALL SKYSOL (A,W,MSKY,EPS,N,N,1,LPU,LCA,NN,IERR)
         IF (IERR.LT.0)   GO TO 100
      ELSE
         CALL TSKYMB (U,X,X,Y,MSKY,N,1,LCU,22,LPU,IERR)
         IF (IERR.LT.0)   GO TO 100
         CALL SKYSOL (A,Y,MSKY,EPS,N,N,1,LPU,LCA,NN,IERR)
         IF (IERR.LT.0)   GO TO 100
      ENDIF
C
      IF (IFLAG.EQ.3)  CALL TSKYMB (U,Y,Y,X,MSKY,N,1,LCU,11,LPU,IERR)
      IF (IFLAG.EQ.4)  CALL TSKYMB (U,W,X,Y,MSKY,N,1,LCU,11,LPU,IERR)
C
      KSA =-1
C
  100 RETURN
      END
      SUBROUTINE TOLCHK (EVERR,MANPRV,MANCUR,TOL,N,NGDEV,NACEV)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  TOLCHK                  GROUP 5 / PRIVATE
C
C     T A S K :  TO UPDATE THE CURRENT MATRIX OF ACCEPTANCE NUMBERS
C                (MANCUR), COMPARE WITH PREVIOUS NUMBERS AND DETERMINE
C                NGDEV (NO. OF "GOOD" EIGENVALUES)
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-10-04 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           N,NACEV,NGDEV,   MANCUR(N),MANPRV(N)
      DOUBLE PRECISION  TOL,       EVERR(N)
C
      INTEGER           I
C ----------------------------------------------------------------------
      NACEV = 0
      DO 100 I=1,N
         IF (EVERR(I).LT.TOL) THEN
            NACEV = NACEV + 1
            IF (MANCUR(I).EQ.0) THEN
               MANCUR(I) = N
            ENDIF
         ELSE
            MANCUR(I) = 0
         ENDIF
  100 CONTINUE
C
      NGDEV = 0
      DO 200 I=1,N
         IF (MANCUR(I).EQ.0)          GO TO 300
         IF (MANCUR(I).NE.MANPRV(I))  GO TO 300
         NGDEV = I
  200 CONTINUE
C
  300 RETURN
      END
      SUBROUTINE VRORT (G,V,TOL,N,NV,MAXIT,IPSW,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  VRORT                   GROUP 5 / PRIVATE
C
C     T A S K :  TO ORTHONORMALIZE A UNIT VECTOR  V  AGAINST THE NV
C                COLUMNVECTORS OF MATRIX  G  WHICH ARE ALREADY ORTHO-
C                NORMAL TO ONE ANOTHER.
C     A MULTI-PASS (MAXIT.LT.0) OR AN ITERATIVE (MAXIT.GT.0) GRAM-
C     SCHMIDT PROCEDURE IS USED.
C     PRINT :  IF IPSW.GT.4  SOME INFORMATION ABOUT THE PROCESS IS
C              PRINTED ON UNIT LPU
C
C     ROUTINES CALLED/REFERENCED :  ORTHO1   (SAM-3)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-07-17 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IPSW,LPU,MAXIT,N,NV
      DOUBLE PRECISION  TOL, G(N,*),V(N)
C
      INTEGER           IFLAG,IT,ITV,KOUNT,NIT
C
      EXTERNAL          ORTHO1
C ----------------------------------------------------------------------
      IERR  = 0
      IF (MAXIT.EQ.0)   GO TO 100
      IF (NV.LT.1)      GO TO 100
      ITV   = NV+1
      IFLAG = 1
C
      IF (MAXIT.LT.0)   THEN
C                                               ** MULTI-PASS  G-S
         NIT =-MAXIT
         IF (IPSW.GT.4) THEN
            WRITE(LPU,610) ITV
            WRITE(LPU,620) NIT
         ENDIF
         DO 20 IT=1,NIT
            CALL ORTHO1 (G,V,TOL,N,NV,IFLAG,IPSW,LPU,KOUNT,IERR)
            IF (IERR.LT.0) GO TO 90
   20    CONTINUE
C
      ELSE
C                                               ** ITERATIVE  G-S
         IF (IPSW.GT.4) THEN
            WRITE(LPU,610) ITV
            WRITE(LPU,630) TOL,MAXIT
         ENDIF
C                                               ** FIRST PASS,
C                                                  NO ITERATION
C
         CALL ORTHO1 (G,V,TOL,N,NV,IFLAG,IPSW,LPU,KOUNT,IERR)
         IF (IERR.LT.0)  GO TO  90
         IF (MAXIT.EQ.1) GO TO 100
         NIT   = MAXIT-1
         IFLAG = 2
         DO 40 IT=1,NIT
            IF (IPSW.GT.4) WRITE(LPU,640) IT
            CALL ORTHO1 (G,V,TOL,N,NV,IFLAG,IPSW,LPU,KOUNT,IERR)
            IF (IERR.LT.0)  GO TO 90
            IF (IPSW.GT.4)  WRITE(LPU,650) KOUNT
            IF (KOUNT.EQ.0) GO TO 100
   40    CONTINUE
         IERR = 1
         IF (LPU.GT.0) WRITE(LPU,660) ITV,MAXIT,TOL
C
      ENDIF
      GO TO 100
C ------------------------------------------------ ERROR EXIT
   90 IF (LPU.GT.0)  WRITE(LPU,690) ITV
C ------------------------------------------------
  100 RETURN
C ------------------------------------------------ FORMATS
C
  610 FORMAT(///5X,'REORTHOGONALIZING VECTOR NO.',I4)
  620 FORMAT(5X,'USING A MULTI-PASS (=',I3,') G-S PROCEDURE'/)
  630 FORMAT(5X,'USING AN ITERATIVE G-S PROCEDURE'/
     +       5X,'WITH TOLERANCE         =',1PE11.3,'  AND'/
     +       5X,'MAX. NO. OF ITERATIONS =',I3 //
     +       5X,'FIRST PASS :  NO ITERATION' / )
  640 FORMAT(5X,'ITERATION NO.',I5)
  650 FORMAT(5X,'NO. OF VECTORS CONTRIBUTING TO THE PURGING :',I5)
  660 FORMAT(///3X,'* WARNING FROM  S A M  LIBRARY ROUTINE VRORT'/
     +          5X,'REORTHOGONALIZATION OF VECTOR NO.',I5,/
     +          5X,'NOT ACCOMPLISHED IN',I4,'  ITERATIONS'/
     +          5X,'(WITH TOL =',1PE11.3,')' )
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE VRORT'/
     +       5X,'VECTOR NO',I5,' CANNOT BE MADE ORTHOGONAL' )
C
      END
      SUBROUTINE CHBDG (B,MSKY,N,KSB,NZDB,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  CHBDG                  GROUP 5 / PRIVATE
C
C     T A S K :  TO DETERMINE THE NUMBER OF ZERO DIAGONAL ELEMENTS IN
C                MATRIX B AND TO CHECK FOR NEGATIVE DEFINITENESS
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-10-11 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KSB,N,NZDB,            MSKY(*)
      DOUBLE PRECISION  B(*)
C
      INTEGER           I,II
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
      NZDB = 0
      DO 10 I=1,N
         II = I
         IF (KSB.EQ.1)  II = MSKY(I)
         IF (B(II).GT.ZERO)  GO TO 10
         IF (B(II).EQ.ZERO)  NZDB = NZDB+1
         IF (B(II).LT.ZERO)  IERR = IERR-1
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE SSSSVC (A,B,EVEC,V,MSKY,RAN,N,NIV,KSB)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SSSSVC                 GROUP 5 / PRIVATE
C
C     T A S K :  TO GENERATE NIV START VECTORS FOR SUBSPACE ITERATION
C                BASED ON THE DIAGONAL ELEMENTS OF MATRICES  A  AND  B
C
C     ROUTINES CALLED/REFERENCED :  RMINT (SAM-0)  AND   RANVEC (SAM-8)
C                                   ABS   (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-10-10 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           KSB,N,NIV,    MSKY(*)
      DOUBLE PRECISION  RAN,          A(*),B(*),EVEC(N,NIV),V(N)
C
      INTEGER           I,II,J,K
      DOUBLE PRECISION  G,ONE,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C
      EXTERNAL          RANVEC,RMINT
C ----------------------------------------------------------------------
      G = ZERO
C
      IF (NIV.GT.1) THEN
         CALL RMINT (EVEC,N,NIV,ZERO)
         DO 20 I=1,N
            II = MSKY(I)
            IF (A(II).EQ.ZERO) THEN
               V(I) = G
            ELSEIF (KSB.EQ.1) THEN
               V(I) = ABS(B(II)/A(II))
            ELSE
               V(I) = ABS(B(I)/A(II))
            ENDIF
            IF (V(I).GT.G)  G = V(I)
            EVEC(I,1) = ONE
   20    CONTINUE
C
         IF (NIV.GT.2) THEN
            DO 50 J=2,NIV-1
               G = ZERO
               K = J
               DO 40 I=1,N
                  IF (V(I).GT.G) THEN
                     G = V(I)
                     K = I
                  ENDIF
   40          CONTINUE
               EVEC(K,J) = ONE
               V(K)      = ZERO
   50       CONTINUE
         ENDIF
      ENDIF
C                                               ** LAST VECTOR IS RANDOM
      CALL RANVEC (EVEC(1,NIV),RAN,N)
C
      RETURN
      END
      SUBROUTINE EVARR (EVAL,EVEC,N,NEV,NV,LORD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EVARR                   GROUP 5 / PUBLIC
C
C     T A S K : TO REARRANGE THE FIRST NV OF NEV EIGENVALUES, STORED IN
C               EVAL(1), EVAL(2),...., EVAL(NEV), AND THE CORRESPONDING
C               EIGENVECTORS, STORED IN COLUMNS 1, 2,...., NEV OF EVEC.
C               REQUIREMENTS:  NEV.LE.N   AND   NV.LE.NEV
C     THE ORDERING OF THE EIGENVALUES IS DETERMINED BY THE VALUE OF
C     ARGUMENT LORD :
C          LORD = 1 :  INCREASING ALGEBRAIC VALUE (ASCENDING ORDER)
C          LORD = 2 :  INCREASING NUMERICAL VALUE
C          LORD = 3 :  DECREASING ALGEBRAIC VALUE (DESCENDING ORDER)
C          LORD = 4 :  DECREASING NUMERICAL VALUE
C
C     ROUTINES CALLED/REFERENCED :   MINDX     (SAM-8)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-03-22 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           LORD,N,NEV,NV
      DOUBLE PRECISION  EVAL(*),EVEC(N,*)
C
      INTEGER           I,J,K,NVL,        MINDX
      DOUBLE PRECISION  E
C
      EXTERNAL          MINDX
C ----------------------------------------------------------------------
      IF (NEV.LT.2)     GO TO 100
      IF (NEV.GT.N)     GO TO 100
      NVL = NV
      IF (NV.GT.NEV)    NVL = NEV
      IF (LORD.LT.1)    GO TO 100
      IF (LORD.GT.4)    GO TO 100
C
      DO 50 J=1,NVL
         K = MINDX(EVAL,J,NEV,LORD)
         IF (K.GT.J) THEN
            E       = EVAL(J)
            EVAL(J) = EVAL(K)
            EVAL(K) = E
            DO 25 I=1,N
               E         = EVEC(I,J)
               EVEC(I,J) = EVEC(I,K)
               EVEC(I,K) = E
   25       CONTINUE
         ENDIF
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE MDFCTA (A,B,MSKY,EPS,SHIFT,N,KSB,LPU,NND,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MDFCTA                  GROUP 5 / PRIVATE
C
C     T A S K :  TO SUBTRACT (IF SHIFT.NE.ZERO)  SHIFT*B  FROM THE
C                SYMMETRIC SKYLINE MATRIX  A  AND TO FACTORIZE THE
C                RESULTING MATRIX INTO  LA*DA*LAT.
C                THE NUMBER OF NEGATIVE ELEMENTS IN  DA  IS RETURNED IN
C                ARGUMENT NND.
C     MATRIX  B  MAY BE A SYMMETRIC SKYLINE MATRIX (KSB=1) OR A DIAGONAL
C     MATRIX (KSB=2)
C
C     ROUTINES CALLED/REFERENCED :  SKYSOL   (SAM-4)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-07-13 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KSB,LPU,N,NND,       MSKY(N)
      DOUBLE PRECISION  A(*),B(*),EPS,SHIFT
C
      INTEGER           I,J
      DOUBLE PRECISION  S,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          SKYSOL
C
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (SHIFT.NE.ZERO) THEN
C ----------------------------------------------------------------------
C  MODIFY MATRIX  A  BY SUBTRACTING  SHIFT*B
C ----------------------------------------------------------------------
         S = SHIFT
         IF (KSB.EQ.1) THEN
C                                               ** SKYLINE B
            J = MSKY(N)
            DO 25 I=1,J
               A(I) = A(I) - S*B(I)
   25       CONTINUE
C
         ELSEIF (KSB.EQ.2) THEN
C                                               ** DIAGONAL B
            DO 50 J=1,N
               I    = MSKY(J)
               A(I) = A(I) - S*B(J)
   50       CONTINUE
C
         ELSE
C                                               ** UNKNOWN B
            IERR =-5
            GO TO 100
         ENDIF
      ENDIF
C ----------------------------------------------------------------------
C  FACTORIZE  A, WHETHER MODIFIED OR NOT
C ----------------------------------------------------------------------
C
      CALL SKYSOL (A,B,MSKY,EPS,N,N,N,LPU,1,NND,IERR)
C
C
  100 RETURN
      END
      SUBROUTINE MDSHFT (A,B,MSKY,SHIFT,N,KSB)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MDSHFT                  GROUP 5 / PRIVATE
C
C     T A S K :  To subtract (if SHIFT.NE.ZERO)  SHIFT*B  from the
C                symmetric "skyline" matrix  A.
C     Matrix  B  may be a symmetric "skyline" matrix (KSB=1) or a diagonal
C     matrix (KSB=2)
C
C     ROUTINES CALLED/REFERENCED :  None.
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   96-12-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           KSB,N,       MSKY(N)
      DOUBLE PRECISION  A(*),B(*),SHIFT
C
      INTEGER           I,J
      DOUBLE PRECISION  S,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
C
      IF (SHIFT .NE. ZERO) THEN
C
C --- Modify matrix  A  by subtracting  SHIFT*B
C
         S = SHIFT
         IF (KSB .EQ. 1) THEN
C                                               ** skyline B
            J = MSKY(N)
            DO 25 I=1,J
               A(I) = A(I) - S*B(I)
   25       CONTINUE
C
         ELSEIF (KSB .EQ. 2) THEN
C                                               ** diagonal B
            DO 50 J=1,N
               I    = MSKY(J)
               A(I) = A(I) - S*B(J)
   50       CONTINUE
C
         ENDIF
      ENDIF
C
      RETURN
      END
