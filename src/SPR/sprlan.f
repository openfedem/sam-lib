C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRLAN (A,B,TOL,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                   MIP,MOP,MSING,EVL,EVERR,V,R,RWA,LWA,
     &                   SHIFT,N,NEVAL,NEVEC,MAXLAN,LPU,IPSW,IERR,
     &                   LMVMULT,LMMMULT,CMVMULT,CMMMULT)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRLAN                  GROUP 9 / PUBLIC
C
C     TASK :  To determine the NEVAL eigenvalues Lambda_i of the
C             generalized, symmetric eigenproblem
C                  (A - Lambda_i*B)*q_i = 0
C             that are closest to a specified SHIFT of origin.
C             The routine also returns NLV (=MOP(3)) ortogonal Lanczos
C             vectors (in V) or, optionally, approximations to NEVEC
C             eigenvectors (corresponding to the NEVEC first eigen-
C             values).  The eigenvectors, which are B-orthonormal, are
C             returned as the first NEVEC vectors in V.
C     A is a symmetric matrix stored in sparse (compressed) format.
C     B is either a symmetric sparse matrix or a diagonal matrix.
C     For both A and B the original matrices or their appropriate
C     factors are accepted as input (in arrays A and B).
C     A truncated Lanczos method, assuming all matrices in primary
C     storage, is used to transform the problem to tridiagonal, special
C     form.  The reduced, tridiagonal eigenproblem is solved by implicit
C     QL transformation.
C     Various reorthogonalization schemes may be specified.
C
C     This is the sparse version of LANCZ2.
C
C
C     ROUTINES CALLED/REFERENCED :  IMP, SCALE                  (SAM-0)
C                                   BIGMUL                      (SAM-3)
C                                   QLS3D, QLS3DV, EVARR        (SAM-5)
C                                   MGSOR1, ABRLER, REORL2      (SAM-5)
C                                   RANVEC, REARRV              (SAM-8)
C                                   SPRFC1, SPRFC2, SPRBS2      (SAM-9)
C                                   SPRDG1, SPRSHA, SPRSMM      (SAM-9)
C                                   SPRLP1, SPRLP2, SPRLP3      (SAM-9)
C                                   SPRLER, I7COPY              (SAM-9)
C                                   DCOPY, DSCAL, DDOT, DAXPY   (BLAS)
C                                   ABS, SQRT         (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   96-12-27 / 1.0
*
*     ACD - 99-01-06 - MSPAR(50) -> MSPAR(*),
*                      added superelement option, and modification in
*                      order to account for evaluation of internal
*                      stiffness only. The evaluation of the internal
*                      stiffness is implemented by setting all vectors
*                      in the computations to zero for the degrees of
*                      freedom corresponding to external dofs.
C
C     KMO - 03-02-25 - Added a separate sparse structure for matrix B.
C                      Replaced SPRZRN by SPRZR2.
C
C     KMO - 05-08-27 - Removed all calls to SPRZR2. The internal DOFs
C                      only option is now handled directly by the
C                      underlying solve routines when N = NEQ1 < NEQ.
C
C     KMO - 11-03-16 - Corrected the special case for N=1.
C
C     KMO - 14-06-18 - Using ICOPY and SCOPY/DCOPY from BLAS library.
C
C     KMO - 17-05-12 - More use of BLAS where possible.
C                      Use I7COPY instead of ICOPY, and some related
C                      changes to facilitate using INTEGER*8 as default.
C                      Added subroutines SPRLCK and SPRLTD.
C
C     KMO - 17-11-24 - Removed conditionals (single/double precision).
C
C     KMO - 21-08-17 - Replaced SPRMDA by SPRSHA which also can handle
C                      two matrices with different sparsity patterns.
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,IPSW,LPU,MAXLAN,N,NEVAL,NEVEC,
     &                  LWA(*),MIP(7),MOP(10),MSING(*),
     &                  MPARA(*),MTREEA(*),MSIFA(*),
     &                  MPARB(*),MTREEB(*),MSIFB(*)
      DOUBLE PRECISION  SHIFT,A(*),B(*),EVERR(MAXLAN),EVL(MAXLAN),
     &                  R(*),RWA(*),TOL(3),V(N,MAXLAN)
C
      include 'integer_kind.inc'
C
      INTEGER           I,IFLAG,IP,IPA,IPB,IPD,IPE,IPK,IPKM1,IPKP1,IPS,
     &                  J,JORT,JP,K,KEIG,KEVEX,KEX,KM1,KORT,KPD,KRST,
     &                  KSA,KSB,KSR,LERR,LLP,LPC,LPP,LPW,NACEV,NEVEXT,
     &                  NGDEV,NLV,NSDIG,NVC
      INTEGER(kind=i4)  JERR,NEQ,NN
      DOUBLE PRECISION  ALFA,AUX,BETA,ONE,RAN,RLPR,
     &                  SIGMA,SMALL,TEN,TRSH,VPMAX,ZERO
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 , TEN = 10.0D0 )
C
      DOUBLE PRECISION  DDOT
      INTEGER(kind=i4)  IMP
C
      EXTERNAL          DDOT
      EXTERNAL          IMP
      EXTERNAL          LMVMULT,LMMMULT, CMVMULT,CMMMULT
C ----------------------------------------------------------------------
C
      IF (IPSW.GE.0) THEN
         LLP = LPU
      ELSE
         LLP = -1
      ENDIF
C
#if FT_DEBUG > 1
      IF (LLP.GT.0) WRITE(LLP,"(/A)") 'ENTERING SPRLAN'
#endif
C
      IF (IPSW.GT.0 .AND. NEVAL.GT.0) THEN
         CALL SPRLP1 (MPARA(8)+MPARA(16),MIP,SHIFT,
     &                N,NEVAL,NEVEC,MAXLAN,LPU)
      ENDIF
C
C --- Check input parameters -------------------------------------------
C
      IERR = 0
      KSA  = MIP(1)
      KSB  = MIP(2)
      KSR  = ABS(MIP(3))
      KORT = MIP(4)
      KEX  = MIP(5)
      KPD  = MIP(6)
      CALL I7COPY (10,0,0,MOP)
      CALL SPRLCK (KSA,KSB,MIP(3),KORT,KEX,KPD,N,MPARA(8),MPARB(8),
     &             NEVAL,NEVEC,MAXLAN,LLP,IERR)
      IF (IERR .LT. 0) GOTO 1000
C
      IF (N .EQ. 1) THEN
C
C ------ Special case:  N = 1 ------------------------------------------
C
         ALFA = A(1+MPARA(8))
         IF (KSB .EQ. 1) THEN
            BETA = B(1+MPARB(8))
         ELSE
            BETA = B(1)
         ENDIF
         IF (BETA .EQ. ZERO) THEN
            CALL SPRLER (23,I,I,I,LLP,IERR)
            IF (LLP .GT. 0) write(LLP,*) '    A =',ALFA,'    B =',BETA
            GOTO 1000
         ENDIF
         EVL(1) = ALFA / BETA
         IF (NEVEC .GT. 0) THEN
            IF (BETA .GT. ZERO) THEN
               V(1,1) = ONE / SQRT(BETA)
            ELSE
               V(1,1) = ONE
            ENDIF
         ENDIF
         MOP(1) = 1
         IF (EVL(1) .LT. SHIFT)  MOP(2) = 1
         GOTO 1000
      ENDIF
C
C --- Set/initiate parameters, variables and pointers ------------------
C
      NEVEXT = 0
      NACEV  = 0
      NGDEV  = 0
      NLV    = 0
      NEQ    = int(N,i4)
      NSDIG  = int(IMP(3))
      RLPR   = TEN**(-NSDIG)
      TRSH   = SQRT(RLPR)
      I      = NSDIG - NSDIG/5
      SMALL  = TEN**(-I)
      RAN    = 0.1234567891234567D0
C                                                ! default tolerances
      I    = 2*NSDIG/3
      AUX  = TEN**(-I)
      IF (TOL(1) .LE. ZERO) TOL(1) = AUX
      IF (TOL(2) .LE. ZERO) TOL(2) = AUX
      IF (TOL(3) .LT. ONE)  TOL(3) = TEN
C                                                ! start vector
      IF (KSR .EQ. 1) THEN
         CALL RANVEC (R,RAN,NEQ)
      ELSEIF (KSR .EQ. 3) THEN
         CALL DCOPY (NEQ,ONE,0_i4,R,1_i4)
      ELSEIF (KSR .EQ. 2) THEN
         IF (KSB .EQ. 1) THEN
            CALL SPRDG1 (MPARB,MTREEB,MSIFB,LWA,LLP,LERR)
            IF (LERR .LT. 0) THEN
               CALL SPRLER (20,I,I,I,LLP,IERR)
               GOTO 1000
            ENDIF
            DO 10 I = 1, N
               R(I) = B(LWA(I))
   10       CONTINUE
         ELSE
            DO 20 I = 1, N
               R(I) = B(I)
   20       CONTINUE
         ENDIF
      ENDIF
C
C --- If necessary, modify and factorize matrix A,
C     if KPD=2 factorize matrix B
C
      IF (ABS(SHIFT) .LT. SMALL*SMALL) THEN
         SHIFT = ZERO
      ELSE
         IF (KPD .EQ. 1) THEN
            CALL SPRLER (2,I,I,I,LLP,IERR)
         ENDIF
         IF (KSA .GT. 1 .OR. KSB .LT. 0) THEN
            CALL SPRLER (19,KSA,KSB,I,LLP,IERR)
         ELSE
            CALL SPRSHA (MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,LWA,
     &                   A,B,SHIFT,N,KSB,LLP,IERR)
         ENDIF
         IF (IERR .LT. 0) GOTO 1000
      ENDIF
C
      IF (MIP(7) .LE. 0) THEN
         CALL SPRLF1 (A,B,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                KSA,KSB,KPD,N,MOP(6),NEVAL,TOL(2),
     &                MSING,RWA,LWA,LLP,IERR,
     &                LMVMULT,LMMMULT,CMVMULT,CMMMULT)
      ELSE
         CALL SPRLF2 (A,B,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                KSA,KSB,KPD,MIP(7),N,MOP(6),NEVAL,TOL(2),TOL(3),
     &                MSING,V,RWA,LWA,LLP,IERR,
     &                LMVMULT,LMMMULT,CMVMULT,CMMMULT)
      ENDIF
      IF (IERR .LT. 0) GOTO 1000
C
      IF (MIP(1) .EQ. 1) THEN
         IF (MIP(6) .EQ. 2 .OR. NEVAL .EQ. 0) THEN
            MOP(4) = MPARA(26)
C                                                ! exit if NEVAL=0
            IF (NEVAL.EQ.0) GOTO 1000
         ENDIF
      ENDIF
C
C                                                ! pointers in RWA
      IPA   = 1
      IPB   = IPA + MAXLAN
      IPD   = IPB + MAXLAN
      IPS   = IPD + MAXLAN
      IPE   = IPS + MAXLAN
      IPK   = IPE + MAXLAN*MAXLAN
      IPKM1 = IPK + MAXLAN
      IPKP1 = IPKM1 + MAXLAN
C                                                ! pointers in LWA
      LPP = 1
      LPC = LPP + MAXLAN
      LPW = LPC + MAXLAN
C                                                ! initialize (to zero)
      NN  = int(MAXLAN*(MAXLAN+7),i4)
      CALL DCOPY (NN,ZERO,0_i4,RWA,1_i4)
      NN  = int(MAXLAN,i4)
      CALL DCOPY (NN,ZERO,0_i4,EVL,1_i4)
      CALL I7COPY (2*MAXLAN,0,0,LWA)
C
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
C
         CALL SPRSMM (A,B,R(1),R(N+1),R(N+1),LWA(LPW),
     &                MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                N,KSA,KSB,IFLAG,LLP,LERR)
         IF (LERR .LT. 0) THEN
            CALL SPRLER (24,K,K,K,LLP,IERR)
            GOTO 1000
         ENDIF
C
      ENDIF
C                                                ! length of R
      AUX = DDOT(NEQ,R,1_i4,R,1_i4)
      IF (AUX .LT. SMALL) THEN
         CALL RANVEC (R,RAN,NEQ)
         AUX = DDOT(NEQ,R,1_i4,R,1_i4)
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
#if FT_DEBUG > 1
      IF (LLP.GT.0) then
         WRITE(LLP,"(/' Start Lanczos iteration',3I4)") K,KEVEX,MAXLAN
      ENDIF
#endif
C
C                                                ! form Lanczos vector
      AUX = ONE/BETA
      CALL SCALE (R,V(1,K),NEQ,1_i4,AUX)
C
      IF (KRST.EQ.0 .AND. K.GT.1) THEN
C
C --- Reorthogonalization - if necessary and/or if specified -----------
C
         KM1 = K-1
         IF (KORT .GT. 0) THEN
            MOP(8) = MOP(8)+1
            DO 150 I=1,KORT
               CALL MGSOR1 (V,V(1,K),NEQ,int(KM1,i4))
               IF (IPSW .GT. 1) WRITE (LPU,6010) K
               MOP(9) = MOP(9)+KM1
  150       CONTINUE
         ELSEIF (KORT.EQ.-1 .AND. K.GT.2) THEN
            IF (JORT .EQ. 1) THEN
               CALL MGSOR1 (V,V(1,K),NEQ,int(KM1,i4))
               IF (IPSW .GT. 1) WRITE (LPU,6010) K
               MOP(8) = MOP(8)+1
               MOP(9) = MOP(9)+KM1
               NN = int(K,i4)
               CALL DCOPY (NN,RLPR,0_i4,RWA(IPK),1_i4)
               CALL DCOPY (NN,RLPR,0_i4,RWA(IPKM1),1_i4)
               JORT = 0
            ELSE
               CALL REORL2 (RWA(IPA),RWA(IPB),RWA(IPK),
     &                      RWA(IPKM1),RWA(IPKP1),BETA,
     &                      RLPR,TRSH,int(KM1,i4),VPMAX)
               IF (IPSW .GT. 1) WRITE (LPU,6020) K,VPMAX
               IF (VPMAX .GT. TRSH) THEN
                  CALL MGSOR1 (V,V(1,K),NEQ,int(KM1,i4))
                  IF (IPSW .GT. 1) WRITE (LPU,6010) K
                  MOP(8) = MOP(8)+1
                  MOP(9) = MOP(9)+KM1
                  JORT   = 1
               ENDIF
            ENDIF
         ELSEIF (KORT .EQ. -2) THEN
            JORT = JORT+1
            IF (JORT .GT. 0) THEN
               CALL MGSOR1 (V,V(1,K),NEQ,int(KM1,i4))
               IF (IPSW .GT. 1) WRITE (LPU,6010) K
               MOP(8) = MOP(8)+1
               MOP(9) = MOP(9)+KM1
            ENDIF
            IF (JORT .EQ. 2)  JORT = -2
         ENDIF
C
      ENDIF
C                                                ! prepare new R
      CALL SPRSMM (A,B,V(1,K),R(1),R(N+1),LWA(LPW),
     &             MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &             N,KSA,KSB,IFLAG,LLP,LERR)
      IF (LERR.LT.0) THEN
         CALL SPRLER (24,K,K,K,LLP,IERR)
         GOTO 1000
      ENDIF
C
      IF (KRST.EQ.0 .AND. K.GT.1) THEN
         AUX = -BETA
         CALL DAXPY (NEQ,AUX,V(1,K-1),1_i4,R,1_i4)
      ENDIF
      KRST = 0
C                                                ! element alpha_k of
C                                                  tridiag. form
      ALFA = DDOT(NEQ,V(1,K),1_i4,R,1_i4)
      RWA(IPA+K-1) = ALFA
C
      NLV = K
C
      IF (KEVEX.EQ.0 .OR. K.EQ.MAXLAN) THEN
C
C --- Solve reduced, tridiagonal eigenproblem --------------------------
C
         IF (NEVEXT .GT. 0) CALL I7COPY (K,LWA(LPC),1,LWA(LPP))
C
         CALL SPRLTD (int(NLV,i4),LWA(LPC),LWA(LPP),
     &                EVERR,EVL,RWA(IPA),RWA(IPB),RWA(IPD),RWA(IPS),
     &                TOL(1),NGDEV,NACEV,IPSW,LLP,IERR)
         IF (IERR .LT. 0) GOTO 1000
C
C                                                ! next eigenvalue
C                                                  extraction
         NEVEXT = NEVEXT+1
         IF (NEVEXT .EQ. 1 .AND. NEVAL .GT. 1) THEN
            KEVEX = NEVAL/2
         ELSE
            KEVEX = KEX
         ENDIF
#if FT_DEBUG > 1
         IF (LLP.GT.0) then
            WRITE(LLP,"(' Solve tridiagonal eigenproblem',4I4)")
     &      NEVEXT,KEVEX,NGDEV,NACEV
         ENDIF
#endif
C                                                ! loop exit
C
         IF (NGDEV.GE.NEVAL .OR. K.EQ.MAXLAN)  GOTO 300
      ENDIF
C
C --- Next (unscaled) Lanczos vector -----------------------------------
      AUX = -ALFA
      CALL DAXPY (NEQ,AUX,V(1,K),1_i4,R,1_i4)
C                                                ! length of R
      AUX  = DDOT(NEQ,R,1_i4,R,1_i4)
      BETA = SQRT(AUX)
      AUX  = ABS(RWA(IPA))
C
      IF (BETA .LT. RLPR*AUX) THEN
C
C ------ Breakdown; generate new (random) start vector and make it
C        ortonormal to previous Lanczos vectors, and start again
C
         CALL RANVEC (R,RAN,NEQ)
         AUX = DDOT(NEQ,R,1_i4,R,1_i4)
         AUX = ONE/SQRT(AUX)
         CALL SCALE (R,V(1,K),NEQ,1_i4,AUX)
C
         CALL MGSOR1 (V,R,NEQ,int(K,i4))
         BETA = ONE
         RWA(IPB+K) = ZERO
C
         MOP(10) = MOP(10) + 1
         IF (MOP(10) .GT. 10) THEN
            CALL SPRLER (31,K,K,K,LLP,LERR)
         ELSE
            KRST = 1
         ENDIF
C
      ELSE
         RWA(IPB+K) = BETA
      ENDIF
C
      GOTO 100
C
C ================================================ End of Lanczos loop
C
  300 SMALL = RLPR*ABS(EVL(1))
C
      IF (NGDEV .EQ. 0) THEN
         CALL SPRLER (32,MAXLAN,K,K,LLP,IERR)
         GOTO 1000
      ENDIF
      IF (NGDEV .LT. NEVAL) THEN
         CALL SPRLER (1,NGDEV,MAXLAN,NEVAL,LLP,LERR)
      ENDIF
C                                                ! restore eigenvalues
C                                                  Lambda_i
      DO 350 I=1,NLV
         IF (ABS(EVL(I)) .GT. SMALL) THEN
            EVL(I) = SHIFT + ONE/EVL(I)
         ELSE
            IF (I-1 .LT. NGDEV) NGDEV=I-1
            IF (I-1 .LT. NEVAL) THEN
               CALL SPRLER (3,I-1,I,I,LLP,LERR)
               GOTO 400
            ENDIF
         ENDIF
  350 CONTINUE
C
  400 MOP(1) = NGDEV
      MOP(2) = NACEV
      MOP(3) = NLV
      MOP(5) = MIN(NEVEC,NGDEV)
      MOP(7) = NEVEXT
C
C                                                ! Print out eigenvalues
      IF (IPSW .GT. 0) CALL SPRLP2 (EVL,EVERR,TOL,LWA(LPC),MOP,NLV,LPU)
C
      IF (NEVEC .LE. 0) GOTO 1000
C
C --- Determine (B-orthonormal) approx. to the NVC first eigenvectors
C
      NVC = MOP(5)
C                                                ! re-solve small
C                                                  eigenproblem
      NN = int(NLV,i4)
      CALL QLS3DV (RWA(IPA),RWA(IPB),RWA(IPE),NN,int(LLP,i4),JERR)
      IF (JERR .LT. 0_i4) THEN
         CALL SPRLER (25,NLV,NLV,NLV,LLP,IERR)
         GOTO 1000
      ENDIF
C                                                ! rearrange ordering
      CALL EVARR (RWA(IPA),RWA(IPE),NN,NN,int(NVC,i4),4_i4)
C
C --- Restore eigenvectors of special problem -----------------------
C
      IF (KPD .EQ. 1) THEN
C                                                ! scale eigenvectors
         SIGMA = TRSH
         IPB   = IPA+NVC-1
         DO 410 IP=IPA,IPB
            IF (RWA(IP) .LT. SIGMA)  SIGMA = RWA(IP)
  410    CONTINUE
         IF (SIGMA .LT. TRSH) THEN
            SIGMA = ABS(SIGMA) + TRSH
            DO 420 IP=IPA,IPB
               RWA(IP) = RWA(IP) + SIGMA
  420       CONTINUE
         ELSE
            SIGMA = ZERO
         ENDIF
         JP = IPE
         DO 430 IP=IPA,IPB
            AUX = ONE / SQRT(RWA(IP))
            CALL DSCAL (NN,AUX,RWA(JP),1_i4)
            RWA(IP) = RWA(IP) - SIGMA
            JP = JP+NLV
  430    CONTINUE
      ENDIF
C
      CALL BIGMUL (V,RWA(IPE),V,R,NEQ,NN,int(NVC,i4),0_i4)
C
C
C --- Restore desired eigenvectors ----------------------------------
C
      IF (KPD .EQ. 1) THEN
C
         CALL SPRBS2 (MPARA,MTREEA,MSIFA,A,V,N,NVC,RWA,LLP,LERR)
         IF (LERR .LT. 0) THEN
            CALL SPRLER (27,NVC,NVC,NVC,LLP,IERR)
            GOTO 1000
         ENDIF
C
      ELSEIF (KSB .EQ. -1) THEN
C                                                ! sparse B
C
         CALL SPRBS2 (MPARB,MTREEB,MSIFB,B,V,N,NVC,RWA,LLP,LERR)
         IF (LERR .LT. 0) THEN
            CALL SPRLER (27,NVC,NVC,NVC,LLP,IERR)
            GOTO 1000
         ENDIF
C
      ELSEIF (MOP(6) .EQ. 0) THEN
C                                                ! diagonal B
         DO 450 J=1,NVC
            DO 440 I=1,N
               V(I,J) = V(I,J) / B(I)
  440       CONTINUE
  450    CONTINUE
C
      ELSE
C
         DO 460 J=1,NVC
C
            CALL SPRSMM (A,B,V(1,J),R,R,LWA(LPW),
     &                   MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                   N,KSA,KSB,1,LLP,LERR)
            IF (LERR .LT. 0) THEN
               CALL SPRLER (26,J,J,J,LLP,IERR)
               GOTO 1000
            ENDIF
C
            AUX = ONE / RWA(IPA-1+J)
            CALL DSCAL (NEQ,AUX,V(1,J),1_i4)
C
  460    CONTINUE
C
      ENDIF
C
C ------------------------------------------------ Exit
C
 1000 CONTINUE
#if FT_DEBUG > 1
      IF (LLP.GT.0) WRITE(LLP,"(A/)") 'LEAVING  SPRLAN'
#endif
      RETURN
C
 6010 FORMAT (5X,'Lanczos step',I5,' :  reorthogonalize')
 6020 FORMAT (5X,'Lanczos step',I5,' :  max. inner product =',1PE12.4)
C
      END


      SUBROUTINE SPRLCK (KSA,KSB,KSR,KORT,KEX,KPD,N,NEQA,NEQB,
     &                   NEVAL,NEVEC,MAXLAN,LPU,IERR)
C
      IMPLICIT           NONE
C
      INTEGER            KSA,KSB,KSR,KORT,KEX,KPD,N,NEQA,NEQB,
     &                   NEVAL,NEVEC,MAXLAN,LPU,IERR
C
C     SPRLCK
C
C     THIS SUBROUTINE CHECKS INPUT PARAMETER CONSISTENCY FOR SPRLAN.
C
C     CREATED   : MAY  19, 2017 (KMO)
C
      INTEGER            I

      IERR = 0
      IF (KSA.LT.1 .OR. KSA.GT.3) CALL SPRLER (11,KSA,KSB,I,LPU,IERR)
      I = ABS(KSB)
      IF (I.LT.1   .OR.   I.GT.3) CALL SPRLER (11,KSA,KSB,I,LPU,IERR)
      IF (ABS(KORT) .GT. 2)       CALL SPRLER (12,KORT,I,I,LPU,IERR)
      IF (KEX.LT.1 .OR. KEX.GT.5) CALL SPRLER (13,KEX,I,I,LPU,IERR)
      IF (KPD.LT.1 .OR. KPD.GT.2) CALL SPRLER (14,KPD,I,I,LPU,IERR)
      IF (KPD .EQ. 1) THEN
         IF (KSA.EQ.2 .OR. KSB.LT.0) THEN
            CALL SPRLER (15,KPD,KSA,KSB,LPU,IERR)
         ENDIF
      ENDIF
      IF (KSB.LT.0 .AND. ABS(KSR).EQ.2) THEN
         CALL SPRLER (16,KSB,KSR,I,LPU,IERR)
      ENDIF
      IF (NEVAL.EQ.0 .AND. KSA.GT.1) THEN
         CALL SPRLER (17,KSA,NEVAL,NEVAL,LPU,IERR)
      ENDIF
      IF (NEVAL .GT. MAXLAN) THEN
         CALL SPRLER (18,NEVAL,NEVEC,MAXLAN,LPU,IERR)
      ENDIF
      IF (MAXLAN.LT.1 .AND. NEVAL.GT.0) THEN
         CALL SPRLER (18,NEVAL,NEVEC,MAXLAN,LPU,IERR)
      ENDIF
      IF (MAXLAN .GT. N)    CALL SPRLER (18,NEVAL,NEVEC,MAXLAN,LPU,IERR)
      IF (NEVEC .GT. NEVAL) CALL SPRLER (18,NEVAL,NEVEC,MAXLAN,LPU,IERR)
      IF (I.EQ.1 .AND. NEQA.NE.NEQB) THEN
         CALL SPRLER (33,NEQA,NEQB,I,LPU,IERR)
      ENDIF

      RETURN
      END


      SUBROUTINE SPRLF1 (A,B,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                   KSA,KSB,KPD,N,NN,NEVAL,EPS,
     &                   MSING,RWA,LWA,LPU,IERR,
     &                   LMVMULT,LMMMULT,CMVMULT,CMMMULT)
C
      IMPLICIT           NONE
C
      INTEGER            KSA,KSB,KPD,N,NN,NEVAL,LWA(*),
     &                   MPARA(*),MTREEA(*),MSIFA(*),
     &                   MPARB(*),MTREEB(*),MSIFB(*),
     &                   MSING(*),LPU,IERR
      DOUBLE PRECISION   EPS,A(*),B(*),RWA(*)
C
C     SPRLF1
C
C     THIS SUBROUTINE ADMINISTERS THE MATRIX FACTORIZATIONS NEEDED
C     FOR THE LANCZOS EIGENSOLVER. WHEN A SINGULARITY OCCUR, THE
C     FACTORIZATION IS ABORTED.
C
C     CREATED   : JUL. 04, 2003 (KMO)
C
      INTEGER            I,LERR
      DOUBLE PRECISION   ZERO,TPR(3)
      PARAMETER        ( ZERO = 0.0D0 )
C
      EXTERNAL           LMVMULT,LMMMULT,CMVMULT,CMMMULT
C
      NN     = 0
      IERR   = 0
      TPR(1) = EPS
      MSING(1) = 0
      MPARA(27) = 0
      MPARB(27) = 0
C
C                                                ! matrix A
      IF (KSA .EQ. 1) THEN
         IF (KPD .EQ. 1 .AND. NEVAL .GT. 0) THEN
C
            KSA = 3
            CALL SPRFC2 (MPARA,MTREEA,MSIFA,A,TPR,LWA,RWA,LPU,IERR,
     &                   CMVMULT,CMMMULT)
C
         ELSE
C
            KSA = 2
            CALL SPRFC1 (MPARA,MTREEA,MSIFA,A,TPR,LWA,RWA,LPU,IERR,
     &                   LMVMULT,LMMMULT)
C
         ENDIF
         IF (IERR .LT. 0) THEN
            CALL SPRLER (21,I,I,I,LPU,IERR)
            IF (MPARA(27) .GT. 0) THEN
               MSING(1) = MPARA(27)
               MSING(2) = -1
            ENDIF
            RETURN
         ENDIF
      ENDIF
C
C                                                ! matrix B
      IF (KSB .GT. 0 .AND. KPD .EQ. 2) THEN
         IF (KSB .EQ. 1) THEN
C
            KSB = -1
            CALL SPRFC2 (MPARB,MTREEB,MSIFB,B,TPR,LWA,RWA,LPU,LERR,
     &                   CMVMULT,CMMMULT)
            IF (LERR .LT. 0) THEN
               IERR = LERR
               CALL SPRLER (22,KSB,I,I,LPU,IERR)
               IF (MPARB(27) .GT. 0) THEN
                  MSING(1) = MPARA(27)
                  MSING(2) = -1
               ENDIF
            ENDIF
C
         ELSE
C                                                ! B is diagonal
            KSB = -2
            DO 100 I = 1, N
               IF (B(I) .GT. ZERO) THEN
                  B(I) = SQRT(B(I))
               ELSEIF (B(I) .EQ. ZERO) THEN
                  NN = NN + 1
               ELSE
                  IERR = -3
                  MPARB(27) = I
                  MSING(1) = I
                  MSING(2) = -1
                  CALL SPRLER (22,KSB,I,I,LPU,IERR)
                  RETURN
               ENDIF
  100       CONTINUE
C
         ENDIF
      ELSEIF (ABS(KSB) .EQ. 2) THEN
C
         DO 200 I = 1, N
            IF (B(I) .EQ. ZERO) NN = NN + 1
  200    CONTINUE
C
      ENDIF
C
      RETURN
      END


      SUBROUTINE SPRLF2 (A,B,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                   KSA,KSB,KPD,KPS,N,NN,NEVAL,EPS,SCALE,
     &                   MSING,DCOR,RWA,LWA,LPU,IERR,
     &                   LMVMULT,LMMMULT,CMVMULT,CMMMULT)
C
      IMPLICIT           NONE
C
      INTEGER            KSA,KSB,KPD,KPS,N,NN,NEVAL,LWA(*),
     &                   MPARA(*),MTREEA(*),MSIFA(*),
     &                   MPARB(*),MTREEB(*),MSIFB(*),
     &                   MSING(*),LPU,IERR
      DOUBLE PRECISION   EPS,SCALE,A(*),B(*),RWA(*),DCOR(2,N)
C
C     SPRLF2
C
C     THIS SUBROUTINE ADMINISTERS THE MATRIX FACTORIZATIONS NEEDED
C     FOR THE LANCZOS EIGENSOLVER. WHEN A SINGULARITY OCCUR, THE
C     AFFECTED DIAGONAL ELEMENT IS MODIFIED SUCH THAT THE FACTORIZATION
C     CAN BE COMPLETED AND ALL SINGULARITIES CAN BE FOUND.
C     THE EIGENSOLVER IS THEN ABORTED ONLY IF KPS=1, OR A NEGATIVE
C     DEFINITE MATRIX IS DETECTED.
C
C     CREATED   : JUL. 04, 2003 (KMO)
C
      INTEGER            I,LERR,NSING,NRING
      DOUBLE PRECISION   ZERO,TPR(4)
      PARAMETER        ( ZERO = 0.0D0 )
C
      EXTERNAL           LMVMULT,LMMMULT,CMVMULT,CMMMULT
C
      NN     = 0
      IERR   = 0
      NSING  = 0
      NRING  = 0
      TPR(1) = EPS
      TPR(4) = SCALE
      MPARA(27) = 0
      MPARB(27) = 0
C
C                                                ! matrix A
      IF (KSA .EQ. 1) THEN
         IF (KPD .EQ. 1 .AND. NEVAL .GT. 0) THEN
C
            KSA = 3
            CALL SPRFC4 (MPARA,MTREEA,MSIFA,A,TPR,LWA,RWA,DCOR,LPU,IERR,
     &                   CMVMULT,CMMMULT)
C
         ELSE
C
            KSA = 2
            CALL SPRFC3 (MPARA,MTREEA,MSIFA,A,TPR,LWA,RWA,DCOR,LPU,IERR,
     &                   LMVMULT,LMMMULT)
C
         ENDIF
C
         IF (IERR .LT. 0) CALL SPRLER (21,I,I,I,LPU,IERR)
         IF (MPARA(27) .NE. 0) CALL SPRF3S (N,NSING,NRING,MSING,DCOR)
         IF (MPARA(27) .GT. 0 .AND. IERR .EQ. -2) RETURN
C
      ENDIF
C
C                                                ! matrix B
      IF (KSB .GT. 0 .AND. KPD .EQ. 2) THEN
         IF (KSB .EQ. 1) THEN
C
            KSB = -1
            CALL SPRFC4 (MPARB,MTREEB,MSIFB,B,TPR,LWA,RWA,DCOR,LPU,LERR,
     &                   CMVMULT,CMMMULT)
            IF (LERR .LT. 0) THEN
               CALL SPRLER (22,KSB,I,I,LPU,LERR)
               IF (LERR .NE. -3 .OR. IERR .GE. 0) IERR = LERR
            ENDIF
C
            IF (MPARB(27) .NE. 0) CALL SPRF3S (N,NSING,NRING,MSING,DCOR)
C
         ELSE
C                                                ! B is diagonal
            KSB = -2
            DO 100 I = 1, N
               IF (B(I) .GT. ZERO) THEN
                  B(I) = SQRT(B(I))
                  DCOR(2,I) = ZERO
               ELSEIF (B(I) .EQ. ZERO) THEN
                  DCOR(2,I) = ZERO
                  NN = NN + 1
               ELSE
                  IERR = -3
                  MPARB(27) = I
                  DCOR(2,I) = B(I)
                  CALL SPRLER (22,KSB,I,I,LPU,IERR)
                  CALL SPRF3S (I,NSING,NRING,MSING,DCOR)
                  RETURN
               ENDIF
  100       CONTINUE
C
         ENDIF
      ELSEIF (ABS(KSB) .EQ. 2) THEN
C
         DO 200 I = 1, N
            IF (B(I) .EQ. ZERO) NN = NN + 1
  200    CONTINUE
C
      ENDIF
C
      IF (KPS .GE. 2 .AND. MPARA(27) .LE. 0 .AND. MPARB(27) .LE. 0) THEN
C        Ignore the singularities encountered
         IF (IERR .GE. -3 .AND. IERR .LE. -2) THEN
            IF (KPS .EQ. 3 .OR. NSING .EQ. NRING) IERR = NSING
         ENDIF
      ENDIF
C
      RETURN
      END


      SUBROUTINE SPRLTD (NLV,LPC,LPP,EVERR,EVL,RWA,RWB,RWD,RWS,EPS,
     &                   NGDEV,NACEV,IPSW,LPU,IERR)
C
      IMPLICIT NONE
C
      include 'integer_kind.inc'
C
      INTEGER(kind=i4)  NLV
      INTEGER           LPC(NLV),LPP(NLV),NGDEV,NACEV,IPSW,LPU,IERR
      DOUBLE PRECISION  EVERR(NLV),EVL(NLV),
     &                  RWA(NLV),RWB(NLV),RWD(NLV),RWS(NLV),EPS
C
C     SPRLTD
C
C     THIS SUBROUTINE SOLVES THE REDUCED, TRIDIAGONAL EIGENPROBLEM.
C
C     CREATED   : MAY  19, 2017 (KMO)
C
      INTEGER(kind=i4)  LERR
      INTEGER           I, NLV_I
C
      NLV_I = int(NLV)
      NGDEV = 0
      NACEV = 0
      IERR  = 0
C
      CALL DCOPY (NLV,RWA(1),1_i4,RWD,1_i4)
      CALL DCOPY (NLV,RWB(1),1_i4,RWS,1_i4)
C
      CALL QLS3D (RWD,RWS,NLV,int(LPU,i4),LERR)
      IF (LERR .LT. 0_i4) THEN
         CALL SPRLER (25,NLV_I,0,0,LPU,IERR)
         RETURN
      ENDIF
      CALL REARRV (RWD,NLV,4_i4)
C
C                                                ! relative error
      CALL ABRLER (EVL,RWD,EVERR,NLV)
      CALL DCOPY (NLV,RWD(1),1_i4,EVL,1_i4)
C
C                                                ! check tolerance
      DO 100 I = 1, NLV_I
         IF (EVERR(I) .LT. EPS) THEN
            NACEV = NACEV + 1
            IF (LPC(I) .EQ. 0) THEN
               LPC(I) = NLV_I
            ENDIF
            IF (NGDEV .LE. 0) THEN
               IF (LPC(I) .NE. 0 .AND. LPC(I) .EQ. LPP(I)) THEN
                  NGDEV = -I
               ELSE
                  NGDEV = -NGDEV
               ENDIF
            ENDIF
         ELSE
            LPC(I) = 0
            NGDEV = ABS(NGDEV)
         ENDIF
 100  CONTINUE
      NGDEV = ABS(NGDEV)
C
      IF (IPSW .GT. 2) THEN
         CALL SPRLP3 (RWA,RWB,EVL,EVERR,LPC,EPS,NLV_I,NGDEV,LPU,IPSW-2)
      ENDIF
C
      RETURN
      END
