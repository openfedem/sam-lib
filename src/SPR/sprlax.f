C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRLAX (A,B,TOL,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &                   IOP,EVL,V,RWORK,IWORK,
     &                   SHIFT,N,NEVAL,MAXLAN,LPU,IPSW,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRLAX                  GROUP 9 / PUBLIC
C
C     TASK :  To determine the NEVAL eigenvalues Lambda_i of the
C             generalized, symmetric eigenproblem
C                  (A - Lambda_i*B)*q_i = 0
C             that are closest to a specified SHIFT of origin.
C             The routine also returns NLV (=MOP(3)) orthogonal Lanczos
C             vectors (in V) or, optionally, approximations to NEVAL
C             eigenvectors (corresponding to the NEVAL first eigen-
C             values).  The eigenvectors, which are B-orthonormal,
C             are returned as the first NEVAL vectors in V.
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,IOP,IPSW,LPU,MAXLAN,N,NEVAL,IWORK(*),
     &                  MPARA(*),MTREEA(*),MSIFA(*),
     &                  MPARB(*),MTREEB(*),MSIFB(*)
      DOUBLE PRECISION  SHIFT,A(*),B(*),EVL(MAXLAN),
     &                  RWORK(*),TOL(3),V(N,MAXLAN)
C
      INTEGER           MIP(7),MOP(10),MSING(2)
C
      EXTERNAL          LMVY18,LMMY18,CMVY18,CMMY18
C
      MIP(1) =  1
      MIP(2) =  1
      MIP(3) = -1
      MIP(4) =  2
      MIP(5) =  2
      MIP(6) =  IOP
      MIP(7) =  0
      CALL SPRLAN (A,B,TOL,MPARA,MTREEA,MSIFA,MPARB,MTREEB,MSIFB,
     &             MIP,MOP,MSING,EVL,RWORK,V,RWORK(1+MAXLAN),
     &             RWORK(1+MAXLAN+N+N),IWORK,
     &             SHIFT,N,NEVAL,NEVAL,MAXLAN,LPU,IPSW,IERR,
     &             LMVY18,LMMY18,CMVY18,CMMY18)
C
      RETURN
      END
