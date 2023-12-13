C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE INSUB (B,S,MB,NB,MS,NS,I,J,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  INSUB                   GROUP 3 / PUBLIC
C
C     INSUB  INSERTS A SUBMATRIX (S) INTO A SPECIFIED POSITION OF A
C     RECTANGULAR MATRIX (B) OR IT ADDS THE SUBMATRIX TO THE CURRENT
C     CONTENT :
C               B(I,J) = S(1,1)             ETC., FOR  IFLAG=1
C               B(I,J) = B(I,J) + S(1,1)    ETC., FOR  IFLAG=2
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           I,IFLAG,J,MB,MS,NB,NS
      DOUBLE PRECISION  B(*),S(*)
C
      INTEGER           IPB,IPS,K,L
C ----------------------------------------------------------------------
      IF (MB.LT.1.OR.NB.LT.1)  GO TO 100
      IF (MS.LT.1.OR.NS.LT.1)  GO TO 100
      IF (I .LT.1.OR.J .LT.1)  GO TO 100
      IF ((I-1+MS).GT.MB)      GO TO 100
      IF ((J-1+NS).GT.NB)      GO TO 100
C
      IPB = MB*(J-1) + I - 1
      IPS = 0
      DO 50 K=1,NS
         DO 25 L=1,MS
            IPB = IPB+1
            IPS = IPS+1
            IF (IFLAG.EQ.1)  B(IPB) = S(IPS)
            IF (IFLAG.EQ.2)  B(IPB) = B(IPB) + S(IPS)
   25    CONTINUE
         IPB = IPB+MB-MS
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE EXTSUB (A,S,MA,NA,MS,NS,I,J)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  EXTSUB                  GROUP 3 / PUBLIC
C
C     EXTSUB  EXTRACTS A SUBMATRIX (S) FROM A SPECIFIED POSITION IN
C     A RECTANGULAR MATRIX (A)   -   S(1,1) = A(I,J)
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-20 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           I,J,MA,MS,NA,NS
      DOUBLE PRECISION  A(*),S(*)
C
      INTEGER           IPA,IPS,K,L
C ----------------------------------------------------------------------
      IF (MA.LT.1.OR.NA.LT.1)  GO TO 100
      IF (MS.LT.1.OR.NS.LT.1)  GO TO 100
      IF (I .LT.1.OR.J .LT.1)  GO TO 100
      IF ((I-1+MS).GT.MA)      GO TO 100
      IF ((J-1+NS).GT.NA)      GO TO 100
C
      IPA = MA*(J-1) + I - 1
      IPS = 0
      DO 50 K=1,NS
         DO 25 L=1,MS
            IPA    = IPA+1
            IPS    = IPS+1
            S(IPS) = A(IPA)
   25    CONTINUE
         IPA = IPA+MA-MS
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE MAB (A,B,C,L,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MAB                     GROUP 3 / PUBLIC
C
C     MAB  MULTIPLIES TWO RECTANGULAR MATRICES  A  AND  B.
C     THE RESULTING MATRIX IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0)
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,L,M,N
      DOUBLE PRECISION  A(L,M),B(M,N),C(L,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (L.LT.1.OR.M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)  GO TO 200
      IF (IFLAG.GT.0)  GO TO 400
C                                               **  C = A*B
      DO 100 J=1,N
         DO 50 I=1,L
            C(I,J) = PRACC(A(I,1),B(1,J),L,1,M)
   50    CONTINUE
  100 CONTINUE
      GO TO 1000
C                                               **  C = C - A*B
  200 DO 300 J=1,N
         DO 250 I=1,L
            C(I,J) = C(I,J) - PRACC(A(I,1),B(1,J),L,1,M)
  250    CONTINUE
  300 CONTINUE
      GO TO 1000
C                                               **  C = C + A*B
  400 DO 500 J=1,N
         DO 450 I=1,L
            C(I,J) = C(I,J) + PRACC(A(I,1),B(1,J),L,1,M)
  450    CONTINUE
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE BIGMUL (A,B,C,WA,L,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  BIGMUL                  GROUP 3 / PUBLIC
C
C     T A S K :  TO MULTIPLY TWO RECTANGULAR MATRICES A AND B.
C                THE RESULTING MATRIX  A*B  IS STORED IN C (IFLAG=0),
C                ADDED TO THE (INPUT) CONTENT OF C  (IFLAG.GT.0) OR
C                SUBTRACTED FROM THE CONTENT OF C (IFLAG.LT.0)
C     BIGMUL PERFORMS THE SAME OPERATION AS SUBROUTINE MAB, THE ONLY
C     DIFFERENCE BEING THE SCRATCH ARRAY  WA
C
C     ROUTINES CALLED/REFERENCED :  PRACC  (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-09-27 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,L,M,N
      DOUBLE PRECISION  A(L,M),B(M,N),C(L,N),WA(M)
C
      INTEGER           I,J,K
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (L.LT.1.OR.M.LT.1.OR.N.LT.1)  GO TO 100
C
      DO 50 I=1,L
         DO 10 K=1,M
            WA(K) = A(I,K)
   10    CONTINUE
         IF (IFLAG.EQ.0) THEN
C                                               ** C = A*B
            DO 20 J=1,N
               C(I,J) = PRACC(WA,B(1,J),1,1,M)
   20       CONTINUE
         ELSEIF (IFLAG.LT.0) THEN
C                                               ** C = C - A*B
            DO 30 J=1,N
               C(I,J) = C(I,J) - PRACC(WA,B(1,J),1,1,M)
   30       CONTINUE
         ELSE
C                                               ** C = C + A*B
            DO 40 J=1,N
               C(I,J) = C(I,J) + PRACC(WA,B(1,J),1,1,M)
   40       CONTINUE
         ENDIF
   50 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE MATB (A,B,C,L,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MATB                    GROUP 3 / PUBLIC
C
C     MATB  MULTIPLIES TWO RECTANGULAR MATRICES,  A-TRANSPOSE  AND  B.
C     THE RESULTING MATRIX IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0)
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,L,M,N
      DOUBLE PRECISION  A(M,L),B(M,N),C(L,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (L.LT.1.OR.M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)  GO TO 200
      IF (IFLAG.GT.0)  GO TO 400
C                                               **  C = AT*B
      DO 100 J=1,N
         DO 50 I=1,L
            C(I,J) = PRACC(A(1,I),B(1,J),1,1,M)
   50    CONTINUE
  100 CONTINUE
      GO TO 1000
C                                               **  C = C - AT*B
  200 DO 300 J=1,N
         DO 250 I=1,L
            C(I,J) = C(I,J) - PRACC(A(1,I),B(1,J),1,1,M)
  250    CONTINUE
  300 CONTINUE
      GO TO 1000
C                                               **  C = C + AT*B
  400 DO 500 J=1,N
         DO 450 I=1,L
            C(I,J) = C(I,J) + PRACC(A(1,I),B(1,J),1,1,M)
  450    CONTINUE
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE MABT (A,B,C,L,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MABT                    GROUP 3 / PUBLIC
C
C     MABT  MULTIPLIES TWO RECTANGULAR MATRICES,  A  AND  B-TRANSPOSE.
C     THE RESULTING MATRIX IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0)
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,L,M,N
      DOUBLE PRECISION  A(L,M),B(N,M),C(L,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (L.LT.1.OR.M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)  GO TO 200
      IF (IFLAG.GT.0)  GO TO 400
C                                               **  C = A*BT
      DO 100 J=1,N
         DO 50 I=1,L
            C(I,J) = PRACC(A(I,1),B(J,1),L,N,M)
   50    CONTINUE
  100 CONTINUE
      GO TO 1000
C                                               **  C = C - A*BT
  200 DO 300 J=1,N
         DO 250 I=1,L
            C(I,J) = C(I,J) - PRACC(A(I,1),B(J,1),L,N,M)
  250    CONTINUE
  300 CONTINUE
      GO TO 1000
C                                               **  C = C + A*BT
  400 DO 500 J=1,N
         DO 450 I=1,L
            C(I,J) = C(I,J) + PRACC(A(I,1),B(J,1),L,N,M)
  450    CONTINUE
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE MATBT (A,B,C,L,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MATBT                   GROUP 3 / PUBLIC
C
C     MATBT  MULTIPLIES TWO RECTANGULAR MATRICES,  A-TRANSPOSE  AND
C     B-TRANSPOSE
C     THE RESULTING MATRIX IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0)
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,L,M,N
      DOUBLE PRECISION  A(M,L),B(N,M),C(L,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (L.LT.1.OR.M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)  GO TO 200
      IF (IFLAG.GT.0)  GO TO 400
C                                               **  C = AT*BT
      DO 100 J=1,N
         DO 50 I=1,L
            C(I,J) = PRACC(A(1,I),B(J,1),1,N,M)
   50    CONTINUE
  100 CONTINUE
      GO TO 1000
C                                               **  C = C - AT*BT
  200 DO 300 J=1,N
         DO 250 I=1,L
            C(I,J) = C(I,J) - PRACC(A(1,I),B(J,1),1,N,M)
  250    CONTINUE
  300 CONTINUE
      GO TO 1000
C                                               **  C = C + AT*BT
  400 DO 500 J=1,N
         DO 450 I=1,L
            C(I,J) = C(I,J) + PRACC(A(1,I),B(J,1),1,N,M)
  450    CONTINUE
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE MATA (A,C,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MATA                    GROUP 3 / PUBLIC
C
C     MATA  MULTIPLIES THE TRANSPOSE OF A MATRIX  WITH THE MATRIX IT-
C     SELF, I.E.,  AT*A
C     THE RESULTING MATRIX IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0)
C     IN ALL THREE CASES  C  RETURNS AS A SYMMETRIC MATRIX
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,M,N
      DOUBLE PRECISION  A(M,N),C(N,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)        GO TO 200
      IF (IFLAG.GT.0)        GO TO 400
C                                               **  C = AT*A
      DO 100 J=1,N
         DO 50 I=J,N
            C(I,J) = PRACC(A(1,I),A(1,J),1,1,M)
            C(J,I) = C(I,J)
   50    CONTINUE
  100 CONTINUE
      GO TO 1000
C                                               **  C = C - AT*A
  200 DO 300 J=1,N
         DO 250 I=J,N
            C(I,J) = C(I,J) - PRACC(A(1,I),A(1,J),1,1,M)
            C(J,I) = C(I,J)
  250    CONTINUE
  300 CONTINUE
      GO TO 1000
C                                               **  C = C + AT*A
  400 DO 500 J=1,N
         DO 450 I=J,N
            C(I,J) = C(I,J) + PRACC(A(1,I),A(1,J),1,1,M)
            C(J,I) = C(I,J)
  450    CONTINUE
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE MAAT (A,C,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MAAT                    GROUP 3 / PUBLIC
C
C     MAAT  MULTIPLIES A MATRIX  A  WITH ITS OWN TRANSPOSE, I.E.,
C     A*AT
C     THE RESULTING MATRIX IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0)
C     IN ALL THREE CASES  C  RETURNS AS A SYMMETRIC MATRIX
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,M,N
      DOUBLE PRECISION  A(M,N),C(M,M)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)        GO TO 200
      IF (IFLAG.GT.0)        GO TO 400
C                                               **  C = A*AT
      DO 100 J=1,M
         DO 50 I=J,M
            C(I,J) = PRACC(A(I,1),A(J,1),M,M,N)
            C(J,I) = C(I,J)
   50    CONTINUE
  100 CONTINUE
      GO TO 1000
C                                               **  C = C - A*AT
  200 DO 300 J=1,M
         DO 250 I=J,M
            C(I,J) = C(I,J) - PRACC(A(I,1),A(J,1),M,M,N)
            C(J,I) = C(I,J)
  250    CONTINUE
  300 CONTINUE
      GO TO 1000
C                                               **  C = C + A*AT
  400 DO 500 J=1,M
         DO 450 I=J,M
            C(I,J) = C(I,J) + PRACC(A(I,1),A(J,1),M,M,N)
            C(J,I) = C(I,J)
  450    CONTINUE
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE SMATB (A,B,C,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SMATB                   GROUP 3 / PUBLIC
C
C     SMATB  MULTIPLIES TWO RECTANGULAR MATRICES,  A-TRANSPOSE  AND  B,
C     ASSUMING THE PRODUCT TO BE A SYMMETRIC MATRIX.
C     THE PRODUCT  AT*B  IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0).
C     IN ALL CASES  C  IS A SYMMETRIC MATRIX.
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-08-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,M,N
      DOUBLE PRECISION  A(M,N),B(M,N),C(N,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)        GO TO 200
      IF (IFLAG.GT.0)        GO TO 400
C                                               **  C = AT*B
      DO 100 J=1,N
         DO 50 I=J,N
            C(I,J) = PRACC(A(1,I),B(1,J),1,1,M)
   50    CONTINUE
  100 CONTINUE
      GO TO 700
C                                               **  C = C - AT*B
  200 DO 300 J=1,N
         DO 250 I=J,N
            C(I,J) = C(I,J) - PRACC(A(1,I),B(1,J),1,1,M)
  250    CONTINUE
  300 CONTINUE
      GO TO 700
C                                               **  C = C + AT*B
  400 DO 500 J=1,N
         DO 450 I=J,N
            C(I,J) = C(I,J) + PRACC(A(1,I),B(1,J),1,1,M)
  450    CONTINUE
  500 CONTINUE
C                                               **  COMPLETE  C
  700 DO 800 J=1,N
         DO 750 I=1,J
            C(I,J) = C(J,I)
  750    CONTINUE
  800 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE SMABT (A,B,C,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SMABT                   GROUP 3 / PUBLIC
C
C     SMABT  MULTIPLIES TWO RECTANGULAR MATRICES,  A  AND  B-TRANSPOSE,
C     ASSUMING THE PRODUCT IS A SYMMETRIC MATRIX.
C     THE PRODUCT  A*BT  IS STORED IN  C  (IFLAG.EQ.0), ADDED TO THE
C     (INPUT) CONTENT OF  C  (IFLAG.GT.0)  OR SUBTRACTED FROM THE
C     (INPUT) CONTENT OF  C  (IFLAG.LT.0).
C     IN ALL CASES  C  IS A SYMMETRIC MATRIX.
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   87-08-29 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,M,N
      DOUBLE PRECISION  A(N,M),B(N,M),C(N,N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)  GO TO 1000
C
      IF (IFLAG.LT.0)        GO TO 200
      IF (IFLAG.GT.0)        GO TO 400
C                                               **  C = A*BT
      DO 100 J=1,N
         DO 50 I=J,N
            C(I,J) = PRACC(A(I,1),B(J,1),N,N,M)
   50    CONTINUE
  100 CONTINUE
      GO TO 700
C                                               **  C = C - A*BT
  200 DO 300 J=1,N
         DO 250 I=J,N
            C(I,J) = C(I,J) - PRACC(A(I,1),B(J,1),N,N,M)
  250    CONTINUE
  300 CONTINUE
      GO TO 700
C                                               **  C = C + A*BT
  400 DO 500 J=1,N
         DO 450 I=J,N
            C(I,J) = C(I,J) + PRACC(A(I,1),B(J,1),N,N,M)
  450    CONTINUE
  500 CONTINUE
C                                               ** COMPLETE  C
  700 DO 800 J=1,N
         DO 750 I=1,J
            C(I,J) = C(J,I)
  750    CONTINUE
  800 CONTINUE
 1000 RETURN
      END
      SUBROUTINE PREMUL (B,C,WA,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PREMUL                  GROUP 3 / PUBLIC
C
C     PREMUL  PREMULTIPLIES MATRIX  C  BY A SQUARE MATRIX  B  OR ITS
C     TRANSPOSE  -  THE RESULT MAY (SUBJECT TO IFLAG) BE  ADDED TO /
C     SUBTRACTED FROM THE (INPUT) CONTENT OF  C
C
C     ROUTINES CALLED/REFERENCED :  PRACC      (SAM-0)
C                                   ABS        (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-22 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,M,N
      DOUBLE PRECISION  B(M,M),C(M,N),WA(M)
C
      INTEGER           I,IOP,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)         GO TO 1000
      IOP = ABS(IFLAG)
      IF (IOP.LT.1.OR.IOP.GT.4)     GO TO 1000
C
      DO 500 J=1,N
         IF (IOP.EQ.2.OR.IOP.EQ.4)  GO TO 150
         DO 100 I=1,M
            WA(I) = PRACC(B(I,1),C(1,J),M,1,M)
  100    CONTINUE
         GO TO 250
  150    DO 200 I=1,M
            WA(I) = PRACC(B(1,I),C(1,J),1,1,M)
  200    CONTINUE
C
  250    IF (IOP.GT.2)              GO TO 400
C
         IF (IFLAG.LT.0)            GO TO 350
         DO 325 I=1,M
            C(I,J) = WA(I)
  325    CONTINUE
         GO TO 500
  350    DO 375 I=1,M
            C(I,J) =-WA(I)
  375    CONTINUE
         GO TO 500
C
  400    IF (IFLAG.LT.0)            GO TO 450
         DO 425 I=1,M
            C(I,J) = C(I,J) + WA(I)
  425    CONTINUE
         GO TO 500
  450    DO 475 I=1,M
            C(I,J) = C(I,J) - WA(I)
  475    CONTINUE
C
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE PSTMUL (B,C,WA,M,N,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PSTMUL                  GROUP 3 / PUBLIC
C
C     PSTMUL  POSTMULTIPLIES MATRIX  C  BY A SQUARE MATRIX  B  OR ITS
C     TRANSPOSE  -  THE RESULT MAY (SUBJECT TO IFLAG) BE ADDED TO /
C     SUBTRACTED FROM THE (INPUT) CONTENT OF  C
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C                                   ABS       (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-22 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IFLAG,M,N
      DOUBLE PRECISION  B(N,N),C(M,N),WA(N)
C
      INTEGER           I,IOP,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)         GO TO 1000
      IOP = ABS(IFLAG)
      IF (IOP.LT.1.OR.IOP.GT.4)     GO TO 1000
C
      DO 500 I=1,M
         IF (IOP.EQ.2.OR.IOP.EQ.4)  GO TO 150
         DO 100 J=1,N
            WA(J) = PRACC(C(I,1),B(1,J),M,1,N)
  100    CONTINUE
         GO TO 250
  150    DO 200 J=1,N
            WA(J) = PRACC(C(I,1),B(J,1),M,N,N)
  200    CONTINUE
C
  250    IF (IOP.GT.2)              GO TO 400
C
         IF (IFLAG.LT.0)            GO TO 350
         DO 325 J=1,N
            C(I,J) = WA(J)
  325    CONTINUE
         GO TO 500
  350    DO 375 J=1,N
            C(I,J) =-WA(J)
  375    CONTINUE
         GO TO 500
C
  400    IF (IFLAG.LT.0)            GO TO 450
         DO 425 J=1,N
            C(I,J) = C(I,J) + WA(J)
  425    CONTINUE
         GO TO 500
  450    DO 475 J=1,N
            C(I,J) = C(I,J) - WA(J)
  475    CONTINUE
C
  500 CONTINUE
C
 1000 RETURN
      END
      SUBROUTINE PRESKY (A,B,C,WA,MSKY,
     +                   M,N,KSA,IFLAG,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  PRESKY                  GROUP 3 / PUBLIC
C
C     PRESKY  PERFORMS ONE OF THE FOLLOWING MATRIX MULTIPLICATIONS :
C             C = AB      (IFLAG = 1)
C             C =-AB      (IFLAG =-1)
C             C = C + AB  (IFLAG = 2)
C             C = C - AB  (IFLAG =-2)
C     A  IS A SYMMETRIC 'SKYLINE' MATRIX (KSA.GT.0) OR A DIAGONAL
C     MATRIX (KSA.LE.0),  WHEREAS  B  AND  C  ARE RECTANGULAR MATRICES
C
C     ROUTINES CALLED/REFERENCED :  PRACC       (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C                       86-06-04 / 1.1    K.BELL
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,KSA,LPU,M,N,     MSKY(*)
      DOUBLE PRECISION  A(*),B(M,N),C(M,N),WA(M)
C
      INTEGER           I,IE,II,IP,IS,J,NC
      DOUBLE PRECISION  BIJ
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1.OR.M.LT.1)  GO TO 920
      IF (KSA.LT.1)          GO TO 100
      IF (M.EQ.1)            GO TO 100
C                                               ** CHECK MSKY
      DO 50 I=2,M
         IF (MSKY(I).LT.I)   GO TO 910
         NC = MSKY(I)-MSKY(I-1)
         IF (NC.GT.I)        GO TO 910
   50 CONTINUE
C ------------------------------------------------ MULTIPLICATION
  100 DO 800 J=1,N
         IF (M.EQ.1)         GO TO 120
         IF (KSA.GT.0)       GO TO 200
C                                               ** A  IS DIAGONAL
  120    DO 140 I=1,M
            WA(I) = A(I)*B(I,J)
  140    CONTINUE
         GO TO 350
C                                               ** A  IS A SYMMETRIC
C                                                  SKYLINE MATRIX
  200    WA(1) = A(1)*B(1,J)
         DO 300 I=2,M
            NC    = MSKY(I)-MSKY(I-1)
            IP    = MSKY(I-1)+1
            IS    = I-NC+1
            WA(I) = PRACC(A(IP),B(IS,J),1,1,NC)
            IF (NC.EQ.1)     GO TO 300
            IE    = I-1
            BIJ   = B(I,J)
            DO 250 II=IS,IE
               WA(II) = WA(II) + A(IP)*BIJ
               IP     = IP+1
  250       CONTINUE
  300    CONTINUE
C
  350    IF (IFLAG.EQ.1)              GO TO 650
         IF (IFLAG.EQ.(-1))           GO TO 550
         IF (IFLAG.EQ.2)              GO TO 450
         IF (IFLAG.EQ.(-2))           GO TO 375
         IF (IFLAG.EQ.88.AND.N.EQ.1)  GO TO 1000
         GO TO 920
C                                               **  C = C - A*B
  375    DO 400 I=1,M
            C(I,J) = C(I,J)-WA(I)
  400    CONTINUE
         GO TO 800
C                                               **  C = C + A*B
  450    DO 500 I=1,M
            C(I,J) = C(I,J)+WA(I)
  500    CONTINUE
         GO TO 800
C                                               **  C =-A*B
  550    DO 600 I=1,M
            C(I,J) =-WA(I)
  600    CONTINUE
         GO TO 800
C                                               **  C = A*B
  650    DO 700 I=1,M
            C(I,J) = WA(I)
  700    CONTINUE
C
  800 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
  920 IERR = IERR-1
      IF (LPU.LT.1)      GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))  WRITE (LPU,6910) M,N,IFLAG
      IF (IERR.EQ.(-2))  WRITE (LPU,6920) I,MSKY(I-1),MSKY(I)
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE PRESKY')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  M, N, IFLAG =',3I8)
 6920 FORMAT(/5X,'MSKY IN ERROR :  I, MSKY(I-1), MSKY(I) =',3I8)
C
      END
      SUBROUTINE TSKYMB (U,B,C,WA,MSKY,
     +                   M,N,KSU,IFLAG,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE : TSKYMB                   GROUP 3 / PUBLIC
C
C     T A S K : TO PREMULTIPLY A MATRIX -B- BY AN UPPER TRIANGULAR (OR
C               DIAGONAL) MATRIX -U- OR ITS TRANSPOSE :
C               C = U*B  (IFLAG = 1)  /  WA = U*B  (IFLAG = 11 AND N=1)
C               C = UT*B (IFLAG = 2)  /  WA = UT*B (IFLAG = 22 AND N=1)
C     U  IS A SKYLINE MATRIX (KSU.GT.0) OR A DIAGONAL MATRIX (KSU.LE.0)
C     C AND B ARE RECTANGULAR (M*N) MATRICES  (THEIR ACTUAL ARGUMENTS
C     MAY BE THE SAME)
C
C     ROUTINES CALLED/REFERENCED :   PRACC AND RMINT      (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-05-25/ 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,KSU,LPU,M,N,      MSKY(M)
      DOUBLE PRECISION  B(*),C(*),U(*),WA(M)
C
      INTEGER           I,II,IP,IS,J,JP,KODE,NC
      DOUBLE PRECISION  BIJ,ZERO,       PRACC
C
      PARAMETER         ( ZERO = 0.0D0 )
C
      EXTERNAL          PRACC,RMINT
C ----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1.OR.M.LT.1)     GO TO 920
      IF (KSU.GT.0.AND.M.GT.1)  THEN
C                                               ** CHECK MSKY
         DO 50 I=2,M
            IF (MSKY(I).LT.I)   GO TO 910
            NC = MSKY(I)-MSKY(I-1)
            IF (NC.GT.I)        GO TO 910
   50    CONTINUE
      ENDIF
C
      KODE = 0
      IF (IFLAG.EQ.1.OR.IFLAG.EQ.11)  KODE = 1
      IF (IFLAG.EQ.2.OR.IFLAG.EQ.22)  KODE = 2
      IF (KODE.EQ.0)                  GO TO 920
C ------------------------------------------------ MULTIPLICATION
      DO 500 J=1,N
         JP = (J-1)*M
         IF (KSU.LE.0.OR.M.EQ.1) THEN
C                                               ** U  IS DIAGONAL
            IF (IFLAG.GT.10.AND.N.EQ.1) THEN
               DO 75 I=1,M
                  JP = JP+1
                  WA(JP) = U(I)*B(JP)
   75          CONTINUE
            ELSE
               DO 100 I=1,M
                  JP    = JP+1
                  C(JP) = U(I)*B(JP)
  100          CONTINUE
            ENDIF
         ELSE
C                                               ** U  IS TRIANGULAR
C
            IF (KODE.EQ.1)  CALL RMINT (WA,M,1,ZERO)
            WA(1) = U(1)*B(JP+1)
            DO 300 I=2,M
               NC = MSKY(I)-MSKY(I-1)
               IP = MSKY(I-1) + 1
               IS = I-NC+1
               IF (KODE.EQ.2) THEN
C                                               ** UT*B
C
                  WA(I) = PRACC(U(IP),B(JP+IS),1,1,NC)
               ELSE
C                                               ** U*B
                  BIJ = B(JP+I)
                  DO 200 II=IS,I
                     WA(II) = WA(II) + U(IP)*BIJ
                     IP     = IP+1
  200             CONTINUE
               ENDIF
  300       CONTINUE
            IF (IFLAG.GT.10.AND.N.EQ.1)  GO TO 1000
            DO 400 I=1,M
               JP    = JP+1
               C(JP) = WA(I)
  400       CONTINUE
         ENDIF
  500 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
  920 IERR =IERR-1
      IF (LPU.LT.1)     GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1)) WRITE (LPU,6910) M,N,IFLAG
      IF (IERR.EQ.(-2)) WRITE (LPU,6920) I,MSKY(I-1),MSKY(I)
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE TSKYMB')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  M, N, IFLAG =',3I8)
 6920 FORMAT(/5X,'MSKY IN ERROR :  I, MSKY(I-1), MSKY(I) =',3I8)
C
      END
      SUBROUTINE CONTRA (B,C,WA,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  CONTRA                  GROUP 3 / PUBLIC
C
C     CONTRA  PERFORMES THE CONGRUENCE TRANSFORMATION
C          C = BT*C*B   (BT = TRANSPOSE OF B)
C     WHERE  C  IS A SYMMETRIC MATRIX AND  B  IS A SQUARE MATRIX
C
C     ROUTINES CALLED/REFERENCED :  PRACC    (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C
      INTEGER           N
      DOUBLE PRECISION  B(N,N),C(N,N),WA(N)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (N.LT.1)       GO TO 100
C
      DO 30 J=1,N
         DO 10 I=1,N
            WA(I) = PRACC(B(1,I),C(1,J),1,1,N)
   10    CONTINUE
         DO 20 I=1,N
            C(I,J) = WA(I)
   20    CONTINUE
   30 CONTINUE
C
      DO 60 I=1,N
         DO 40 J=I,N
            WA(J) = PRACC(C(I,1),B(1,J),N,1,N)
   40    CONTINUE
         DO 50 J=I,N
            C(I,J) = WA(J)
   50    CONTINUE
   60 CONTINUE
C                                               ** RESTORE SYMMETRY
      DO 80 I=1,N
         DO 70 J=I,N
            C(J,I) = C(I,J)
   70    CONTINUE
   80 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE SYMTRA (A,B,C,WA,M,N)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SYMTRA                  GROUP 3 / PUBLIC
C
C     SYMTRA  PERFORMS THE MATRIX MULTIPLICATION
C          C = BT*A*B   (BT = TRANSPOSE OF B)
C     WHERE  A  IS A SYMMRTIC MATRIX, AND  B  IS A RECTANGULAR MATRIX
C
C     ROUTINES CALLED/REFERENCED :  PRACC    (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           M,N
      DOUBLE PRECISION  A(M,M),B(M,N),C(N,N),WA(M)
C
      INTEGER           I,J
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.1)  GO TO 100
C
      DO 75 I=1,N
         DO 25 J=1,M
            WA(J) = PRACC(B(1,I),A(1,J),1,1,M)
   25    CONTINUE
         DO 50 J=I,N
            C(I,J) = PRACC(WA(1),B(1,J),1,1,M)
            C(J,I) = C(I,J)
   50    CONTINUE
   75 CONTINUE
C
  100 RETURN
      END
      SUBROUTINE SKYTRA (A,B,C,WA,MSKY,M,N,KSA,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SKYTRA                  GROUP 3 / PUBLIC
C
C     SKYTRA  PERFORMS THE MATRIX MULTIPLICATION
C          C = BT*A*B    (BT = TRANSPOSE OF B)
C     WHERE  A  IS A SYMMETRIC 'SKYLINE' MATRIX (KSA.GT.0) OR A
C     DIAGONAL MATRIX (KSA.LE.0), AND  B  IS A RECTANGULAR MATRIX
C
C     ROUTINES CALLED/REFERENCED :  PRACC     (SAM-0)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-23 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KSA,LPU,M,N,    MSKY(*)
      DOUBLE PRECISION  A(*),B(M,N),C(N,N),WA(M)
C
      INTEGER           I,IE,II,IP,IS,J,NC
      DOUBLE PRECISION  BIJ
      DOUBLE PRECISION  PRACC
C
      EXTERNAL          PRACC
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (N.LT.1.OR.M.LT.1)  GO TO 920
      IF (KSA.LT.1)          GO TO 100
      IF (M.EQ.1)            GO TO 100
C                                               ** CHECK  MSKY
      DO 50 I=2,M
         IF (MSKY(I).LT.I)   GO TO 910
         NC   = MSKY(I)-MSKY(I-1)
         IF (NC.GT.I)        GO TO 910
   50 CONTINUE
C ------------------------------------------------ MULTIPLICATION
  100 DO 500 J=1,N
         IF (M.EQ.1)         GO TO 120
         IF (KSA.GT.0)       GO TO 200
C                                               ** A  IS DIAGONAL
  120    DO 140 I=1,M
            WA(I) = A(I)*B(I,J)
  140    CONTINUE
         GO TO 350
C                                               ** A  IS A SYMMETRIC
C                                                  SKYLINE MATRIX
  200    WA(1) = A(1)*B(1,J)
         DO 300 I=2,M
            NC    = MSKY(I)-MSKY(I-1)
            IP    = MSKY(I-1)+1
            IS    = I-NC+1
            WA(I) = PRACC(A(IP),B(IS,J),1,1,NC)
            IF (NC.EQ.1)     GO TO 300
            IE    = I-1
            BIJ   = B(I,J)
            DO 250 II=IS,IE
               WA(II) = WA(II) + A(IP)*BIJ
               IP     = IP+1
  250       CONTINUE
  300    CONTINUE
C
  350    DO 400 I=J,N
            C(I,J) = PRACC(B(1,I),WA(1),1,1,M)
            C(J,I) = C(I,J)
  400    CONTINUE
  500 CONTINUE
      GO TO 1000
C ----------------------------------------------------- ERROR EXIT
  910 IERR =-1
  920 IERR = IERR-1
      IF (LPU.LT.1)      GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))  WRITE (LPU,6910)  M,N
      IF (IERR.EQ.(-2))  WRITE (LPU,6920)  I,MSKY(I-1),MSKY(I)
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE SKYTRA')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  M, N =',2I8)
 6920 FORMAT(/5X,'MSKY IN ERROR :  I, MSKY(I-1), MSKY(I) =',3I8)
C
      END
      SUBROUTINE VECTRA (T,V,WA,M,N,K,IFLAG,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  VECTRA                  GROUP 3 / PUBLIC
C
C     VECTRA  PREMULTIPLIES A VECTOR  V  BY A TRANSFORMATION MATRIX
C     (IFLAG.EQ.1) OR ITS TRANSPOSE (IFLAG.EQ.2).
C     THE TRANSFORMATION MATRIX IS A UNIT MATRIX WITH AN  N BY N  SUB-
C     MATRIX  T  ON THE DIAGONAL.
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-25 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,K,LPU,M,N
      DOUBLE PRECISION  T(N,N),V(M),WA(N)
C
      INTEGER           I,II,J
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (M.LT.1.OR.N.LT.1)  GO TO 90
      IF (K.LT.1)            GO TO 90
      IF ((K+N-1).GT.M)      GO TO 90
C
      IF (IFLAG.EQ.1)        GO TO 10
      IF (IFLAG.EQ.2)        GO TO 40
      GO TO 100
C                                               **  V = T*V
   10 DO 30 I=1,N
         II    = K
         WA(I) = ZERO
         DO 20 J=1,N
            WA(I) = WA(I) + T(I,J)*V(II)
            II    = II+1
   20    CONTINUE
   30 CONTINUE
      GO TO 70
C                                               **  V = TT*V
   40 DO 60 I=1,N
         II = K
         WA(I) = ZERO
         DO 50 J=1,N
            WA(I) = WA(I) + T(J,I)*V(II)
            II    = II+1
   50    CONTINUE
   60 CONTINUE
C
   70 II = K
      DO 80 I=1,N
         V(II) = WA(I)
         II    = II+1
   80 CONTINUE
      GO TO 100
C ----------------------------------------------- ERROR EXIT
   90 IERR =-1
      IF (LPU.GT.0)     WRITE (LPU,690) M,N,K
C ----------------------------------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE VECTRA'
     +     //5X,'ILLEGAL MATRIX DIMENSION(S) :  M, N  =',2I8,
     +     / 5X,'AND/OR SUBMATRIX POSITION            =',I8  )
C
      END
      SUBROUTINE MATTRA (T,A,WA,M,N,K,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  MATTRA                  GROUP 3 / PUBLIC
C
C     MATTRA  PERFORMS THE MATRIX MULTIPLICATION
C          A = TT*A*T   (TT = TRANSPOSE OF T)
C     WHERE  A  IS A FULL, SYMMETRIC MATRIX AND  T  IS A UNIT MATRIX
C     WITH AN N BY N SUBMATRIX (ARGUMENT T) ON THE DIAGONAL
C
C     ROUTINES CALLED/REFERENCED :  NONE
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-25 / 1.0
C                       86-10-09 / 1.1    K.BELL
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,K,LPU,M,N
      DOUBLE PRECISION  A(M,M),T(N,N),WA(N)
C
      INTEGER           I,II,J,JJ,JS,KN,L
      DOUBLE PRECISION  ZERO
C
      PARAMETER         ( ZERO = 0.0D0 )
C ----------------------------------------------------------------------
      IERR = 0
C
      IF (M.LT.1.OR.N.LT.1)  GO TO 900
      IF (K.LT.1)            GO TO 900
      KN   = K+N-1
      IF (KN.GT.M)           GO TO 900
C
      DO 100 JJ=K,M
         DO 50 I=1,N
            II    = K
            WA(I) = ZERO
            DO 25 L=1,N
               WA(I) = WA(I) + T(L,I)*A(II,JJ)
               II    = II+1
   25       CONTINUE
   50    CONTINUE
         II = K
         DO 75 I=1,N
            A(II,JJ) = WA(I)
            IF (JJ.GT.KN)    A(JJ,II) = WA(I)
            II       = II+1
   75    CONTINUE
  100 CONTINUE
C
      DO 200 II=1,KN
         JS = 1
         IF (II.GT.K)  JS = II-K+1
         DO 150 J=JS,N
            JJ    = K
            WA(J) = ZERO
            DO 125 L=1,N
               WA(J) = WA(J) + A(II,JJ)*T(L,J)
               JJ    = JJ+1
  125       CONTINUE
  150    CONTINUE
         JJ = K
         IF (II.GT.K)  JJ=II
         DO 175 J=JS,N
            A(II,JJ) = WA(J)
            IF (II.LT.K)  A(JJ,II) = WA(J)
            JJ       = JJ+1
  175    CONTINUE
  200 CONTINUE
      DO 250 II=K,KN
         DO 225 JJ=II,KN
            A(JJ,II) = A(II,JJ)
  225    CONTINUE
  250 CONTINUE
C
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  900 IERR =-1
      IF (LPU.GT.0)     WRITE (LPU,6900) M,N,K
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE MATTRA'
     +      //5X,'ILLEGAL MATRIX DIMENSION(S) :  M, N =',2I8,
     +      / 5X,'AND/OR SUBMATRIX POSITION           =',I8  )
C
      END
      SUBROUTINE NORVEC (V,N,IFLAG,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NORVEC                  GROUP 3 / PUBLIC
C
C     NORVEC  NORMALIZES VECTOR  V  SUCH THAT ITS LENGTH BECOMES UNITY
C     (FOR IFLAG.GT.0) OR ITS NUMERICALLY LARGEST ELEMENT BECOMES UNITY
C     (FOR IFLAG.LE.0)
C
C     ROUTINES CALLED/REFERENCED :  PRACC         (SAM-0)
C                                   SCALE         (SAM-3)
C                                   ABS AND SQRT  (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-25 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,N
      DOUBLE PRECISION  V(N)
C
      INTEGER           I
      DOUBLE PRECISION  T,ONE,ZERO
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C
      EXTERNAL          PRACC,SCALE
C ----------------------------------------------------------------------
      IERR = 0
      IF (N.LT.1)       GO TO 90
      IF (IFLAG.GT.0)   GO TO 50
C                                               ** FIND LARGEST ELEMENT
      T = ZERO
      DO 25 I=1,N
         IF (ABS(V(I)) .GT. ABS(T))  T = V(I)
   25 CONTINUE
      GO TO 75
C                                               ** VECTOR LENGTH
   50 T = PRACC(V,V,1,1,N)
      T = SQRT(T)
C
   75 IF (T.EQ.ZERO)    GO TO 90
C                                               ** SCALE VECTOR
      T = ONE/T
      CALL SCALE (V,V,N,1,T)
      GO TO 100
C ------------------------------
   90 IERR =-1
C ------------------------------
  100 RETURN
      END
      SUBROUTINE NORSKY (A,V,WA,MSKY,N,KSA,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  NORSKY                  GROUP 3 / PUBLIC
C
C     NORSKY  NORMALIZES VECTOR  V  WITH RESPECT TO MATRIX  A  TO GIVE
C     VT*A*V = 1
C     MATRIX  A  IS A SYMMETRIC 'SKYLINE' MATRIX  (KSA.GT.0)  OR A
C     DIAGONAL MATRIX  (KSA.LE.0)
C
C     ROUTINES CALLED/REFERENCED :  SKYTRA AND SCALE  (SAM-3)
C                                   SQRT              (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-24 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,KSA,LPU,N,   MSKY(*)
      DOUBLE PRECISION  A(*),V(N),WA(N)
C
      DOUBLE PRECISION  S,T(1)
C
      EXTERNAL          SCALE,SKYTRA
C ----------------------------------------------------------------------
      IERR = 0
      CALL SKYTRA (A,V,T,WA,MSKY,N,1,KSA,LPU,IERR)
      IF (IERR.LT.0)    GO TO 90
      IF (T(1).LE.0.0D0)GO TO 90
C                                               ** SCALE VECTOR
      S = SQRT(T(1))
      CALL SCALE (V,V,N,1,S)
      GO TO 100
C ------------------------------------------------ ERROR EXIT
   90 IERR =-1
      IF (LPU.GT.0)     WRITE (LPU,690)
C ----------------------------------------------------------------------
  100 RETURN
C ----------------------------------------------------------------------
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE NORSKY'
     + //5X,'VECTOR  V  CANNOT BE NORMALIZED WITH RESPECT TO MATRIX  A')
C
      END
      SUBROUTINE ORTHO1 (G,V,TOL,N,NV,IFLAG,IPSW,LPU,KOUNT,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  ORTHO1                  GROUP 3 / PUBLIC
C
C     T A S K :  TO ORTHONORMALIZE A UNIT LENGTH VECTOR  V  AGAINST
C                THE NV COLUMNVECTORS OF MATRIX  G  WHICH ARE ALREADY
C                ORTHONORMAL TO ONE ANOTHER
C     A GRAM-SCHMIDT PROCEDURE, DEPENDING ON IFLAG, IS USED :
C     IFLAG = 1 :  V IS ORTHONORMALIZED AGAINST ALL VECTORS IN G
C     IFLAG = 2 :  V IS ORTHONORMALIZED ONLY AGAINST THOSE (KOUNT)
C                  VECTORS IN G WHOSE INNER-PRODUCT WITH V EXCEEDS TOL
C     PRINT :
C      IF IPSW.GT.5 :  ALL INNER-PRODUCTS, BEFORE PURGING, ARE PRINTED
C                      ON UNIT LPU
C      IF IPSW.GT.6 :  ALL INNER-PRODUCTS, AFTER PURGING, ARE PRINTED
C                      ON UNIT LPU
C
C     ROUTINES CALLED/REFERENCED :  PRACC, SCALE AND VASCA     (SAM-0)
C                                   ABS AND SQRT     (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   86-07-17 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,IFLAG,IPSW,KOUNT,LPU,N,NV
      DOUBLE PRECISION  TOL, G(N,*),V(N)
C
      INTEGER           J
      DOUBLE PRECISION  ONE,S,   PRACC
C
      PARAMETER         ( ONE = 1.0D0 )
C
      EXTERNAL          PRACC,SCALE,VASCA
C ----------------------------------------------------------------------
      IERR  = 0
      KOUNT = 0
      IF (NV.LT.1)      GO TO 100
C
      IF (IFLAG.EQ.1)   THEN
C                                               ** NO TESTING
         DO 20 J=1,NV
            S =-PRACC(V,G(1,J),1,1,N)
            CALL VASCA (V,G(1,J),1,1,S,N)
            IF (IPSW.GT.5) WRITE(LPU,610) J,S
            IF (IPSW.GT.6) THEN
               S =-PRACC(V,G(1,J),1,1,N)
               WRITE(LPU,620) J,S
            ENDIF
   20    CONTINUE
C
      ELSEIF (IFLAG.EQ.2) THEN
C                                               ** TESTING
         DO 40 J=1,NV
            S =-PRACC(V,G(1,J),1,1,N)
            IF (IPSW.GT.5) WRITE(LPU,610) J,S
            IF (ABS(S).GT.TOL) THEN
               CALL VASCA (V,G(1,J),1,1,S,N)
               KOUNT = KOUNT+1
               IF (IPSW.GT.6) WRITE(LPU,620) J,S
            ENDIF
   40    CONTINUE
C
      ELSE
C
         GO TO 100
C
      ENDIF
C
      S = PRACC(V,V,1,1,N)
      S = SQRT(S)
      IF (S.LT.TOL)  GO TO 90
      S = ONE/S
      CALL SCALE (V,V,N,1,S)
      GO TO 100
C ------------------------------------------------ ERROR EXIT
   90 IERR =-1
      IF (LPU.GT.0) WRITE(LPU,690)
C ------------------------------------------------
  100 RETURN
C ------------------------------------------------ FORMATS
C
  610 FORMAT(5X,'INNER-PRODUCT WITH VECTOR',I5,'  (BEFORE PURGING) = ',
     +       1PE11.3)
  620 FORMAT(5X,'INNER-PRODUCT WITH VECTOR',I5,'  (AFTER PURGING)  = ',
     +       1PE11.3)
  690 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE ORTHO1'
     +      /5X,'VECTOR IS LINEAR DEPENDENT')
C
      END
      SUBROUTINE GRAMS (G,EPS,N,NV,K,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  GRAMS                   GROUP 3 / PUBLIC
C
C     GRAMS  ORTHONORMALIZES THE COLUMN-VECTORS OF  G(N,NV),
C     BY A MODIFIED GRAM-SCHMIDT PROCEDURE,  ASSUMING THAT THE  K-1
C     FIRST COLUMNS ALREADY ARE ORTHONORMAL TO ONE ANOTHER
C
C     ROUTINES CALLED/REFERENCED :  PRACC                     (SAM-0)
C                                   NORVEC, SCADD AND SCALE   (SAM-3)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-24 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,K,LPU,N,NV
      DOUBLE PRECISION  EPS,   G(N,NV)
C
      INTEGER           J,JJ,JJM1
      DOUBLE PRECISION  S,ONE
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ONE = 1.0D0 )
C
      EXTERNAL          NORVEC,PRACC,SCADD,SCALE
C ----------------------------------------------------------------------
      IF (N.LT.1.OR.NV.LT.1)  GO TO 910
      IF (NV.GT.N)            GO TO 910
      IF (K.LT.1.OR.K.GT.NV)  GO TO 910
C
      IERR  = 0
C
C                                               ** NORMALIZE VECTORS
      DO 50 JJ=K,NV
         CALL NORVEC (G(1,JJ),N,1,IERR)
         IF (IERR.LT.0)       GO TO 920
   50 CONTINUE
C
      DO 200 JJ=K,NV
         IF (JJ.EQ.1)         GO TO 200
         JJM1 = JJ-1
         DO 100 J=1,JJM1
            S =-PRACC(G(1,J),G(1,JJ),1,1,N)
            CALL SCADD (G(1,JJ),G(1,J),G(1,JJ),N,1,S)
  100    CONTINUE
C
         S = PRACC(G(1,JJ),G(1,JJ),1,1,N)
         S = SQRT(S)
C                                               ** LINEAR DEP. ?
         IF (S.LT.EPS)        GO TO 930
C                                               ** NORMALIZE JJ'TH COL.
         S = ONE/S
         CALL SCALE (G(1,JJ),G(1,JJ),N,1,S)
  200 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      GO TO 950
  930 IERR =-3
C
  950 IF (LPU.LE.0)           GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))       WRITE (LPU,6910) N,NV,K
      IF (IERR.EQ.(-2))       WRITE (LPU,6920) JJ
      IF (IERR.EQ.(-3))       WRITE (LPU,6930) JJ
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE GRAMS')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) :  N, NV, K =',3I8)
 6920 FORMAT(/5X,'COLUMN NO.',I6,' OF  G  IS A NULL VECTOR')
 6930 FORMAT(/5X,'COLUMN NO.',I6,' OF  G  IS LINEAR DEPENDENT')
C
      END
      SUBROUTINE GORSKY (A,G,WA,MSKY,
     +                   EPS,N,NV,K,KSA,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  GORSKY                  GROUP 3 / PUBLIC
C
C     GORSKY  ORTHONORMALIZES THE COLUMN-VECTORS OF  G(N,NV) WITH
C     RESPECT TO MATRIX  A  TO GIVE    GT*A*G = I ,  BY USE OF A
C     MODIFIED GRAM-SCHMIDT PROCEDURE.
C     THE FIRST  K-1  VECTORS OF  G  ARE ASSUMED TO SATISFY
C     A-ORTHONORMALITY ALREADY.
C     A  IS A SYMMETRIC 'SKYLINE' MATRIX  (KSA.GT.0)  OR A DIAGONAL
C     MATRIX  (KSA.LE.0).
C
C     ROUTINES CALLED/REFERENCED :  PRACC                   (SAM-0)
C                                   NORSKY, PRESKY, SCADD,
C                                   SKYTRA AND SCALE        (SAM-3)
C                                   SQRT          (FORTRAN LIBRARY)
C
C     PROGRAMMED BY :   KOLBEIN BELL
C     DATE/VERSION  :   84-05-24 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      INTEGER           IERR,K,KSA,LPU,N,NV,    MSKY(*)
      DOUBLE PRECISION  EPS,   A(*),G(N,NV),WA(N)
C
      INTEGER           J,JJ,JJM1
      DOUBLE PRECISION  S,ONE,ZERO,   C(1)
      DOUBLE PRECISION  PRACC
C
      PARAMETER         ( ZERO = 0.0D0 , ONE = 1.0D0 )
C
      EXTERNAL          NORSKY,PRACC,PRESKY,SCADD,SCALE,SKYTRA
C ----------------------------------------------------------------------
      IF (N.LT.1.OR.NV.LT.1)  GO TO 910
      IF (NV.GT.N)            GO TO 910
      IF (K.LT.1.OR.K.GT.NV)  GO TO 910
C
      IERR  = 0
C                                               ** NORMALIZE VECTORS
C                                                  WITH RESPECT TO  A
      DO 100 JJ=K,NV
         CALL NORSKY (A,G(1,JJ),WA,MSKY,
     +                N,KSA,LPU,IERR    )
         IF (IERR.LT.0)       GO TO 920
  100 CONTINUE
C
      DO 200 JJ=K,NV
         IF (JJ.EQ.1)         GO TO 200
         JJM1 = JJ-1
         CALL PRESKY (A,G(1,JJ),WA,WA,MSKY,N,1,KSA,88,LPU,IERR)
         IF (IERR.LT.0)       GO TO 920
         DO 150 J=1,JJM1
            S =-PRACC(G(1,J),WA,1,1,N)
            CALL SCADD (G(1,JJ),G(1,J),G(1,JJ),N,1,S)
  150    CONTINUE
C
         CALL SKYTRA (A,G(1,JJ),C,WA,MSKY,N,1,KSA,LPU,IERR)
         IF (IERR.LT.0)       GO TO 930
C
         IF (C(1).LE.ZERO)    GO TO 930
         S = SQRT (C(1))
C                                               ** LINEAR DEP. ?
         IF (S.LT.EPS)        GO TO 940
C
         S = ONE/S
         CALL SCALE (G(1,JJ),G(1,JJ),N,1,S)
  200 CONTINUE
      GO TO 1000
C ------------------------------------------------ ERROR EXIT
  910 IERR =-1
      GO TO 950
  920 IERR =-2
      GO TO 950
  930 IERR =-3
      GO TO 950
  940 IERR =-4
C
  950 IF (LPU.LE.0)           GO TO 1000
      WRITE (LPU,6900)
      IF (IERR.EQ.(-1))       WRITE (LPU,6910) N,NV,K
      IF (IERR.EQ.(-2))       WRITE (LPU,6920) JJ
      IF (IERR.EQ.(-3))       WRITE (LPU,6930) JJ
      IF (IERR.EQ.(-4))       WRITE (LPU,6940) JJ
C ----------------------------------------------------------------------
 1000 RETURN
C ----------------------------------------------------------------------
 6900 FORMAT(///' *** ERROR RETURN FROM  S A M  LIBRARY ROUTINE GORSKY')
 6910 FORMAT(/5X,'ILLEGAL PARAMETER(S) : N, NV, K = ',3I8)
 6920 FORMAT(/5X,'ERROR OCCURED FOR COLUMN NO.',I6,' OF MATRIX G')
 6930 FORMAT(/5X,'COLUMN NO.',I5,' OF G CANNOT BE NORMALIZED' /
     +        5X,'WITH RESPECT TO MATRIX  A' )
 6940 FORMAT(/5X,'COLUMN NO.',I5,' OF  G  IS LINEAR DEPENDENT')
C
      END
