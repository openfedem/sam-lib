C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRSHA (MPARA,MTREEA,MSIFA,
     $                   MPARB,MTREEB,MSIFB,MDLOC,
     $                   A,B,SHIFT,N,KSB,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRSHA                 GROUP 9 / PRIVATE
C
C     T A S K :  To subtract SHIFT*B  from the
C                symmetric sparse matrix  A.
C     Matrix  B  may be a symmetric sparse (compressed) matrix (KSB=1)
C     or a diagonal matrix (KSB=2)
C
C     ROUTINES CALLED/REFERENCED :  SPRSH2, SPRDG1              (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   96-12-22 / 1.0
C                       21-08-16 / 2.0   K.M.Okstad
C
C **********************************************************************
C
      IMPLICIT           NONE
C
      INTEGER            N,KSB,LPU,IERR
      INTEGER            MPARA(50),MTREEA(*),MSIFA(*),
     $                   MPARB(50),MTREEB(*),MSIFB(*),MDLOC(N)
      DOUBLE PRECISION   A(*),B(*),SHIFT
C
      INTEGER            NSUPRA,NOFSUA,XSUPRA,XLINDA,LINDXA,XLNZA,
     $                   NSUPRB,NOFSUB,XSUPRB,XLINDB,LINDXB,XLNZB,
     $                   NZEROL,PERE2I,LNZA,LNZB,I,J
C
C ----------------------------------------------------------------------
C
      IF (KSB .EQ. 1) THEN
C                                               ** sparse B
         NZEROL = MPARB(16)
         LINDXB = MPARB(47)
         IF (MPARA(16) .EQ. NZEROL) THEN
C
            J = MPARA(8)
            DO 10 I = 1, NZEROL
               A(J+I) = A(J+I) - SHIFT*B(I)
   10       CONTINUE
C
         ELSEIF (LINDXB .GT. 1) THEN
C                                               ** different pattern
            NSUPRA = MPARA(11)
            NOFSUA = MPARA(15)
            XSUPRA = MPARA(45)
            XLINDA = MPARA(46)
            LINDXA = MPARA(47)
            XLNZA  = MPARA(48)
C
            NSUPRB = MPARB(11)
            NOFSUB = MPARB(15)
            XSUPRB = MPARB(45)
            XLINDB = MPARB(46)
            XLNZB  = MPARB(48)
C
            LNZA   = 1 + MPARA(8)
            LNZB   = 1 + MPARB(8)
            PERE2I = 1 + MPARB(8)
C
            CALL SPRSH2 (NSUPRA, NOFSUA,
     $                   MTREEA(XSUPRA), MSIFA(XLNZA),
     $                   MTREEA(XLINDA), MSIFA(LINDXA), A(LNZA),
     $                   NSUPRB, NOFSUB, MPARB(8),
     $                   MTREEB(XSUPRB), MSIFB(XLNZB) ,
     $                   MTREEB(XLINDB), MSIFB(LINDXB), B(LNZB),
     $                   MSIFB(PERE2I) , SHIFT, IERR )
C
         ELSE
            IERR = -1
         ENDIF
         IF (IERR .LT. 0) THEN
            IF (LPU .GT. 0) WRITE(LPU,600) IERR
            MPARA(1) = -66
            RETURN
         ENDIF
C
      ELSEIF (KSB .EQ. 2) THEN
C                                               ** diagonal B
         CALL SPRDG1 (MPARA,MTREEA,MSIFA,MDLOC,LPU,IERR)
         IF (IERR .LT. 0) RETURN
C
         DO 20 J = 1, N
            I    = MDLOC(J)
            A(I) = A(I) - SHIFT*B(J)
   20    CONTINUE
C
      ELSE
         RETURN
      ENDIF
C
      MPARA(1) = 66
C
      RETURN
 600  FORMAT(//'*** ERROR RETURN FROM SAM ROUTINE SPRSHA, IERR =',I8)
      END


      SUBROUTINE SPRSH2 (NSUPRA, NOFSUA, XSUPRA, XLNZA ,
     $                   XLINDA, LINDXA, LNZA  ,
     $                   NSUPRB, NOFSUB, NEQ   , XSUPRB, XLNZB ,
     $                   XLINDB, LINDXB, LNZB  , PERE2I, SHIFT , IERR  )
C
      IMPLICIT           NONE
C
      INTEGER            NSUPRA, NOFSUA, NSUPRB, NOFSUB, NEQ   , IERR
      INTEGER            XSUPRA(NSUPRA+1)      , XLNZA(NSUPRA+1)       ,
     $                   XLINDA(NSUPRA+1)      , LINDXA(NOFSUA)        ,
     $                   XSUPRB(NSUPRB+1)      , XLNZB(NSUPRB+1)       ,
     $                   XLINDB(NSUPRB+1)      , LINDXB(NOFSUB)        ,
     $                   PERE2I(NEQ)
      DOUBLE PRECISION   LNZA(*), LNZB(*), SHIFT
C
C     ---------------------------------------------------------------
C --- To subtract  SHIFT*B  from  A  where both A and B are
C     symmetric sparse matrices but with different sparsity pattern.
C     It is assumed that all non-zero elements in matrix B has a
C     corresponding non-zero element in matrix A, but not vice versa.
C     ---------------------------------------------------------------
C
      INTEGER            FSTVAR, I     , IEQ   , IPLNZA, IPLNZB, ISTRT ,
     $                   ISTOP , J     , JEQ   , JS    , LSTVAR
C
      IPLNZA = -1
      DO 400 JS = 1, NSUPRB
         FSTVAR = XSUPRB(JS)
         LSTVAR = XSUPRB(JS+1) - 1
         ISTRT  = XLINDB(JS)
         ISTOP  = XLINDB(JS+1) - 1
         IPLNZB = XLNZB(JS)
         DO 300 J = FSTVAR, LSTVAR
            JEQ = PERE2I(J)
            CALL SPRSH3 (NSUPRA, NOFSUA, XSUPRA, XLNZA ,
     $                   XLINDA, LINDXA, IPLNZA, JEQ   , JEQ   )
            IF (IPLNZA .LT. 0) GOTO 900
C
            LNZA(IPLNZA) = LNZA(IPLNZA) - SHIFT*LNZB(IPLNZB)
            IPLNZB = IPLNZB + 1
C
            DO 200 I = ISTRT+1, ISTOP
               IEQ = PERE2I(LINDXB(I))
               CALL SPRSH3 (NSUPRA, NOFSUA, XSUPRA, XLNZA ,
     $                      XLINDA, LINDXA, IPLNZA, IEQ   , JEQ   )
               IF (IPLNZA .LT. 0) GOTO 900
C
               LNZA(IPLNZA) = LNZA(IPLNZA) - SHIFT*LNZB(IPLNZB)
               IPLNZB = IPLNZB + 1
  200       CONTINUE
C
            ISTRT = ISTRT + 1
  300    CONTINUE
  400 CONTINUE
C
  900 CONTINUE
      IERR = MIN(0,IPLNZA)
      RETURN
      END


      SUBROUTINE SPRSH3 (NSUPER, NOFSUB, XSUPER, XLNZ  ,
     $                   XLINDX, LINDX , IPLNZ , I, J  )
C
      IMPLICIT           NONE
C
      INTEGER            NSUPER, NOFSUB, IPLNZ , I, J
      INTEGER            XSUPER(NSUPER+1)      , XLNZ(NSUPER+1)        ,
     $                   XLINDX(NSUPER+1)      , LINDX(NOFSUB)
C
C     ------------------------------------------------------------------
C --- TO RETRIEVE THE STOARGE INDEX IN  LNZ  OF THE MATRIX ELEMENT (I,J)
C     ------------------------------------------------------------------
C
      INTEGER            FSTVAR, IEQ   , II    , ISTRT ,
     $                   ISTOP , JEQ   , JJ    , JS    , LSTVAR
C
      IF (J .GT. I) THEN
         IEQ = J
         JEQ = I
      ELSE
         IEQ = I
         JEQ = J
      ENDIF
C
      DO 300 JS = 1, NSUPER
         FSTVAR = XSUPER(JS)
         LSTVAR = XSUPER(JS+1) - 1
         IF (JEQ .GE. FSTVAR .AND. JEQ .LE. LSTVAR) THEN
            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1) - 1
            IPLNZ = XLNZ(JS)
            DO 100 JJ = FSTVAR, JEQ-1
               IPLNZ = IPLNZ + 1+ISTOP-ISTRT
               ISTRT = ISTRT + 1
  100       CONTINUE
            IF (IEQ .EQ. JEQ) RETURN
            DO 200 II = ISTRT+1, ISTOP
               IF (IEQ .EQ. LINDX(II)) THEN
                  IPLNZ = IPLNZ + II-ISTRT
                  RETURN
               ENDIF
  200       CONTINUE
C           Corrupt matrix data structure,
C           equation pair (IEQ,JEQ) not found
            IPLNZ = -IEQ
            RETURN
         ENDIF
  300 CONTINUE
C
C     Corrupt matrix data structure, equation JEQ not found
      IPLNZ = -XSUPER(NSUPER+1)
      RETURN
      END
