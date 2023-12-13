C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE   SPRER1   ( N     , NAME  , I1    , I2    , I3    ,
     $                        LPU   , IERR                            )

C     ------------------------------------------------------------------

      INTEGER           N     , I1    , I2    , I3    , LPU   , IERR
      CHARACTER*(*)     NAME

C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRER1                 GROUP 9 / PRIVATE
C
C     T A S K :  To write (on LPU) warning/error messages for SPR-
C                routines.
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-01-11 / 1.0
C
C --- JEG HAR TATT MINE FEILMELDINGER INN I DENNE RUTINEN, OENSKER DU
C     AA ENDRE LAY-OUT ER DET HELT OPP TIL DEG.
C
C **********************************************************************

      IF (LPU .LT. 0 .OR. IERR .EQ. 0) GO TO 100

C     --------------
C     Write heading.
C     --------------

      IF ( N .LE. 20 ) THEN
         IF ( IERR .GT. 0 ) THEN
            WRITE(LPU,400) NAME
         ELSE
            WRITE(LPU,600) NAME
         ENDIF
      ELSEIF ( N .LT. 60 ) THEN
         WRITE ( LPU, '(///A,A/,4X,A,I4/)' )
     $   ' *** ERROR RETURN FROM  S A M / SPR  LIBRARY ROUTINE : ',
     $   NAME,
     $   ' WITH ERROR FLAG : ', IERR
      ELSE
         WRITE ( LPU, '(///A,A/)' )
     $   '  ** WARNING FROM  S A M / SPR  LIBRARY ROUTINE : ', NAME
      ENDIF

C     ==================================================================

C     --------------------
C     Error/Warnings BELL.
C     --------------------
      IF (N .EQ.  1) THEN
         WRITE (LPU,601) I1
         WRITE (LPU,701) I2,I3
      ELSEIF (N .EQ. 2) THEN
         WRITE (LPU,701) I1,I2
         WRITE (LPU,602) I3
      ELSEIF (N .EQ. 3) THEN
         WRITE (LPU,603) I1
      ELSEIF (N .EQ. 4) THEN
         WRITE (LPU,701) I1,I2
         WRITE (LPU,604)
      ELSEIF (N .EQ. 5) THEN
         WRITE (LPU,605)
      ELSEIF (N .EQ. 6) THEN
         WRITE (LPU,701) I1,I2
         WRITE (LPU,606)
      ELSEIF (N .EQ. 7) THEN
         WRITE (LPU,607) I1
      ENDIF

C     ==================================================================

C     --------------
C     INCONCISTENCY.
C     --------------
      IF ( N .EQ. 21 ) THEN
         WRITE ( LPU, '(4X,A)' )
     $   ' INCONSISTENCY IN CONTROL DATA, CHECK MPAR AND MSPAR'
      ELSEIF ( N .EQ. 22 ) THEN
         WRITE ( LPU, '(4X,A)' )
     $   ' INCONSISTENCY IN CONTROL DATA, CHECK SAM CONTROL STRUCTURES'
      ELSEIF ( N .EQ. 23 ) THEN
         WRITE ( LPU, '(4X,A)' )
     $   ' INCONSISTENCY IN CONTROL DATA, CHECK SPR CONTROL STRUCTURES'
      ELSEIF ( N .EQ. 24 ) THEN
         WRITE ( LPU, '(4X,A/,4X,A,I12/,4X,A,I12/,A)' )
     $   ' INCONSISTENCY IN CONTROL DATA',
     $   ' ORDER OF B, C IS GIVEN AS : ', I3,
     $   ' ORDER OF A IS GIVEN AS    : ', I2,
     $   ' ORDER OF B,C MUST BE GREATER THAN OR EQUAL TO ORDER OF A'
      ELSEIF ( N .EQ. 25 ) THEN
         WRITE ( LPU, '(4X,A/,4X,A,I12/,4X,A,I12)' )
     $   ' OPTION IFLAG IS SET TO AN ILLEGAL VALUE',
     $   ' THE VALUE IS : ', I1,
     $   ' NUMBER OF COLUMNS IN RECTANGULAR MATRIX :', I2
      ELSEIF ( N .EQ. 26 ) THEN
         WRITE ( LPU, '(4X,A/,4X,A,I12)' )
     $   ' OPTION KSA IS SET TO AN ILLEGAL VALUE',
     $   ' THE VALUE IS : ', I1
      ENDIF
      IF ( N .GT. 20 .AND. N .LT. 24 ) THEN
         IF ( I1 .GT. 0 ) THEN
            WRITE(LPU,610) I1, I2, I3
         ELSE IF ( I1 .LT. 0 ) THEN
            WRITE(LPU,611) I2, I3
         ELSE
            WRITE(LPU,"('')")
         ENDIF
      ENDIF

C     ------------------
C     ERROR NOT CLEARED.
C     ------------------
      IF ( N .EQ. 30 ) THEN
         WRITE ( LPU, '(4X,A/)' )
     $   ' ERROR FROM PRECEEDING  S A M / SPR  ROUTINE IS NOT CLEARED'
      ENDIF

C     ------------
C     ARRAY SIZES.
C     ------------
      IF ( N .EQ. 31 ) THEN
         WRITE(LPU,801) 'MSICA, AND/OR IWORK'
      ELSEIF ( N .EQ. 32 ) THEN
         WRITE(LPU,801) 'MSICA'
         WRITE(LPU,'(5X,A,A)') 'PROBABLY CAUSED BY CORRUPTED ',
     $        'SAM CONTROL STRUCTURES'
      ELSEIF ( N .EQ. 33 ) THEN
         WRITE(LPU,801) 'MTREES, AND/OR IWORK'
      ELSEIF ( N .EQ. 34 ) THEN
         WRITE(LPU,801) 'IWORK'
      ELSEIF ( N .EQ. 35 ) THEN
         WRITE(LPU,801) 'MSICA/MTREES/MSIFA, AND/OR IWORK'
      ELSEIF ( N .EQ. 36 ) THEN
         WRITE(LPU,801) 'MTREES/MSIFA/SM, AND/OR IWORK/RWORK'
      ELSEIF ( N .EQ. 37 ) THEN
         WRITE(LPU,801) 'MTREES/MSIFA/SM/B, AND/OR IWORK/RWORK'
      ELSEIF ( N .EQ. 38 ) THEN
         WRITE(LPU,801) 'MTREES'
      ELSEIF ( N .EQ. 39 ) THEN
         WRITE(LPU,801) 'MSIFA'
      ENDIF
      IF ( N .GT. 30 .AND. N .LT. 40 .AND. I1 .GT. I2 ) THEN
         WRITE(LPU,802) I1, I2
      ENDIF

C     -----------------
C     NUMERICAL ERRORS.
C     -----------------
      IF ( N .EQ. 40 ) THEN
         IF     ( I3 .EQ. -10 ) THEN
            WRITE ( LPU, '(4X,A)' )
     $      ' INSUFFICIENT REAL WORKSPACE, CHECK SPR CONTROL STRUCTURES'
         ELSEIF ( I3 .EQ. -20 ) THEN
            WRITE ( LPU, '(4X,A,A)' )
     $      ' INSUFFICIENT SPACE FOR FACTORS,',
     $      ' CHECK SPR CONTROL STRUCTURES'
         ELSEIF ( I3 .EQ. -30 .OR. I3 .EQ. -31 ) THEN
            WRITE ( LPU, '(4X,A)' )
     $      ' ABSOLUTE PIVOT VALUE REDUCED TO SIZE LESS THAN',
     $      ' OR EQUAL TO THE THRESHOLD VALUE GIVEN BY TOL(1)'
         ELSEIF ( I3 .EQ. -32 .OR. I3 .EQ. -33  ) THEN
            WRITE ( LPU, '(4X,A)' )
     $      ' PIVOT VALUE REDUCED TO SIZE LESS THAN OR',
     $      ' EQUAL TO THE THRESHOLD VALUE GIVEN BY TOL(1)'
         ELSEIF ( I3 .EQ. -35 ) THEN
            WRITE ( LPU, '(4X,A,I6,A)' )
     $      ' SINGULAR EQUATION SYSTEM, THERE ARE', I1, ' ZERO PIVOT(S)'
         ELSEIF ( I3 .EQ. -40 ) THEN
            WRITE ( LPU, '(4X,A)' )
     $      ' INSUFFICIENT REAL WORKSPACE, CHECK SPR CONTROL STRUCTURES'
         ENDIF
         IF ( I3 .EQ. -31 .OR. I3 .EQ. -33 ) THEN
            WRITE ( LPU, '(/4X,A,I6,2X,A/)' )
     $      ' OCCURED IN', I2, 'EQUATION(S)'
         ELSE
            WRITE ( LPU, '(/4X,A,I12/)' )
     $      ' OCCURED IN EQUATION:', I2
         ENDIF
         IF ( I1 .GT. 0 .AND. I3 .EQ. -30 ) THEN
            WRITE ( LPU, '(4X,A,I6,A/)' )
     $      ' THERE ARE ALSO',I1,' NEGATIVE PIVOTS BEFORE THIS EQUATION'
         ELSE IF ( I1 .GT. 0 .AND. I3 .EQ. -31 ) THEN
            WRITE ( LPU, '(4X,A,I6,A/)' )
     $      ' THERE ARE ALSO IN ADDITION',I1,' NEGATIVE PIVOTS'
         ENDIF
      ENDIF

C     ----------------
C     ASSEMBLY ERRORS.
C     ----------------
      IF ( N .EQ. 57 ) THEN
         WRITE ( LPU, '(4X,A,I6,A,I6,A/,4X,A,I10)' )
     $   ' INCORRECT VALUE OF ARGUMENT NSV (=',I2,', NPDOF =',I3,
     $   ')',' FOR ELEMENT', I1
      ENDIF

C     -------------------
C     NUMERICAL WARNINGS.
C     -------------------
      IF ( N .EQ. 60 ) THEN
         WRITE ( LPU, '(4X,A/,4X,A,I12,4X,A)' )
     $   ' THE STIFFNESS MATRIX IS NOT POSITIVE DEFINITE',
     $   ' THERE ARE', I3, 'NEGATIVE PIVOT(S)'
      ENDIF
      IF ( N .EQ. 61 ) THEN
         WRITE ( LPU, '(4X,A/,4X,A,I12,4X,A,I12)' )
     $   ' THE SIZE OF THE SPR CONNECTIVITY IS LESS THAN ANTICIPATED',
     $   ' ESTIMATED SIZE =', I1, 'ACTUAL SIZE =', I2
      ENDIF

C     ==================================================================

#ifdef FT_DEBUG
      IF (N .EQ. 40 .AND. I3/10 .EQ. -3) GOTO 100
      IF (N .GE. 60) GOTO 100
      IERR = IERR/0 ! Force a core dump for easy debugging
#endif

  100 RETURN

  400 FORMAT (///' *** WARNING from SAM library routine ',A)
  600 FORMAT (///' *** ERROR return from SAM library routine ',A)
  601 FORMAT (5X,'Illegal status code',I8,'  encountered for')
  602 FORMAT (5X,'is coupled to too many (=',I4,') master DOFs')
  603 FORMAT (5X,'Illegal DOF (number',I8,') is constrained')
  604 FORMAT (5X,'is constrained but not specified (in MSC)')
  605 FORMAT (5X,'Control matrix MPMCEQ is in error (corrupt?)')
  606 FORMAT (5X,'is coupled to a specified DOF')
  607 FORMAT (5X,'Illegal or inconsistent number of DOFs (=',I8,')')
  610 FORMAT (5X,'DETECTED FOR ELEMENT',I8,' :',2I8/)
  611 FORMAT (5X,'DISCREPANCY IN THE NUMBER OF SPR-NODES :',2I8/)
  701 FORMAT (5X,'DOF number',I4,' of node',I8)
  801 FORMAT (5X,A,': NOT GREAT ENOUGH')
  802 FORMAT (5X,'NEED',I8,'  WORDS, HAVE ONLY',I8,'  AVAILABLE'/)
C     ------------------------------------------------------------------
      END
