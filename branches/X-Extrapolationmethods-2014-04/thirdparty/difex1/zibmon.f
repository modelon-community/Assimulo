      SUBROUTINE MON
C
C* Begin Prologue MON
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Time monitor to discover cpu time consumption of interesting parts
C  of a code.
C
C* Written by        U. Nowak, U. Poehle, L. Weimann
C* Purpose           Cpu time measuring
C* Category          ...
C* File              mon.f
C* Version           2.1
C* Latest Change     96/03/04 (1.3)
C* Library           CodeLib
C* Code              Fortran 77
C                    Single precision
C
C  ---------------------------------------------------------------------
C
C* Summary
C  -------
C 
C  MON provides different entries 
C   - to initialize the time monitor, 
C   - to start and to stop the measuring, 
C   - to start and to stop a number of stop-watches which may be nested,
C   - to print a summary table,
C   - to define text labels for the summary table, and
C   - to store the vector of average cpu times
C
C  Errors detected during the measuring will disable the time monitor
C  but will not affect the calling program.
C
C  ---------------------------------------------------------------------
C
C* Parameter list description
C  --------------------------
C
C* Input parameters
C  ----------------
C
C  TEXTH      Char    Text to identify the measuring (up to 31
C                     characters will be printed in the headline of the
C                     summary table)
C                     Used by MONINI
C
C  IOUNIT     Int     Output unit for error messages and summary table
C                     Used by MONINI
C
C  INDX       Int     Number of a range of source code lines (a distinct
C                     part of the algorithm normally) for which the
C                     cpu time will be measured by a kind of stop-watch
C                     Used by MONDEF, MONON, and MONOFF
C
C  NAMEH      Char    Text to label the summary lines (up to 17
C                     characters will be printed)
C                     Used by MONDEF
C
C* Output parameters
C  -----------------
C
C  IRET       Int     Return code of ZIBSEC
C                     Used by MONSTR
C
C  AVER       Dble    Vector of average cpu time for all measurings
C                     Used by MONGET
C
C
C* End Prologue
C  ------------
C
C
C* Constants
C  ---------
C
C  MAXTAB     Int     Maximum number of stop-watches
C  MNEST      Int     Maximum number of nested measurings
C
C
C* Local Variables
C  ---------------
C
C  ASEC       Real    Array of averages per call
C  INDACT     Int     Array for indices of nested measurings
C  NCALL      Int     Array of call counts
C  PC1        Real    Array of per cent values with respect to all
C                     measurings
C  PC2        Real    Array of per cent values with respect to the sum
C                     of all stop-watches
C  QDISAB     Log     Time monitor disabled
C  QON        Log     Array reflecting if an index is active
C  QSTART     Log     Time monitor started
C  SEC        Real    Array of cpu time measurings
C
C
C
      INTEGER MAXTAB, MNEST
C
      PARAMETER (MAXTAB = 25, MNEST=20)
C
      INTEGER INDACT(MNEST), NCALL(0:MAXTAB)
      REAL ASEC(0:MAXTAB), PC1(0:MAXTAB), PC2(0:MAXTAB), SEC(0:MAXTAB)
      DOUBLE PRECISION AVER(0:MAXTAB)
      LOGICAL QDISAB, QON(0:MAXTAB), QSTART
      CHARACTER NAME(0:MAXTAB)*17, TEXT*30 
      CHARACTER*(*) NAMEH, TEXTH
C
      SAVE ASEC, CPTIM, INDACT, IONCNT, MAXIND, MONI, NAME, NCALL, PC1,
     $     PC2, QDISAB, QON, QSTART, SEC, TEXT
C
C
      DATA CPTIM /0.0E0/, MONI /6/, INFO /1/, IONCNT /-1/,
     $     QDISAB /.FALSE./, QSTART /.FALSE./
C
      RETURN
C
C
      ENTRY MONINI (TEXTH, IOUNIT)
C     Initialize monitor.
C     Has to be called first. May be called again after MONHLT.
C
      TEXT = TEXTH
      MONI = IOUNIT
C
      IF (IONCNT .GT. 0 .AND. .NOT. QDISAB) GOTO 1070
C
      MAXIND = 0
      IONCNT = 0
      QDISAB = .FALSE.
C
      DO 1000 I = 0,MAXTAB
         SEC(I) = 0.
         ASEC(I) = 0.
         NCALL(I) = 0
         QON(I) = .FALSE.
         WRITE (NAME(I), '(A, I2)') 'Part ', I
 1000 CONTINUE
      NAME(0) = 'General'
C
      DO 1010 I = 1,MNEST
         INDACT(I) = 0
 1010 CONTINUE
C
      RETURN
C
C
      ENTRY MONDEF(INDX, NAMEH)
C     Define one monitor entry.
C     May be called at any time before MONPRT.
C
      IF (QDISAB) RETURN
      IF (INDX.LT.0 .OR. INDX.GT.MAXTAB) GOTO 1080
C
      NAME(INDX) = NAMEH
C
      RETURN
C
C
      ENTRY MONSTR (IRET)
C     Start monitor measurings.
C     Has to be called once after initialization.
C
      IF (QDISAB) RETURN
      IF (IONCNT .LT. 0) GOTO 1090
      IF (IONCNT .GT. 0) GOTO 1100
      IF (QON(0)) GOTO 1110
C
      IFAIL = 0
      CALL ZIBSEC (CPTIM, IFAIL)
C
      IF (IFAIL .EQ. 0) THEN
C        Switch on general stop-watch
         SEC(0) = -CPTIM
         QON(0) = .TRUE.
         IONCNT = 1
         QSTART = .TRUE.
      ENDIF
      IRET = IFAIL
C
      RETURN
C
C
      ENTRY MONON (INDX)
C     Start one measuring.
C     A running stop-watch will be deactivated until the new stop-watch
C     stops.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
      IF (IONCNT .LT. 1) GOTO 1120
      IF (INDX .GT. MAXTAB .OR. INDX .LE. 0) GOTO 1130
      IF (QON(INDX)) GOTO 1140
C
      MAXIND = MAX(MAXIND, INDX)
      CALL ZIBSEC (CPTIM, IFAIL)
C
C     Hold actual stop-watch
      SEC(INDACT(IONCNT)) = SEC(INDACT(IONCNT)) + CPTIM
C
C     Switch on new stop-watch
      IF (INFO .GT. 1) WRITE (MONI,*) ' Enter ', NAME(INDX), SEC(INDX)
C
      IONCNT = IONCNT + 1
      IF (IONCNT .GT. MNEST) GOTO 1150
C
      INDACT(IONCNT) = INDX
      SEC(INDX) = SEC(INDX) - CPTIM
      QON(INDX) = .TRUE.
C
      RETURN
C
C
      ENTRY MONOFF (INDX)
C     Stop one measuring.
C     May be called for the stop-watch started most recently.
C     The stop-watch deactivated most recently will be activated.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
      IF (INDX .GT. MAXTAB .OR. INDX .LE. 0) GOTO 1160
      IF (INDACT(IONCNT) .NE. INDX) GOTO 1170
C
      CALL ZIBSEC (CPTIM, IFAIL)
C
C     Switch off actual stop-watch
      QON(INDX) = .FALSE.
      SEC(INDX) = SEC(INDX) + CPTIM
      NCALL(INDX) = NCALL(INDX) + 1
      IONCNT = IONCNT - 1
      IF (INFO .GT. 1) WRITE (MONI,*) ' Exit ', NAME(INDX), SEC(INDX)
C
C     Continue previous stop-watch
      SEC(INDACT(IONCNT)) = SEC(INDACT(IONCNT)) - CPTIM
C
      RETURN
C
C
      ENTRY MONHLT
C     Terminate monitor.
C     Stops all active stop-watches.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
C
      CALL ZIBSEC (CPTIM, IFAIL)
C
      DO 1020 I=IONCNT,1,-1
         QON(INDACT(I)) = .FALSE.
         SEC(INDACT(I)) = SEC(INDACT(I)) + CPTIM
         NCALL(INDACT(I)) = NCALL(INDACT(I)) + 1
 1020 CONTINUE
C
      IONCNT = 0
C     This means that the time monitor has to be started by MONSTR
C     before another measuring.
C
      RETURN
C
C
C
      ENTRY MONPRT
C     Print statistics.
C     May be called after MONHLT only.
C
      IF (IONCNT .GT. 0) GOTO 1180
      IF (.NOT. QSTART) GOTO 1190
C
      SUM = 1.E - 10
      DO 1030 I = 1,MAXIND
         SUM = SUM + SEC(I)
         IF (NCALL(I) .GT. 0) ASEC(I) = SEC(I)/FLOAT(NCALL(I))
 1030 CONTINUE
      SUM0 = SUM + SEC(0)
      IF (NCALL(0) .GT. 0) ASEC(0) = SEC(0)/FLOAT(NCALL(0))
C
      DO 1040 I = 1,MAXIND
         PC1(I) = 100.*SEC(I)/SUM0
         PC2(I) = 100.*SEC(I)/SUM
 1040 CONTINUE
      PC1(0) = 100.*SEC(0)/SUM0
C
      WRITE (MONI,9000)
      WRITE (MONI,9010)
      WRITE (MONI,9020) ' '
 9000 FORMAT (///)
 9010 FORMAT (1X, 77('#'))
 9020 FORMAT (' #   ', A, T75,  '   #')
C
      IF (QDISAB) THEN
         WRITE (MONI,9020) ' '
         WRITE (MONI,9020)
     $        'Warning  The following results may be misleading'
         WRITE (MONI,9020)
     $        'because an error occured and disabled the time monitor'
      ENDIF
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9030) 'Results from time monitor program for:', TEXT
 9030 FORMAT (' #   ', A40, 1X, A31, '#')
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9040) 'Total time:', SUM0, 'Sum of parts:', SUM
 9040 FORMAT (' #   ', A11, F11.3, 5X, A13, F11.3, 21X, '#')
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9050)
 9050 FORMAT (' #   ', 2X, 'Name', 14X, 'Calls', 7X, 'Time', 4X,
     $     'Av-time', 4X, '% Total', 6X, '% Sum   #')
C
      WRITE (MONI,9060) NAME(0), NCALL(0), SEC(0), ASEC(0), PC1(0)
 9060 FORMAT (' #   ', A17, I8, F11.3, F11.4, F11.2, 14X, '#')
C
      DO 1050 I = 1,MAXIND
         WRITE (MONI,9070) NAME(I), NCALL(I), SEC(I), ASEC(I), PC1(I),
     $        PC2(I)
 1050 CONTINUE
C
 9070 FORMAT (' #   ', A17, I8, F11.3, F11.4, F11.2, F11.2, 3X, '#')
C
C
      WRITE (MONI,9020) ' '
      WRITE (MONI,9010)
      WRITE (MONI,9000)
C
      RETURN
C
C
      ENTRY MONGET(AVER)
C     Store average cpu times per call.
C     May be called at any time after MONPRT.
C
      IF (.NOT. QSTART) RETURN
      IF (QDISAB) RETURN
C
      DO 1060 I=0,MAXIND
         AVER(I) = DBLE(ASEC(I))
 1060 CONTINUE
C
      RETURN
C
C  Error exits
C
 1070 CONTINUE
      WRITE (MONI,9080) 'MONINI', 'Time monitor is running already.'
      GOTO 1200
C
 1080 CONTINUE
      WRITE (MONI,9090) 'MONDEF', 'Index out of range', INDX
      GOTO 1200
C
 1090 CONTINUE
      WRITE (MONI,9080) 'MONSTR',
     $     'Time monitor has to be initialized by MONINI first.'
      GOTO 1200
C
 1100 CONTINUE
      WRITE (MONI,9080) 'MONSTR', 'Time monitor is running already.'
      GOTO 1200
C
 1110 CONTINUE
      WRITE (MONI,9080) 'MONSTR',
     $     'Time monitor has been started already.'
      GOTO 1200
C
 1120 CONTINUE
      WRITE (MONI,9080) 'MONON', 'Time monitor is not yet started.'
      GOTO 1200
C
 1130 CONTINUE
      WRITE (MONI,9090) 'MONON', 'Index out of range', INDX
      GOTO 1200
C
 1140 CONTINUE
      WRITE (MONI,9090) 'MONON',
     $     'Measuring is running already for this INDX', INDX
      GOTO 1200
C
 1150 CONTINUE
      WRITE (MONI,9100) 'MONON', 'Nesting is too deep.',
     $     'The following indices are active', (INDACT(I),I=0,IONCNT)
      GOTO 1200
C
 1160 CONTINUE
      WRITE (MONI,9090) 'MONOFF', 'Index out of range', INDX
      GOTO 1200
C
 1170 CONTINUE
      WRITE (MONI,9110) 'MONOFF', 'Measuring ', INDX,
     $     'cannot be stopped.',
     $     'The following indices are active', (INDACT(I),I=0,IONCNT)
      GOTO 1200
C
 1180 CONTINUE
      WRITE (MONI,9080) 'MONPRT', 'Time monitor is still running.'
      GOTO 1200
C
 1190 CONTINUE
      WRITE (MONI,9080) 'MONPRT', 'Time monitor was not started.'
      GOTO 1200
C
 1200 QDISAB = .TRUE.
      RETURN
C
 9080 FORMAT (/, ' ++ Error in subroutine ', A, ' ++',/, 4X, A)
C
 9090 FORMAT (/, ' ++ Error in subroutine ', A, ' ++',/, 4X, A,
     $     ' (', I6, ').')
C
 9100 FORMAT (/, ' ++ Error in subroutine ', A, ' ++', 4X, A,/, 4X, A,
     $     (10I4))
C
 9110 FORMAT (/, ' ++ Error in subroutine ', A, ' ++', 4X, A, I3, 1X, A,
     $     /, 4X, A, (10I4))
C
C  End subroutine monitor
C
      END
