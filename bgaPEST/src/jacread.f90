	module jacread

        use utilities
        
        
        public readJCO
        
        contains 
        
        subroutine readJCO(jacfle, X)
        IMPLICIT NONE

! -- Program JACWRIT translates a binary Jacobian "jco" file to an ascii file.

        INTEGER NESPAR,NXROW,IERR,J,I,IPP,IOBS,NBLC,NEWFLAG,IROW,IES,ICOUNT
        INTEGER NBLNK
        integer jacunit
        DOUBLE PRECISION DTEMP
        CHARACTER*12 AVERSION
        CHARACTER*(*) JACFLE
        CHARACTER*200 COMLIN


        DOUBLE PRECISION, pointer :: X(:,:)
        CHARACTER*12 APAR(:)
        CHARACTER*12 AOBS1(:)
        CHARACTER*20 AOBS(:)

        ALLOCATABLE::APAR,AOBS,AOBS1


        NEWFLAG=0
        jacunit = utl_nextunit()
        OPEN(UNIT=jacunit,FILE=JACFLE,FORM='binary', STATUS='OLD',IOSTAT=IERR)
        IF(IERR.NE.0)THEN
          WRITE(*,125) JACFLE(1:NBLNK(JACFLE))
125       FORMAT(/,' *** Cannot open file ',A,' ***',/)
          GO TO 9990
        END IF
         READ(jacunit,ERR=9000,END=9100)NESPAR,NXROW
        IF(NESPAR.LT.0)THEN
          NESPAR=-NESPAR
          NXROW=-NXROW
          NEWFLAG=1
        ELSE
          NEWFLAG=0
        END IF

        IF(NEWFLAG.EQ.1)THEN
          ALLOCATE(X(NXROW,NESPAR),AOBS(NXROW),APAR(NESPAR), STAT=IERR)
        ELSE
          ALLOCATE(X(NXROW,NESPAR),AOBS1(NXROW),APAR(NESPAR),AOBS(NXROW),STAT=IERR)
        END IF
        IF(IERR.NE.0)THEN
          WRITE(*,200)
200       FORMAT(/,' *** Cannot allocate sufficient memory to ','continue execution ***',/)
          GO TO 9990
        END IF

        IF(NEWFLAG.EQ.1)THEN
          X=0.0D0                ! AN ARRAY
          READ(jacunit,ERR=9000,END=9100)ICOUNT
          DO I=1,ICOUNT
            READ(jacunit,ERR=9000,END=9100) J,DTEMP
            IES=(J-1)/NXROW+1
            IROW=J-(IES-1)*NXROW
            X(IROW,IES)=DTEMP
          END DO
        ELSE
          READ(jacunit,ERR=9000,END=9100) ((X(J,I),J=1,NXROW),I=1,NESPAR)
        END IF
        DO IPP=1,NESPAR
          READ(jacunit,ERR=9000,END=9100) APAR(IPP)
        END DO
        IF(NEWFLAG.EQ.1)THEN
          DO IOBS=1,NXROW
            READ(jacunit,ERR=9000,END=9100) AOBS(IOBS)
          END DO
        ELSE
          DO IOBS=1,NXROW
            READ(jacunit,ERR=9000,END=9100) AOBS1(IOBS)
            AOBS(IOBS)=AOBS1(IOBS)
          END DO
        END IF
        CLOSE(UNIT=jacunit)



        GO TO 9995
        
9000    WRITE(*,9010) JACFLE(1:NBLNK(JACFLE))
9010    FORMAT(/,' *** Error encountered in reading file ',A,' ***',/)
        GO TO 9990
9100    WRITE(*,9110) JACFLE(1:NBLNK(JACFLE))
9110    FORMAT(/,' *** Unexpected end encountered to file ',A, ' ***',/)
        GO TO 9990


9900    WRITE(*,9910,ERR=9000)
9910    FORMAT(' JACWRIT is run using the command:',/)
        WRITE(*,9920,ERR=9000)
9920    FORMAT('    jacwrit jacfile1 jacfile2',/,/,' where',/)
        WRITE(*,9930,ERR=9000)
9930    FORMAT('    "jacfile1" is a PEST binary Jacobian file ','(ext ".jco"), and')
        WRITE(*,9940,ERR=9000)
9940    FORMAT('    "jacfile2" is a text Jacobian file to be ','written by JACWRIT.')

        GO TO 9990


9990    CONTINUE

        DEALLOCATE(X,APAR,AOBS,STAT=IERR)
        IF(NEWFLAG.EQ.0)THEN
          DEALLOCATE(AOBS1,STAT=IERR)
        END IF


9995    CONTINUE

        DEALLOCATE(APAR,AOBS,STAT=IERR)

        END subroutine



        SUBROUTINE UPCAS(ASTRNG)

! -- SUBROUTINE UPCAS CONVERTS A STRING TO UPPER CASE

        INTEGER NBLNK
        INTEGER I,J
        CHARACTER*(*) ASTRNG

        DO 10 I=1,NBLNK(ASTRNG)
        J=ICHAR(ASTRNG(I:I))
        IF((J.GE.97).AND.(J.LE.122)) ASTRNG(I:I)=CHAR(J-32)
10      CONTINUE
        RETURN
        END subroutine



        SUBROUTINE SHIFTL(AA)

! -- SUBROUTINE SHIFTL REMOVES LEADING BLANKS FROM A STRING

        INTEGER L,I,J,II
        CHARACTER*(*) AA

        L=LEN(AA)
        DO 10 I=1,L
        IF((AA(I:I).NE.' ').AND.(ICHAR(AA(I:I)).NE.9)) GO TO 50
10      CONTINUE
        RETURN
50      IF(I.EQ.1) RETURN
        II=I-1
        DO 100 J=I,L
100     AA(J-II:J-II)=AA(J:J)
        DO 110 J=1,II
110     AA(L-J+1:L-J+1)=' '
        RETURN
        END subroutine
end module jacread
! -- function to determine the number of blanks at the end of a string
        INTEGER FUNCTION NBLNK(ASTRNG)
        CHARACTER*(*) ASTRNG 

        NBLNK=LEN_TRIM(ASTRNG)
        RETURN
        END FUNCTION NBLNK
        



