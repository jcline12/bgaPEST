	module jacread

        
        use utilities
        use bayes_pest_control
        
        implicit none
        type dmatrix
             integer                         :: nrow,ncol,icode
             double precision, pointer       :: vector(:)
             character*20, pointer           :: arow(:)
             character*20, pointer           :: acol(:)
        end type dmatrix
        
        public readJCO, readJAC
        
        contains 
        
        subroutine readJAC(jacfle,X)
        implicit none
! -- Program readJCO reads a ascii Jacobian "jac" file
! -- Code adapted from mat_read subroutine in matman.f code distributed with PEST
 
        integer jacunit
        integer               :: ifail       ! return as zero if an error
        character*(*)         :: jacfle      ! the jacobian matrix file
        character(500)        :: amessage    ! string to write error message to
        character(500)        :: cline       ! a character work string 
        type (dmatrix)        :: mat         ! the matrix to be read
        DOUBLE PRECISION, pointer :: X(:,:)  ! read the Jacobian into X rather than mat%array
        integer        :: ierr,ncol,nrow,icode,irow,icol
        integer        :: lw(3),rw(3)
        character*6 aarow
        character*200  :: afile

! -- Initialisation
       ifail=0
       jacunit = utl_nextunit()
       call addquote(jacfle,afile)
! -- The matrix file is opened.

       open(unit=jacunit,file=jacfle,status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,10) trim(afile)
10       format(' Cannot open jacobian file ',a,'.')
         go to 9800
       end if
       
! -- The header line is read.

       read(jacunit,'(a)',err=9000,end=9000) cline
       call linspl(ifail,3,lw,rw,cline)
       if(ifail.ne.0)then
         write(amessage,40) trim(afile)
40       format(' Three integers are expected on first line of file ', a,'.')
         go to 9800
       end if
       call intread(ifail,cline(lw(1):rw(1)),nrow)
       if(ifail.ne.0) go to 9000
       call intread(ifail,cline(lw(2):rw(2)),ncol)
       if(ifail.ne.0) go to 9000
       call intread(ifail,cline(lw(3):rw(3)),icode)
       if(ifail.ne.0) go to 9000
       if((ncol.le.0).or.(nrow.le.0))then
         write(amessage,50) trim(afile)
50       format(' NCOL or NROW is less than or equal to zero at ', &
         'first line of file ',a,'.')
         go to 9800
       end if
! -- for this case, ICODE must equal 2 --- no other options supported
       if(icode.ne.2)then
         write(amessage,70) trim(afile)
70       format(' ICODE must be "2" on first line of ', &
        'file ',a,'.')
         go to 9800
       end if
     
! -- Arrays in the matrix structure are dimensioned.

       mat%nrow=nrow
       mat%ncol=ncol
       mat%icode=icode
       allocate(X(nrow,ncol),stat=ierr)
       if(ierr.ne.0) go to 9400
       allocate(mat%arow(nrow),mat%acol(ncol),stat=ierr)
       if(ierr.ne.0) go to 9400
    
       
      
! -- The matrix is read.

       do irow=1,nrow
           read(jacunit,*,err=9100,end=9200) (X(irow,icol),icol=1,ncol)
         end do
         
! -- The row and column labels are read. NOTE - not currently used for anything
       
       read(jacunit,'(a)',err=9300,end=9300) cline
       call UTL_CASETRANS(cline, 'lo')
       if(index(cline,'* row names').eq.0)then
         write(amessage,130) trim(afile)
130      format(' "* row names" header expected immediately ', &
        'folowing matrix in file ',a,'.')
         go to 9800
       end if
! -- read row names        
       do irow=1,nrow
131      read(jacunit,*,err=9300,end=9300) mat%arow(irow)
         if(mat%arow(irow).eq.' ') go to 131
         mat%arow(irow)=adjustl(mat%arow(irow))
         call UTL_CASETRANS(mat%arow(irow),'lo')
       end do
! -- read column names
         read(jacunit,'(a)',err=9500,end=9500) cline
         call UTL_CASETRANS(cline,'lo')
         if(index(cline,'* column names').eq.0) go to 9500
         do icol=1,ncol
132        read(jacunit,*,err=9300,end=9300) mat%acol(icol)
           if(mat%acol(icol).eq.' ') go to 132
           mat%acol(icol)=adjustl(mat%acol(icol))
           call UTL_CASETRANS(mat%acol(icol),'lo')
         end do

       close(unit=jacunit)
       return       
       
9000   write(amessage,9010) trim(afile)
9010   format(' Error reading integer matrix header line from first ', &
      'line of file ',a,'.')
       go to 9800
9100   write(aarow,'(i6)') irow
       aarow=adjustl(aarow)
       write(amessage,9110) trim(aarow),trim(afile)
9110   format(' Error reading matrix row number ',a,' from file ',a,'.')
       go to 9800
9200   write(amessage,9210) trim(afile)
9210   format(' Unexpected end encountered to file ',a,' while ', &
       'reading matrix.')
       go to 9800
9300   write(amessage,9310) trim(afile)
9310   format(' Error reading row and/or column names from matrix ', &
       'file ',a,'.')
       go to 9800
9400   write(amessage,9410) trim(afile)
9410   format(' Cannot allocate sufficient memory to hold matrix ', &
       'located in file ',a,'.')
       go to 9800
9500   write(amessage,9510) trim(afile)
9510   format(' "* column names" header expected immediately ', &
       'following row names in file ',a,'.')
       go to 9800

9800   ifail=1
       close(unit=jacunit,iostat=ierr)
       return
        
        
        end subroutine readJAC
        
        
        subroutine readJCO(jacfle, X)
        IMPLICIT NONE

! -- Program readJCO reads a binary Jacobian "jco" file 

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
        

        subroutine linspl(ifail,num,lw,rw,cline)

! -- Subroutine LINSPL splits a line into whitespace-separated substrings.

        integer ifail,nw,nblc,j,i
        integer num
        integer lw(num),rw(num)
        character*(*) cline

        ifail=0
        nw=0
        nblc=len_trim(cline)
        if(nblc.eq.0) then
          ifail=1
          return
        end if
        j=0
5       if(nw.eq.num) return
        do 10 i=j+1,nblc
        if((cline(i:i).ne.' ').and.(cline(i:i).ne.',') &
        .and.(ichar(cline(i:i)).ne.9)) go to 20
10      continue
        ifail=1
        return
20      nw=nw+1
        lw(nw)=i
        do 30 i=lw(nw)+1,nblc
        if((cline(i:i).eq.' ').or.(cline(i:i).eq.',') &
        .or.(ichar(cline(i:i)).eq.9)) go to 40
30      continue
        rw(nw)=nblc
        if(nw.lt.num) ifail=1 
        return
40      rw(nw)=i-1
        j=rw(nw)
        go to 5
 
        end subroutine LINSPL


       subroutine addquote(afile,aqfile)

! -- Subroutine ADDQUOTE adds quotes to a filename if it has a space in it.

         character (len=*), intent(in)  :: afile
         character (len=*), intent(out) :: aqfile
         integer nbb

         if(index(trim(afile),' ').eq.0)then
           aqfile=afile
         else
           aqfile(1:1)='"'
           aqfile(2:)=trim(afile)
           nbb=len_trim(aqfile)+1
           aqfile(nbb:nbb)='"'
         end if

         return

       end subroutine addquote