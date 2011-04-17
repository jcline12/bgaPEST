C -- THE SENSITIVITY MATRIX IS NEXT CALCULATED

        J=0
        ISENREUSEFLAG=0
        DO 145 IPP=1,NPAR
        IF(ITRANS(IPP).LT.0) GO TO 145
        IF((ICOVOBS.EQ.0).AND.(JSTK(IPP).LT.0)) GO TO 145
        IF(SCALE(IPP).LT.-1.0D35) GO TO 145
        IF(ISENREUSE.NE.0)THEN
          IF(SCREUSE(IPP).GE.0.0D0) THEN
            ISENREUSEFLAG=ISENREUSEFLAG+1
            GO TO 145
          END IF
        END IF
        I=IPARGP(IPP)
        IF(FORCEN(I).EQ.2) THEN
          J=J+1
        ELSE IF(FORCEN(I).EQ.3) THEN
          J=J+2
        ELSE IF(FORCEN(I).EQ.1) THEN
          J=J+I2OR3-1
        ELSE IF(FORCEN(I).EQ.4)THEN
          IF(I2OR3.EQ.2)THEN
            J=J+1
          ELSE IF(I2OR3.EQ.3)THEN
            J=J+4
          END IF
        ELSE IF(FORCEN(I).EQ.5)THEN
          J=J+4
        END IF
        CALL WRITINT(ATEMP1,J)
145     CONTINUE
        IF((SVDA_SUPDERCALC.EQ.1).AND.(IOPT.EQ.1))THEN
          WRITE(6,1512)
1512      FORMAT  (/,'    Calculating Jacobian matrix from ',
     +    'base parameter sensitivities.....')
        ELSE
          IF(JACFILE.EQ.0)THEN
            CALL WRITINT(ATEMP1,J)
            WRITE(6,151) TRIM(ATEMP1)
151         FORMAT(/,'    Calculating Jacobian matrix: running ',
     +      'model ',A,' times .....')
            IF(ISENREUSEFLAG.NE.0)THEN
              WRITE(6,1515)
1515          FORMAT('    (Sensitivity re-use operational this ',
     +        'iteration.)')
              WRITE(IREC,15151)
15151         FORMAT(/,4X,'Sensitivity re-use is operational during ',
     +        'this iteration.')
              WRITE(IREC,15152) TRIM(ATEMP1)
15152         FORMAT(4X,'Model is run ',a,' times for derivatives ',
     +        'calculation.')
            END IF
          ELSE
            WRITE(6,1511)
1511        FORMAT(/,'    Calculating Jacobian matrix .....')
          END IF
        END IF
        JTIMES=J
#ifdef PARALLEL
#ifndef MPEST
        WRITE(IRMR,151) TRIM(ATEMP1)
#endif
#endif
        NRUN=J
        IF((IIRST.EQ.2).OR.(IRST1JAC.EQ.1))THEN
          IF(IIRST.EQ.2)THEN
            WRITE(6,152)
152         FORMAT('    Jacobian already calculated.')
#ifdef PARALLEL
#ifndef MPEST
            WRITE(IRMR,152)
#endif
#endif
          END IF
          FLENME=CASEFL(:LEN_TRIM(CASEFL))//'.jac'
          CALL FFOPEN(JFAIL,-IRSF,'r',' ',131,CLINE)
          IF(JFAIL.NE.0) GO TO 9891
          READ(IRSF,ERR=9750,END=9750) ITEMP,JJPRED,JJREG
          IF(JJPRED.NE.IPRED) THEN
            IF(JJPRED.EQ.0) THEN
              CALL STPERR(139,0,' ',0,' ',CLINE)
              GO TO 9891
            END IF
            IF(JJPRED.NE.0) THEN
              CALL STPERR(138,0,' ',0,' ',CLINE)
              GO TO 9891
            END IF
          END IF
          IF(JJREG.NE.IREG)THEN
            IF(JJREG.EQ.0) THEN
              CALL STPERR(149,0,' ',0,' ',CLINE)
              GO TO 9891
            END IF
            IF(JJREG.NE.0) THEN
              CALL STPERR(148,0,' ',0,' ',CLINE)
              GO TO 9891
            END IF
          END IF
          IF(IIRST.EQ.2)THEN
            IF(ITEMP.NE.IOPT)THEN
              CALL STPERR(127,1,' ',0,' ',CLINE)
              GO TO 9891
            END IF
          END IF
          IF(MAXCOMPDIM.LE.1)THEN
            X=0.0D0                     ! AN ARRAY
            READ(IRSF,ERR=9750,END=9750) ICOUNT
            IF(ICOUNT.NE.0)THEN
              DO I=1,ICOUNT
                READ(IRSF,ERR=9750,END=9750) IROW,IES,X(IROW,IES)
              END DO
            END IF
          ELSE
            READ(IRSF,ERR=9750,END=9750)NCOMPDIM
            IF(NCOMPDIM.GE.MAXCOMPDIM-3) GO TO 9970
            IF(NCOMPDIM.NE.0)THEN
              DO I=1,NCOMPDIM
                READ(IRSF,ERR=9750,END=9750) IXC(I),XC(I)
              END DO
            END IF
          END IF
          READ(IRSF,ERR=9750,END=9750) WFSOL
          CLOSE(UNIT=IRSF)

C -- THE SENSITIVITY MATRIX IS NOW ROTATED.
C    HOWEVER THIS IS ONLY DONE IF OBSERVATION COVARIANCE MATRICES ARE PROVIDED.

          IF(NUMCOV.NE.0)THEN
            MATDIM1=NXROW
            NM=MXOBSCOV
            IPSTART=1
            DO 6913 JCOV=1,NUMCOV
              IGROUP=COVGP(JCOV)
              IF(NPRIOR.NE.0)THEN
                DO 6914 I=NOBS+1,NXROW
                  IF(NOBGNM(I).EQ.IGROUP) GO TO 6913
6914            CONTINUE
              END IF
              DO 6915 I=IPSTART,NXROW
                IF(NOBGNM(I).EQ.IGROUP)THEN
                  CALL ROTATE(2,NOBSCOV(JCOV),NOBSCOV(JCOV),MATDIM1,
     +            NESPAR,IPSTART,COVAR(JCOV)%DVAL,W1,X,W2,IGROUP,
     +            NXROW,NOBGNM)
                  GO TO 6913
                END IF
6915          CONTINUE
6913        CONTINUE
          END IF

          JIRST=1
          IRST1JAC=0
          IF(IIRST.EQ.2) THEN
            IIRST=0
            IF(IREG.EQ.1)THEN
C            IF((IREG.EQ.1).AND.(IREGADJ.NE.0))THEN
              GO TO 6907
            ELSE
              GO TO 302
            END IF
          END IF
        END IF
C -- Note that the above code is repeated below.

C -- If we are using SVD assist, a new parcalc template file is written.

        IF((SVDA.EQ.1).AND.(IOPT.GT.1))THEN
          IF(SVDA_MULBPA.EQ.0)THEN
            CALL BASE_PARDEF(IFAIL,IREC,0,NREGADJPAR,-1,
     +      SVDA_SCALADJ)
          ELSE
            CALL BASE_PARDEF(IFAIL,IREC,0,NREGADJPAR,IOPT-1,
     +      SVDA_SCALADJ)
          END IF
          IF(IFAIL.NE.0) GO TO 9890
        END IF

C -- Now derivatives are calculated.

        IES=0
        IRUN=0
        IIRUN=0
        ISTART_C=1
        DERFLAG=0

C -- If we start with the "/i" switch, we now read the nominated JCO file.

        IF((RESTART.EQ.5).AND.(IOPT.EQ.1))THEN
          CALL READ_JCO_FIRST_FILE(IFAIL)
          IF(IFAIL.NE.0) GO TO 9890
          GO TO 2013
        END IF

C -- If we are using SVDA and super parameter derivatives can be calculated
C    from base parameter derivatives, this is now done.

        IF(SVDA.NE.0)THEN
          IF(SVDA_SUPDERCALC.NE.0)THEN
            IF(IOPT.EQ.1)THEN
              CALL BASE_SUPDERCALC(IFAIL,SVDA_SCALADJ,1)
              IF(IFAIL.EQ.1) THEN
                 GO TO 9970
              ELSE IF(IFAIL.NE.0)THEN
                GO TO 9890
              END IF
              GO TO 2013
            END IF
          END IF
          IF(JACFILE.NE.0)THEN
            IF(IOPT.GT.1)THEN
#ifndef PARALLEL
              NXDIM=NXROW
              TEMPJACDIM=MAX(NXROW,NPAR)
              IF(JACFILE.EQ.2)THEN
                 DIM_DI_PAR=DI_NPAR
                 DIM_DI_OBS=DI_NDEP
                 DIM_DI_IES=DI_NPAR
                 DIM_DI_PARENT=DI_NPAR
              ELSE
                 DIM_DI_PAR=1
                 DIM_DI_OBS=1
                 DIM_DI_IES=1
                 DIM_DI_PARENT=1
              END IF
              IEXTJACMODE=1
              CALL EXTJAC(JFAIL,NXDIM,NOBS,NESPAR,NPAR,
     +        DIM_DI_PAR,DIM_DI_OBS,DIM_DI_IES,DIM_DI_PARENT,
     +        PVAL,ORGVAL,ITRANS,
     +        SCALE,W2,X,COMJAC,EXTDERFLE,CLINE,APAR,AOBS,DERCOM,PRECIS,
     +        NOPNT,NTPLFLE,NW,OFFSET,PARDEL,PWORD,INFLE,TEMPFLE,
     +        MESSFILE,
     +        DI_NSKIP,DI_NDEP,DI_NPAR,DI_ORIENTATION,DI_DERFORMAT,
     +        DI_PAR,
     +        DI_IES,DI_PARENT,DI_OBS,JACFILE,TEMPJACDIM,NUMTIED,W1,
     +        IEXTJACMODE)
              IF(JFAIL.EQ.1) GO TO 9891
              IF(JFAIL.EQ.2) GO TO 9890
#endif

              CALL PROCESS_BASE_PAR_DERIV(IFAIL)
              IF(IFAIL.NE.0) GO TO 9890
              CALL BASE_SUPDERCALC(IFAIL,SVDA_SCALADJ,2)
              IF(IFAIL.EQ.1) THEN
                GO TO 9970
              ELSE IF(IFAIL.NE.0)THEN
                GO TO 9890
              END IF
              GO TO 2013
            END IF
          END IF
        END IF

C -- If necessary the derivatives restart file is opened.

        IDCOUNT=0
        IF(DRST.EQ.1)THEN
          DIRESTART=2
        ELSE
          DIRESTART=1
        END IF
        IF(RSTFLE.EQ.0) DIRESTART=0
        IF(IPRLL.NE.0) DIRESTART=0
        IF(DIRESTART.EQ.2)THEN
          INQUIRE(FILE=DTFILE,EXIST=LEXIST)
          IF(.NOT.LEXIST)THEN
            DIRESTART=1
            DRST=0
          ELSE
            JPERR=1
            CALL FFOPEN(JFAIL,-DTUNIT,'r',DTFILE,22,CLINE)
            JPERR=0
            IF(IPERR.NE.0)THEN
              DIRESTART=1
              DRST=0
            ELSE
              ITEMP=0
              READ(DTUNIT,IOSTAT=IERR) ITEMP
              IF(ITEMP.NE.IOPT)THEN
                CLOSE(UNIT=DTUNIT,IOSTAT=IERR)
                DIRESTART=1
                DRST=0
              END IF
            END IF
          END IF
        END IF
        IF(DIRESTART.EQ.1)THEN
          CALL DELFILE1(JFAIL,DTFILE,CLINE)
          IF(JFAIL.NE.0) GO TO 9891
          JPERR=1
          CALL FFOPEN(JFAIL,-DTUNIT,'w',DTFILE,25,CLINE)
          JPERR=0
          IF(IPERR.NE.0)THEN
            DRST=0
            DIRESTART=0
            GO TO 217
          END IF
          WRITE(DTUNIT,IOSTAT=IERR) IOPT
#ifdef FLUSHFILE
          CALL FLUSH(DTUNIT)
#endif
217       CONTINUE
        END IF

C -- Now derivatives are calculated.

C -- But first, if we are reading external derivatives, we zero all existing derivatives
C    except for prior information.
        IF(JACFILE.NE.0)THEN
          IF(MAXCOMPDIM.GT.1)THEN
            CALL ZERO_MATRIX(IFAIL,NCOMPDIM,XC,IXC,NOBS)
          ELSE
            DO IES=1,NESPAR
              DO IROW=1,NOBS
                X(IROW,IES)=0.0D0
              END DO
            END DO
          END IF
        END IF

        IES=0
        DO 200 IPP=1,NPAR
        IF(ITRANS(IPP).LT.0) GO TO 200
        IES=IES+1
        IF((JACFILE.NE.0).AND.(DERCOM(IPP).EQ.0)) GO TO 200
        IF((ICOVOBS.EQ.0).AND.(JSTK(IPP).LT.0)) GO TO 200
        IF(SCALE(IPP).LT.-1.0D35)THEN
          ITEMP=NINT(OFFSET(IPP))
          DO IOBS=1,NXROW
C Note - we are assuming that no
C iw_ variable is written to a template file (and hence affects any observations)
C and thus that it is only cited in prior information.
            IF(NOBGNM(IOBS).EQ.ITEMP)THEN
              IF(ITRANS(IPP).EQ.0)THEN
                IF(MAXCOMPDIM.LE.1)THEN
                  X(IOBS,IES)=-1.0/PVAL(IPP)/PVAL(IPP)*
     +            (REFOBS(IOBS)-OVAL(IOBS))
                ELSE
                  RTEMP=-1.0/PVAL(IPP)/PVAL(IPP)*
     +            (REFOBS(IOBS)-OVAL(IOBS))
                  CALL STORE_VALUE(IFAIL,NCOMPDIM,XC,IXC,RTEMP,
     +            IOBS,IES)
                  IF(IFAIL.NE.0) GO TO 9970
                  ISTART_C=IFOUND_C
                END IF
              ELSE
                IF(MAXCOMPDIM.LE.1)THEN
                  X(IOBS,IES)=-2.30259/PVAL(IPP)*
     +            (REFOBS(IOBS)-OVAL(IOBS))
                ELSE
                  RTEMP=-2.30259/PVAL(IPP)*
     +            (REFOBS(IOBS)-OVAL(IOBS))
                  CALL STORE_VALUE(IFAIL,NCOMPDIM,XC,IXC,RTEMP,
     +            IOBS,IES)
                  IF(IFAIL.NE.0) GO TO 9970
                  ISTART_C=IFOUND_C
                END IF
              END IF
            ELSE
              IF(MAXCOMPDIM.LE.1)THEN
                IF(IOBS.LE.NOBS)THEN
                  X(IOBS,IES)=0.0D0
                END IF
              END IF
            END IF
          END DO
          GO TO 200
        END IF

C -- Now we check whether these derivatives are available in a derivatives
C    restart file.

        IF(DIRESTART.EQ.2)THEN
194       CONTINUE
          READ(DTUNIT,END=197,ERR=197) ITEMP,ITEMP1,ITEMP2
          IF(ITEMP.GT.IES) GO TO 197
          IF(MAXCOMPDIM.LE.1)THEN
            READ(DTUNIT,END=197,ERR=197) (X(IOBS,IES),IOBS=1,NOBS)
            IF(ITEMP.LT.IES) GO TO 194
            IDCOUNT=IDCOUNT+1
            NCALL=ITEMP1
            IIRUN=ITEMP2
            GO TO 200
          ELSE
            READ(DTUNIT,END=197,ERR=197) (WORK_C(IOBS),IOBS=1,NOBS)
            IF(ITEMP.LT.IES) GO TO 194
            ISTART_C=1
            DO IOBS=1,NOBS
              CALL STORE_VALUE(IFAIL,NCOMPDIM,XC,IXC,WORK_C(IOBS),
     +        IOBS,IES)
              IF(IFAIL.NE.0) GO TO 9970
              ISTART_C=IFOUND_C
            END DO
            IDCOUNT=IDCOUNT+1
            NCALL=ITEMP1
            IIRUN=ITEMP2
            GO TO 200
          END IF
        END IF

196     CONTINUE

        IF(ISENREUSE.NE.0)THEN
          IF(SCREUSE(IPP).GE.0.0D0) GO TO 1961
        END IF

        IF(MAXCOMPDIM.LE.1)THEN
          CALL DERCLC(JFAIL,DERFLAG,NSCALL,MMCALL,NSCALLP,MMCALLP,
     +    SOPDIM,I2OR3,IPP,NPAR,NOBS,NPARGP,ASIZE,NINSTR,NTPLFLE,
     +    NINSFLE,NUML,NBLBMX,PRECIS,NOPNT,PVAL,PARLBND,PARUBND,ITRANS,
     +    SCALE,OFFSET,IPARGP,W2,NW,PWORD,APAR,AOBS,DERINC,IDBND,
     +    INCTYP,DERINCLB,FORCEN,DERMTHD,DERINCMUL,REFOBS,
     +    X(1,IES),W1,LCINS,LL,OBSN1,OBSN2,IIOBS,INSFLE,TEMPFLE,INFLE,
     +    OUTFLE,PARDEL,MRKDEL,A,CLINE,BUF,COMLIN(DERCOM(IPP)),
     +    NRUN,IRUN,1,INCPAR,
     +    MESSFILE,IIRUN,NORMRETURN,WV5DIM,WORKVEC5,SPLITFLAG)
          IF(JFAIL.NE.0) GO TO 9891
        ELSE
          IF((SPLITFLAG.NE.0).AND.(JSPLIT.NE.0))THEN
            IF((SPLITTHRESH(IPARGP(IPP)).NE.0.0D0).AND.
     +         (SPLITACTION(IPARGP(IPP)).EQ.3))THEN
              IF(NCOMPDIM.EQ.0)THEN
                DO I=1,NOBS
                  WORK_C(I)=0.0D0
                END DO
              ELSE
                ISTART_C=1
                CALL GET_VECTOR(NCOMPDIM,NXROW,XC,IXC,WORK_C,IES)
                ISTART_C=IFOUND_C
              END IF
            END IF
          END IF
          CALL DERCLC(JFAIL,DERFLAG,NSCALL,MMCALL,NSCALLP,MMCALLP,
     +    SOPDIM,I2OR3,IPP,NPAR,NOBS,NPARGP,ASIZE,NINSTR,NTPLFLE,
     +    NINSFLE,NUML,NBLBMX,PRECIS,NOPNT,PVAL,PARLBND,PARUBND,ITRANS,
     +    SCALE,OFFSET,IPARGP,W2,NW,PWORD,APAR,AOBS,DERINC,IDBND,
     +    INCTYP,DERINCLB,FORCEN,DERMTHD,DERINCMUL,REFOBS,
     +    WORK_C,W1,LCINS,LL,OBSN1,OBSN2,IIOBS,INSFLE,TEMPFLE,INFLE,
     +    OUTFLE,PARDEL,MRKDEL,A,CLINE,BUF,COMLIN(DERCOM(IPP)),
     +    NRUN,IRUN,1,INCPAR,
     +    MESSFILE,IIRUN,NORMRETURN,WV5DIM,WORKVEC5,SPLITFLAG)
          IF(JFAIL.NE.0) GO TO 9891
          ISTART_C=1
          DO IOBS=1,NOBS
            CALL STORE_VALUE(IFAIL,NCOMPDIM,XC,IXC,WORK_C(IOBS),
     +      IOBS,IES)
            IF(IFAIL.NE.0) GO TO 9970
            ISTART_C=IFOUND_C
          END DO
        END IF

1961    CONTINUE

C -- If necessary, values are stored in the mid-run restart file.

        IF(DIRESTART.EQ.1)THEN
          IF(NORMRETURN.EQ.1)THEN
            WRITE(DTUNIT,IOSTAT=IERR) IES,NCALL,IIRUN
            IF(MAXCOMPDIM.LE.1)THEN
              WRITE(DTUNIT,IOSTAT=IERR) (X(IOBS,IES),IOBS=1,NOBS)
            ELSE
              WRITE(DTUNIT,IOSTAT=IERR) (WORK_C(IOBS),IOBS=1,NOBS)
            END IF
#ifdef FLUSHFILE
            CALL FLUSH(DTUNIT)
#endif
          END IF
        END IF

C -- Has a stop condition been encountered?

        IF(ISTOP.EQ.2)THEN
          IF(DIRESTART.NE.0) CLOSE(UNIT=DTUNIT,IOSTAT=IERR)
          IFIN=10
          GO TO 6000
        ELSE IF(ISTOP.EQ.1) THEN
          IF(DIRESTART.NE.0) CLOSE(UNIT=DTUNIT,IOSTAT=IERR)
          IPFAIL=-1
          GO TO 9891
        END IF
        GO TO 200

197     CONTINUE
        CALL WRITINT(ATEMP3,IDCOUNT)
        WRITE(6,198) TRIM(ATEMP3)
198     FORMAT('    Derivatives for ',A,' parameters found ',
     + 'in restart file.')
        NOWRITE=0
        IF(IDCOUNT.EQ.NESPAR) NOWRITE=1
#ifdef INTEL
        IF((IDCOUNT.NE.0).AND.(IDCOUNT.NE.NESPAR))THEN
          WRITE(6,2561)
2561      FORMAT('    - number of runs completed...')
          WRITE(6,'(a)',ADVANCE='NO') '     '
        END IF
#endif
        CLOSE(UNIT=DTUNIT,STATUS='DELETE',IOSTAT=IERR)
        IF(IERR.NE.0)THEN
          DIRESTART=0
          DRST=0
        ELSE
          DIRESTART=1
          JPERR=1
          CALL FFOPEN(JFAIL,-DTUNIT,'w',DTFILE,25,CLINE)
          JPERR=0
          IF(IPERR.NE.0) THEN
            DIRESTART=0
            DRST=0
            GO TO 218
          END IF
          WRITE(DTUNIT,ERR=195) IOPT
          IF(IES.GT.1)THEN
            DO I=1,IES-1
              WRITE(DTUNIT,ERR=195) I,NCALL,IIRUN
              IF(MAXCOMPDIM.LE.1)THEN
                WRITE(DTUNIT,ERR=195) (X(IOBS,I),IOBS=1,NOBS)
              ELSE
                CALL GET_VECTOR(NCOMPDIM,NXROW,XC,IXC,WORK_C,I)
                WRITE(DTUNIT,ERR=195) (WORK_C(IOBS),IOBS=1,NOBS)
              END IF
            END DO
          END IF
#ifdef FLUSHFILE
          CALL FLUSH(DTUNIT)
#endif
218       CONTINUE
          GO TO 196
        END IF
195     CLOSE(UNIT=DTUNIT,IOSTAT=IERR)
        DIRESTART=0
        DRST=0
        GO TO 196

200     CONTINUE

C -- Close the derivatives restart file if necessary.

        IF(DIRESTART.NE.0) THEN
          INQUIRE(UNIT=DTUNIT,OPENED=LOPENED)
          IF(LOPENED) THEN
            CLOSE(UNIT=DTUNIT,IOSTAT=IERR)
          END IF
        END IF
        IF(DIRESTART.EQ.2)THEN
          CALL WRITINT(ATEMP3,IDCOUNT)
          WRITE(6,198) TRIM(ATEMP3)
          NOWRITE=0
          IF(IDCOUNT.EQ.NESPAR) NOWRITE=1
#ifdef INTEL
        IF((IDCOUNT.NE.0).AND.(IDCOUNT.NE.NESPAR))THEN
          WRITE(6,2561)
          WRITE(6,'(a)',ADVANCE='NO') '     '
        END IF
#endif
        END IF
        DRST=0

#ifdef PARALLEL
        IF(SRST.EQ.1)THEN
          PIRESTART=2
        ELSE IF(SRST.EQ.0)THEN
          PIRESTART=1
        END IF
        IF(RSTFLE.EQ.0) PIRESTART=0
#ifdef BEO
       if (BEOMASTER) then
         call RUNMASTER(PARREG,OBSREG,NRUN,JFAIL,                               !jd
     +   pirestart,ptunit,iopt,ptfile,workvec7dim,workvec7)                     !jd
         if (JFAIL.ne.0) goto 9891
       else
#endif
#ifdef MPEST
        FLENME=TRIM(CASEFL)//'.jacobian_runs'
        CALL FFOPEN(JFAIL,IRMR,'w',' ',6,CLINE)
        IF(JFAIL.NE.0) GO TO 9891
        WRITE(IRMR,'(I5)') IOPT
        RTEMPLAM=LAMLST
        IF(LAMFLAG.EQ.0) RTEMPLAM=-1.0D0
        WRITE(IRMR,'(1X,1PG14.7)') RTEMPLAM
        WRITE(IRMR,'(1X,I5,2X,I5)') SWITCHFLAG,SBACKFLAG
        CLOSE(UNIT=IRMR)
        IF(SWITCHFLAG.EQ.1)SWITCHFLAG=0
        IF(SBACKFLAG.EQ.1)SBACKFLAG=0
        FLENME=TRIM(CASEFL)//'.lambda_runs'
        CALL DELFILE1(JFAIL,FLENME,CLINE)
        CALL DORUNS_M(JFAIL,NSLAVE,NRUN,
     +  ASLDIR,NINSTR,NINSFLE,ASIZE,NUML,NOBS,NBLBMX,LCINS,LL,OBSN1,
     +  OBSN2,IIOBS,AOBS,A,MRKDEL,CLINE,BUF,AFILE,
     +  NPAR,PRECIS,NOPNT,NTPLFLE,NW,SCALE,
     +  OFFSET,PARDEL,PWORD,TEMPFLE,APAR,
     +  OUTFLE,INSFLE,LDSIN,LDSOU,SINFLE,
     +  SOUFLE,ITRIAL,MESSFILE,INCPAR,ITRANS,MSINFLE,MSOUFLE,
     +  MODFLE,WORKVEC7DIM,WORKVEC7)
        IF(JFAIL.EQ.2) GO TO 9890
        IF(JFAIL.NE.0) GO TO 9891
#else
        CALL DORUNS(JFAIL,NSLAVE,NRUN,ISTATS,ISTATR,
     +  ASLDIR,NINSTR,NINSFLE,ASIZE,NUML,NOBS,NBLBMX,LCINS,LL,OBSN1,
     +  OBSN2,IIOBS,AOBS,A,MRKDEL,CLINE,BUF,AFILE,ISTRTME,
     +  IRUNTME,JRUN,NPAR,PRECIS,NOPNT,NTPLFLE,NW,SCALE,
     +  OFFSET,PARDEL,PWORD,TEMPFLE,APAR,OREADFLE,PREADFLE,
     +  MANFLE,ASLAVE,NNRUN,OUTFLE,INSFLE,LDSIN,LDSOU,SINFLE,
     +  SOUFLE,ITRIAL,MESSFILE,INCPAR,ITRANS,IDET,SREADFLE,MREADFLE,
     +  SCOM,PIRESTART,PTUNIT,PTFILE,IOPT,REPEATRUN,0,SLAVEGROUP,
     +  WORKVEC7DIM,WORKVEC7)
        IF(JFAIL.NE.0) GO TO 9891
#endif
#ifdef BEO
       endif
#endif
        IF(ISTOP.EQ.2) THEN
          IFIN=10
          GO TO 6000
        ELSE IF(ISTOP.EQ.1) THEN
          IPFAIL=-1
          GO TO 9891
        END IF
        SRST=0
        IES=0
        IRUN=0
        DO 201 IPP=1,NPAR
        IF(ITRANS(IPP).LT.0) GO TO 201
        IF(SCALE(IPP).LT.-1.0D35) GO TO 201
        IES=IES+1
        IF((ICOVOBS.EQ.0).AND.(JSTK(IPP).LT.0)) GO TO 201
        IF(ISENREUSE.NE.0)THEN
          IF(SCREUSE(IPP).GE.0.0D0) GO TO 201
        END IF
        IF(MAXCOMPDIM.LE.1)THEN
          CALL DERCLC(JFAIL,DERFLAG,NSCALL,MMCALL,NSCALLP,MMCALLP,
     +    SOPDIM,I2OR3,IPP,NPAR,NOBS,NPARGP,ASIZE,NINSTR,NTPLFLE,
     +    NINSFLE,NUML,NBLBMX,PRECIS,NOPNT,PVAL,PARLBND,PARUBND,
     +    ITRANS,SCALE,OFFSET,IPARGP,W2,NW,PWORD,APAR,AOBS,DERINC,
     +    IDBND,INCTYP,DERINCLB,FORCEN,DERMTHD,DERINCMUL,REFOBS,
     +    X(1,IES),W1,LCINS,LL,OBSN1,OBSN2,IIOBS,INSFLE,TEMPFLE,INFLE,
     +    OUTFLE,PARDEL,MRKDEL,A,CLINE,BUF,MODFLE,
     +    NRUN,IRUN,2,INCPAR,
     +    MESSFILE,IIRUN,NORMRETURN,WV5DIM,WORKVEC5,SPLITFLAG)
          IF(JFAIL.NE.0) GO TO 9891
        ELSE
          IF((SPLITFLAG.NE.0).AND.(JSPLIT.NE.0))THEN
            IF((SPLITTHRESH(IPARGP(IPP)).NE.0.0D0).AND.
     +         (SPLITACTION(IPARGP(IPP)).EQ.3))THEN
              IF(NCOMPDIM.EQ.0)THEN
                DO I=1,NOBS
                  WORK_C(I)=0.0D0
                END DO
              ELSE
                ISTART_C=1
                CALL GET_VECTOR(NCOMPDIM,NXROW,XC,IXC,WORK_C,IES)
                ISTART_C=IFOUND_C
              END IF
            END IF
          END IF
          CALL DERCLC(JFAIL,DERFLAG,NSCALL,MMCALL,NSCALLP,MMCALLP,
     +    SOPDIM,I2OR3,IPP,NPAR,NOBS,NPARGP,ASIZE,NINSTR,NTPLFLE,
     +    NINSFLE,NUML,NBLBMX,PRECIS,NOPNT,PVAL,PARLBND,PARUBND,
     +    ITRANS,SCALE,OFFSET,IPARGP,W2,NW,PWORD,APAR,AOBS,DERINC,
     +    IDBND,INCTYP,DERINCLB,FORCEN,DERMTHD,DERINCMUL,REFOBS,
     +    WORK_C,W1,LCINS,LL,OBSN1,OBSN2,IIOBS,INSFLE,TEMPFLE,INFLE,
     +    OUTFLE,PARDEL,MRKDEL,A,CLINE,BUF,MODFLE,
     +    NRUN,IRUN,2,INCPAR,
     +    MESSFILE,IIRUN,NORMRETURN,WV5DIM,WORKVEC5,SPLITFLAG)
          IF(JFAIL.NE.0) GO TO 9891
          ISTART_C=1
          DO IOBS=1,NOBS
            CALL STORE_VALUE(IFAIL,NCOMPDIM,XC,IXC,WORK_C(IOBS),
     +      IOBS,IES)
            IF(IFAIL.NE.0) GO TO 9970
            ISTART_C=IFOUND_C
          END DO
        END IF
201     CONTINUE
#else

#ifdef LAHEY
        IF(NOWRITE.EQ.0) WRITE(6,*)
#endif

#endif

#ifndef PARALLEL
        IF(JACFILE.NE.0)THEN
          NXDIM=NXROW
          TEMPJACDIM=MAX(NXROW,NPAR)
          IF(JACFILE.EQ.2)THEN
             DIM_DI_PAR=DI_NPAR
             DIM_DI_OBS=DI_NDEP
             DIM_DI_IES=DI_NPAR
             DIM_DI_PARENT=DI_NPAR
           ELSE
             DIM_DI_PAR=1
             DIM_DI_OBS=1
             DIM_DI_IES=1
             DIM_DI_PARENT=1
           END IF
          IEXTJACMODE=0
          CALL EXTJAC(JFAIL,NXDIM,NOBS,NESPAR,NPAR,
     +    DIM_DI_PAR,DIM_DI_OBS,DIM_DI_IES,DIM_DI_PARENT,
     +    PVAL,ORGVAL,ITRANS,
     +    SCALE,W2,X,COMJAC,EXTDERFLE,CLINE,APAR,AOBS,DERCOM,PRECIS,
     +    NOPNT,NTPLFLE,NW,OFFSET,PARDEL,PWORD,INFLE,TEMPFLE,MESSFILE,
     +    DI_NSKIP,DI_NDEP,DI_NPAR,DI_ORIENTATION,DI_DERFORMAT,DI_PAR,
     +    DI_IES,DI_PARENT,DI_OBS,JACFILE,TEMPJACDIM,NUMTIED,W1,
     +    IEXTJACMODE)
          IF(JFAIL.EQ.1) GO TO 9891
          IF(JFAIL.EQ.2) GO TO 9890
        END IF
c        note, in the above, well have to give the "jacobian model" information
c        on which parameters are tied and fixed.

#endif

! -- Set a flag in the PESTDATA module saying that sensitivities are now available.

2013    CONTINUE

! -- Do we need to cease or pause execution?

        CALL STOPRESS(0)
        IF(ISTOP.EQ.2)THEN
          IFIN=10
          GO TO 6000
        ELSE IF(ISTOP.EQ.1) THEN
          IPFAIL=-1
          GO TO 9891
        END IF

        SENFLAG=1

! -- Save restart information.

        IF(RSTFLE.NE.0)THEN
          FLENME=CASEFL(:LEN_TRIM(CASEFL))//'.jac'
          CALL FFOPEN(JFAIL,-IRSF,'w',' ',25,CLINE)
          IF(JFAIL.NE.0) GO TO 9891
          WRITE(IRSF) IOPT,IPRED,IREG
          IF(MAXCOMPDIM.LE.1)THEN
            ICOUNT=0
            DO IES=1,NESPAR
              DO IROW=1,NXROW
                IF(X(IROW,IES).NE.0.0D0) ICOUNT=ICOUNT+1
              END DO
            END DO
            WRITE(IRSF) ICOUNT
            IF(ICOUNT.GT.0)THEN
              DO IES=1,NESPAR
                DO IROW=1,NXROW
                  IF(X(IROW,IES).NE.0.0D0) WRITE(IRSF) IROW,IES,
     +               X(IROW,IES)
                END DO
              END DO
            END IF
          ELSE
            WRITE(IRSF) NCOMPDIM
            DO I=1,NCOMPDIM
              WRITE(IRSF) IXC(I),XC(I)
            END DO
          END IF
          WRITE(IRSF) WFSOL
        END IF

C -- The JCO file is written if this is the first optimisation iteration.

        IF((NOPTMAX.EQ.-1).OR.(NOPTMAX.EQ.-2).OR.(IOPT.EQ.1))THEN

C -- Prior information is un-rotated.

          IF(NUMCOV.NE.0)THEN
            IF(NPRIOR.NE.0)THEN
              MATDIM1=NXROW
              NM=MXOBSCOV
              IPSTART=NOBS+1
              DO 6899 JCOV=1,NUMCOV
                IGROUP=COVGP(JCOV)
                DO I=IPSTART,NXROW
                  IF(NOBGNM(I).EQ.IGROUP)THEN
                    CALL ROTATE(-2,NOBSCOV(JCOV),NOBSCOV(JCOV),MATDIM1,
     +              NESPAR,IPSTART,COVAR(JCOV)%DVAL,W1,X,W2,IGROUP,
     +              NXROW,NOBGNM)
                    GO TO 6899
                  END IF
                END DO
6899          CONTINUE
            END IF
          END IF
          IFLAG_X=1
          FLENME=CASEFL(:LEN_TRIM(CASEFL))//'.jco'
          CALL FFOPEN(JFAIL,-IRSF,'w',' ',25,CLINE)
          IF(JFAIL.NE.0) GO TO 9891
          CALL JCOWRITE(NPAR,NESPAR,NXROW,IRSF,X,ITRANS,APAR,
     +    AOBS,REFOBS,OVAL,PVAL,SCALE,OFFSET)
          CLOSE(UNIT=IRSF)
          IF(NOPTMAX.EQ.-2) THEN
            IFIN=11
            GO TO 6000
          END IF
          IF(NUMCOV.NE.0)THEN
            IF(NPRIOR.NE.0)THEN
              MATDIM1=NXROW
              NM=MXOBSCOV
              IPSTART=NOBS+1
              DO 6898 JCOV=1,NUMCOV
                IGROUP=COVGP(JCOV)
                DO I=IPSTART,NXROW
                  IF(NOBGNM(I).EQ.IGROUP)THEN
                    CALL ROTATE(2,NOBSCOV(JCOV),NOBSCOV(JCOV),MATDIM1,
     +              NESPAR,IPSTART,COVAR(JCOV)%DVAL,W1,X,W2,IGROUP,
     +              NXROW,NOBGNM)
                    GO TO 6898
                  END IF
                END DO
6898          CONTINUE
            END IF
          END IF
        END IF

C -- THE SENSITIVITY MATRIX IS NOW ROTATED.
C    HOWEVER THIS IS ONLY DONE IF OBSERVATION COVARIANCE MATRICES ARE PROVIDED.

        IF(NUMCOV.NE.0)THEN
          MATDIM1=NXROW
          NM=MXOBSCOV
          IPSTART=1
          DO 6900 JCOV=1,NUMCOV
            IGROUP=COVGP(JCOV)
            IF(NPRIOR.NE.0)THEN
              DO 6890 I=NOBS+1,NXROW
                IF(NOBGNM(I).EQ.IGROUP) GO TO 6900
6890          CONTINUE
            END IF
            DO 6920 I=IPSTART,NXROW
              IF(NOBGNM(I).EQ.IGROUP)THEN
                CALL ROTATE(2,NOBSCOV(JCOV),NOBSCOV(JCOV),MATDIM1,
     +          NESPAR,IPSTART,COVAR(JCOV)%DVAL,W1,X,W2,IGROUP,NXROW,
     +          NOBGNM)
                GO TO 6900
              END IF
6920        CONTINUE
6900      CONTINUE
        END IF

C -- ALSO, NOW THAT THE SENSITIVITY MATRIX HAS BEEN CALCULATED, WE ROTATE
C    THE REFERENCE OBSERVATIONS BACK AGAIN.

        IF(NUMCOV.GT.0)THEN
          DO 6905 JCOV=1,NUMCOV
            IF(NOBSCOV(JCOV).LE.0) GO TO 6905
            IGROUP=COVGP(JCOV)
            I=0
            DO 6903 J=1,NXROW
              IF(NOBGNM(J).EQ.IGROUP)THEN
                I=I+1
                W1(I)=REFOBS(J)
              END IF
6903        CONTINUE
            CALL ROTATE(1,NOBSCOV(JCOV),NOBSCOV(JCOV),1,1,1,
     +      COVAR(JCOV)%DVAL,W1,X,W2,IGROUP,NXROW,NOBGNM)
            I=0
            DO 6904 J=1,NXROW
              IF(NOBGNM(J).EQ.IGROUP)THEN
                I=I+1
                REFOBS(J)=W1(I)
              END IF
6904        CONTINUE
6905      CONTINUE
        END IF


6907    CONTINUE

C -- If subspace-enhanced Tikhonov is requested, that is done here.

        IF(IREG.NE.0)THEN
          IF(((IREGADJ.EQ.4).OR.(IREGADJ.EQ.5)).AND.
     +       (((IOPT-1)/NOPTREGADJ)*NOPTREGADJ.EQ.(IOPT-1)))THEN
              WRITE(IREC,6916)
              WRITE(6,6916)
6916          FORMAT(/,'    Undertaking Tikhonov regularisation sub',
     +        'space enhancment....')
#ifdef FLUSHFILE
              CALL FLUSH(IREC)
#endif

C -- First all regularisation weights are set to zero.

            TEMPDIM1=0
            DO IROW=1,NXROW
              ITEMP=NOBGNM(IROW)
              IF(IRGP(ITEMP).NE.0) THEN
                OWGHT(IROW)=0
              ELSE
                TEMPDIM1=TEMPDIM1+1
              END IF
            END DO

C -- SVD is undertaken on Q1/2X

            TEMPDIM2=NESPAR
            LWORKTEMP=MAX(3*MIN(TEMPDIM1,NESPAR)+MAX(TEMPDIM1,NESPAR),
     +                5*MIN(TEMPDIM1,NESPAR))
            LWORKTEMP=NINT(FLOAT(LWORKTEMP)*1.3)       ! ARBITRARY
            ALLOCATE(TEMPSVD(TEMPDIM1,TEMPDIM2),
     +      WORKVEC6(LWORKTEMP),STAT=IERR)
            IF(IERR.NE.0)THEN
              WRITE(ERRMSG,6911)
6911          FORMAT('Cannot allocate temporary memory required for ',
     +        'Tikhonov regularisation subspace enhancement.')
              GO TO 9890
            END IF
            ISTART_C=1
            IES=0
            IESS=0
            DO IPP=1,NPAR
              IF(ITRANS(IPP).LT.0)THEN
                IF(ITRANS(IPP).LT.-100001)IESS=IESS+1
                CYCLE
              END IF
              IES=IES+1
              IESS=IESS+1
              IF(MAXCOMPDIM.LE.1)THEN
                IRR=0
                DO IROW=1,NXROW
                  ITEMP=NOBGNM(IROW)
                  IF(IRGP(ITEMP).EQ.0) THEN
                    IRR=IRR+1
                    TEMPSVD(IRR,IES)=X(IROW,IESS)*SQRT(OWGHT(IROW))
                  END IF
                END DO
              ELSE
                CALL GET_VECTOR(NCOMPDIM,NXROW,XC,IXC,W1,IESS)
                ISTART_C=IFOUND_C
                IRR=0
                DO IROW=1,NXROW
                  ITEMP=NOBGNM(IROW)
                  IF(IRGP(ITEMP).EQ.0) THEN
                    IRR=IRR+1
                    TEMPSVD(IRR,IES)=W1(IROW)*SQRT(OWGHT(IROW))
                  END IF
                END DO
              END IF
            END DO
            IESTOT=IES
            NRMDIM=1
            NRMLODIM=1
            CALL DGESVD('N','O',TEMPDIM1,IESTOT,TEMPSVD,TEMPDIM1,W1,NRM,
     +      NRMDIM,NRMLO,NRMLODIM,WORKVEC6,LWORKTEMP,INFO)


C -- Before we do anything else we calculate the maximum magnitude of any regularisation row
C    of the Jacobian matrix.

            REGSENMAX=-1.0D200
            DO IROW=1,NXROW
              ITEMP=NOBGNM(IROW)
              IF(IRGP(ITEMP).NE.0)THEN
                IES=0
                IESS=0
                RSUM2=0.0D0
                DO IPP=1,NPAR
                  IF(ITRANS(IPP).LT.0)THEN
                    IF(ITRANS(IPP).LT.-100001)IESS=IESS+1
                    CYCLE
                  ELSE
                    IES=IES+1
                    IESS=IESS+1
                    IF(MAXCOMPDIM.LE.1)THEN
                      RSUM2=RSUM2+X(IROW,IESS)*X(IROW,IESS)
                    ELSE
                      CALL GET_VALUE(NCOMPDIM,XC,IXC,RTEMP,IROW,IESS)
                      RSUM2=RSUM2+RTEMP*RTEMP
                      ISTART_C=IFOUND_C
                    END IF
                  END IF
                END DO
                IF(RSUM2.GT.REGSENMAX)REGSENMAX=RSUM2
              END IF
            END DO

C -- New weights are now calculated.

            WEIGHTMAX=-1.0D200
            WEIGHTMIN=1.0D200
            NNEIG=MIN(IESTOT,TEMPDIM1)
            KREG=0
            DO IROW=1,NXROW
              ITEMP=NOBGNM(IROW)
              IF(IRGP(ITEMP).NE.0)THEN
                KREG=KREG+1
                PROJTOT=0.0D0
                WEIGHTPROJTOT=0.0D0
                RSUM2=0.0D0
                DO I=1,NNEIG
                  RSUM1=0.0D0
                  IES=0
                  IESS=0
                  DO IPP=1,NPAR
                    IF(ITRANS(IPP).LT.0)THEN
                      IF(ITRANS(IPP).LT.-100001)IESS=IESS+1
                      CYCLE
                    ELSE
                      IES=IES+1
                      IESS=IESS+1
                      IF(MAXCOMPDIM.LE.1)THEN
                        RSUM1=RSUM1+TEMPSVD(I,IES)*X(IROW,IESS)
                        IF(I.EQ.1)
     +                  RSUM2=RSUM2+X(IROW,IESS)*X(IROW,IESS)
                      ELSE
                        CALL GET_VALUE(NCOMPDIM,XC,IXC,RTEMP,IROW,IESS)
                        IF(RTEMP.NE.0.0D0) THEN
                          RSUM1=RSUM1+TEMPSVD(I,IES)*RTEMP
                          IF(I.EQ.1)
     +                    RSUM2=RSUM2+RTEMP*RTEMP
                        END IF
                        ISTART_C=IFOUND_C
                      END IF
                    END IF
                  END DO
                  IF(RSUM2.LT.1.0D-30*REGSENMAX)THEN   ! Arbitrary
                    IF(OWGHT(IROW).GT.0.0D0)THEN
                      OWGHT(IROW)=-1.1D300               ! A marker value
                    END IF
                    GO TO 6901
                  END IF
                  RTEMP1=W1(I)/W1(1)
                  IF(ABS(RTEMP1).LT.1.0D-200) GO TO 6909
                  RTEMP1=1/RTEMP1
                  IF(IREGADJ.EQ.4)THEN
                    IF(RTEMP1.GE.REGWEIGHTNUL) GO TO 6909
                  ELSE
                    IF(RTEMP1.GE.REGSINGTHRESH) GO TO 6909
                  END IF
                  RSUM1=RSUM1/SQRT(RSUM2)
                  PROJTOT=PROJTOT+RSUM1*RSUM1
                  RSUM1=RSUM1*RTEMP1
                  WEIGHTPROJTOT=WEIGHTPROJTOT+RSUM1*RSUM1
                END DO
6909            CONTINUE
                IF(IREGADJ.EQ.5)THEN
                  IF(OWGHTKP(KREG).EQ.0.0D0)THEN
                    OWGHT(IROW)=0.0D0
                  ELSE
                    IF(PROJTOT.GE.0.5D0)THEN
                      OWGHT(IROW)=1.0D0
                    ELSE
                      OWGHT(IROW)=REGWEIGHTRAT*REGWEIGHTRAT
                    END IF
                  END IF
                ELSE
                  IF(OWGHTKP(KREG).EQ.0.0D0)THEN
                    OWGHT(IROW)=0.0D0
                  ELSE
                    RTEMP1=1.0D0-PROJTOT
                    IF(RTEMP1.GT.0.0)WEIGHTPROJTOT=WEIGHTPROJTOT+
     +              RTEMP1*RTEMP1*REGWEIGHTNUL*REGWEIGHTNUL
                    OWGHT(IROW)=WEIGHTPROJTOT
                  END IF
                END IF
                IF(OWGHT(IROW).GT.WEIGHTMAX)WEIGHTMAX=OWGHT(IROW)
                IF(OWGHT(IROW).GT.0.0D0)THEN
                  IF(OWGHT(IROW).LT.WEIGHTMIN)WEIGHTMIN=OWGHT(IROW)
                END IF
              END IF
6901          CONTINUE
            END DO

! -- Regularisation observations which have proved to be totally insensitive are given a minimum weight.

            DO IROW=1,NXROW
              IF(OWGHT(IROW).LT.-1.0D300) OWGHT(IROW)=WEIGHTMAX
            END DO

C -- The spread of weights is now reduced in accordance with user settings.

            TOTWEIGHTSUM=0.0D0
            WEIGHTMAX=SQRT(WEIGHTMAX)
            WEIGHTMIN=SQRT(WEIGHTMIN)
            IF(ABS(WEIGHTMAX/WEIGHTMIN-1.0D0).LT.1.0D-3)THEN
              DO IROW=1,NXROW
                ITEMP=NOBGNM(IROW)
                IF(IRGP(ITEMP).GT.0)THEN
                  IF(OWGHT(IROW).GT.0.0D0)THEN
                    OWGHT(IROW)=1.0D0
                    TOTWEIGHTSUM=TOTWEIGHTSUM+OWGHT(IROW)
                  END IF
                END IF
              END DO
            ELSE
              IF(REGWEIGHTRAT.GT.0.0D0)THEN
                WEIGHTRANGE=WEIGHTMAX-WEIGHTMIN
                WEIGHTRATIO=(REGWEIGHTRAT*WEIGHTMIN-WEIGHTMIN)/
     +          WEIGHTRANGE
                DO IROW=1,NXROW
                  ITEMP=NOBGNM(IROW)
                  IF(IRGP(ITEMP).GT.0)THEN
                    IF(OWGHT(IROW).GT.0.0D0)THEN
                      OWGHT(IROW)=
     +                WEIGHTMIN+(SQRT(OWGHT(IROW))-
     +                WEIGHTMIN)*WEIGHTRATIO
                      OWGHT(IROW)=OWGHT(IROW)*OWGHT(IROW)
                      TOTWEIGHTSUM=TOTWEIGHTSUM+OWGHT(IROW)
                    END IF
                  END IF
                END DO
              ELSE
                RTEMP=WEIGHTMIN*ABS(REGWEIGHTRAT)
                WEIGHTMIN=LOG(WEIGHTMIN)
                WEIGHTMAX=LOG(WEIGHTMAX)
                WEIGHTRANGE=WEIGHTMAX-WEIGHTMIN
                WEIGHTRATIO=(LOG(RTEMP)-WEIGHTMIN)/WEIGHTRANGE
                DO IROW=1,NXROW
                  ITEMP=NOBGNM(IROW)
                  IF(IRGP(ITEMP).GT.0)THEN
                    IF(OWGHT(IROW).GT.0.0D0)THEN
                      OWGHT(IROW)=WEIGHTMIN+
     +                (0.5*LOG(OWGHT(IROW))-WEIGHTMIN)*WEIGHTRATIO
                      OWGHT(IROW)=EXP(OWGHT(IROW))
                      OWGHT(IROW)=OWGHT(IROW)*OWGHT(IROW)
                      TOTWEIGHTSUM=TOTWEIGHTSUM+OWGHT(IROW)
                    END IF
                  END IF
                END DO
              END IF
            END IF

C -- Weights are now normalised (so as not to be too high, this resulting in the
C    possible need for overloaded weight factors).

            IF(TOTWEIGHTSUM.NE.0.0D0)THEN
              TOTWEIGHTSUM=1.0D0/TOTWEIGHTSUM*WFLAST
              DO IROW=1,NXROW
                ITEMP=NOBGNM(IROW)
                IF(IRGP(ITEMP).NE.0)THEN
                  OWGHT(IROW)=OWGHT(IROW)*TOTWEIGHTSUM
                END IF
              END DO
            END IF

            DEALLOCATE(TEMPSVD,WORKVEC6,STAT=IERR)
            IF(IERR.NE.0)THEN
              WRITE(ERRMSG,6912)
6912          FORMAT('Cannot deallocate temporary memory required for ',
     +        'Tikhonov subspace enhancement.')
              GO TO 9890
            END IF
            WRITE(IREC,6917)
            WRITE(6,6917)
6917        FORMAT('    Regularisation weights have been ',
     +      're-calculated.')
            WRITE(6,*)
#ifdef FLUSHFILE
            CALL FLUSH(IREC)
#endif
          END IF

          WFLAST=WF
          WF=1.0
          IF(WF.LT.WFMIN)WF=WFMIN
          IF(WF.GT.WFMAX)WF=WFMAX

        END IF

C -- IF WE ARE DOING ADAPTIVE REGULARISATION, PERTINENT OBSERVATION WEIGHTS ARE WORKED
C    OUT BASED ON PERTINENT IW_ PARAMETER VALUES.

        IF(NREGADJPAR.GT.0)THEN
          ICOUNT=0
          DO IPP=1,NPAR
            IF(SCALE(IPP).LT.-1.0D35)THEN
              ITEMP=NINT(OFFSET(IPP))
              DO IOBS=1,NXROW
                IF(NOBGNM(IOBS).EQ.ITEMP)THEN
                  ICOUNT=ICOUNT+1
                  OWGHT(IOBS)=ORIGWGHT(ICOUNT)/PVAL(IPP)/PVAL(IPP)
                END IF
              END DO
            END IF
          END DO
          WF=1.0D0
          IF(WF.LT.WFMIN)WF=WFMIN
          IF(WF.GT.WFMAX)WF=WFMAX

C -- WE CALCULATE THE TOTAL CONTRIBUTION TO THE OBJECTIVE FUNCTION BY ALL
C    REGULARISATION GROUPS

          IRCOUNT=0
          TOTPHIREG=0.0D0
          DO IG=1,NOBSGP
            IF(IRGP(IG).NE.0)THEN
              IRCOUNT=IRCOUNT+1
              TOTPHIREG=TOTPHIREG+PSISUB(IG)
            END IF
          END DO

C -- WE ALSO ALTER THE OBSERVED VALUE OF THE PRIOR INFORMATION EQUATIONS WHICH
C    CITE THE INVERSE WEIGHT PARAMETERS. NOTE THAT THIS ASSUMES LOGARITHMIC
C    TRANSFORMATION OF THE IW PARAMETERS.

          ICOUNT=0
          DO IPP=1,NPAR
            IF(ITRANS(IPP).GE.0)THEN
              IF(SCALE(IPP).LT.-1.0D35)THEN
                ICOUNT=ICOUNT+1
                IOBS=PRIORPAR(ICOUNT)
                IF((IOBS.EQ.0).OR.(IOBS.GT.NXROW))THEN
                  WRITE(ERRMSG,6908)
6908              FORMAT('Incorrect adaptive regularisation ',
     +            'information supplied in PEST control file; ',
     +            'check this file with PESTCHEK.')
                  GO TO 9890
                END IF
                IF(OVAL(IOBS).GE.REFOBS(IOBS))THEN
                  OVAL(IOBS)=NINT(REFOBS(IOBS))
                  IF(OVAL(IOBS).GE.REFOBS(IOBS))
     +            OVAL(IOBS)=OVAL(IOBS)-1.0D0
                ELSE
                  IF(IOPT.GT.1)THEN
                    IF(IWSENS(ICOUNT).GT.GEOMAVSENS*2.0D0)THEN
                      RTEMP1=REFOBS(IOBS)-OVAL(IOBS)
                      RTEMP1=RTEMP1*0.667
                      OVAL(IOBS)=REFOBS(IOBS)-MAX(RTEMP1,0.3D0)           ! 0.3 is arbitrary
                    ELSE IF(IWSENS(ICOUNT).LT.GEOMAVSENS*0.5)THEN
                      IG=NOBGNM(IOBS)
                      IF(PSISUB(IG).LT.TOTPHIREG*0.95)THEN
                        RTEMP1=REFOBS(IOBS)-OVAL(IOBS)
                        RTEMP1=RTEMP1*1.3
                        OVAL(IOBS)=REFOBS(IOBS)-RTEMP1
                      END IF
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END DO

        END IF

C -- IF WE ARE WORKING IN REGULARISATION MODE, HAVE MORE THAN ONE REGULARISATION
C    OBSERVATION GROUP, AND HAVE RELATIVE WEIGHTS ADJUSTMENT ACTIVATED, THE
C    RELATIVE WEIGHTS ARE NOW WORKED OUT.
C -- BUT FIRST WE FIND THE REGULARISATION OBSERVATION GROUP WITH THE HIGHEST SENSITIVITY;
C    THIS WILL BE THE "PIVOT" OBSERVATION GROUP TO WHICH THE WEIGHT FACTOR DIRECTLY APPLIES
C    AND FOR WHICH THE RELATIVE ADJUSTMENT IS 1.0.


        IF((IREG.EQ.1).AND.((IREGADJ.GT.0).AND.(IREGADJ.LT.4)))THEN
          RMAXTOT=0.0D0
          IMAXTOT=0
          DO 5710 I=1,NOBSGP
            IF(IRGP(I).EQ.0)THEN
              SEOGP(I)=-1.0D0
              GO TO 5710
            ELSE
              RSUM2=0.0D0
              RMAX2=0.0D0
              RCOUNT=0.0D0
              ISTART_C=1
              DO 5720 K=1,NXROW
                IF(NOBGNM(K).NE.I) GO TO 5720
                IF(OWGHT(K).LE.0.0D0) GO TO 5720
                RSUM1=0.0D0
                IES=0
                IF(MAXCOMPDIM.LE.1)THEN
                  DO 5730 IES=1,NESPAR
                     RSUM1=RSUM1+X(K,IES)*X(K,IES)
5730              CONTINUE
                ELSE
                  DO IES=1,NESPAR
                    CALL GET_VALUE(NCOMPDIM,XC,IXC,RTEMP,K,IES)
                    IF(RTEMP.NE.0.0D0) RSUM1=RSUM1+RTEMP*RTEMP
                    ISTART_C=IFOUND_C
                  END DO
                END IF
                RSUM2=RSUM2+SQRT(RSUM1*OWGHT(K))
                IF(RSUM1.GT.RMAX2)RMAX2=RSUM1
                RCOUNT=RCOUNT+SQRT(OWGHT(K))
5720          CONTINUE
              IF((IREGADJ.EQ.1).OR.(IREGADJ.EQ.3))THEN
                SEOGP(I)=RSUM2
              ELSE IF(IREGADJ.EQ.2)THEN
                SEOGP(I)=RCOUNT
              END IF
              IF(RMAXTOT.LT.SEOGP(I))THEN
                RMAXTOT=SEOGP(I)
                IMAXTOT=I
              END IF
            END IF
5710      CONTINUE
          IPIVOTREG=IMAXTOT
          IF(IPIVOTREG.EQ.0)IPIVOTREG=1
          SEOGPPIV=SEOGP(IPIVOTREG)
          DO 5740 I=1,NOBSGP
            IF(IRGP(I).NE.0)THEN
              IF(SEOGP(I).GT.0.0D0)THEN
                SEOGP(I)=SEOGP(I)/SEOGPPIV
                IF(SEOGP(I).GT.0.0D0)THEN
                  IF(SEOGP(I).GT.1.0D-10)THEN  !ARBITRARY
                    SEOGP(I)=(1.0/(SEOGP(I)*SEOGP(I)))
                  ELSE
                    SEOGP(I)=0.0D0
                  END IF
                  IF(IOPT.GT.1)THEN
                    IF(SEOGP(I).GT.10.0D0)SEOGP(I)=10.0D0  ! This is supposed to combat massive adjustment - seawat problem.
                  END IF
                END IF
              END IF
            END IF
5740      CONTINUE
          KREG=0
          DO 5750 I=1,NXROW
            IGPNM=NOBGNM(I)
            IF(IRGP(IGPNM).NE.0)THEN
              KREG=KREG+1
              IF(IGPNM.NE.IPIVOTREG)THEN
                IF(SEOGP(IGPNM).GT.0.0D0)
     +          OWGHT(I)=OWGHT(I)*SEOGP(IGPNM)
              END IF
              IF(IREGADJ.EQ.3)THEN
                OWGHT(I)=OWGHT(I)*OWGHTKP(KREG)
              END IF
            END IF
5750      CONTINUE
          WRITE(IREC,5760)
5760      FORMAT(/,'    Regularisation weights adjustment....')
          WRITE(IREC,5762)
5762      FORMAT('    Group Name    Weights adjustment ratio from ',
     +    'previous iteration')
          DO 5765 I=1,NOBSGP
            IF(IRGP(I).NE.0)THEN
              IF(SEOGP(I).GE.0.0)THEN
                WRITE(IREC,5763) OBGNME(I)(1:LEN_TRIM(OBGNME(I))),
     +          SQRT(SEOGP(I))
5763            FORMAT(T6,A,T20,1PG13.6)
              END IF
            END IF
5765      CONTINUE
#ifdef FLUSHFILE
          CALL FLUSH(IREC)
#endif
        END IF

C -- If individual group target objective functions are supplied, weights are adjusted.

        IF(IGTARG.NE.0)THEN
          DO IG=1,NOBSGP
            IF(IRGP(IG).EQ.0)THEN

C -- First we compute contributions to the objective function on the basis of
C    original weights and compute our first cut at weight factor adjustment.

              TPSISUB=PSISUB(IG)/GFAC(IG)
              DTEMP=TPSISUB/GTARG(IG)
              IF(DTEMP.GT.2.0D0)DTEMP=2.0D0
              IF(DTEMP.LT.0.5D0)DTEMP=0.5D0
              GFAC(IG)=GFAC(IG)*DTEMP
              IF(GFAC(IG).GT.128.0) GFAC(IG)=128.0
              IF(GFAC(IG).LT.1.0D0/128.0D0)GFAC(IG)=1.0D0/128.0D0
            END IF
          END DO

C -- Now we adjust the adjustment so that the overall target measurement
C -- objective function is the same.

          TOTGPHIMLIM=0.0D0
          GPHIMLIM=0.0D0
          DO IG=1,NOBSGP
            IF(IRGP(IG).EQ.0)THEN
              TOTGPHIMLIM=TOTGPHIMLIM+GTARG(IG)*GFAC(IG)
              GPHIMLIM=GPHIMLIM+GTARG(IG)
            END IF
          END DO
          FACADJUST=GPHIMLIM/TOTGPHIMLIM
          DO IG=1,NOBSGP
            IF(IRGP(IG).EQ.0)THEN
              GFAC(IG)=GFAC(IG)*FACADJUST
            END IF
          END DO
          PHIMLIM=GPHIMLIM
          PHIMACCEPT=PHIMLIM*1.05
          PHIMLIMKP=PHIMLIM
          PD1RFAC=1.05

C -- Weights are now adjusted.

          DO IG=1,NOBSGP
            IF(IRGP(IG).EQ.0)THEN
              DTEMP=GFAC(IG)/OLDGFAC(IG)
              OLDGFAC(IG)=GFAC(IG)
              DO IROW=1,NXROW
                IF(NOBGNM(IROW).EQ.IG) OWGHT(IROW)=OWGHT(IROW)*DTEMP
              END DO
            END IF
          END DO
          WRITE(6,5766)
          WRITE(IREC,5766)
5766      FORMAT(/,'    Measurement weights adjusted for ',
     +    'fitting individual group phi targets.')
        END IF

C -- REGULARISATION WEIGHT FACTOR IS CALCULATED

        IF(IREG.EQ.1)THEN
          IF((NOPTMAX.NE.-1).AND.(NOPTMAX.NE.-2))THEN
          IF(FRACPHIM.GT.0.0)THEN
            PHIMLIM=PHIMLO*FRACPHIM
            IF(PHIMLIM.LE.PHIMLIMKP)PHIMLIM=PHIMLIMKP
            PHIMACCEPT=PHIMLIM*PD1RFAC
            WRITE(6,7454) PHIMLIM
            WRITE(IREC,7454) PHIMLIM
7454        FORMAT('    FRACPHIM-adjusted target measurement ',
     +      'objective function',T62,': ',1PG12.5)
#ifdef FLUSHFILE
            CALL FLUSH(IREC)
#endif
          END IF
7546      CONTINUE
          WRITE(6,5250)
5250      FORMAT('    Calculating optimal regularisation weight ',
     +    'factor .....')
          NSP4=NESPAR
          NXDIM=NXROW
          IF((MEMSAVE.NE.0).OR.(LSQRMODE.NE.0))THEN
          IF(CONJGRAD.EQ.1)THEN
            LHSD1=1
            LHSD2=1
            LHSVD=LHSVDIM
          ELSE IF(LSQRMODE.NE.0)THEN
            LHSD1=1
            LHSD2=1
            LHSVD=1
          ELSE IF(SVDMODE.EQ.2)THEN
            LHSD1=NXROW
            LHSD2=NESPAR
            LHSVD=1
          ELSE
            LHSD1=NSP4
            LHSD2=NESPAR
            LHSVD=1
          END IF
          IF(LSQRMODE.NE.0)THEN
            IF(PR_INDEX.NE.0)THEN
              IES=0
              IESS=0
              DO IPP=1,NPAR
                IF(ITRANS(IPP).LT.0)THEN
                  IF(ITRANS(IPP).LT.-100001)IESS=IESS+1
                  CYCLE
                END IF
                IES=IES+1
                IESS=IESS+1
                IESTRANS(IESS)=IES
              END DO
            END IF
          END IF
          CALL OPTWT_SL(IFAIL,NXDIM,NSP4,NXROW,NPRIOR,NOBS,NPAR,NESPAR,
     +    W1DIM,RHSDIM,REGITN,WF,WFFAC,PHIMLIM,WFMAX,WFMIN,WFSOL,WFTOL,
     +    ITRANS,REFOBS,OVAL,OWGHT,NOBGNM,W2,TMPOBS,SC,W1,X,RHS,LHS,
     +    INFOCOUNT,IUNSTABLE,LAMBDA,RLAMFAC,ICOVOBS,JSTK,NOBSGP,IRGP,
     +    CONJGRAD,LHSD1,LHSD2,npcgdim,cgr1,cgr2,cgv,cgw,cgy,
     +    cgprecon,cgshift,0,cgitnlim,cgeps,cgrtol,cgistop,cgitn,
     +    cganorm,cgacond,cgrnorm,cgxnorm,lhsvdim,lhsvec,
     +    WORKVEC1,WORKVEC2,WORKVEC3,
     +    SVDMODE,MAXSING,LWORK,EIGTHRESH,WORKVEC4)
          IF(IFAIL.EQ.1) GO TO 9970
          ELSE
          CALL OPTWT(IFAIL,NXDIM,NSP4,NXROW,NPRIOR,NOBS,NPAR,NESPAR,
     +    W1DIM,REGITN,WF,WFFAC,PHIMLIM,WFMAX,WFMIN,WFSOL,WFTOL,ITRANS,
     +    REFOBS,OVAL,OWGHT,NOBGNM,W2,TMPOBS,SC,W1,X,RHS,LHS,INFOCOUNT,
     +    IUNSTABLE,LAMBDA,RLAMFAC,ICOVOBS,JSTK,NORM,NRM,
     +    WORKVEC1,WORKVEC2,WORKVEC3,WORKVEC4,NOBSGP,IRGP,LHSFLAG,
     +    LASTOBSROW,FIRSTREGROW,INOCOV,IOPTCALL,LINREG,WFSTART,
     +    SVDMODE,MAXSING,LWORK,EIGTHRESH,WORKVEC4)
          IF(IFAIL.EQ.1) GO TO 9970
          END IF

          IF(ISTOP.EQ.2)THEN
            IFIN=10
            GO TO 6000
          ELSE IF(ISTOP.EQ.1) THEN
            IPFAIL=-1
            GO TO 9891
          END IF

          WRITE(IREC,*,ERR=9350)
          WRITE(IREC,5220,ERR=9350) WFSOL
          WRITE(6,5220) WFSOL
5220      FORMAT(T5,'Re-calculated regularisation weight factor',T62,
     +    ': ',1PG12.5)
          IF(ABS(WFSOL-WFMAX).LT.1.0D-5*WFMAX)THEN
            WRITE(6,5221)
            WRITE(IREC,5221,ERR=9350)
5221        FORMAT(T5,'***Warning: weight factor at its upper ',
     +      'bound***')
            IF(INFOCOUNT.NE.0)THEN
              WRITE(6,5222)
              WRITE(IREC,5222,ERR=9350)
5222          FORMAT(T5,'***Warning: near singular normal matrix***')
              WRITE(6,5223)
              WRITE(IREC,5223,ERR=9350)
5223          FORMAT(T5,'***Regularisation strategy may need ',
     +        'improvement***')
            ELSE
              IF(FRACPHIM.GT.0.0)THEN
                WRITE(6,5227)
                WRITE(IREC,5227,ERR=9350)
5227            FORMAT(T5,'***Consider lowering FRACPHIM***')
              END IF
            END IF
          ELSE IF(ABS(WFSOL-WFMIN).LT.1.0D-5*WFMIN)THEN
            WRITE(6,5224)
            WRITE(IREC,5224,ERR=9350)
5224        FORMAT(T5,'***Warning: weight factor at its lower ',
     +      'bound***')
          END IF
          IF(IUNSTABLE.EQ.1)THEN
            IF((FRACPHIM.GT.0.0).AND.(PHIMLIM.GT.PHIMLIMKP))THEN
                WRITE(6,7456)
                WRITE(IREC,7456)
7456            FORMAT('    Instability in weight factor calc. - ',
     +          'reducing target meas. obj. fn.')
                PHIMLIM=PHIMLIM/2.0
                IF(PHIMLIM.LT.PHIMLIMKP)PHIMLIM=PHIMLIMKP
                PHIMACCEPT=PHIMLIM*PD1RFAC
                WRITE(6,7454) PHIMLIM
                WRITE(IREC,7454) PHIMLIM
                GO TO 7546
            END IF
          END IF
        END IF
#ifdef FLUSHFILE
        CALL FLUSH(IREC)
#endif
        END IF

        IF(RSTFLE.NE.0)THEN
          INQUIRE(UNIT=IRSF,OPENED=LOPENED)
          IF(LOPENED) THEN
            WRITE(IRSF) WFSOL
            CLOSE(UNIT=IRSF)
          END IF
        END IF
302     CONTINUE
        IF(IAUI.EQ.0)THEN
          IF(JIRST.EQ.1)THEN
            JIRST=0
            IF(IREG.EQ.1)THEN
              IF((NOPTMAX.NE.-1).AND.(NOPTMAX.NE.-2))THEN
                IF(FRACPHIM.GT.0.0)THEN
                  PHIMLIM=PHIMLO*FRACPHIM
                  IF(PHIMLIM.LE.PHIMLIMKP)PHIMLIM=PHIMLIMKP
                  PHIMACCEPT=PHIMLIM*PD1RFAC
                  WRITE(6,7454) PHIMLIM
                  WRITE(IREC,7454) PHIMLIM
                  WRITE(IREC,*,ERR=9350)
                  WRITE(IREC,5220,ERR=9350) WFSOL
                  WRITE(6,5220) WFSOL
#ifdef FLUSHFILE
                  CALL FLUSH(IREC)
#endif
                END IF
              END IF
            END IF
          END IF
        END IF

        IF(HOLDFLAG.NE.0)THEN
          IF(IAUI.EQ.0)THEN
            CALL HLDREAD(JFAIL,NPAR,NPARGP,NESPAR,NLAMBDA,NRELPMX,
     +      NFACPMX,APAR,PARGNME,PHOLD,GHOLD,EHOLD,CLINE,IREG,NVECBND,
     +      IREC,NRELPREDSTP,NABSPREDSTP,NINITSCHFAC,NMULSCHFAC,
     +      NNSEARCH,IPRED)
            IF(JFAIL.NE.0) GO TO 9891
#ifdef FLUSHFILE
            CALL FLUSH(IREC)
#endif
          END IF
          IF(NLAMBDA.GT.0.0D0)LAMBDA=NLAMBDA
          IF(NRELPMX.GT.0.0D0) RELPARMAX=NRELPMX
          IF(NRELPREDSTP.GE.0.0D0) RELPREDSTP=NRELPREDSTP
          IF(NABSPREDSTP.GE.0.0D0) ABSPREDSTP=NABSPREDSTP
          IF(NINITSCHFAC.GT.0.0D0) INITSCHFAC=NINITSCHFAC
          IF(NMULSCHFAC.GT.1.0D0) MULSCHFAC=NMULSCHFAC
          IF(NNSEARCH.GE.1) NSEARCH=NNSEARCH
          IF(NSEARCH.GT.MAXSEARCH)NSEARCH=MAXSEARCH
          IF(NFACPMX.GT.0.0D0)THEN
            DMAX2=NFACPMX
            FACPARMAX=LOG10(NFACPMX)
          END IF
          IF(NVECBND.GE.0)UPVECBEND=NVECBND
        END IF

        IF(IAUI.EQ.0)THEN
          IF(IIIRST.NE.2) THEN
            WRITE(ISNS,128) IOPT
128         FORMAT(/,/,' OPTIMISATION ITERATION NO.',I3,' ----->')
            DO 123 IPP=1,NPAR
              IF(JSTK(IPP).LT.0) GO TO 124
123         CONTINUE
            GO TO 125
124         WRITE(ISNS,327)
327         FORMAT(1X,'Note that sensitivities for parameters glued ',
     +      'to their bounds may be in error.')
125         CONTINUE
            WRITE(ISNS,126)
126         FORMAT(1X,'Parameter name',t20,'Group',t35,
     +      'Current value',t52,'Sensitivity')
#ifdef FLUSHFILE
          CALL FLUSH(ISNS)
#endif
          END IF

          IF(IREG.EQ.1)THEN
            IF((NOPTMAX.NE.-1).AND.(NOPTMAX.NE.-2))THEN
            WFTEMP=WFSOL*WFSOL/WF/WF
            DO 5240 IROW=1,NXROW
              IGPNM=NOBGNM(IROW)
              IF(IRGP(IGPNM).NE.0)OWGHT(IROW)=OWGHT(IROW)*WFTEMP
5240        CONTINUE
            WF=WFSOL
            CALL OBJCLC(1,PSI,NXROW,NPRIOR,NOBS,REFOBS,OVAL,OWGHT,
     +      NOBGNM,NOBSGP,IRGP,SUM1)
            IF(PREDNOISE.GT.0) PSI=PSI+RES_PRED
            PSIL=PSI
            WRITE(IREC,5245,ERR=9350) PSI
            WRITE(6,5245) PSI
5245        FORMAT(T5,'New starting objective function for this itn. ',
     +      '(ie. phi)',T62,': ',1PG12.5)
            DO 5292 I=1,NOBSGP
              CALL OBJCLC(-I,PSISUB(I),NXROW,NPRIOR,NOBS,REFOBS,OVAL,
     +        OWGHT,NOBGNM,NOBSGP,IRGP,SUM1)
              WRITE(IREC,5293,ERR=9350)
     +        OBGNME(I)(1:LEN_TRIM(OBGNME(I))),PSISUB(I)
              WRITE(6,5293) OBGNME(I)(1:LEN_TRIM(OBGNME(I))),PSISUB(I)
5293          FORMAT('    Contribution to phi from observation group "',
     +        A,'"',T62,': ',1PG12.5)
5292        CONTINUE
            IF(NPRIOR.NE.0)THEN
              CALL OBJCLC(-999,PSISUB(NOBSGP+1),NXROW,NPRIOR,NOBS,
     +        REFOBS,OVAL,OWGHT,NOBGNM,NOBSGP,IRGP,SUM1)
              IF(PSISUB(NOBSGP+1).GT.-1.0D299)THEN
                 WRITE(IREC,5294,ERR=9350) PSISUB(NOBSGP+1)
                 WRITE(6,5294) PSISUB(NOBSGP+1)
5294             FORMAT('    Contribution to phi from ungrouped prior ',
     +           'information',T62,': ',1PG12.5)
              END IF
            END IF
            END IF
            IF(IREGADJ.NE.0)THEN
              IF(OWGHTLO(1).LT.-1.0E35)THEN
                DO 5751 I=1,NXROW
                  OWGHTLO(I)=OWGHT(I)
5751            CONTINUE
              END IF
            END IF
          END IF

#ifdef FLUSHFILE
            CALL FLUSH(IREC)
#endif

#ifdef PARALLEL
#ifndef MPEST
          WRITE(IRMR,2360)
2360      FORMAT(/,' Testing parameter upgrades .....')
#endif
#endif
        END IF

        NOWRLM=-1
        ILAMPL=1
        NOPAR=0
215     ITN=0
        NREV=0
        IMM=0
        IPASS1=0
        IF(DOAUI.EQ.'aui')THEN
          IF(LAMBDA.GT.1.0D8)LAMBDA=1.0D8
        END IF
        IF(RLAMFAC_ADJUST.EQ.1)THEN
          IF(LAMBDA.LT.1.0D-8) LAMBDA=1.0D-8
        END IF
        LAMSTR=LAMBDA
        NOWRLM=NOWRLM+1
        JACUPDATEFLAG=1

C -- NEXT THE GRADIENT VECTOR IS CALCULATED

        IES=0
        IESS=0
        ISTART_C=1
        DO 150 IPP=1,NPAR
        IF(ITRANS(IPP).LT.0)THEN
          IF(ITRANS(IPP).LT.-100001) IESS=IESS+1
          GO TO 150
        END IF
        IES=IES+1
        IESS=IESS+1
        RTEMP=0.0D0
        IF(MAXCOMPDIM.LE.1)THEN
        DO 170 IROW=1,NXROW
170     RTEMP=RTEMP+X(IROW,IESS)*OWGHT(IROW)*(OVAL(IROW)-REFOBS(IROW))
        ELSE
          CALL SINGLE_VECTOR_MUL1(IFAIL,NCOMPDIM,1,NXROW,XC,IXC,
     +    IESS,OWGHT,OVAL,REFOBS,RTEMP)
          IF(IFAIL.NE.0) GO TO 9970
          ISTART_C=IFOUND_C
        END IF
        GRAD(IES)=-2.0D0*RTEMP
150     CONTINUE
        NESTMP=IES
        IF(NESTMP.EQ.0) THEN
          WRITE(IREC,305,ERR=9350)
          WRITE(6,305)
305       FORMAT(T8,'all parameters frozen')
          IF(JMM.NE.0)THEN
            NOAUI=1
            GO TO 1200
          END IF
          IF((I2OR3.EQ.2).AND.(ISWTCH.EQ.1)) GO TO 1250
          IFIN=4
          GO TO 6000
        END IF
        CALL OBJCLC1(GRDNRM,NESTMP,GRAD)
        GRDNRM=SQRT(GRDNRM)
        IF(GRDNRM.LE.0.0D0) THEN
          IF(NOWRLM.NE.0)THEN
            WRITE(IREC,310,ERR=9350)
            WRITE(6,310)
310         FORMAT(T8,'phi gradient zero in non-frozen parameter ',
     +      'space')
          ELSE
            WRITE(IREC,311,ERR=9350)
            WRITE(6,311)
311         FORMAT(T5,'phi gradient zero in non-frozen parameter ',
     +      'space')
          END IF
          IF(SENSTO.EQ.1)THEN
            WRITE(ISNS,312)
312         FORMAT(' Phi gradient zero:- ',
     +      '  all parameter sensitivities are zero or near-zero.')
#ifdef FLUSHFILE
          CALL FLUSH(ISNS)
#endif
          END IF
          IF(JMM.NE.0) THEN
            NOAUI=0
            GO TO 1200
          END IF
          IF((I2OR3.EQ.2).AND.(ISWTCH.EQ.1)) GO TO 1250
          IF(PESTMODE.EQ.4) GO TO 12501
          IFIN=5
          GO TO 6000
        END IF