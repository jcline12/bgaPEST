        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:33 2010
        MODULE DSPMV__genmod
          INTERFACE 
            SUBROUTINE DSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: AP(*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: Y(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DSPMV
          END INTERFACE 
        END MODULE DSPMV__genmod
