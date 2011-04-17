        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:33 2010
        MODULE DTRMV__genmod
          INTERFACE 
            SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: TRANS
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: X(*)
              INTEGER(KIND=4) :: INCX
            END SUBROUTINE DTRMV
          END INTERFACE 
        END MODULE DTRMV__genmod
