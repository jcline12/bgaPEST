        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:34 2010
        MODULE DSYTRF__genmod
          INTERFACE 
            SUBROUTINE DSYTRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSYTRF
          END INTERFACE 
        END MODULE DSYTRF__genmod
