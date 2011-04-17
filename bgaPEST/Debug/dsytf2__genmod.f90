        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:34 2010
        MODULE DSYTF2__genmod
          INTERFACE 
            SUBROUTINE DSYTF2(UPLO,N,A,LDA,IPIV,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSYTF2
          END INTERFACE 
        END MODULE DSYTF2__genmod
