        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:34 2010
        MODULE DSYSV__genmod
          INTERFACE 
            SUBROUTINE DSYSV(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,   &
     &INFO)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DSYSV
          END INTERFACE 
        END MODULE DSYSV__genmod
