        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:34 2010
        MODULE DLASYF__genmod
          INTERFACE 
            SUBROUTINE DLASYF(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
              INTEGER(KIND=4) :: LDW
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NB
              INTEGER(KIND=4) :: KB
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: W(LDW,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLASYF
          END INTERFACE 
        END MODULE DLASYF__genmod
