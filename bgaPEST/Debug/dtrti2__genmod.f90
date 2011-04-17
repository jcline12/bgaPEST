        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:33 2010
        MODULE DTRTI2__genmod
          INTERFACE 
            SUBROUTINE DTRTI2(UPLO,DIAG,N,A,LDA,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DTRTI2
          END INTERFACE 
        END MODULE DTRTI2__genmod
