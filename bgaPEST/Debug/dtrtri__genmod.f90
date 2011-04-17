        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:33 2010
        MODULE DTRTRI__genmod
          INTERFACE 
            SUBROUTINE DTRTRI(UPLO,DIAG,N,A,LDA,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              CHARACTER(LEN=1) :: DIAG
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DTRTRI
          END INTERFACE 
        END MODULE DTRTRI__genmod
