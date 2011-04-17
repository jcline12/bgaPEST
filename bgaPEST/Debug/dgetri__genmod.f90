        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:33 2010
        MODULE DGETRI__genmod
          INTERFACE 
            SUBROUTINE DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGETRI
          END INTERFACE 
        END MODULE DGETRI__genmod
