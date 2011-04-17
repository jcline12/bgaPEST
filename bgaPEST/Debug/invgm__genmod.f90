        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 09:10:33 2010
        MODULE INVGM__genmod
          INTERFACE 
            SUBROUTINE INVGM(N,A,INVA)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: A(N,N)
              REAL(KIND=8), INTENT(INOUT) :: INVA(N,N)
            END SUBROUTINE INVGM
          END INTERFACE 
        END MODULE INVGM__genmod
