        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec 31 13:51:51 2010
        MODULE JMIN__genmod
          INTERFACE 
            FUNCTION JMIN(RHO,D_XQR,D_S,CV_PAR,D_A,D_PAR,D_PM,CV_OBS,   &
     &D_OBS,CV_PM,D_MOD,MIOSTRUC,ERRSTRUC) RESULT(JMIN_0)
              USE ERROR_MESSAGE
              USE MODEL_INPUT_OUTPUT
              USE BAYES_PEST_CONTROL
              REAL(KIND=8) :: RHO
              TYPE (KERNEL_XQR), INTENT(IN) :: D_XQR
              TYPE (D_STRUCT), INTENT(INOUT) :: D_S
              TYPE (CV_PARAM), INTENT(IN) :: CV_PAR
              TYPE (D_ALGORITHMIC), INTENT(INOUT) :: D_A
              TYPE (D_PARAM), INTENT(INOUT) :: D_PAR
              TYPE (D_PRIOR_MEAN), INTENT(IN) :: D_PM
              TYPE (CV_OBSERV), INTENT(IN) :: CV_OBS
              TYPE (D_OBSERV), INTENT(IN) :: D_OBS
              TYPE (CV_PRIOR_MEAN), INTENT(IN) :: CV_PM
              TYPE (D_COMLIN) :: D_MOD
              TYPE (MIO_STRUC) :: MIOSTRUC
              TYPE (ERR_FAILURE_STRUC) :: ERRSTRUC
              REAL(KIND=8) :: JMIN_0
            END FUNCTION JMIN
          END INTERFACE 
        END MODULE JMIN__genmod
