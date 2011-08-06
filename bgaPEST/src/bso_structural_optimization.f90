module struct_param_optimization

   use jupiter_input_data_support
   use bayes_pest_control
   use model_input_output
   use bayes_pest_mio_setup
   use bayes_pest_model
   use bayes_pest_reader
   use bayes_pest_finalize
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use param_trans
   use nelder_mead_simplex
 contains
 
 
 subroutine marginal_structural_parameter_optim ()
 
 end subroutine marginal_structural_parameter_optim

end module struct_param_optimization

real (kind = 8) function SP_min(theta, sig, d_XQR, Q0_all,cv_OBS, d_OBS, cv_A, d_A, d_PAR, cv_S,d_PM, cv_PAR)

   use bayes_pest_control
   use bayes_matrix_operations
   implicit none
   type(kernel_XQR),     intent(in)     :: d_XQR
   type(cv_observ),      intent(in)     :: cv_OBS
   type(d_observ),       intent(in)     :: d_OBS
   type(d_algorithmic),  intent(inout)  :: d_A
   type(cv_algorithmic), intent(inout)  :: cv_A
   type(d_prior_mean),   intent(in)     :: d_PM
   type(cv_param),       intent(in)     :: cv_PAR 
   type(d_param),        intent(inout)  :: d_PAR        
   type(cv_struct),      intent(in)     :: cv_S 
   type(Q0_compr),       intent(in)     :: Q0_All(:)
   double precision,     intent(in)     :: theta(:,:)
   double precision,     intent(in)     :: sig
   double precision                     :: z(cv_OBS%nobs) 
   double precision                     :: HXB(cv_OBS%nobs)
   double precision                     :: HXQbb(cv_OBS%nobs,cv_PAR%p)
   double precision                     :: OMEGA(cv_OBS%nobs,cv_OBS%nobs)
   double precision                     :: Gyy(cv_OBS%nobs,cv_OBS%nobs)

   !------------------------@@@@@@@@@@@@@@@@@@@@----------------------@@@@@@@@@@@@@@@@@@@@@@--------------
   !
   !
   !    VITAL ISSUE HERE! -- We need to bring in theta and sig as a single vector but also parse it accordingly
   !
   !
   !------------------------@@@@@@@@@@@@@@@@@@@@----------------------@@@@@@@@@@@@@@@@@@@@@@--------------
   


   !----------------------------- Form the linearization-corrected residuals --------------------------
   HXB = UNINIT_REAL ! -- matrix (nobs)
   call DGEMV('n',cv_OBS%nobs, cv_PAR%p, 1.D0, d_A%HX, cv_OBS%nobs, &
            d_PM%beta_0, 1, 0.D0, HXB,1)
              
   z = d_OBS%obs - d_OBS%h + d_A%Hsold - HXB
   
   !----------------------------- Form Gyy with the current values of theta and sigma------------------
   call bmo_form_Qss_Qsy(d_XQR, theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All)
   call bmo_form_HQsy_Qyy(d_XQR, sig, cv_PAR, cv_OBS, d_A)

   HXQbb = UNINIT_REAL ! -- matrix (nobs x p)
   
   ! -- form HXQbb
   call dgemm('n', 'n', cv_OBS%nobs, cv_PAR%p,  cv_PAR%p, 1.D0, d_A%HX, &
              cv_OBS%nobs, d_PM%Qbb, cv_PAR%p, 0.D0, HXQbb, cv_OBS%nobs)
   
   ! -- now multiply HXQbb x (HX)' to form OMEGA

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MNF --- need to sort out all these dimensions
   call dgemm('n', 't', cv_OBS%nobs, cv_OBS%nobs,  cv_PAR%p, 1.D0, HXQbb, &
              cv_OBS%nobs, d_A%HX, cv_PAR%p, 0.D0, OMEGA, cv_OBS%nobs)
! MNF fix the above
!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! -- add OMEGA + Qyy to form Gyy
   Gyy = OMEGA + d_A%Qyy

   !----------------------------- Calculate the determinant term --------------------------------------
   
   !-- First perform LU decomposition on Gyy

   ! -- WHEN WE GET TO FORMING z*'Gyy(^-1)*z WE SHOULD CREATE A PSEUDO LINEAR SYSTEM FOR Gyy^(-1)*z

SP_min = 0.0 !theta objective function here!
return
end function SP_min