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

real (kind = 8) function SP_min(theta, cv_OBS, d_OBS, d_A, d_PM, cv_PAR)

   use bayes_pest_control

   implicit none
   type(cv_observ),     intent(in)     :: cv_OBS
   type(d_observ),      intent(in)     :: d_OBS
   type(d_algorithmic), intent(in)     :: d_A
   type(d_prior_mean),  intent(in)     :: d_PM
   type(cv_param),      intent(in)     :: cv_PAR
   double precision, pointer           :: theta(:)
   double precision                    :: z(cv_OBS%nobs) 
   double precision                    :: HXB(cv_OBS%nobs)


   !----------------------------- Form the linearization-corrected residuals --------------------------
   HXB = UNINIT_REAL ! -- matrix (nobs x p)
   call DGEMV('n',cv_OBS%nobs, cv_PAR%p, 1.D0, d_A%HX, cv_OBS%nobs, &
            d_PM%beta_0, 1, 0.D0, HXB,1)
              
   z = d_OBS%obs - d_OBS%h + d_A%Hsold - HXB
   
   !----------------------------- 

SP_min = 0.0!theta objective function here!
return
end function SP_min