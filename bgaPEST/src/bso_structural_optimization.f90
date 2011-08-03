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

real (kind = 8) function SP_min(theta, cv_OBS, d_OBS, d_A.......)

   type(cv_observ),     intent(in)     :: cv_OBS
   type(d_observ),      intent(in)     :: d_OBS
   type(d_algorithmic),  intent(in)    :: d_A
   double precision, pointer           :: theta(:)
   double precision, pointer           :: z(cv_OBS%nobs) 



   z = d_OBS%obs - d_OBS%h + d_A%Hsold - ! NEED TO CALCULATE HXB here


SP_min = !theta objective function here!
return
end function SP_min