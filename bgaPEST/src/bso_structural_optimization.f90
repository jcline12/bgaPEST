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
   use nelder_mead_simplex_linesearch
 contains
 
 
 subroutine marginal_structural_parameter_optim (cv_S,d_S, cv_PAR)
  
        use bayes_pest_control
        use jupiter_input_data_support
        use utilities  
        implicit none
        type (cv_struct),intent(inout)     :: cv_S
        type (cv_param), intent(in)        :: cv_PAR
        type(d_struct), intent(inout)      :: d_S
        integer                            :: i,j,k ! counters
        
 

!-------------------------------------------------------------------------
!------------------ CALL NELDER MEAD HERE --------------------------------
!-------------------------------------------------------------------------


 end subroutine marginal_structural_parameter_optim



end module struct_param_optimization

real (kind = 8) function SP_min(str_par_opt_vec, d_XQR, Q0_all,cv_OBS, d_OBS, cv_A, d_A, d_PAR, cv_S, d_S, d_PM, cv_PAR,cv_PM)

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
   type(d_struct),       intent(inout)  :: d_S
   type(Q0_compr),       intent(in)     :: Q0_All(:)
   type(cv_prior_mean),  intent(in)     :: cv_PM
   
   integer                              :: pivot(cv_OBS%nobs), errcode, i, j, k
   double precision,     intent(inout)  :: str_par_opt_vec(:) !This must be the pars vector to be optimized for. Must be the first argument in SP_min (NelMead requires this) 
   double precision                     :: z(cv_OBS%nobs) 
   double precision                     :: lndetGyy, ztiGyyz, dthQttdth
   double precision                     :: UinvGyy(cv_OBS%nobs,cv_OBS%nobs) ! used as both U and InvGyy
   double precision                     :: Gyy(cv_OBS%nobs,cv_OBS%nobs),  dtheta(cv_S%num_theta_opt)
   double precision, pointer            :: HXB(:), TMPV(:)  
   double precision, pointer            :: HXQbb(:,:)
   double precision, pointer            :: OMEGA(:,:)
   
   !----- intitialize variables
     errcode   = UNINIT_INT
     lndetGyy  = 0.D0
     ztiGyyz   = 0.D0
     dthQttdth = 0.D0
     UinvGyy   = UNINIT_REAL   ! matrix
     Gyy       = UNINIT_REAL   ! matrix
     dtheta    = UNINIT_REAL   ! matrix
   
   !First we need to reform d_S%theta and d_S%sig overwriting the elements that must be optimized for and  
   !leaving unchanged the others. These elements are into str_par_opt_vec (passed by NelMead). 
   !d_S%theta and d_S%sig are used during the matrix operations.(Qss Qsy HQsy Qyy)
   
   !******************************************************************************
   !******************************************************************************
   k = 0
   if ((maxval(cv_S%struct_par_opt).eq.1)) then   
     do i = 1, cv_PAR%p
       if (cv_S%struct_par_opt(i).eq.1) then
         do j = 1,cv_S%num_theta_type (i)
           k = k + 1
           d_S%theta(i,j) = str_par_opt_vec(k)
         end do
       end if
     end do
   end if
   
   if (d_S%sig_opt.eq.1) then 
     k = k + 1
     d_S%sig = str_par_opt_vec(k)
   endif
   !*****************************************************************************
   !*****************************************************************************
   
   !*************************************************************************************************************
   ! At this point d_S%theta, d_S%sig and str_par_opt_vec are ready to be used in the calculation of the obj func.
   
   !We need to recalculate Qss Qsy HQsy only if at least one theta optimization is required. Otherwise, if the optimization
   !is just for sig, we need to recalculate just Qyy
   if ((maxval(cv_S%struct_par_opt).eq.1)) then            !Only if theta opimization is required.  
      call bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All) 
   endif
   !Here Qss Qsy HQsy are ready, recalculated if necessary. We need to recalculate Qyy.
   call bmo_form_Qyy(d_XQR, d_S%sig, cv_OBS, d_A)
   
   ! At this point we need HXB, HXQbb, OMEGA only if we have prior information about beta. 
   if (cv_PM%betas_flag.ne.0) then   !------> we have prior information about beta   
      !Calculate HXB
      allocate(HXB(cv_OBS%nobs))
      HXB = UNINIT_REAL ! -- matrix (nobs)
      call DGEMV('n',cv_OBS%nobs, cv_PAR%p, 1.D0, d_A%HX, cv_OBS%nobs, &
              d_PM%beta_0, 1, 0.D0, HXB,1)
      !Form the linearization-corrected residuals and deallocate HXB no more necessary.
      z = d_OBS%obs - d_OBS%h + d_A%Hsold - HXB
      if (associated(HXB))  deallocate(HXB)
      !Form HXQbb
      call dgemm('n', 'n', cv_OBS%nobs, cv_PAR%p,  cv_PAR%p, 1.D0, d_A%HX, &
              cv_OBS%nobs, d_PM%Qbb, cv_PAR%p, 0.D0, HXQbb, cv_OBS%nobs)
      !Form OMEGA = HXQbb x (HX)' and deallocate HXQbb no more necessary.
      allocate(OMEGA(cv_OBS%nobs,cv_OBS%nobs))
      call dgemm('n', 't', cv_OBS%nobs, cv_OBS%nobs,  cv_PAR%p, 1.D0, HXQbb, &
              cv_OBS%nobs, d_A%HX,  cv_OBS%nobs, 0.D0, OMEGA, cv_OBS%nobs)
      if (associated(HXQbb))  deallocate(HXQbb)
      !Form Gyy = Qyy + OMEGA and deallocate OMEGA no more necessary.
      Gyy = OMEGA + d_A%Qyy
      if (associated(OMEGA))  deallocate(OMEGA)
   else !------> we don't have prior information about beta
     !Form the linearization residuals z  
     z = d_OBS%obs - d_OBS%h + d_A%Hsold
     !Form Gyy = Qyy 
     Gyy = d_A%Qyy
   endif  
   !***********************************************************************************************************
   
   !*******************************************************************************   
   !--- Calculate the determinant term of the str. pars objective function --------
   !----------------------------- 0.5ln(det(Gyy)) ---------------------------------
   !-- First perform LU decomposition on Gyy 
   UinvGyy = Gyy !-nobs x nobs --- note that this is used as U in this context
   call dgetrf(cv_OBS%nobs, cv_OBS%nobs, UinvGyy, cv_OBS%nobs, pivot, errcode)
   do i = 1,cv_OBS%nobs                            !Note that using DGETRF, the sign of the determinant should be not correct due to the pivoting
     lndetGyy = lndetGyy + dlog(abs(UinvGyy(i,i))) !And even if the determinant is positive, some element on the diagonal of U should be negative.
   end do                                          !If we are sure that the determinant is always positive the ABS before DLOG is enough.   
   lndetGyy = 0.5 * lndetGyy
   !*******************************************************************************
   
   !*******************************************************************************   
   !------------ Calculate the misfit term of the objective function --------------
   !------------------------ 0.5(z' x invGyy x z) ---------------------------------
   !Calculate the inverse of Gyy
   UinvGyy = Gyy   ! nobs x nobs, re-use UinvGyy, now as InvGyy
   call INVGM(cv_OBS%nobs,UinvGyy)
   !Form inv(Gyy)*z
   allocate(TMPV(cv_OBS%nobs))
     call DGEMV('n',cv_OBS%nobs, cv_OBS%nobs, 1.D0, UinvGyy, cv_OBS%nobs, &
            z, 1, 0.D0, TMPV,1) !On exit TMPV is invGyy*z
   !Multiply z' * TMPV and 0.5 
   call DGEMV('t',cv_OBS%nobs, 1, 5.0D-1, z, cv_OBS%nobs, &
          TMPV, 1, 0.D0, ztiGyyz,1)
   if (associated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
   !*******************************************************************************
   
   !**************************************************************************************   
   !--------- Calculate the prior theta/sig term of the objective function --------
   !-------------------- 0.5(dtheta x invQtheta x dtheta) -------------------------
   !Form dtheta
   dtheta = str_par_opt_vec - d_S%struct_par_opt_vec_0
   !Form invQtt*dtheta
   allocate(TMPV(cv_S%num_theta_opt))
   call DGEMV('n',cv_S%num_theta_opt, cv_S%num_theta_opt, 1.D0, d_S%invQtheta, &
           cv_S%num_theta_opt, dtheta, 1, 0.D0, TMPV,1) !On exit TMPV is invQtt*dtheta
   !Multiply dtheta' * TMPV and 0.5
   call DGEMV('t',cv_S%num_theta_opt, 1, 0.5D-1, dtheta, cv_S%num_theta_opt, &
           TMPV, 1, 0.D0, dthQttdth,1)
   if (associated(TMPV)) deallocate(TMPV) !Deallocate TMPV no more necessary here
   !**************************************************************************************
   
   !****************************************************************************** 
   !----------------- OBJECTIVE FUNCTION FOR STRUCTURAL PARAMETERS ---------------
   SP_min = lndetGyy + ztiGyyz + dthQttdth
   !*****************************************************************************
   
return
end function SP_min