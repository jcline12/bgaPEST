module param_trans

! ************************************************************************************************************
! Module with subroutines useful to transform sensitivity and parameters in estimation space and backtransform 
! ****************************************** M.D. 14/06/2010  ************************************************
     
 use bayes_pest_control
 use utilities  

contains

  subroutine sen_par_trans(cv_PAR, cv_OBS, d_PAR, d_A, d_PM)
  !Subroutine to transform the sensitivity and the parameters in
  !the estimation space.
  !Only if required and for the parameters that require the LOG 
  !transformation. 
  !First we do H=H*s (s is in the physical space)
  !Second s_old = log(s_old) !After s_old is in the estimation space
       implicit none
       
       integer                      :: i, j
       type (d_algorithmic)         :: d_A
       type (cv_param)              :: cv_PAR
       type (d_param)               :: d_PAR
       type (cv_observ)             :: cv_OBS
       type (d_prior_mean)          :: d_PM
       
        
       do i = 1, cv_PAR%p 
          if (d_PM%Partrans(i).eq.1) then
            do j = 1, cv_OBS%nobs
               where (d_PAR%BetaAssoc.eq.i)
                 d_A%H(j,:)=d_A%H(j,:)*d_PAR%pars
               end where
            enddo
            where (d_PAR%BetaAssoc.eq.i)
              d_PAR%pars_old = log(d_PAR%pars)   !MD At the beginning pars is the vector of the initial values of the parameters 
                 ! as read in the file. Then became the best estimate. Here we transform in the estimation space if required
            end where
          endif
        enddo

  end subroutine sen_par_trans
  
   subroutine par_back_trans(cv_PAR, d_PAR, d_PM)
  !Subroutine to back-transform the parameters in the physical spac/e.
  !Only if required and for the parameters that were LOG transformed. 
  !We do s = exp(s) After this s is again in the physical space
  
       implicit none
       
       integer                      :: i
       type (cv_param)              :: cv_PAR
       type (d_param)               :: d_PAR
       type (d_prior_mean)          :: d_PM
               
       do i = 1, cv_PAR%p
          if (d_PM%Partrans(i).eq.1) then
            where (d_PAR%BetaAssoc.eq.i)
              d_PAR%pars = exp(d_PAR%pars) !Back-transform the parameters in the physical space
            end where
          endif
       enddo

  end subroutine par_back_trans
  
  subroutine par_back_trans_lns(cv_PAR, d_PAR, d_PM) 
  !Subroutine to back-transform the parameters in the physical space.
  !Only if required and for the parameters that were LOG transformed. 
  !We do s = exp(s) After this s is again in the physical space
  
       implicit none
       
       integer                      :: i
       type (cv_param)              :: cv_PAR
       type (d_param)               :: d_PAR
       type (d_prior_mean)          :: d_PM
        
       do i = 1, cv_PAR%p
          if (d_PM%Partrans(i).eq.1) then
            where (d_PAR%BetaAssoc.eq.i)
              d_PAR%pars_lns = exp(d_PAR%pars_lns) !Back-transform the parameters in the physical space
            end where
           endif
       enddo

  end subroutine par_back_trans_lns

end module param_trans