module make_kernels
       
      contains
      subroutine bxq_make_X0_Q0_R0_InvQbb(d_PAR,cv_S,d_S,cv_PAR,d_XQR,cv_A,d_OBS,nobs,d_PM,Q0_All,cv_PM)
!--  Subroutine to create covariance matrix Q0, X matrix, R0 matrix, and if necessary Qbb^-1 and beta0*Qbb^-1      
!
!--  The covariance matrix can be:  (The flag (cv_A%Q_compression_flag) controls this)
!    - the full Q0 matrix ------------------------->   (cv_A%Q_compression_flag) = 0
!    - compressed form (in block for each beta) --->   (cv_A%Q_compression_flag) = 1
!--  In case of compressed form we have 3 different types:  (The flag Q0_All(x)%Toep_flag control this)
!    - full matrix for the specified beta --------------------->  (Q0_All(x)%Toep_flag) = 0
!    - just a vector for the specified beta  ------------------>  (Q0_All(x)%Toep_flag) = 1
!    - 1 value in case of nugget for the specified beta ------->  (either if the Toep_flag is [1] or [0])
!    
! Initial version 9/9/08  a m!ke@usgs joint
! Modified by M.D. 23/9/09
!
!
        use bayes_pest_control
        use jupiter_input_data_support
        use utilities      
        
        implicit none
        ! declarations
        type(d_param),intent(in)           :: d_PAR
        type(cv_struct),intent(inout)      :: cv_S
        type(d_struct), intent(in)         :: d_S
        type(cv_algorithmic), intent(in)   :: cv_A
        type(cv_param), intent(in)         :: cv_PAR
        type (cv_prior_mean), intent(in)   :: cv_PM
        type(Q0_compr), intent(inout)      :: Q0_All(:)
        type(d_observ), intent(in)         :: d_OBS
        type(kernel_XQR),intent(inout)     :: d_XQR
        type(d_prior_mean), intent(inout)  :: d_PM
        integer, intent(in)                :: nobs
        integer, pointer                   :: cnp(:)
        integer                            :: i,j,k,p ! local counters
        double precision                   :: ltmp ! Temporary value of Lmax
        character (len=ERRORWIDTH)         :: retmsg
        
        ! Allocate memory for  X and initialize to 0
        allocate(d_XQR%X(cv_PAR%npar,cv_PAR%p))
        d_XQR%X = 0.   ! matrix
        
select case (cv_A%Q_compression_flag)  !Select the compressed or not form of Q0 matrix    
     
  case(0) !Calculate full Q0 matrix    
        
        !*******************************************************************************************
        !************************** Make full Q0 and X matrix **************************************
        !*******************************************************************************************
        
           !Allocate memory for cnp **** Cnp is a counter to verify that in case of not nugget 
           !                             variogram, the number of parameters is gt 1
           allocate(cnp(cv_PAR%p))
           cnp=0 ! 2 or more means no problem, 1 found only one parameter, 0 no parameters found 
             
           select case (cv_A%store_Q)        
            case (.TRUE.)
            ! Allocate memory for  Q and initialize to 0
            allocate(d_XQR%Q0(cv_PAR%npar,cv_PAR%npar))
            d_XQR%Q0=0. ! matrix
                         
            !*** Start to fill the Q0 matrix based on the variogram type and to fill X matrix with 1 
            !*** to associate the correct beta to each parameter with a loop over all the parameters
            !*** In the same loop are 3 controls. 
            !*** 1. To check that each beta has at least one parameter defined
            !*** 2. To verify that in case of not nugget variogram, the number of parameters is gt 1
            !*** 3. To avoid that the same beta corresponds to parameters of different type.
            
            do i=1, cv_PAR%npar      !Loop over all the parameters
                if (cv_S%var_type(d_PAR%BetaAssoc(i))==0) then  ! put 1 on diagonal for nugget
                d_XQR%Q0(i,i)=1.
                cnp(d_PAR%BetaAssoc(i))= 2 
                do j=i+1, cv_PAR%npar                                    !*** This control avoid that the same
                  if (d_PAR%BetaAssoc(i).eq.d_PAR%BetaAssoc(j)) then     !*** beta corresponds to parameters
                    if (d_PAR%Group_type(i).ne.d_PAR%Group_type(j)) then !*** of different type.
                       write(retmsg,30) d_PAR%BetaAssoc(i), i, j
30                      format('Error: Beta association value ',i6, ' corresponds to different parameter types.' &
                        ' Check rows ',i6,' and ',i6,' of the parameter table. Excecution stopped.')
                        call utl_writmess(6,retmsg)
                       stop
                    endif
                  endif
                enddo  !Finished control
               else  ! all other variograms require distances
                cnp(d_PAR%BetaAssoc(i))= cnp(d_PAR%BetaAssoc(i))+1
                do j=i+1, cv_PAR%npar
                  if (d_PAR%BetaAssoc(i).eq.d_PAR%BetaAssoc(j)) then !Search in the parameters list the associated parameters and calculate the distances  
                    if (d_PAR%Group_type(i).ne.d_PAR%Group_type(j)) then !*** This control avoid that the same beta corresponds
                      write(retmsg,40) d_PAR%BetaAssoc(i), i, j         !*** to parameters of different type 
40                    format('Error: Beta association value ',i6, ' corresponds to different parameter types.' &
                      ' Check rows ',i6,' and ',i6,' of the parameter table. Excecution stopped.')
                      call utl_writmess(6,retmsg)
                      stop
                    endif   ! Finished control
                    do k = 1,cv_PAR%ndim            ! calculate the distances                   
                     d_XQR%Q0(i,j) = d_XQR%Q0(i,j) + (d_PAR%lox(i,k) - d_PAR%lox(j,k))**2 !Here the squared distance
                    enddo
                    d_XQR%Q0(i,j) = sqrt(d_XQR%Q0(i,j)) ! Now calculate the sqrt 
                    d_XQR%Q0(j,i)=d_XQR%Q0(i,j) ! Because the Q0 matrix is symmetric
                  endif
                enddo
               endif
              
              d_XQR%X(i,d_PAR%BetaAssoc(i))= 1. !Fill the X matrix to associate the correct beta to each parameter
             
             enddo
             if (minval(cnp).eq.0) then
               write(retmsg,10) minloc(cnp)
10              format('Error: No parameters correspond to beta association value',i6, &
                '. Excecution stopped.')
                call utl_writmess(6,retmsg)
                stop
             elseif (minval(cnp).eq.1) then
               write(retmsg,20) minloc(cnp)
20              format('Error: Found only one parameter that corresponds to beta association value',i6, &
                '. Variogram type must be nugget. Excecution stopped.')
                call utl_writmess(6,retmsg)
                stop             
             endif
       
            d_XQR%L = 10 * maxval(d_XQR%Q0)             ! Define L as 10 times the maximum distance
       
          case (.FALSE.) ! We need to address this option
            allocate(d_XQR%Q0(1,1))
            d_XQR%Q0 = UNINIT_REAL
          end select   ! d_A%store_Q
       
 
  case(1)!Calculate compressed form of Q0 matrix 
         !(in block for each beta or vector for each beta if toepl_flag is 1 and 1 value for nugget)
  
     !******************************************************************************************************
     !******* Make the Q0_C matrix or vector in case of Toeplitz or single value in case of nugget *********
     !******* and X matrix *********************************************************************************
     !******************************************************************************************************     
        
     d_XQR%L = 0. !Initialize the 10 times maximum distance in the Q0_C matrices

     select case (cv_A%store_Q)        
       case (.TRUE.)
        do p = 1, cv_PAR%p  !Loop for each beta that correspond to each different Q0_C 
          if (cv_S%var_type(Q0_All(p)%BetaAss)==0) then  ! Q0_C is just a single 1 for nugget
            allocate (Q0_All(p)%Q0_C(1,1)) !Allocation Just a value
            Q0_All(p)%Q0_C(1,1) = 1.
          else
           select case (Q0_All(p)%Toep_flag) ! Toep_flag [0] full Q0 matrix for that beta [1] just vector with distances for that beta
            case(0) !full matrix for this beta ---> allocate the matrix [npar * npar] for the p-th beta    
              allocate (Q0_All(p)%Q0_C(Q0_All(p)%npar,Q0_All(p)%npar)) !Allocation
              Q0_All(p)%Q0_C = 0. !Initialization
              do i =1, Q0_All(p)%npar 
                do j=i+1, Q0_All(p)%npar 
                  do k = 1,cv_PAR%ndim            ! calculate the distances                   
                    Q0_All(p)%Q0_C(i,j) = Q0_All(p)%Q0_C(i,j) + &
                    & (d_PAR%lox(Q0_All(p)%Beta_Start+i-1,k) - d_PAR%lox(Q0_All(p)%Beta_Start+j-1,k))**2 
                    !Here the squared distance. 
                    !Beta_Start identify where in the parameter list, starts the value with the p-th beta association
                  enddo
                  Q0_All(p)%Q0_C(i,j) = sqrt(Q0_All(p)%Q0_C(i,j)) ! Now calculate the sqrt 
                  Q0_All(p)%Q0_C(j,i) = Q0_All(p)%Q0_C(i,j) ! Because the Q0 matrix is symmetric
                enddo
              enddo
            case(1) !just a vector for this beta ---> allocate the matrix [npar * 1] for the p-th beta  
              allocate (Q0_All(p)%Q0_C(Q0_All(p)%npar,1)) !Allocation a vector
              Q0_All(p)%Q0_C = 0. !Initialization
                do j=2, Q0_All(p)%npar 
                  do k = 1,cv_PAR%ndim            ! calculate the distances                   
                    Q0_All(p)%Q0_C(j,1) = Q0_All(p)%Q0_C(j,1) + &
                    & (d_PAR%lox(Q0_All(p)%Beta_Start,k) - d_PAR%lox(Q0_All(p)%Beta_Start+j-1,k))**2 
                    !Here the squared distance. 
                    !Beta_Start identifies where in the parameter list, starts the value with the p-th beta association
                  enddo
                  Q0_All(p)%Q0_C(j,1) = sqrt(Q0_All(p)%Q0_C(j,1)) ! Now calculate the sqrt 
                enddo
                          
           end select  !Q0_All(p)%Toep_flag   
          
           ltmp = maxval(Q0_All(p)%Q0_C) !Temporary value of maximum distance in the p-th Q0_C matrix
           if (ltmp.gt.d_XQR%L)  d_XQR%L = ltmp
         
          endif ! cv_S%var_type(Q0_All(p)%BetaAss)==0
        enddo   !p = 1, cv_PAR%p 
       
        d_XQR%L = 10 * d_XQR%L !before here L was just the maximum distance in all the Q0_C matrices 
      
       case (.FALSE.) ! We need to address this option
         allocate(d_XQR%Q0(1,1))
         d_XQR%Q0 = UNINIT_REAL
     end select   ! d_A%store_Q
     
    ! Make the X matrix. The 1 values are associated in the order of the parameter list ***************
    do i=1, cv_PAR%npar      
      d_XQR%X(i,d_PAR%BetaAssoc(i))= 1. !Fill the X matrix to associate the correct beta to each parameter
    enddo
    !******************************************************************************************************
    !******* End Make the Q0_C matrix
    !******************************************************************************************************

end select ! (cv_A%Q_compression_flag) 


!*********************************************************************************************************
!*********************************************************************************************************
!********* The next lines are valid for both the full and compressed form of Q0 cases ********************
!*********************************************************************************************************
!*********************************************************************************************************

                

!******************************************************************************************************
!*************************************** Make the R0 matrix. ******************************************
!******************************************************************************************************
  allocate(d_XQR%R0(nobs,nobs))
    d_XQR%R0 = 0.  ! array
    !We use just a loop because R0 is diagonal *** Must change if allow full R0 matrix
    do i=1,nobs
     d_XQR%R0(i,i) = 1./(d_OBS%weight(i)**2)
    end do
!******************************************************************************************************
!********************************** End  Make the R0 matrix. ******************************************
!******************************************************************************************************
        

!******************************************************************************************************
!********* Calculate the inverse of the Qbb matrix and InvQbb * beta0 *********************************
!***************** We do that just if the betas_flag is not 0 *****************************************
!****************************************************************************************************** 
 if (cv_PM%betas_flag .ne. 0) then 
   !Calculate the inverse of the Qbb matrix
   d_PM%InvQbb=d_PM%Qbb
   call INVGM(cv_PAR%p,d_PM%InvQbb) !-- N.B. -> InvQbb must be Qbb on the way in and is returned as InvQbb
   !Calculate the product of inverse of the Qbb matrix and beta_0
   d_PM%InvQbbB0=matmul(d_PM%InvQbb,d_PM%beta_0)
 endif   
!******************************************************************************************************
!********* End Calculate the inverse of the Qbb matrix and InvQbb * beta0 *****************************
!******************************************************************************************************
end subroutine bxq_make_X0_Q0_R0_InvQbb


!******************************************************************************************************
!********* Determine the number of theta parameter to be optimized ************************************
!******************************************************************************************************

   
subroutine bxq_theta_cov_calcs(cv_PAR,cv_S,d_S,cv_PM,cv_A)
        
   use bayes_pest_control
        
   implicit none
   type (cv_prior_mean), intent(in)    :: cv_PM
   type (cv_param),      intent(in)    :: cv_par
   type (cv_struct),     intent(inout) :: cv_S
   type (d_struct),      intent(inout) :: d_S
   type (cv_algorithmic),intent(in)    :: cv_A
   integer                             :: i, j, k , m ! counters

   !----------------------------- Determine the number of theta parameter to be optimized ----------------
   ! -- first consider the prior means and add to num_th_opt based on the number of theta parameters for each 
   ! -- beta association for which theta parameters are meant 
   if ((maxval(cv_S%struct_par_opt).eq.1)) then  !MD This if statement is not strictly necessary 
     cv_S%num_theta_opt =sum(cv_S%struct_par_opt*cv_S%num_theta_type)
   endif
   !-- now add 1 if sigma is to be optimized for
   if (d_S%sig_opt .eq. 1) then 
     cv_S%num_theta_opt = cv_S%num_theta_opt + 1
   endif
   
   !----------------------------- Form struct_par_opt_0 --> local variable with starting values for all stuctural parameters - may include sigma
   allocate(d_S%struct_par_opt_vec_0(cv_S%num_theta_opt))
   allocate(d_S%struct_par_opt_vec(cv_S%num_theta_opt))
   
   k = 0
   if ((maxval(cv_S%struct_par_opt).eq.1)) then  !MD This if statement is not strictly necessary 
     do i = 1, cv_PAR%p
       if (cv_S%struct_par_opt(i).eq.1) then
         do j = 1,cv_S%num_theta_type (i)
           k = k + 1
           d_S%struct_par_opt_vec_0(k) = d_S%theta_0(i,j)
         end do
       end if
     end do
   end if
   
   if (d_S%sig_opt.eq.1) then 
     k = k + 1
     d_S%struct_par_opt_vec_0(k) = d_S%sig_0
   endif
      
   !----------------------------- Form Qtheta --> the prior covariance matrix for theta (can include sigma) --
   ! -- Finally, calculate the theta prior term
   allocate(d_S%invQtheta(cv_S%num_theta_opt,cv_S%num_theta_opt)) ! allocate  
     d_S%invQtheta = 0.D0 ! matrix, initialize to zero  MD important in case cv_A%theta_cov_form = 0 (no cov matrix)
        
   !       j = 0 ! counter up to cv_S%num_theta_opt
   !       do i = 1, cv_PAR%p                                           !****
   !          if (cv_S%struct_par_opt(i).eq.1) then                     !****  MD  This loop was not correct for different reasons
   !             select case (cv_PM%Qbb_form)                           !****    1. we need to select based on cv_A%theta_cov_form and not Qbb_form 
   !                 case (1)! diagonal                                 !****    2. in case of diagonal, when we read from the input file we have just a vector 
   !                     do k = 1, cv_S%num_theta_type(i)               !****       (dimension sum(num_theta_type)) with the covariance values and not a px2 matrix  
   !                        j = j+1                                     !****       even if we have,for some BetaAssoc, num_theta_type=2 
   !                        d_S%invQtheta(j,j) = d_S%theta_cov(i,k)     !****                           
   !                    end do                                          !****   So I rewrote the loop to be consistent with the input instruction, we need to change  
   !                 case (2) ! full matrix                             !****   here and the input reading if we don't like this way.
                    ! MNF --- need to add this eventually               !****
   !              end select                                            !****
   !          end if                                                    !****
   !       end do                                                       !****
 
    k=0 ! counter up to invQtheta
    m=0 ! counter up to vector that contains the covariance values d_S%theta_cov
    select case (cv_A%theta_cov_form)    
      case (1)! diagonal 
        do i = 1, cv_PAR%p 
          if (cv_S%struct_par_opt(i).eq.1) then
            do j = 1, cv_S%num_theta_type(i)              
              m=m+1
              k = k+1                                    
              d_S%invQtheta(k,k) = d_S%theta_cov(m,1)  ! In case of diagonal d_S%theta_cov is a vector                  
             end do                                            
          else
            m=m+cv_S%num_theta_type(i)
          end if  
        enddo  
      case (2) ! full matrix                        
        ! MNF --- need to add this eventually              
    end select                  

    if (d_S%sig_opt .eq. 1) then 
        d_S%invQtheta(cv_S%num_theta_opt,cv_S%num_theta_opt) = d_S%sig_p_var
    endif
   
!----------------------------- Finally peform the inversion of Qtheta to invQtheta ---------------------
   select case (cv_PM%Qbb_form)
    case (1)! diagonal
     do i = 1, cv_S%num_theta_opt
       d_S%invQtheta(i,i) =  1.D0/d_S%invQtheta(i,i)                            
     end do
       case (2) ! full matrix
         ! MNF --- need to add this eventually
   end select   
            
end subroutine bxq_theta_cov_calcs


end module make_kernels 