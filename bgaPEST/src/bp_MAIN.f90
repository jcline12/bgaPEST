program bp_main

! *********************************************
! ***     MAIN BAYES MODULE PEST PROGRAM    ***
! ***      Implementation of: Bayesian      ***
! *** Geostatistical inverse method in PEST.***
! ***                                       ***
! ***             a m!ke@usgs joint         ***
! ***          Michael N. Fienen, PhD       ***
! ***    UNITED STATES GEOLOGICAL SURVEY    ***
! ***             mnfienen@usgs.gov         ***
! ***              March 19, 2008           ***
! ***      Modified by M.D. 1/10/2009       ***
! ***  Further Modifications MNF/MD 11/2010 ***
! ***  Further Modifications   MD   09/2011 ***
! *********************************************
   
   use jupiter_input_data_support
   use bayes_pest_control
   use model_input_output
   use bayes_pest_mio_setup
   use bayes_pest_model
   use bayes_pest_reader
   use bayes_pest_finalize
   use bayes_output_control
   use error_message
   use utilities
   use make_kernels
   use bayes_matrix_operations
   use param_trans
   use linesearch
   use extern_derivs
   use bayes_output_control
   use struct_param_optimization
   use posterior_cov_operations
   implicit none


!--  Main Data Arrays for OBS and PARS and ALGORITHM   
       integer                      :: n1, i, j
       integer                      :: forward_flag_der
       integer                      :: s_ind, p_ind, b_ind  !Indices for structural parameters, quasi-linear and bga method loops
       type (mio_struc)             :: miostruc
       type (err_failure_struc)     :: errstruc  
       integer                      :: ifail, restart, outunit, bprunit, cparunit, cobsunit
       type (cv_algorithmic)        :: cv_A
       type (d_algorithmic)         :: d_A
       type (cv_prior_mean)         :: cv_PM
       type (d_prior_mean)          :: d_PM  
       type (cv_struct)             :: cv_S  
       type (d_struct)              :: d_S 
       type (cv_param)              :: cv_PAR
       type (Q0_compr), pointer     :: Q0_All(:)
       type (d_param)               :: d_PAR
       type (cv_observ)             :: cv_OBS
       type (d_observ)              :: d_OBS
       type (d_comlin)              :: d_MOD
       type (cv_minout)             :: cv_MIO
       type (d_minout)              :: d_MIO
       type (kernel_XQR)            :: d_XQR
       character (len=ERRORWIDTH)   :: retmsg
       character (len=100)          :: command_line, curr_par_file, curr_resid_file
       character (len=FILEWIDTH)    :: ctlfile
       character (len=FILEWIDTH)    :: casename
       character (len=FILEWIDTH)    :: atemp
       character (len=20)           :: inner_iter ! aka p_ind - this is the temporary holder for printing out the inner iteration number
       character (len=20)           :: outer_iter ! aka b_ind - this is the temporary holder for printing out the outer iteration number 
       double precision,dimension(1) :: curr_structural_conv, curr_phi_conv !Current iteration convergence values for structural parameters and quasi linear objective function
       double precision,dimension(1) :: curr_phi !Current value for quasi linear objective function
       double precision, pointer    :: curr_struct_vec(:) !Current vector of theta and sigma values to be optimized for
       double precision, pointer    :: VV(:,:), V(:) !VV is the posterior covariance matrix, V is only the diagonal of VV 
       double precision             :: huge_val=huge(huge_val) !Largest machine number

! -- PRINT OUT THE BANNER INFORMATION
   call bpo_write_banner()

!-- READ AND PARSE THE COMMAND LINE TO OBTAIN THE CONTROL FILE NAME       
  call UTL_GET_COMMAND_LINE(COMMAND_LINE)
  call UTL_PARSE_COMMAND_LINE(IFAIL,COMMAND_LINE,CTLFILE,RESTART)

!-- handle the case where no control file was indicated
  IF (IFAIL.NE.0) THEN
    call bpo_write_nocmd()
    stop
  END IF
  
  
! -- An extension of ".bgp" is added to the bgaPEST control file if necessary.
      i=LEN_TRIM(CTLFILE)
      IF(i.GE.5) THEN
        ATEMP=CTLFILE(I-3:I)
        CALL UTL_CASETRANS(ctlfile,'lo')
        IF(ATEMP.NE.'.bgp') THEN
          CASENAME = CTLFILE
          CTLFILE(I+1:)='.bgp'
        ELSE
          CASENAME = CTLFILE(1:I-4)
        END IF
      END IF
   
!--  INITIALIZE MIO STRUCTURE SO IT CAN BE PASSED    
    if(mio_initialise(errstruc,miostruc).ne.0) then  !MD mio_initialise is an integer
      call utl_bomb_out(errstruc)
      n1=mio_finalise(errstruc,miostruc)
      stop 
    end if  


! open up the main output record file 
       bprunit = utl_nextunit()          
       call bpc_openfile(bprunit,trim(trim(casename) // '.bpr'),1) ![1] at end indicates open with write access

!--  READ INPUT FILE AND PERFORM ASSOCIATED ALLOCATIONS AND PARSING     
    call bpr_read(errstruc,CTLFILE,cv_A, d_A, cv_PM, d_PM, cv_S, d_S, cv_PAR,Q0_All, d_PAR, &
                        cv_OBS, d_OBS, d_MOD, cv_MIO, d_MIO, miostruc)

!--  SETUP THE MODEL INPUT AND OUTPUT INFORMATION (TPL/INPUT AND INS/OUTPUT PAIRINGS)
    call bpm_setup_mio(errstruc, cv_MIO, d_MIO,  cv_PAR%npargp, &
                       cv_OBS%nobsgp, cv_PAR%npar, cv_OBS%nobs, miostruc)
                    
!--  INITIALIZE THE INVERSE MODEL
    !Make Q0, R0, X0, and InvQbb if necessary   
    call bxq_make_X0_Q0_R0_InvQbb(d_PAR,cv_S,d_S,cv_PAR,d_XQR,cv_A,d_OBS,cv_OBS%nobs,d_PM,Q0_All,cv_PM)
   
    allocate(d_OBS%h(cv_OBS%nobs)) ! Allocate the current model output [y]
    
!-- IF STRUCTURAL PARAMETERS WILL BE OPTIMIZED FOR, SET UP REQUIRED INFORMATION
    if ((maxval(cv_S%struct_par_opt).eq.1).or.(d_S%sig_opt.eq.1)) then
       call bxq_theta_cov_calcs(cv_PAR,cv_S,d_S,cv_PM,cv_A)
       allocate(curr_struct_vec(cv_S%num_theta_opt))
    end if
    
!-- CALL THE SETUP OF EXTERNAL DERIVATIVES FILES (IF REQUIRED).  THIS HAPPENS ONLY ONCE FOR ALL BUT PARAMETERS FILE
    if (cv_A%deriv_mode .eq. 0) then
        call bxd_write_ext_PEST_files(d_MOD, cv_MIO, d_MIO, cv_OBS, cv_PAR, d_OBS)
    end if
    
!-- WRITE THE HEADER INFORMATION TO THE REC FILE
    call bpo_write_bpr_header(bprunit,casename,cv_PAR,cv_OBS,d_MOD, cv_A, &
                cv_MIO, d_MIO,Q0_all,cv_PM,d_PM,cv_S,d_S,d_PAR)

    do b_ind = 1, cv_A%it_max_bga  !*********************************************************************** (more external loop)
    
    !***************************************************************************************************************************  
    !****************************** FROM HERE THE QUASI-LINEAR PARAMETER ESTIMATION LOOP ***************************************
    !***************************************************************************************************************************
    
          curr_phi_conv = huge_val !Initialize current quasi linear objective function convergence
          curr_phi      = huge_val !Initialize current quasi-linear objective function value
          do p_ind = 1, cv_A%it_max_phi !************************************************************* (first intermediate loop)
                                        !********** quasi-liner parameter estimation for given structural parameters ***********
      
             !-- RUN THE FORWARD MODEL (INCLUDES DELETING OLD OUTPUT, WRITING NEW INPUT, RUNNING MODEL, AND READING NEW OUTPUT)
             select case(cv_A%deriv_mode)
                case (0)
                    forward_flag_der = 1
                case (1)
                    forward_flag_der = 2
                    
             end select
             call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS, cv_A,  d_OBS%h, d_A%H, forward_flag_der, miostruc)
  
  
   
            d_PAR%pars_old = d_PAR%pars   !MD At the beginning pars is the vector of the initial values of the parameters 
                                          !as read in the parameter file. Then became the best estimate. 
            
            !-- CONVERT OR NOT SENSITIVITY AND PARAMETERS IN THE ESTIMATION SPACE
            if (maxval(d_PM%Partrans).eq.1) then  !If yes, the parameters transformation is required  
               call sen_par_trans(cv_PAR, cv_OBS, d_PAR, d_A, d_PM) !Converting sensitivity and parameters in the estimation space
            endif
          
            !-- SOLVE THE BAYESIAN LINEAR SYSTEM AND CALCULATE THE OBJECTIVE FUNCTIONS           
            call  bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All)
            call  bmo_form_Qyy(d_XQR, d_S%sig, cv_OBS, d_A)
            call  bmo_H_only_operations(d_XQR, d_A,cv_OBS,d_PAR,cv_PAR)
            call  bmo_solve_linear_system(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS, d_A, d_PAR,cv_PM)
            
            
            !-- PERFORMING LINESEARCH IF REQUIRED
            if (cv_A%lns_flag.eq.1) then  !If yes, we perform the linesearch procedure  
               call lns_proc(d_XQR,d_S,cv_PAR,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,cv_A,p_ind,miostruc,errstruc)
            endif
            
            !-- BACK-TRANSFORM OR NOT PARAMETERS IN THE PHYSICAL SPACE
            if (maxval(d_PM%Partrans).eq.1) then  !If yes, we need to back-transform the parameters in the physical space  
               call par_back_trans(cv_PAR, d_PAR, d_PM)
            endif 
            
            !-- set temporary string version of iteration numbers and phi to write out
            curr_phi_conv = abs(curr_phi - d_PAR%phi_T) 
            curr_phi = d_PAR%phi_T
            call UTL_INT2CHAR(p_ind,inner_iter)
            call UTL_INT2CHAR(b_ind,outer_iter)  
            curr_par_file = trim(casename) // '.bpp.' // trim(outer_iter) // '_' // trim(inner_iter)
            curr_resid_file = trim(casename) // '.bre.' // trim(outer_iter) // '_' // trim(inner_iter)
            !-- Write intermediate values out to BPR record file
            call bpo_write_bpr_intermed(bprunit,p_ind,b_ind,curr_par_file,curr_resid_file,d_PAR) 
            ! --Write the intermediate parameter and residuals files
      
            cparunit = utl_nextunit()
            call bpc_openfile(cparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
            call bpo_write_allpars(cv_PAR,d_PAR,d_PAR%pars,cparunit)
            close(cparunit)
            cobsunit = utl_nextunit()
            call bpc_openfile(cobsunit,trim(curr_resid_file),1) ![1] at end indicates open with write access
            call bpo_write_residuals(cv_OBS,d_OBS,d_OBS%h,cobsunit)
            close(cobsunit) 
            !-- check for convergence - exit if convergence has been achieved  
            if (curr_phi_conv(1).le.cv_A%phi_conv) exit

          enddo  !(first intermediate loop) quasi-linear method  --> p_ind
          
    !***************************************************************************************************************************  
    !************************************* END OF QUASI-LINEAR PARAMETER ESTIMATION LOOP ***************************************
    !***************************************************************************************************************************
    
    
    !***************************************************************************************************************************  
    !********************** FROM HERE THE STRUCTURAL PARAMETER ESTIMATION LOOP  (ONLY IF REQUIRED) *****************************
    !************************************************************************** *************************************************
       if ((maxval(cv_S%struct_par_opt).eq.1).or.(d_S%sig_opt.eq.1)) then !Enter the structural pars estimation loop only if required
         curr_struct_vec = d_S%struct_par_opt_vec !Current struct pars vector. At the first bga loop d_S%struct_par_opt_vec is d_S%struct_par_opt_vec_0  
         call marginal_struct_param_optim(d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,b_ind,cv_S%num_theta_opt)
         !Here d_S%struct_par_opt_vec is the vector of the optimized theta and sigma values
         curr_structural_conv(1) = sqrt(sum((curr_struct_vec - d_S%struct_par_opt_vec)**2)) !Calculate norm of difference between actual and previous vectors
         if (curr_structural_conv(1).le.cv_A%structural_conv) then !If yes, structural parameters have converged, the structural parameters estimation
           cv_S%struct_par_opt = 0                     !loop is no more required. Set to zero cv_S%struct_par_opt and d_S%sig_opt so the structural 
           d_S%sig_opt = 0                             !parameters estimation loop is no more entered. The optimized struct_par_opt_vec is used to run the 
                                                       !quasi-linear loop that should be the last one. Probably we need some output before exit.
           if (associated(curr_struct_vec)) deallocate(curr_struct_vec) 
         endif
       else
         exit !If the structural pars optimization is not required or structural pars have converged (run the last quasi_linear), exit the bga_loop
       endif
    !***************************************************************************************************************************  
    !*************************** END OF STRUCTURAL PARAMETER ESTIMATION LOOP  (ONLY IF REQUIRED) *******************************
    !***************************************************************************************************************************  
       
    enddo      !(more external loop) --> b_ind
    
    !*************************************************************************************************************************
    !******** FROM HERE THE EVALUATION OF THE POSTERIOR COVARIANCE (ONLY IF REQUIRED --> cv_A%post_cov_flag = 1 **************
    !*********** The posterior covariance is the full matrix (matrix VV) in case of no compression of Q, *********************
    !********************* it is only the diagonal (vector V) in case of compression of Q ************************************
    !*************************************************************************************************************************
     if (cv_A%post_cov_flag.eq.1) then
      call form_post_covariance(d_XQR, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR,Q0_All,cv_PM,d_PM,d_S,VV,V)
     end if
    !*************************************************************************************************************************
    !*********** END OF THE EVALUATION OF THE POSTERIOR COVARIANCE (ONLY IF REQUIRED --> cv_A%post_cov_flag = 1 **************
    !*************************************************************************************************************************
  
  
  write(*,*) 'Parameter estimation is complete!'
  do i = 1,cv_PAR%npar
    write(*,*) d_PAR%pars(i)
  enddo
  
!-- FINALIZE and CLEANUP - deallocate all the structures
    call bpd_finalize(d_PM, d_S, cv_PAR, d_PAR, cv_OBS, d_OBS, d_MOD, d_MIO, d_XQR)


end program bp_main                   


 


