!*****  subroutine bpf_model_run performs the following tasks
!       delete last set of model output files
!       write new model input files
!       run the forward model(s)
!       read and return the output data
module bayes_pest_model

contains
   subroutine bpf_model_run(errstruc, d_MOD, cv_PAR, d_PAR, cv_OBS, cv_A, obsval, H,forward_flag, miostruc)
    use utilities
    use bayes_pest_control
    use model_input_output
    use error_message
    use extern_derivs
    use jacread

 
    implicit none
!--  Main Data Arrays for OBS and PARS
    type (mio_struc)         :: miostruc
    type (err_failure_struc) :: errstruc
    type (d_comlin)                     :: d_MOD
    type (d_param)                      :: d_PAR
    type (cv_param)                     :: cv_PAR
    type (cv_observ)                    :: cv_OBS
    type (cv_algorithmic)               :: cv_A
    double precision, pointer           :: obsval(:)
    integer, intent(in)                 :: forward_flag ! 0, 1, 2, or 3
                                        ! 0 is forward run
                                        ! 1 is external PEST-style Jacobian
                                        ! 2 is dercom alternative Jacobian
                                        ! 3 is forward run for linesearch using d_PAR%pars_lns
    double precision, pointer           :: H(:,:)  ! Jacobian matrix  
    integer                             :: i, ifail
    character (len=100)                 :: adjfle

!-- MIO delete the last set of output files
    if(mio_delete_model_output_files(errstruc,miostruc).ne.0) then
      call utl_bomb_out(errstruc)
    end if   
 
!-- MIO write the model input files
    select case (forward_flag)
        case (3)
            if(mio_write_model_input_files(errstruc,miostruc, d_PAR%pars_lns).ne.0) then
                call utl_bomb_out(errstruc)
            end if 
        case default
            if(mio_write_model_input_files(errstruc,miostruc, d_PAR%pars).ne.0) then
              call utl_bomb_out(errstruc)
            end if   
    end select


!-- RUN THE MODEL IN THE MODE DESIRED
    select case (forward_flag)
        case (0) ! single forward run
            call system(d_MOD%com)
            !-- MIO read the ouput file results and update 
            if(mio_read_model_output_files(errstruc,miostruc, obsval).ne.0) then
              call utl_bomb_out(errstruc)
            end if 
        case (1) ! external PEST-style Jacobian
            call system(d_MOD%com) ! run the model once, forward, to have outputs with current parameters 
            !-- MIO read the ouput file results and update 
            if(mio_read_model_output_files(errstruc,miostruc, obsval).ne.0) then
              call utl_bomb_out(errstruc) 
            end if 
            !-- create PEST input files and run PEST          
            call bxd_write_param_file(cv_PAR,d_PAR) ! write the parameter file
            call system('pst_generator.exe')        ! create the necessary PEST control file
            call system('run_pest_scratch.bat')     ! run PEST externally for derivatives
            call readJCO('scratch.jco', H)
        case (2) ! dercom alternative Jacobian
            call system(d_MOD%dercom) 
            ! for now, assume this alternative is MODFLOW_ADJOINT
            select case (cv_A%jacobian_format)
              case ('binary')
                call readJCO(cv_A%jacfle,H)
              case ('ascii')
                call readJAC(cv_A%jacfle,H)
            end select
            !-- MIO read the ouput file results and update 
            if(mio_read_model_output_files(errstruc,miostruc, obsval).ne.0) then
              call utl_bomb_out(errstruc)
            end if 
        case (3) ! same as case 0, but linesearch parameters written as indicated above
            call system(d_MOD%com)
            !-- MIO read the ouput file results and update 
            if(mio_read_model_output_files(errstruc,miostruc, obsval).ne.0) then
              call utl_bomb_out(errstruc)
            end if 
    end select
    
  

   end subroutine bpf_model_run
   
end module bayes_pest_model