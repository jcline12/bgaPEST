module bayes_output_control
    ! ** DECLARATIONS **
       use jupiter_input_data_support
       use bayes_pest_control
       use bdp_data_parsers
       use bpi_initializations      
       use utilities
       use model_input_output
       use error_message
    
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine bpo_write_banner()
      write(6,*)
      write(6,*) '*******************************************'
      write(6,*) '*          Welcome to bgaPEST             *'
      write(6,*) '*******************************************' 
      write(6,*) '           Brought to you by:             ' 
      write(6,*) '     United States Geological Survey      '
      write(6,*) '      Watermark Numerical Computing       '
      write(6,*) '            University of Parma           '
      write(6,*) '*******************************************'
      write(6,*)
      
   end subroutine bpo_write_banner
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_nocmd()     
      write(6,*) 'bgaPEST is run using the command :-'
      write(6,*) 'bgaPEST bgapestfile '
      write(6,*) 'where bgapestfile is a bgaPEST control file'
   end subroutine bpo_write_nocmd
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine bpo_write_allpars(cv_PAR,d_PAR,bprunit,iternum)
      type (cv_param)             :: cv_PAR
      type (d_param)              :: d_PAR
      integer                     :: iternum
   end subroutine bpo_write_allpars
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine bpo_write_bpr_header(bprunit,casename,cv_PAR,cv_OBS, &
                d_MOD,cv_A,cv_MIO,d_MIO,Q0_all,cv_PM,cv_S)
   
    integer, intent(in)                 :: bprunit
    type(cv_param), intent(in)          :: cv_PAR
    type(cv_observ), intent(in)         :: cv_OBS
    type (cv_struct)                    :: cv_S  
    type(cv_prior_mean), intent(in)     :: cv_PM  
    type (Q0_compr), intent(in)         :: Q0_All(:)
    type(d_comlin), intent(in)          :: d_MOD
    type(cv_algorithmic), intent(in)    :: cv_A
    type(d_minout), intent(in)          :: d_MIO
    type(cv_minout), intent(in)         :: cv_MIO
    character (len=4)                   :: indent = '    '
    character (len=FILEWIDTH)           :: casename
    integer i,j ! local counters
    
    write(bprunit,*)
    write(bprunit,*)
    write(bprunit,10) casename
10  format('          bgaPEST RUN RECORD: CASE ', 1A)
!!! dimensions
    write(bprunit,*) 'Case Dimensions :-'
    write(bprunit,15) 'Number of Parameters',cv_PAR%npar   
    write(bprunit,15) 'Number of Parameter Groups',cv_PAR%npargp 
    write(bprunit,15) 'Number of Observations',cv_OBS%nobs    
    write(bprunit,15) 'Number of Observation Groups',cv_OBS%nobsgp    
15 format(A40, ':', I8)
   
!!!  Model calls
   write(bprunit,*)
   write(bprunit,*) 'Model Command Line :-'
   write(bprunit,20)  indent,d_MOD%com
   write(bprunit,*) 'Jacobian Command Line :-'
   if (cv_A%deriv_mode .eq. 1) then
      write(bprunit,20) indent,d_MOD%dercom
   else
      write(bprunit,20) indent,'Not Used'
   end if

!!! model input and output
   write(bprunit,*)
   write(bprunit,*) 'Model Interface Files:-'
   write(bprunit,20) indent,'Template files:'
   do i = 1,cv_MIO%ntplfle
     write(bprunit,25)indent,indent,d_MIO%tpl(i)
   end do
   write(bprunit,20) indent,'for model input files:'
   do i = 1,cv_MIO%ntplfle
     write(bprunit,25)indent,indent,d_MIO%infle(i)
   end do
   write(bprunit,20) indent,'Instruction files:'
   do i = 1,cv_MIO%ninsfle
     write(bprunit,25)indent,indent,d_MIO%ins(i)
   end do
   write(bprunit,20) indent,'for model ouput files:'
   do i = 1,cv_MIO%ntplfle
     write(bprunit,25)indent,indent,d_MIO%outfle(i)
   end do
20 format(2A )     ! single indent and str format
25 format(3A)   ! double indent and str format
30 format(1A, 1ES10.4)     ! single indent and ES format
35 format(3A,1ES10.4 )   ! double indent and ES format
40 format(2A, 1I10)     ! single indent and integer format
45 format(3A,1I10)   ! double indent and integer format
50 format(2A 1L5)     ! single indent and logical format
55 format(3A,1L5)   ! double indent and logical format
   
!!! Algorithmic control variables
    write(bprunit,*)
    write(bprunit,*) 'Algorithmic control variables:-'
    write(bprunit,20) indent,'Structural Paramter Convergence'
    write(bprunit,35) indent,indent,'structural_conv: ',cv_A%structural_conv
    write(bprunit,20) indent,'Objective Function Convergence'
    write(bprunit,35) indent,indent,'phi_conv: ',cv_A%phi_conv
    write(bprunit,20) indent,'Outermost BGA Convergence'
    write(bprunit,35) indent,indent,'bga_conv: ',cv_A%bga_conv       
    write(bprunit,20) indent,'Maximum Number of Structural Paramter Iterations'
    write(bprunit,45) indent,indent,'it_max_structural: ',cv_A%it_max_structural       
    write(bprunit,20) indent,'Maximum Number of Objective Function Iterations'
    write(bprunit,45) indent,indent,'it_max_phi: ',cv_A%it_max_phi       
    write(bprunit,20) indent,'Maximum Number of Outermost BGA Iterations'
    write(bprunit,45) indent,indent,'it_max_bga: ',cv_A%it_max_bga       
    write(bprunit,20) indent,'Linesearch Flag: [0] indicates no linesearh, [1] indicates perform linesarch'
    write(bprunit,45) indent,indent,'lns_flag: ',cv_A%lns_flag       
    write(bprunit,20) indent,'Maximum Number of Linesearch Iterations'
    write(bprunit,45) indent,indent,'it_max_lns: ',cv_A%it_max_lns
    write(bprunit,20) indent,'Logical Flag for storing Q Matrix'
    write(bprunit,55) indent,indent,'store_Q: ',cv_A%store_Q
    write(bprunit,20) indent,'Form of theta covariance: [0] none, [1], diag, [2] full matrix'
    write(bprunit,45) indent,indent,'theta_cov_form: ',cv_A%theta_cov_form
    write(bprunit,20) indent,'Compression of Q0 matrix: [0] none - compute full Q0, [1], Q0 for each beta'
    write(bprunit,45) indent,indent,'Q_compression_flag: ',cv_A%Q_compression_flag
    write(bprunit,20) indent,'Derivatives mode: [0] External PEST Perturbations, [1] specified Jacobian command line'
    write(bprunit,45) indent,indent,'deriv_mode: ',cv_A%deriv_mode
    
!!! Beta associations (facies associations)
    write(bprunit,*)
    select case (cv_PAR%p)
        case (1)
            write(bprunit,56) indent,cv_PAR%p
        case default
            write(bprunit,57) indent,cv_PAR%p
    end select
56  format(1A,I4, ' Beta Association was defined:')
57  format(1A,I4, ' Beta Associations were defined:')

    if (cv_A%Q_compression_flag .ne. 0) then
        ! write out each beta association's details
        do i = 1,cv_PAR%p
            write(bprunit,60) Q0_all%BetaAss
            write(bprunit,20) indent,'Parameter number at which this Beta Association starts'
            write(bprunit,45) indent,indent,'Beta_start: ',Q0_all%Beta_start
            write(bprunit,20) indent,'Toeplitz Flag'
            write(bprunit,45) indent,indent,'Toep_flag: ',Q0_all%Toep_flag
            write(bprunit,20) indent,'Number of Rows'
            write(bprunit,45) indent,indent,'Nrow: ',Q0_all%Nrow
            write(bprunit,20) indent,'Number of Columns'
            write(bprunit,45) indent,indent,'Ncol: ',Q0_all%Ncol
            write(bprunit,20) indent,'Number of Layers'
            write(bprunit,45) indent,indent,'Nlay: ',Q0_all%Nlay
            write(bprunit,20) indent,'Number of Parameters'
            write(bprunit,45) indent,indent,'Npar: ',Q0_all%Npar
        end do
    else
        !indicate no compression specified
            write(bprunit,65) indent,' No Compression Requested.  No further details provided about Beta Associations'
    end if
    
60 format('Variables for Beta Association: ',I3)
65 format(2A)
   
!!! Derivatives Calculations
    write(bprunit,*)
    !!! placeholder here, in case we implement group-specific derivatives  

!!! Prior Means Information if Supplied
    write(bprunit,*)
    if (cv_PM%betas_flag .eq. 1) then
        write(bprunit,*) 'Prior Information on Betas:-'
        write(bprunit,20) indent,'Prior information format on prior menas (betas)'
        write(bprunit,20) indent,'[0] none, [1] diagonal only, [2] full covariance matrix'
        write(bprunit,45) indent,indent,'Qbb_form: ', cv_PM%Qbb_form    
    else
        write(bprunit,20) indent,'No Prior Information Provided for Betas'
    end if
70 format('Prior mean information for beta association ', I2, ' supplied')
!!! Structural Parameter Definitions
    write(bprunit,*) 'Structural parameter control variables:-'
    write(bprunit,20) indent,'Flag determining whether structural parameters are optimized or fixed:'
    write(bprunit,45) indent,indent,'struct_par_opt: ',cv_S%struct_par_opt
    
    !MNF   FILL THIS IN!!!!!
    write(bprunit,*)
    write(bprunit,*) '*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*'
    write(bprunit,*) '                   MNF --- FILL IN THESE VARIABLES!      '        
    write(bprunit,*) '*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*'
    write(bprunit,*)
    
!!! Parameter definitions
    write(bprunit,*)
    write(bprunit,*) 'Initial Parameter Definitions:-'
     
    



   
   
   
   
   
   
   end subroutine bpo_write_bpr_header
end module bayes_output_control