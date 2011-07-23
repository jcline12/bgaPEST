module extern_derivs
 
        use utilities 
        use bayes_pest_control 


contains

subroutine bxd_write_ext_PEST_files(d_MOD, cv_MIO, d_MIO, cv_OBS, cv_PAR, d_OBS)

implicit none


!-- declarations

       type (cv_param)              :: cv_PAR
       type (cv_observ)             :: cv_OBS
       type (d_observ)              :: d_OBS
       type (d_comlin)              :: d_MOD
       type (cv_minout)             :: cv_MIO
       type (d_minout)              :: d_MIO
       integer  i,j,k,iunit
!-- write the command line file
       iunit = utl_nextunit()
       call bpc_openfile(iunit,'bgaPEST.#mc',1)
       write(iunit,*) d_MOD%com
       close(iunit)
       
!-- write the MIO file
       iunit = utl_nextunit()
       call bpc_openfile(iunit,'bgaPEST.#mio',1)
       write(iunit,50) 'MIO_FILE','MOD_FILE'
       do i = 1,cv_MIO%ntplfle
        write(iunit,50) d_MIO%tpl(i),d_MIO%infle(i)
       end do
       do i = 1,cv_MIO%ninsfle
        write(iunit,50) d_MIO%ins(i),d_MIO%outfle(i)
       end do
50     format(A100,A100) 
       close(iunit)
    
!-- write the parameter group file       
       iunit = utl_nextunit()
       call bpc_openfile(iunit,'bgaPEST.#pargp',1)
       write(iunit,52) 'PARGPNME','DERINC','FORCEN'
52     format(A50,A10,A12)
       do i = 1,cv_PAR%npargp
        write(iunit,55) cv_PAR%grp_name(i), 0.01, 'always_2'
       end do       
55      format(A50,F10.5,A12)
       close(iunit)
           
!-- write the observations file       
       iunit = utl_nextunit()
       call bpc_openfile(iunit,'bgaPEST.#obs',1)
       write(iunit,60) 'OBSNME','OBSVAL','OBGNME','WEIGHT'
60     format(A50,' ',A18,' ',A12,' ',A18)       
       do i = 1,cv_OBS%nobs
        write(iunit,65) d_OBS%obsnme(i), d_OBS%obs(i), d_OBS%group(i), d_OBS%weight(i)
       end do       
65      format(A50,' ',ES18.8,' ',A12,' ',ES18.8)
       close(iunit)
       
 end subroutine bxd_write_ext_PEST_files      
      
      
      
subroutine bxd_write_param_file(cv_PAR,d_PAR)
       type (cv_param)             :: cv_PAR
       type (d_param)              :: d_PAR
       integer i,j,k, iunit
       
       iunit = utl_nextunit()
       call bpc_openfile(iunit,'bgaPEST.#par',1)
       write(iunit,70) 'PARNME','PARVAL1','PARGP'
70     format(A50,A20,A50)     
       do i = 1,cv_PAR%npar
            write(iunit,75) d_PAR%parnme(i),d_PAR%pars(i),d_PAR%group(i)    
       end do
75     format(A50,ES20.12,' ',A50)
       close(iunit)


end subroutine bxd_write_param_file
end module extern_derivs