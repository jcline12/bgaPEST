	module adjoint_read

        use utilities
        
        
        public readMF_ADJOINT
        
        contains 
        
        subroutine readMF_ADJOINT(adjfle, nobs,npar, X)
        IMPLICIT NONE
        double precision, pointer     :: X(:,:)
        integer, intent(in)           :: nobs
        integer, intent(in)           :: npar
        integer                       :: i, j, adjunit
        character                     :: header
        character*(*), intent(in)         :: adjfle

 
             !Read the sensitivity from an Adjoint State output file          
             adjunit = utl_nextunit()
             open(adjunit,FILE = adjfle)
             do i=1,nobs+6
               read (adjunit,*) header
             enddo
             do i=1,nobs
               read (adjunit,*) header
               read (adjunit,*) header
               do j=1,npar
                  read (adjunit,*) X(i,j)
              enddo
             enddo
             close(adjunit)
        end subroutine readMF_ADJOINT
     end module adjoint_read