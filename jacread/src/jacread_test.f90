program jacread_test

use jacread

implicit none

character(len=120) :: infile
double precision, pointer :: H(:,:)

infile = 'test.jco'
call readJCO(infile,H)


end program jacread_test