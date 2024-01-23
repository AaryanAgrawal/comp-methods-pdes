! data_type :: name
! always have implicit none

program variables
  implicit none

  integer :: amount
  real :: pi
  complex :: frequency
  character :: initial
  logical :: isOkay
  
  read(*, *) amount
  print *, 'amount = ', amount

end program variables
