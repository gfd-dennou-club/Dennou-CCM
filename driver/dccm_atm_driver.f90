program main
  use dccm_atm_mod
  use mpi
  implicit none 

  logical :: loop_flag = .true.
  integer :: ierr
  
  call atm_init()

  do while(loop_flag)
    call atm_run(loop_flag)
  end do

  call atm_fin()
!!$  stop

end program main
