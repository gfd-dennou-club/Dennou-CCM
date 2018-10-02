program main
  use dccm_ocn_mod
  use mpi
  implicit none 

  logical :: loop_flag = .true.
  integer :: ierr
  
  call ocn_init()

  do while(loop_flag)
    call ocn_run(loop_flag)
  end do

  call ocn_fin()

end program main
