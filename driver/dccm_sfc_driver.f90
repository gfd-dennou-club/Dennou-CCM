program main
  use dccm_sfc_mod
  use mpi
  implicit none 

  logical :: loop_flag = .true.
  integer :: ierr
  
  call sfc_init()

  do while(loop_flag)
    call sfc_run(loop_flag)
  end do

  call sfc_fin()

end program main
