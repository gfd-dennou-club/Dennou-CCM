program main
  use mod_atm
  use mpi
  implicit none 

  logical :: loop_flag = .true.
  integer :: ierr
  
  call atm_init()

  do while(loop_flag)
    call atm_run(loop_flag)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr)  
  call atm_fin()

  stop

end program main
