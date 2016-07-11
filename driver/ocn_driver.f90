program main
  use mod_ocn
  use mpi
  implicit none 

  logical :: loop_flag = .true.
  integer :: ierr
  
  call ocn_init()

  do while(loop_flag)
    call ocn_run(loop_flag)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call ocn_fin()

  stop

end program main
