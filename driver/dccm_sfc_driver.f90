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

  call MPI_Barrier(MPI_COMM_WORLD, ierr)  
  call sfc_fin()

  stop

end program main
