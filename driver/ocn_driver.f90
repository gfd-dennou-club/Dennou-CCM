program main
  use mod_ocn
  implicit none 

  logical :: loop_flag = .true.

  call ocn_init()

  do while(loop_flag)
    call ocn_run(loop_flag)
  end do

  call ocn_fin()

  stop

end program main
