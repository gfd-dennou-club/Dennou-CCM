program main
  use mod_atm
  implicit none 

  logical :: loop_flag = .true.

  call atm_init()

  do while(loop_flag)
    call atm_run(loop_flag)
  end do

  call atm_fin()

  stop

end program main
