!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module interpolation_data_latlon_mod 

  ! モジュール引用; Use statements
  !
  use dc_types

  use field_common, only: ATM, OCN

  use jcup_interface, only: &
       & OPERATION_COEF, SEND_COEF, RECV_COEF, &
       & jcup_send_coef, jcup_recv_coef, &
       & jcup_set_local_coef, &
       & jcup_get_mpi_parameter, jcup_get_component_name, &
       & jcup_get_comp_num_from_name, &
       & jcup_error
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: interpolation_data_latlon_Init, interpolation_data_latlon_Final
  public :: set_operation_index, set_A_to_O_coef, set_O_to_A_coef
  public :: interpolate_data_latlon
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'interpolation_data_latlon_mod' !< Module Name

!  integer, parameter, private :: NO_GRID = -999999999
  integer, parameter, private :: NO_GRID = 0
  integer, parameter, private :: LARGE_VALUE = 999999999

  integer, parameter, private :: AGCM = 1
  integer, parameter, private :: OGCM = 2

  integer, private :: my_model

  type operation_index_type
    integer :: num_of_operation  ! num of my interpolation operation
    integer, pointer :: my_operation_index(:)   ! index of my operation 
    integer, pointer :: send_data_index(:) ! local send data index of each operation
    integer, pointer :: recv_data_index(:) ! local recv data index of each operation
    integer, pointer :: recv_coef_index(:) ! local recv coef index of each operation
    integer, pointer :: send_coef_index(:) ! local send coef index of each operation

    integer :: num_of_recv_coef ! array size of my local recv coef
    integer :: num_of_send_coef ! array size of my local send coef
    integer, pointer :: global_recv_coef_index(:)      ! global index of local recv coef array
    integer, pointer :: global_send_coef_index(:)      ! global index of local send coef array
    integer :: send_model_id
    integer :: index_tag
  end type

  type(operation_index_type), pointer :: operation_index(:,:,:)
  type(operation_index_type), pointer :: coi ! current operation index

  real(DP), allocatable :: coefS(:)

contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolation_data_latlon_Init(num_of_model, num_of_mapping_tag, my_model_id)
  implicit none
  integer, intent(IN) :: num_of_model
  integer, intent(IN) :: num_of_mapping_tag
  integer, intent(IN) :: my_model_id

  allocate(operation_index(num_of_model, num_of_model, num_of_mapping_tag))

  my_model = my_model_id

end subroutine interpolation_data_latlon_Init

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolation_data_latlon_Final()
  implicit none

end subroutine interpolation_data_latlon_Final


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_operation_index(recv_model_name, send_model_name, mapping_tag)
  use jcup_mpi_lib, only : jml_GetMyrankGlobal
  use jcup_interface, only : jcup_get_model_id, jcup_get_local_operation_index, &
                             jcup_get_num_of_send_grid, jcup_get_num_of_recv_grid
  implicit none
  character(len=*), intent(IN) :: recv_model_name, send_model_name
  integer, optional, intent(IN) :: mapping_tag
  integer :: recv_model_id, send_model_id
  integer :: grid_tag

  if (present(mapping_tag)) then
    grid_tag = mapping_tag
  else
    grid_tag = 1
  end if
  
  call jcup_get_model_id(recv_model_name, recv_model_id)
  call jcup_get_model_id(send_model_name, send_model_id)

  coi => operation_index(recv_model_id, send_model_id, grid_tag)

  coi%send_model_id = send_model_id
  coi%index_tag = grid_tag

  call jcup_get_local_operation_index(recv_model_name, send_model_name, grid_tag, &
                                      coi%num_of_operation, &
                                      coi%my_operation_index, &
                                      coi%send_data_index, &
                                      coi%recv_data_index, &
                                      coi%send_coef_index, &
                                      coi%recv_coef_index)

  if ((recv_model_id==2).and.(send_model_id==1)) write(700,*) "set operation index ", recv_model_id, send_model_id, size(coi%send_data_index)

  call jcup_get_num_of_send_grid(recv_model_name, send_model_name, grid_tag, coi%num_of_send_coef)
  call jcup_get_num_of_recv_grid(recv_model_name, send_model_name, grid_tag, coi%num_of_recv_coef)

  
end subroutine set_operation_index

subroutine set_A_to_O_coef(mapping_tag, coefS_global)

  use jcup_grid

  integer, intent(in) :: mapping_tag  
  real(DP), intent(in), optional :: coefS_global(:)

  
  integer :: my_comm, my_group, my_size, my_rank
  character(TOKEN) :: my_comp_name
  integer, pointer :: local_coef_index(:)
  integer :: num_of_index
  
  my_comp_name = jcup_get_component_name(my_model)
  call jcup_get_mpi_parameter(my_comp_name, my_comm, my_group, my_size, my_rank)
  
  coi => operation_index(jcup_get_comp_num_from_name(OCN), jcup_get_comp_num_from_name(ATM),mapping_tag)

  if(my_comp_name == OCN) then

     allocate(coefS(coi%num_of_operation))
     write(*,*) "OCN rank=", my_rank, "nOperation=", coi%num_of_operation, &
          & size(coi%send_data_index), size(coi%recv_data_index), &
          & coi%num_of_send_coef, coi%num_of_recv_coef, &
          & ", mapping_tag=", mapping_tag
     call get_operation_index(OCN, ATM, mapping_tag, num_of_index, local_coef_index)
     write(*,*) "call get_operation_index: nIndex=", num_of_index, "nLcCoefIndex=", &
          & size(local_coef_index)
     
     call jcup_recv_coef(OCN, ATM, mapping_tag, coefS, OPERATION_COEF)
  else if(my_comp_name == ATM .and. my_rank==0) then
     if(.not. present(coefS_global)) &
          & call jcup_error("set_A_to_O_coef", "NO coefS")
     
     call jcup_send_coef(ATM, OCN, coefS_global)
  else
     return
  end if

end subroutine set_A_to_O_coef

subroutine set_O_to_A_coef(mapping_tag, coefS_global)

  integer, intent(in) :: mapping_tag  
  real(DP), intent(in), optional :: coefS_global(:)

  
  integer :: my_comm, my_group, my_size, my_rank
  character(TOKEN) :: my_comp_name
  
  my_comp_name = jcup_get_component_name(my_model)
  call jcup_get_mpi_parameter(my_comp_name, my_comm, my_group, my_size, my_rank)
  
  coi => operation_index(jcup_get_comp_num_from_name(ATM), jcup_get_comp_num_from_name(OCN), mapping_tag)

  if(my_comp_name == ATM) then
     allocate(coefS(coi%num_of_operation))

     call jcup_set_local_coef(ATM, OCN, mapping_tag, coefS_global, coefS, OPERATION_COEF)

  end if
end subroutine set_O_to_A_coef

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data_latlon(recv_model, send_model, send_data, & 
                                   recv_data, num_of_data, grid_num, exchange_tag)
  implicit none
  character(len=*), intent(IN) :: recv_model, send_model
  real(kind=8), intent(IN) :: send_data(:,:)
  real(kind=8), intent(INOUT) :: recv_data(:,:)
  integer, intent(IN) :: num_of_data
  integer, intent(IN) :: grid_num
  integer, intent(IN) :: exchange_tag(:)

  integer :: i, d, j, ri
  integer :: send_point, recv_point, operation_coef, send_coef, recv_coef

  integer :: my_comm, my_group, my_size, my_rank


    coi => operation_index(jcup_get_comp_num_from_name(recv_model), jcup_get_comp_num_from_name(send_model),grid_num)
    !write(0,*) "interpolation data test ", send_model, recv_model, grid_num, num_of_data, &
    !            size(operation_index(recv_model, send_model, grid_num)%send_data_index), & 
    !            send_data(1,1), send_data(2,1)

    recv_data(:,:) = 0d0

    do d = 1, num_of_data
      do i = 1, size(coi%send_data_index)
        send_point = coi%send_data_index(i)
        recv_point = coi%recv_data_index(i) 
        recv_data(recv_point,d) = recv_data(recv_point,d) &
             & + send_data(send_point,d)*CoefS(i)
      end do
    end do

    return

  end subroutine interpolate_data_latlon

!=======+=========+=========+=========+=========+=========+=========+=========+
  
end module interpolation_data_latlon_mod

