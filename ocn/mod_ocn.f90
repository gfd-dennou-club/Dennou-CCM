!-------------------------------------------------------------
! Copyright (c) 2015-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module mod_ocn

  ! モジュール引用; Use statements
  !

  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify
  
  !* Dennou-OGCM
  
  use DOGCM_main_mod, only: &
       & ogcm_main_init => dogcm_main_Init,               &
       & ogcm_main_final => dogcm_main_Final,             &
       & ogcm_setup => dogcm_main_setup,                  &
       & ogcm_shutdown => dogcm_main_shutdown,            &       
       & ogcm_advance_timestep => dogcm_main_advance_timestep

  use DOGCM_Admin_TInteg_mod, only: &
       & TimeSecB, TimeSecN, TimeSecA, &
       & OCN_TLN => TIMELV_ID_N

  use DSIce_Admin_TInteg_mod, only:   &
       & SICE_TLN => TIMELV_ID_N
    
  use DOGCM_Admin_Grid_mod, only: &
       & ISO => IS, IEO => IE, IMO => IM,   &
       & JSO => JS, JEO => JE, JMO => JM,   &
       & KSO => KS, KEO => KE, KMO => KM,   &
       & IAO => IA, JAO => JA, KAO => KA
  
  
  !* DCPCM

  use component_field, only: &
       & component_field_type
  

  use mod_common_params, only: &
       & NUM_DCCM_COMP,             &
       & NUM_OCN_GMAPTAG,           &
       & COMPNAME_OCN, GN25,        &
       & ATM_GRID_2D, OCN_GRID_2D,  &
       & GMAPTAG_ATM2D_OCN2D
       
  use mod_common_compdef, only:     &
       & my_comp   => CompDef_ocn,  &
       & atm_comp  => CompDef_atm,  &
       & sice_comp => CompDef_sice, &
       & GMAPFILENAME_AO,           &
       & GMAPFILENAME_OA,           &
       & AO_COUPLING_CYCLE_SEC
  
  ! 宣言文; Declareration statements
  !    
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !    
  public :: ocn_init
  public :: ocn_run
  public :: ocn_fin

  ! 非公開変数
  ! Private variable
  !
  
  integer :: tstep
  logical :: loop_end_flag
  
  type(component_field_type) :: field
  type(component_field_type) :: field3d

  integer, allocatable :: lnx_field(:)
  integer, allocatable :: lny_field(:)
  
!!$  real(DP), allocatable :: xy_TmpTAvg(:,:)

  character(*), parameter :: module_name = 'mod_ocn'
  
contains

  subroutine ocn_init()

    ! モジュール引用; Use statements
    !

    use field_def, only: &
         & cal_mn

    !* Dennou-OGCM
   
    use DOGCM_Admin_TInteg_mod, only: &
         & InitialDate => InitDate, &
         & RestartTime, DelTime
    
    ! 宣言文; Declareration statements
    !    

    ! 実行文; Executable statement
    !
    
    ! Initialze varibale to manage component  ********************
    !

    call ComponentDef_Init( my_comp,    ,    & ! (inout)
         & COMPNAME_OCN, "Dennou-OGCM" )       ! (in)
    

    !  Intialize a model associated with  ocean component ***
    !
    call MessageNotify( 'M', module_name, &
         & 'Initialize %a .. (MPI_MY_COMM=%d, MY_RANK=%d)', &
         & ca=(/ trim(my_comp%MODELNAME) /), i=(/ my_comp%PRC_comm, my_comp%PRC_rank/))

    call ogcm_main_init( my_comp%PRC_comm, isCoupledRun=.true. )
    call ogcm_setup()

    my_comp%InitTimeInfo(:) = (/ InitialDate%year, InitialDate%month, InitialDate%day,        &
                        &        InitialDate%hour, InitialDate%min, int(InitialDate%sec)  /)
    my_comp%RestartTimeSec  = RestartTime
    my_comp%DelTime         = DelTime
    my_comp%tstep           = 0    
    call MessageNotify( 'M', module_name, &
         & "Inital Time %4d - %2d - %2d:%2d%2d", i=my_comp%InitTimeInfo )
    call MessageNotify( 'M', module_name, &
         & "RestartTiem = %f [sec]", d=(/ my_comp%RestartTimeSec /) )

    !- Initialize information about grid, variables
    
    ! call cal_mn(my_size, DIV_X, DIV_Y)
    my_comp%PRC_NX = 1
    my_comp%GNX    = IMO
    allocate( lnx_field(my_comp%PRC_NX) )
    lnx_field(:)   = IMO
    my_comp%LNX    = IMO
    
    allocate( lny_field(my_comp%PRC_NY) )
    my_comp%PRC_NY = 1  !size(a_jmax)
    my_comp%GNY    = JMO
    lny_field(:)   = JMO !a_jmax(:)
    my_comp%LNY    = JMO !a_jmax(my_comp%PRC_rank)
    
!!$    write(*,*) "ATM: rank=", my_rank, "DIV_X,Y=", DIV_X, DIV_Y, "mysize=", my_size
    

    ! Initialize grid for Jcup
    call init_jcup_grid()

    ! Initialize variables for JCup
    call init_jcup_var()

    !- Initialize an module to interpolate and grid mapping between model components in JCup
    call init_jcup_interpolate()
    
    !
    !

    my_comp%tstep = 1
    
    call jcup_init_time(my_comp%InitTimeInfo)
    
    call output_prepare()
    
    call set_and_put_data(TimeSecN)
    
  end subroutine ocn_init

  subroutine ocn_run(loop_flag)

    ! モジュール引用; Use statements
    !    
    use jcup_interface, only: &
         & jcup_set_time, jcup_inc_time

    !* DCPAM

    ! 宣言文; Declareration statements
    !        
    logical, intent(inout) :: loop_flag

    ! 局所変数
    ! Local variables
    !
!!$    integer, parameter :: end_of_tstep = 2*24*181 + 1  !  6 month
!!$    integer, parameter :: end_of_tstep = 2*24*731 + 1  !  2 year
!!$    integer, parameter :: end_of_tstep = 2*24*1461 + 1 !  4 year
    integer, parameter :: end_of_tstep = 2*24*1826 + 1 !  5 year
!!$    integer, parameter :: end_of_tstep = 2*24*3651 + 1 ! 10 year
!!$    integer, parameter :: end_of_tstep = 2*24*18251 + 1   ! 50 year

    ! 実行文; Executable statement
    !
    
    my_comp%loop_end_flag = .false.
    my_comp%DelTime = TimeSecA - TimeSecN
    
    do while(.not. loop_end_flag)

       call jcup_set_time( my_comp%name,                  & ! (in)
            & my_comp%InitTimeInfo, int(my_comp%DelTime) )  ! (in)

!!$       write(*,*) "* COUPLER Get: atm my_rank=", my_rank, "tstep=", tstep, "time=", tstep*delta_t              
       call get_and_write_data( TimeSecN )

       if (my_comp%PRC_rank==0 .and. mod(my_comp%tstep, 20) == 0) then
          write(*,*) "-> atm my_rank=", my_comp%PRC_rank, "tstep=", my_comp%tstep, "TimeSec=", TimeSecN, &
               & "end_testep=", end_of_tstep
       end if
       
       call agcm_advance_timestep( my_comp%tstep, & ! (in)
            & my_comp%loop_end_flag               & ! (out)
            & )
!!$       write(*,*) "<- atm my_rank=", my_rank, "tstep=", tstep, "time=", tstep*delta_t
!!$       write(*,*) "* COUPLER Put: atm my_rank=", my_rank, "tstep=", tstep, "time=", tstep*delta_t              

       my_comp%tstep = my_comp%tstep + 1       
       call set_and_put_data( TimeSecA ) ! (in)
       call jcup_inc_time( my_comp%name, my_comp%InitTimeInfo ) ! (in)


       if(my_comp%tstep == end_of_tstep) loop_end_flag = .true.
       
    end do
    loop_flag = .false.
    
  end subroutine ocn_run

  subroutine ocn_fin()

    ! モジュール引用; Use statement
    !    
    use jcup_interface, only: &
         & jcup_coupling_end

    ! 実行文; Executable statement
    !
    write(*,*) 'atm fin: my_rank=',my_comp%PRC_rank
    call jcup_coupling_end(my_comp%InitTimeInfo, .false.)
    
    write(*,*) ' = DCPAM fin: my_rank=', my_comp%PRC_rank    
    call agcm_shutdown()
    write(*,*) ' --------- DCPAM fin: my_rank=', my_comp%PRC_rank        

    write(*,*) '-----atm fin: my_rank=', my_comp%PRC_rank
    
  end subroutine ocn_fin

  !- Private subroutines -----------------------------------------------------------------


  !
  !
  subroutine init_jcup_grid()

    ! モジュール引用; Use statement
    !    
    use jcup_interface, only: &
         & jcup_def_grid, jcup_end_grid_def

    use field_def, only: &
         & init_field_def, set_field_def,   &
         & get_local_field, cal_grid_index

    use mod_common_params, only: &
         & ATM_GRID_2D

    use component_field, only: &
         & init_field

    !* Dennou-OGCM


    ! 局所変数
    ! Local variables
    !

    integer, parameter :: ATM_NUM_GRIDTYPE = 1
    integer :: lis
    integer :: lie
    integer :: ljs
    integer :: lje

    ! 宣言文; Declareration statements
    !            
    real(DP), allocatable :: grid_index(:)

    ! 実行文; Executable statement
    !

    call init_field_def(ATM_NUM_GRIDTYPE)

    !-- Define atm_grid_2d
    call set_field_def( my_comp%name, ATM_GRID_2D, & ! (in)
         & my_comp%GNX, my_comp%GNY, 1, 2,         & ! (in) gnx, gny, gnz, halo
         & my_comp%PRC_NX, my_comp%PRC_NY,         & ! (in) 
         & lnx_field, lny_field                    & ! (in)
         & )

    call get_local_field( component_name=my_comp%name, grid_name=ATM_GRID_2D, & ! (in) 
         & local_is=lis, local_ie=lie, local_js=ljs, local_je=lje             & ! (out)
         & )

    write(*,*) "ATM: rank=", my_comp%PRC_rank, "li=", lis, lie, " ,lj=", ljs, lje

    call init_field(field, ATM_GRID_2D, lis, lie, ljs, lje, 1, 1) ! halo=1
    
    ! Skip to define atm_grid_3d ------------
!!$    call set_field_def(ATM, ATM_GRID_3D, GNXA, GNYA, GNZA, 2, DIV_X, DIV_Y)
!!$    call init_field(field3d, ATM_GRID_3D, lis, lie, ljs, lje, 1, GNZA) ! halo=1    


!!$    call gen_grid_index()
    !    call cal_grid_index(ATM, ATM_GRID_2D, field%grid_index)
!    write(*,*) "ATM: grid_index=", field%grid_index

!    call cal_grid_index(ATM, ATM_GRID_3D, field3d%grid_index)

    call jcup_def_grid(field%grid_index, my_comp%name, ATM_GRID_2D, GN25)
!    call jcup_def_grid(field3d%grid_index, ATM, ATM_GRID_3D)

    !---- Finish defining grid for jcup ----------------------------------------
    call jcup_end_grid_def()
    
  end subroutine init_jcup_grid

  !---------------------------------------------------------------------------------------------------
  
  subroutine init_jcup_var()

    ! モジュール引用; Use statement
    !        
    use jcup_interface, only: &
         & jcup_def_varp, jcup_def_varg, jcup_end_var_def

    use component_field, only: &
         & init_field_data

    use mod_common_params
    
    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !


    !
    call init_field_data( field, num_of_25d=GN25,                      & ! (in)
         & num_of_varp=NUM_OCN_PUTVAR2D, num_of_varg=NUM_OCN_GETVAR2D )  ! (in)
    
!!$    call init_field_data(field3d, num_of_25d=1, num_of_varp=1, num_of_varg=1)    

    !- Variable put by my own component

    call jcup_def_varp( field%varp(o2d_SfcTemp_id)%varp_ptr, my_comp%name, o2d_SfcTemp, OCN_GRID_2D )
    call jcup_def_varp( field%varp(o2a_SfcAlbedo_id)%varp_ptr, my_comp%name, o2d_SfcAlbedo, OCN_GRID_2D )
    call jcup_def_varp( field%varp(o2d_SfcSnow_id)%varp_ptr, my_comp%name, o2d_SfcSnow, OCN_GRID_2D )
    
    !- Variable getten from  other components

    call jcup_def_varg( field%varg(a2o_WindStressX_id)%varg_ptr, my_comp%name, 'a2o_WindStressX', OCN_GRID_2D, 1, & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_WindStressX,                                          & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1 )            ! (in)

    call jcup_def_varg( field%varg(a2o_WindStressY_id)%varg_ptr, my_comp%name, 'a2o_WindStressY', OCN_GRID_2D, 1, & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_WindStressY,                                          & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1 )            ! (in)


    call jcup_def_varg( field%varg(a2o_LDwRFlx_id)%varg_ptr, my_comp%name, 'a2o_LDwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_LDwRFlx,                                       & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_SDwRFlx_id)%varg_ptr, my_comp%name, 'a2o_SDwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_SDwRFlx,                                       & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_LUwRFlx_id)%varg_ptr, my_comp%name, 'a2o_LUwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_LUwRFlx,                                       & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_SUwRFlx_id)%varg_ptr, my_comp%name, 'a2o_SUwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_SUwRFlx,                                       & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    
    call jcup_def_varg( field%varg(a2o_LatHFlx_id)%varg_ptr, my_comp%name, 'a2o_LatHFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_LatHFlx,                                       & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=3 )     ! (in)

    call jcup_def_varg( field%varg(a2o_SenHFlx_id)%varg_ptr, my_comp%name, 'a2o_SenHFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_SenHFlx,                                       & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=3 )     ! (in)

    call jcup_def_varg( field%varg(a2o_DSfcHFlxDTs_id)%varg_ptr, my_comp%name, 'a2o_DSfcHFlxDTs', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_LatHFlx,                                               & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=3 )             ! (in)
    

    call jcup_def_varg( field%varg(a2o_SnowFall_id)%varg_ptr, my_comp%name, 'a2o_SnowFall', OCN_GRID_2D, 1,   & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_SnowFall,                                         & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=4 )        ! (in)

    call jcup_def_varg( field%varg(a2o_RainFall_id)%varg_ptr, my_comp%name, 'a2o_RainFall', OCN_GRID_2D, 1,   & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_RainFall,                                         & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=4 )        ! (in)

    call jcup_def_varg( field%varg(a2o_SfcAirTemp_id)%varg_ptr, my_comp%name, 'a2o_SfcAirTemp', OCN_GRID_2D, 1, & ! (in)
         & SEND_MODEL_NAME=my_comp%name, SEND_DATA_NAME=a2d_SfcAirTemp,                                         & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=4 )          ! (in)
    

    !- Finish defining variables for jup ---------------------
    
    call jcup_end_var_def()
    
  end subroutine init_jcup_var

  !---------------------------------------------------------------------------------------------------
  
  subroutine init_jcup_interpolate()

    ! モジュール引用; Use statement
    !            
    use jcup_interface, only: &
         & jcup_set_mapping_table

    use jcup_mpi_lib, only: &
         & jml_finalize

    use field_def, only : cal_mn, &
         & set_grid_mapping, set_grid_mapping_3d

    use interpolation_data_latlon_mod, only: &
         & init_interpolation => interpolation_data_latlon_Init, &
         & set_operation_index, &
         & set_A_to_O_coef, set_O_to_A_coef

    use grid_mapping_util, only: &
         & set_mappingTable_interpCoef

    use mod_common_params
    
    ! 局所変数
    ! Local variables
    !
    integer, allocatable :: send_grid_ao(:)
    integer, allocatable :: send_grid_oa(:)
    integer, allocatable :: recv_grid_ao(:)
    integer, allocatable :: recv_grid_oa(:)
    real(DP), allocatable :: coefS_ao_global(:)
    real(DP), allocatable :: coefS_oa_global(:)
    
    ! 実行文; Executable statement
    !
    
    !
    call init_interpolation(NUM_DCCM_COMP, NUM_ATM_GMAPTAG, my_comp%id)

    
    ! ATM -> OCN grid mapping    **************************x
    if(my_comp%PRC_rank==0) then
       call set_mappingTable_interpCoef( &
            & GMAPFILENAME_AO, my_comp%GNX, my_comp%GNX,          & ! (in)
            & send_grid_ao, recv_grid_ao, coefS_ao_global          & ! (inout)
            & )
       write(*,*) "A20:", size(send_grid_ao), size(recv_grid_ao), size(coefS_ao_global)
    end if
    call jcup_set_mapping_table( my_comp%name,                       & ! (in)
         & my_comp%name, ATM_GRID_2D, my_comp%name, OCN_GRID_2D,    & ! (in) ATM_GRID_2D -> OCN_GRID_2D
         & GMAPTAG_ATM2D_OCN2D, send_grid_ao, recv_grid_ao )           ! (in)

    ! OCN -> ATM grid mapping   *****************************    
    if(my_comp%PRC_rank==0) then
       call set_mappingTable_interpCoef( &
            & GMAPFILENAME_OA, my_comp%GNX, my_comp%GNX,      & ! (in)
            & send_grid_oa, recv_grid_oa, coefS_oa_global      & ! (inout)
            & )
       write(*,*) "send_grid_oa:", send_grid_oa
       write(*,*) "recv_grid_oa:", recv_grid_oa       
    end if
    call jcup_set_mapping_table( my_comp%name,                       & ! (in)
         & my_comp%name, OCN_GRID_2D, my_comp%name, ATM_GRID_2D,    & ! (in) OCN_GRID_2D -> ATM_GRID_2D
         & GMAPTAG_ATM2D_OCN2D, send_grid_oa, recv_grid_oa )           ! (in)

    
    call set_operation_index(my_comp%name, my_comp%name, GMAPTAG_ATM2D_OCN2D)         ! (in)

    if(my_comp%PRC_rank==0) then
       call set_A_to_O_coef(GMAPTAG_ATM2D_OCN2D, coefS_ao_global)   ! (in)
    else
       call set_A_to_O_coef(GMAPTAG_ATM2D_OCN2D)                    ! (in)
    end if

    if(my_comp%PRC_rank==0) then
       call set_O_to_A_coef(GMAPTAG_ATM2D_OCN2D, coefS_oa_global)   ! (in)
    else
       call set_O_to_A_coef(GMAPTAG_ATM2D_OCN2D)                    ! (in)
    end if

  end subroutine init_jcup_interpolate

  !----------------------------------------------------------------------------------
  
  subroutine set_and_put_data( CurrentTimeSec )

    ! モジュール引用; Use statement
    !                

    use mpi
    
    !* Dennou-OGCM


    use DSIce_Admin_Constants_mod, only: &
         & AlbedoOcean, IceMaskMin

    use DSIce_Admin_Variable_mod, only: &
         & xya_SIceCon, xya_SIceSfcTemp, &
         & xya_IceThick
    
    use DSIce_Boundary_vars_mod, only: &
         & xy_SfcAlbedoAI, &
         & xy_SeaSfcTemp

    
    !* DCCM
    
    use jcup_interface, only: &
         & jcup_put_data

    use field_def, only: &
         & set_send_data_2d    

    use mod_common_params


    ! 宣言文; Declareration statements
    !            
    real(DP), intent(in) :: CurrentTimeSec

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xy_SfcAlbedo4Atm(IAO,JAO)
    real(DP) :: xy_SfcTemp4Atm(IAO,JAO)
    real(DP) :: xy_SfcSnow4Atm(IAO,JAO)
    integer :: i
    integer :: j
    real(DP), parameter :: PI = acos(-1d0)

    
    ! 実行文; Executable statement
    !

    !$omp parallel do private(i)
    do j = JSO, JEO
       do i = ISO, IEO
          if ( xya_SIceCon(i,j,SICE_TLN) > IceMaskMin ) then
             xy_SfcAlbedo4Atm(i,j) = xy_SfcAlbedoAI(i,j)
             xy_SfcTemp4Atm(i,j) = xya_SIceSfcTemp(i,j,SICE_TLN)
             xy_SfcSnow4Atm(i,j) = xya_IceThick(i,j,SICE_TLN)
          else
             xy_SfcAlbedo4Atm(i,j) = AlbedoOcean
             xy_SfcTemp4Atm(i,j) = xy_SeaSfcTemp(i,j)
             xy_SfcSnow4Atm(i,j) = 0d0
          end if
       end do
    end do
    
    call ocn_set_send_2d( o2a_SfcTemp_id, xy_SfcTemp4Atm(ISO:IEO,JSO:JEO) )
    call ocn_set_send_2d( o2a_SfcAlbedo_id, xy_SfcAlbedo4Atm(ISO:IEO,JSO:JEO) )
    call ocn_set_send_2d( o2a_SfcSnow_id, xy_SfcSnow4Atm(ISO:IEO,JSO:JEO) )
    
  contains
    subroutine ocn_set_send_2d(varpID, send_data)
      integer, intent(in) :: varpID
      real(DP), intent(in) :: send_data(:,:)
      
      call jcup_put_data(field%varp(varpID)%varp_ptr, pack(send_data, maSK=field%mask2d))
    end subroutine ocn_set_send_2d
    
!!$    subroutine ocn_set_send(varpID, code)
!!$      integer, intent(in) :: varpID, code
!!$      
!!$      call set_send_data_2d(my_comp%name, OCN_GRID_2D, field%send_2d(:,:), step, code)
!!$      call jcup_put_data(field%varp(varpID)%varp_ptr, pack(field%send_2d, Mask=field%mask2d))
!!$    end subroutine ocn_set_send
  end subroutine set_and_put_data

  subroutine get_and_write_data( CurrentTimeSec )

    ! モジュール引用; Use statement
    !                

    !* DCCM
    
    use jcup_interface, only: &
         & jcup_get_data

    use field_def, only: &
         & write_data_2d

    use mod_common_params
    
    ! 宣言文; Declareration statements
    !            
    real(DP), intent(in) :: CurrentTimeSec

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: xy_SfcTemp(IAO,JAO,KAO)
    real(DP) :: xy_SfcAlbedo(IAO,JAO,KAO)
    real(DP) :: xy_SfcSnow(IAO,JAO,KAO)

    ! 実行文; Executable statement
    !
    
    
    ! Get oceanic surface temerature send by OGCM.
    !

    call ocn_get_write( o2a_SfcTemp_id, o2d_SfcTemp,   & ! (in)
         & xy_SfcTemp )                                  ! (out)

    call ocn_get_write( o2a_SfcAlbedo_id, o2d_SfcAlbedo, & ! (in)
         & xy_SfcAlbedo )                                  ! (out)

    call ocn_get_write( o2a_SfcSnow_id, o2d_SfcSnow,   & ! (in)
         & xy_SfcSnow )                                  ! (out)


    !
    !
    if( mod(CurrentTimeSec, dble(AO_COUPLING_CYCLE_SEC)) == 0d0 ) then

!!$    if(my_rank>12) then
!!$       write(*,*) "atm: rank=", my_rank, "SurfTemp=", xy_SurfTemp(0,1:2), "lat=", y_Lat(1:2)/acos(-1d0)*180d0
!!$    end if
       
    end if

  contains
    subroutine ocn_get_write(vargID, vargName, xy_getdata)
      integer, intent(in) :: vargID
      character(*), intent(In) :: vargname
      real(DP), intent(out) :: xy_getdata(IAO,JAO)
      
      field%buffer1d(:) = 0d0
      call jcup_get_data(field%varg(vargID)%varg_ptr, field%buffer1d)

      if( mod(CurrentTimeSec, dble(AO_COUPLING_CYCLE_SEC)) == 0d0 ) then      
         field%recv_2d(:,:) = unpack(field%buffer1d, field%mask2d, field%recv_2d)      
         xy_getdata(ISO:IEO,JSO:JEO) = field%recv_2d
         call output_var( CurrentTimeSec, vargName, xy_getdata )
      end if
      
    end subroutine ocn_get_write
    
  end subroutine get_and_write_data

  !-----------------------------------------------------------------------------
  
  subroutine output_prepare()
    
    use gtool_historyauto, only: &
         & HistoryAutoAddVariable

    ! 局所変数
    ! Local variables
    !    
    character(TOKEN) :: dims_XYT(3)

    ! 実行文; Executable statement
    !
    
    dims_XYT = (/ 'lon ', 'lat ', 'time'  /)
    call HistoryAutoAddVariable('SurfTempOcn', &
         & dims=dims_XYT, longname='surface temperature calculated ocean and sea-ice model', units='K') 

    call HistoryAutoAddVariable('SurfAlbedoOcn', &
         & dims=dims_XYT, longname='surface albedo calculated by ocean and sea-ice model', units='1') 

    call HistoryAutoAddVariable('SurfSnowOcn', &
         & dims=dims_XYT, longname='surface snode depth  calculated by ocean and sea-ice model', units='m') 
    
    call HistoryAutoAddVariable('CheckVar', &
         & dims=dims_XYT, longname='CheckVar', units='1') 

  end subroutine output_prepare

  subroutine output_var(CurrentTime, varname, data2d)

    use gtool_historyauto, only: &
         & HistoryAutoPut
    
    real(DP), intent(in) :: CurrentTime
    character(*), intent(in) :: varName
    real(DP), intent(in) :: data2d(IAO,JAO)

    call HistoryAutoPut(CurrentTime, varname, data2d(ISO:IEO,JSO:JEO))
    
  end subroutine output_var
  !----------------------------------------------------------------------------------------------

!!$  subroutine gen_grid_index()
!!$
!!$    ! モジュール引用; Use statement
!!$    !            
!!$    
!!$    !* DCPAM5
!!$    use gridset, only: a_jmax
!!$
!!$
!!$    ! 局所変数
!!$    ! Local variables
!!$    !    
!!$    integer :: my_rank
!!$    integer :: g_js
!!$    integer :: counter
!!$
!!$    integer :: i
!!$    integer :: j
!!$    
!!$    ! 実行文; Executable statement
!!$    !
!!$
!!$    my_rank = my_comp%PRC_rank
!!$    
!!$    g_js =   sum(a_jmax)/2                      &
!!$         & - sum(a_jmax(0:my_rank))/2
!!$
!!$!    write(*,*) "ATM grid_index: (rank=", my_rank, ") :: jc", jc, "a_jmax=", a_jmax, "lb;", lbound(a_jmax)    write(*,*) "ATM grid_index: (rank=", my_rank, ") :: g_js", g_js
!!$    counter = 0
!!$    do j=1, a_jmax(my_rank)/2
!!$       do i=1, my_comp%GNX
!!$          counter = counter + 1          
!!$          field%grid_index(counter) = i + (g_js + j - 1)*my_comp%GNX
!!$       end do
!!$    end do
!!$
!!$    g_js = sum(a_jmax)/2
!!$    if( my_rank > 0 ) then
!!$       g_js = g_js + sum(a_jmax(0:my_rank-1))/2
!!$    end if
!!$    do j=1, a_jmax(my_rank)/2
!!$       do i=1, my_comp%GNX
!!$          counter = counter + 1          
!!$          field%grid_index(counter) = i + (g_js + j - 1)*my_comp%GNX
!!$       end do
!!$    end do
!!$
!!$    write(*,*) "ATM grid_index: (rank=", my_rank, ") ::", field%grid_index
!!$  end subroutine gen_grid_index
!!$

end module mod_ocn

