!-------------------------------------------------------------
! Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module mod_atm

  ! モジュール引用; Use statements
  !

  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify
  
  !* DCPAM
  
  use dcpam_main_mod, only: &
       & agcm_main_init => dcpam_main_Init,               &
       & agcm_main_final => dcpam_main_Final,             &
       & agcm_setup => MainInit,                          &
       & agcm_shutdown => MainTerminate,                  &       
       & agcm_advance_timestep => dcpam_advance_timestep, &
       & agcm_update_surfprop => dcpam_UpdateSurfaceProperties

  use timeset, only: &
       & TimeSecB => TimeB, &
       & TimeSecN => TimeN, &
       & TimeSecA => TimeA, &
       & EndTimeSec => EndTime
  
  !* DCPCM

  use component_field, only: &
       & component_field_type
  

  use mod_common_params, only: &
       & DEFAULT_DCCM_CONFNAME,     &
       & NUM_DCCM_COMP,             &
       & NUM_ATM_GMAPTAG,           &
       & COMPNAME_ATM, GN25,        &
       & ATM_GRID_2D, OCN_GRID_2D,  &
       & GMAPTAG_ATM2D_OCN2D
       
  use mod_common_compdef, only:     &
       & ComponentDef_Init, ComponentDef_Final, &
       & ComponentDef_Share,                    &
       & common_read_config,                    &
       & my_comp   => CompDef_atm,              &
       & ocn_comp  => CompDef_ocn,              &
       & sice_comp => CompDef_sice,             &
       & GMAPFILENAME_AO,                       &
       & GMAPFILENAME_OA,                       &
       & AO_COUPLING_CYCLE_SEC
  
  ! 宣言文; Declareration statements
  !    
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !    
  public :: atm_init
  public :: atm_run
  public :: atm_fin

  ! 非公開変数
  ! Private variable
  !
  
  integer :: tstep
  
  type(component_field_type) :: field
  type(component_field_type) :: field3d

  integer, allocatable :: lnx_field(:)
  integer, allocatable :: lny_field(:)
  
!!$  real(DP), allocatable :: xy_TmpTAvg(:,:)

  character(*), parameter :: module_name = 'mod_atm'
  
contains

  subroutine atm_init()

    ! モジュール引用; Use statements
    !
    use jcup_interface, only: &
         & jcup_init_time

    use field_def, only: &
         & cal_mn

    !* dcpam5
    
    use gridset, only: &
         & a_jmax, iMax, jMax

    use timeset, only: &
         & InitialDate, DelTime, RestartTime
    
    ! 宣言文; Declareration statements
    !    

    
    ! 実行文; Executable statement
    !
    
    ! Initialze varibale to manage component  ********************
    !

    call common_read_config( DEFAULT_DCCM_CONFNAME )
    
    call ComponentDef_Init( my_comp,         & ! (inout)
         & COMPNAME_ATM, "DCPAM" )             ! (in)
    

    !  Intialize a model associated with  atmosphere component ***
    !
    call MessageNotify('M', 'atm component', &
         & 'Initialize %a .. (MPI_MY_COMM=%d, MY_RANK=%d)', &
         & ca=(/ trim(my_comp%MODELNAME) /),  i=(/ my_comp%PRC_comm, my_comp%PRC_rank/) &
         & )

    call agcm_main_init( my_comp%PRC_comm )
    call agcm_setup()

    my_comp%InitTimeInfo(:) = (/ InitialDate%year, InitialDate%month, InitialDate%day,        &
                        &        InitialDate%hour, InitialDate%min, int(InitialDate%sec)  /)
    my_comp%RestartTimeSec  = RestartTime
    my_comp%DelTime         = DelTime
    my_comp%tstep           = 0    
    call MessageNotify( "M", module_name, &
         & "Inital Time %4d - %2d - %2d:%2d%2d", i=my_comp%InitTimeInfo )
    call MessageNotify( "M", module_name, &
         & "RestartTiem = %f [sec]", d=(/ my_comp%RestartTimeSec /) )


    !- Initialize information about grid, variables
    
    ! call cal_mn(my_size, DIV_X, DIV_Y)
    my_comp%PRC_NX = 1
    my_comp%GNX    = iMax
    allocate( lnx_field(my_comp%PRC_NX) )
    lnx_field(:)   = iMax
    my_comp%LNX    = iMax
    
    my_comp%PRC_NY = size(a_jmax)
    my_comp%GNY    = jMax    
    allocate( lny_field(my_comp%PRC_NY) )
    lny_field(:)   = a_jmax(:)
    my_comp%LNY    = a_jmax(my_comp%PRC_rank)
    
!!$    write(*,*) "ATM: rank=", my_rank, "DIV_X,Y=", DIV_X, DIV_Y, "mysize=", my_size

    call ComponentDef_Share()

    ! Initialize grid for Jcup
    call MessageNotify( 'M', module_name, "Initialize grid for Jcup..")
    call init_jcup_grid()

    ! Initialize variables for JCup
    call MessageNotify( 'M', module_name, "Initialize varibles for Jcup..")    
    call init_jcup_var()

    !- Initialize an module to interpolate and grid mapping between model components in JCup
    call MessageNotify( 'M', module_name, "Initialize interpolation data for Jcup ..")
    call init_jcup_interpolate()

    !
    !


    
    call jcup_init_time(my_comp%InitTimeInfo)

    call MessageNotify( 'M', module_name, "Prepare data output..")    
    call output_prepare()

    my_comp%tstep = 0    
    call set_and_put_data( TimeSecN )
    my_comp%tstep = 1
    
    call MessageNotify( 'M', module_name, "atm_init has been finished.")
    
  end subroutine atm_init

  subroutine atm_run(loop_flag)

    ! モジュール引用; Use statements
    !    
    use jcup_interface, only: &
         & jcup_set_time, jcup_inc_time

    !* Dennou-OGCM
    

    ! 宣言文; Declareration statements
    !        
    logical, intent(inout) :: loop_flag

    ! 局所変数
    ! Local variables
    !

    ! 実行文; Executable statement
    !
    
    my_comp%loop_end_flag = .false.
    my_comp%DelTime = TimeSecA - TimeSecN
    
    do while(.not. my_comp%loop_end_flag)

       call jcup_set_time( my_comp%name,                  & ! (in)
            & my_comp%InitTimeInfo, int(my_comp%DelTime) )  ! (in)

       call get_and_write_data( TimeSecN )

       
       if (my_comp%PRC_rank==0 .and. mod(my_comp%tstep, 100) == 0) then
          call MessageNotify( 'M', module_name,                                &
               & "TimeSecN=%f, EndTimeSec=%f",  d=(/  TimeSecN, EndTimeSec /)  &
               & )
       end if

       !------------------------------------------------------
          
       call agcm_advance_timestep( my_comp%tstep, & ! (in)
            & my_comp%loop_end_flag               & ! (out)
            & )

       !------------------------------------------------------
       
       call set_and_put_data( TimeSecN ) ! (in)
       call jcup_inc_time( my_comp%name, my_comp%InitTimeInfo ) ! (in)
       my_comp%tstep = my_comp%tstep + 1       


       if ( EndTimeSec < TimeSecN ) then
          my_comp%loop_end_flag = .true.
       end if
       
    end do
    loop_flag = .false.
    
  end subroutine atm_run

  subroutine atm_fin()

    ! モジュール引用; Use statement
    !    
    use jcup_interface, only: &
         & jcup_coupling_end

    ! 実行文; Executable statement
    !

    call jcup_coupling_end(my_comp%InitTimeInfo, .false.)

    !------------------------
    call MessageNotify( 'M', module_name, "Shutdown atmospheric model..")
    call agcm_shutdown()
    !-------------------------
    
    call MessageNotify( 'M', module_name, "atm_fin has been finished. (rank=%d)", &
         & i=(/ my_comp%PRC_rank  /) )
    
  end subroutine atm_fin

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
         & ATM_NUM_GRIDTYPE, &
         & ATM_GRID_2D

    use component_field, only: &
         & init_field

    !* DCPAM5
    use gridset, only: a_jmax

    ! 局所変数
    ! Local variables
    !

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


    call gen_grid_index()
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
         & num_of_varp=NUM_ATM_PUTVAR2D, num_of_varg=NUM_ATM_GETVAR2D )  ! (in)
    
!!$    call init_field_data(field3d, num_of_25d=1, num_of_varp=1, num_of_varg=1)    

    !- Variable put by my own component

    call jcup_def_varp( field%varp(a2d_WindStressX_id)%varp_ptr, my_comp%name, a2d_WindStressX, ATM_GRID_2D )
    call jcup_def_varp( field%varp(a2d_WindStressY_id)%varp_ptr, my_comp%name, a2d_WindStressY, ATM_GRID_2D )

    call jcup_def_varp( field%varp(a2d_LDwRFlx_id)%varp_ptr, my_comp%name, a2d_LDwRFlx, ATM_GRID_2D ) 
    call jcup_def_varp( field%varp(a2d_SDwRFlx_id)%varp_ptr, my_comp%name, a2d_SDwRFlx, ATM_GRID_2D )
    call jcup_def_varp( field%varp(a2d_LUwRFlx_id)%varp_ptr, my_comp%name, a2d_LUwRFlx, ATM_GRID_2D )
    call jcup_def_varp( field%varp(a2d_SUwRFlx_id)%varp_ptr, my_comp%name, a2d_SUwRFlx, ATM_GRID_2D )

    call jcup_def_varp( field%varp(a2d_LatHFlx_id)%varp_ptr, my_comp%name, a2d_LatHFlx, ATM_GRID_2D )
    call jcup_def_varp( field%varp(a2d_SenHFlx_id)%varp_ptr, my_comp%name, a2d_SenHFlx, ATM_GRID_2D )
    call jcup_def_varp( field%varp(a2d_DSfcHFlxDTs_id)%varp_ptr, my_comp%name, a2d_DSfcHFlxDTs, ATM_GRID_2D )
    
    call jcup_def_varp( field%varp(a2d_RainFall_id)%varp_ptr, my_comp%name, a2d_RainFall, ATM_GRID_2D )
    call jcup_def_varp( field%varp(a2d_SnowFall_id)%varp_ptr, my_comp%name, a2d_SnowFall, ATM_GRID_2D )

    call jcup_def_varp( field%varp(a2d_SfcAirTemp_id)%varp_ptr, my_comp%name, a2d_SfcAirTemp, ATM_GRID_2D )
    
    
    !- Variable getten from  other components

    call jcup_def_varg( field%varg(o2a_SfcTemp_id)%varg_ptr, my_comp%name, 'o2a_SfcTemp', ATM_GRID_2D, 1,   & ! (in)
         & SEND_MODEL_NAME=ocn_comp%name, SEND_DATA_NAME=o2d_SfcTemp,                                       & ! (in)
         & RECV_MODE="SNP", INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1)       ! (in)

    call jcup_def_varg( field%varg(o2a_SfcAlbedo_id)%varg_ptr, my_comp%name, 'o2a_SfcAlbedo', ATM_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=ocn_comp%name, SEND_DATA_NAME=o2d_SfcAlbedo,                                        & ! (in)
         & RECV_MODE="SNP", INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1)          ! (in)

    call jcup_def_varg(field%varg(o2a_SfcSnow_id)%varg_ptr, my_comp%name, 'o2a_SfcSnow', ATM_GRID_2D, 1,       & ! (in)
         & SEND_MODEL_NAME=ocn_comp%name, SEND_DATA_NAME=o2d_SfcSnow,                                          & ! (in)
         & RECV_MODE="SNP", INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1)          ! (in)


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
            & GMAPFILENAME_AO, my_comp%GNX, ocn_comp%GNX,          & ! (in)
            & send_grid_ao, recv_grid_ao, coefS_ao_global          & ! (inout)
            & )
       write(*,*) "A2O:", size(send_grid_ao), size(recv_grid_ao), size(coefS_ao_global)
    end if
    call jcup_set_mapping_table( my_comp%name,                       & ! (in)
         & my_comp%name, ATM_GRID_2D, ocn_comp%name, OCN_GRID_2D,    & ! (in) ATM_GRID_2D -> OCN_GRID_2D
         & GMAPTAG_ATM2D_OCN2D, send_grid_ao, recv_grid_ao )           ! (in)

    ! OCN -> ATM grid mapping   *****************************    
    if(my_comp%PRC_rank==0) then
       call set_mappingTable_interpCoef( &
            & GMAPFILENAME_OA, ocn_comp%GNX, my_comp%GNX,      & ! (in)
            & send_grid_oa, recv_grid_oa, coefS_oa_global      & ! (inout)
            & )
       write(*,*) "send_grid_oa:", send_grid_oa
       write(*,*) "recv_grid_oa:", recv_grid_oa       
    end if
    call jcup_set_mapping_table( my_comp%name,                       & ! (in)
         & ocn_comp%name, OCN_GRID_2D, my_comp%name, ATM_GRID_2D,    & ! (in) OCN_GRID_2D -> ATM_GRID_2D
         & GMAPTAG_ATM2D_OCN2D, send_grid_oa, recv_grid_oa )           ! (in)

    
    call set_operation_index(my_comp%name, ocn_comp%name, GMAPTAG_ATM2D_OCN2D)         ! (in)
    
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
    
    !* dcpam5
    
    use axesset, only: x_Lon, y_Lat
    use gridset,only: imax, jmax
    use constants, only: RPlanet

    use dcpam_main_mod, only: &
         & xy_TauXAtm, xy_TauYAtm, xy_SensAtm, xy_LatentAtm, &
         & xy_LDWRFlxAtm, xy_LUWRFlxAtm, xy_SDWRFlxAtm, xy_SUWRFlxAtm, &
         & xy_SurfAirTemp, xy_DSurfHFlxDTs, xy_DSurfLatentFlxDTs,      &
         & xy_RainAtm => xy_Rain, xy_SnowAtm => xy_Snow, &
         & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux, xy_SurfTemp
!         & xy_RainAtm, xy_SnowAtm

    use intavr_operate, only: &
         & IntLonLat_xy

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
    integer :: p, j, m, n
    real(DP) :: Tmpavg, TmpAvgGlobal, SurfArea
    integer :: ierr
    real(DP), parameter :: PI = acos(-1d0)

    ! 実行文; Executable statement
    !
    
    call atm_set_send_2d( a2d_WindStressX_id, -xy_TauXAtm )
    call atm_set_send_2d( a2d_WindStressY_id, -xy_TauYAtm )    

    call atm_set_send_2d( a2d_LDwRFlx_id, xy_LDWRFlxAtm )    
    call atm_set_send_2d( a2d_SDwRFlx_id, xy_SDWRFlxAtm )
    call atm_set_send_2d( a2d_LUwRFlx_id, xy_LUWRFlxAtm )    
    call atm_set_send_2d( a2d_SUwRFlx_id, xy_SUWRFlxAtm )
    call atm_set_send_2d( a2d_LatHFlx_id, xy_LatentAtm )    
    call atm_set_send_2d( a2d_SenHFlx_id, xy_SensAtm )
    call atm_set_send_2d( a2d_DSfcHFlxDTs_id, xy_DSurfHFlxDTs )
    
    call atm_set_send_2d( a2d_RainFall_id, xy_RainAtm )
    call atm_set_send_2d( a2d_SnowFall_id, xy_SnowAtm )

    call atm_set_send_2d( a2d_SfcAirTemp_id, xy_SurfAirTemp )
    
  contains
    subroutine atm_set_send_2d(varpID, send_data)
      integer, intent(in) :: varpID
      real(DP), intent(in) :: send_data(:,:)
      
      call jcup_put_data(field%varp(varpID)%varp_ptr, pack(send_data, maSK=field%mask2d))
    end subroutine atm_set_send_2d
    
!!$    subroutine atm_set_send(varpID, code)
!!$      integer, intent(in) :: varpID, code
!!$      
!!$      call set_send_data_2d(my_comp%name, ATM_GRID_2D, field%send_2d(:,:), step, code)
!!$      call jcup_put_data(field%varp(varpID)%varp_ptr, pack(field%send_2d, Mask=field%mask2d))
!!$    end subroutine atm_set_send
  end subroutine set_and_put_data

  subroutine get_and_write_data( CurrentTimeSec )

    ! モジュール引用; Use statement
    !                

    !* dcpam5

    use gridset, only: &
         & iMax, jMax

    use axesset,only: &
         & y_Lat

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
    real(DP) :: xy_SfcTemp(0:iMax-1,jMax)
    real(DP) :: xy_SfcAlbedo(0:iMax-1,jMax)
    real(DP) :: xy_SfcSnow(0:iMax-1,jMax)

    ! 実行文; Executable statement
    !
    
    
    ! Get oceanic surface temerature send by OGCM.
    !

    call atm_get_write( o2a_SfcTemp_id, o2d_SfcTemp,     & ! (in)
         & xy_SfcTemp )                                    ! (out)
    
    call atm_get_write( o2a_SfcAlbedo_id, o2d_SfcAlbedo, & ! (in)
         & xy_SfcAlbedo )                                  ! (out)

    call atm_get_write( o2a_SfcSnow_id, o2d_SfcSnow,     & ! (in)
         & xy_SfcSnow )                                    ! (out)


    !
    !
    if( mod(CurrentTimeSec, dble(AO_COUPLING_CYCLE_SEC)) == 0d0 ) then

!!$    if(my_rank>12) then
!!$       write(*,*) "atm: rank=", my_rank, "SurfTemp=", xy_SurfTemp(0,1:2), "lat=", y_Lat(1:2)/acos(-1d0)*180d0
!!$    end if

!!$       write(*,*) "Check the unit of SfcSnow.."
       
       call agcm_update_surfprop( &
            & xy_SurfTempRecv=xy_SfcTemp, xy_SurfAlbedoRecv=xy_SfcAlbedo,  & ! (in)
            & xy_SurfSnowRecv=1d3*xy_SfcSnow                               & ! (in)
            & )
    end if

  contains
    subroutine atm_get_write(vargID, vargName, xy_getdata)
      integer, intent(in) :: vargID
      character(*), intent(In) :: vargname
      real(DP), intent(inout) :: xy_getdata(:,:)

      
      field%buffer1d(:) = 0d0
      call jcup_get_data(field%varg(vargID)%varg_ptr, field%buffer1d)
      
      field%recv_2d(:,:) = unpack(field%buffer1d, field%mask2d, field%recv_2d)      
      if( mod(CurrentTimeSec, dble(AO_COUPLING_CYCLE_SEC)) == 0d0 ) then
         xy_getdata(:,:) = field%recv_2d
         call output_var( CurrentTimeSec, vargName, xy_getdata )
      end if
    end subroutine atm_get_write
    
  end subroutine get_and_write_data

  !-----------------------------------------------------------------------------
  
  subroutine output_prepare()

    ! モジュール引用; Use statement
    !                    
    use gtool_historyauto, only: &
         & HistoryAutoAddVariable

    use mod_common_params
    
    ! 局所変数
    ! Local variables
    !    
    character(TOKEN) :: dims_XYT(3)

    ! 実行文; Executable statement
    !
    
    dims_XYT = (/ 'lon ', 'lat ', 'time'  /)
    call HistoryAutoAddVariable( o2d_SfcTemp,  &
         & dims=dims_XYT, longname='surface temperature calculated ocean and sea-ice model', units='K') 

    call HistoryAutoAddVariable( o2d_SfcAlbedo, &
         & dims=dims_XYT, longname='surface albedo calculated by ocean and sea-ice model', units='1') 

    call HistoryAutoAddVariable( o2d_SfcSnow, &
         & dims=dims_XYT, longname='surface snode depth  calculated by ocean and sea-ice model', units='m') 
    
    call HistoryAutoAddVariable('CheckVar', &
         & dims=dims_XYT, longname='CheckVar', units='1') 

  end subroutine output_prepare

  subroutine output_var(CurrentTime, varname, data2d)

    use gtool_historyauto, only: &
         & HistoryAutoPut
    
    real(DP), intent(in) :: CurrentTime
    character(*), intent(in) :: varName
    real(DP), intent(in) :: data2d(:,:)

    call HistoryAutoPut(CurrentTime, varname, data2d)
    
  end subroutine output_var
  !----------------------------------------------------------------------------------------------

  subroutine gen_grid_index()

    ! モジュール引用; Use statement
    !            
    
    !* DCPAM5
    use gridset, only: a_jmax


    ! 局所変数
    ! Local variables
    !    
    integer :: my_rank
    integer :: g_js
    integer :: counter

    integer :: i
    integer :: j
    
    ! 実行文; Executable statement
    !

    my_rank = my_comp%PRC_rank
    
    g_js =   sum(a_jmax)/2                      &
         & - sum(a_jmax(0:my_rank))/2

!    write(*,*) "ATM grid_index: (rank=", my_rank, ") :: jc", jc, "a_jmax=", a_jmax, "lb;", lbound(a_jmax)    write(*,*) "ATM grid_index: (rank=", my_rank, ") :: g_js", g_js
    counter = 0
    do j=1, a_jmax(my_rank)/2
       do i=1, my_comp%GNX
          counter = counter + 1          
          field%grid_index(counter) = i + (g_js + j - 1)*my_comp%GNX
       end do
    end do

    g_js = sum(a_jmax)/2
    if( my_rank > 0 ) then
       g_js = g_js + sum(a_jmax(0:my_rank-1))/2
    end if
    do j=1, a_jmax(my_rank)/2
       do i=1, my_comp%GNX
          counter = counter + 1          
          field%grid_index(counter) = i + (g_js + j - 1)*my_comp%GNX
       end do
    end do

    write(*,*) "ATM grid_index: (rank=", my_rank, ") ::", field%grid_index
  end subroutine gen_grid_index

end module mod_atm

