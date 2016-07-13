!-------------------------------------------------------------
! Copyright (c) 2015-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief A module for ocean component
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
       & ogcm_main_init => DOGCM_main_Init,                   &
       & ogcm_main_final => DOGCM_main_Final,                 &
       & ogcm_setup => DOGCM_main_setup,                      &
       & ogcm_shutdown => DOGCM_main_shutdown,                &       
       & ogcm_advance_timestep => DOGCM_main_advance_timestep

  use DSIce_main_mod, only: &
       & sice_main_init => DSIce_main_Init,                   &
       & sice_main_final => DSIce_main_Final,                 &
       & sice_setup => DSIce_main_setup,                      &
       & sice_shutdown => DSIce_main_shutdown,                &       
       & sice_advance_timestep => DSIce_main_advance_timestep
  
  use DOGCM_Admin_TInteg_mod, only: &
       & TimeSecB, TimeSecN, TimeSecA, &
       & EndTimeSec => EndTime,        &
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
       & DEFAULT_DCCM_CONFNAME,     &       
       & NUM_DCCM_COMP,             &
       & NUM_OCN_GMAPTAG,           &
       & COMPNAME_OCN, GN25,        &
       & ATM_GRID_2D, OCN_GRID_2D,  &
       & GMAPTAG_ATM2D_OCN2D
       
  use mod_common_compdef, only:     &
       & ComponentDef_Init, ComponentDef_Final, &
       & ComponentDef_Share,        &
       & common_read_config,                    &
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

    use jcup_interface, only: &
         & jcup_init_time
    
    use field_def, only: &
         & cal_mn

    !* Dennou-OGCM

    use OptionParser_mod, only: &
         & OptionParser_Init,   &
         & OptionParser_Final,  &
         & OptionParser_GetInfo
    
    use DOGCM_Admin_TInteg_mod, only: &
         & InitialDate => InitDate, &
         & RestartTime, DelTime

    use DOGCM_Exp_driver_mod, only: &
         & DOGCM_Exp_driver_SetInitCond
    
    ! 宣言文; Declareration statements
    !    

    character(STRING) :: configNmlFile

    ! 実行文; Executable statement
    !
    
    ! Initialze varibale to manage component  ********************
    !
    
    call common_read_config( DEFAULT_DCCM_CONFNAME )

    call ComponentDef_Init( my_comp,         & ! (inout)
         & COMPNAME_OCN, "Dennou-OGCM" )       ! (in)
    

    !  Intialize a model associated with  ocean component ***
    !
    call MessageNotify( 'M', module_name, &
         & 'Initialize %a .. (MPI_MY_COMM=%d, MY_RANK=%d)', &
         & ca=(/ trim(my_comp%MODELNAME) /), i=(/ my_comp%PRC_comm, my_comp%PRC_rank/))

    
    call ogcm_main_init( my_comp%PRC_comm, isCoupledRun=.true. )
    call sice_main_init( isCoupledRun=.true. )

    call OptionParser_Init()
    call OptionParser_GetInfo( configNmlFile )
    call OptionParser_Final()
        
    call ogcm_setup( configNmlFile )
    call sice_setup( configNmlFile )

    !--
    
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

    my_comp%PRC_NY = 1  !size(a_jmax)
    my_comp%GNY    = JMO    
    allocate( lny_field(my_comp%PRC_NY) )
    lny_field(:)   = JMO !a_jmax(:)
    my_comp%LNY    = JMO !a_jmax(my_comp%PRC_rank)
    
!!$    write(*,*) "ATM: rank=", my_rank, "DIV_X,Y=", DIV_X, DIV_Y, "mysize=", my_size
    
    call ComponentDef_Share()
    

    ! Initialize grid for Jcup
    call MessageNotify( 'M', module_name, "Initialize grid for Jcup..")
    call init_jcup_grid()

    ! Initialize variables for JCup
    call MessageNotify( 'M', module_name, "Initialize varibles for Jcup..")
    call init_jcup_var()

    !- Initialize some data  to interpolate and grid mapping between model components in JCup
    call MessageNotify( 'M', module_name, "Initialize interpolation data for Jcup ..")
    call init_jcup_interpolate()

    
    !
    !
    
    call jcup_init_time(my_comp%InitTimeInfo)

    call MessageNotify( 'M', module_name, "Prepare data output..")
    call output_prepare()

    my_comp%tstep = 0
    my_comp%loop_end_flag = .false.
    call DOGCM_Exp_driver_SetInitCond()
    call sice_advance_timestep(my_comp%tstep, my_comp%loop_end_flag)    
    call ogcm_advance_timestep(my_comp%tstep, my_comp%loop_end_flag)    
    call MessageNotify( 'M', module_name, "Put data ..")
    call set_and_put_data(TimeSecN)
    my_comp%tstep = 1

    call MessageNotify( 'M', module_name, "ocn_init has been finished.")
    
  end subroutine ocn_init

  subroutine ocn_run(loop_flag)

    ! モジュール引用; Use statements
    !    
    use jcup_interface, only: &
         & jcup_set_time, jcup_inc_time

    !* DCPAM

    use DOGCM_Admin_Variable_mod
    
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
       
       !-----------------------------------------------------
       
       if (my_comp%PRC_rank==0 .and. mod(my_comp%tstep, 100) == 0) then
          call MessageNotify( 'M', module_name,                                &
               & "TimeSecN=%f, EndTimeSec=%f",  d=(/  TimeSecN, EndTimeSec /)  &
               & )
       end if

          
       ! Advance sea-ice component
       call sice_advance_timestep( my_comp%tstep, & ! (in)
            & my_comp%loop_end_flag               & ! (out)
            & )

       !- Advance ocean component
       call ogcm_advance_timestep( my_comp%tstep, & ! (in)
            & my_comp%loop_end_flag               & ! (out)
            & )


       !-----------------------------------------------------
       
!!$       write(*,*) "-> COUPLER Put: ocn my_rank=", my_comp%PRC_rank
       call set_and_put_data( TimeSecN ) ! (in)
       call jcup_inc_time( my_comp%name, my_comp%InitTimeInfo ) ! (in)
       my_comp%tstep = my_comp%tstep + 1       

       if ( EndTimeSec < TimeSecN ) then
          my_comp%loop_end_flag = .true.
       end if
       
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

    call jcup_coupling_end(my_comp%InitTimeInfo, .false.)

    !------------------------
    call MessageNotify( 'M', module_name, "Shutdown ocean and sea-ice model..")
    call sice_shutdown()
    call ogcm_shutdown()
    !-------------------------
    
    call MessageNotify( 'M', module_name, "ocn_fin has been finished. (rank=%d)", &
         & i=(/ my_comp%PRC_rank  /) )
    
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
         & OCN_NUM_GRIDTYPE,     &
         & OCN_GRID_2D

    use component_field, only: &
         & init_field

    !* Dennou-OGCM


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

    call init_field_def(OCN_NUM_GRIDTYPE)

    !-- Define atm_grid_2d
    call set_field_def( my_comp%name, OCN_GRID_2D, & ! (in)
         & my_comp%GNX, my_comp%GNY, 1, 2,         & ! (in) gnx, gny, gnz, halo
         & my_comp%PRC_NX, my_comp%PRC_NY          & ! (in) 
!         & , lnx_field, lny_field                    & ! (in)
         & )

    call get_local_field( component_name=my_comp%name, grid_name=OCN_GRID_2D, & ! (in) 
         & local_is=lis, local_ie=lie, local_js=ljs, local_je=lje             & ! (out)
         & )

    write(*,*) "OCN: rank=", my_comp%PRC_rank, "li=", lis, lie, " ,lj=", ljs, lje

    call init_field(field, OCN_GRID_2D, lis, lie, ljs, lje, 1, 1) ! halo=1
    
    ! Skip to define atm_grid_3d ------------
!!$    call set_field_def(ATM, ATM_GRID_3D, GNXA, GNYA, GNZA, 2, DIV_X, DIV_Y)
!!$    call init_field(field3d, ATM_GRID_3D, lis, lie, ljs, lje, 1, GNZA) ! halo=1    


!!$    call gen_grid_index()
    call cal_grid_index(my_comp%name, OCN_GRID_2D, field%grid_index)
!    write(*,*) "ATM: grid_index=", field%grid_index

!    call cal_grid_index(ATM, ATM_GRID_3D, field3d%grid_index)

    call jcup_def_grid(field%grid_index, my_comp%name, OCN_GRID_2D) !, GN25)
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
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_WindStressX,                                         & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1 )            ! (in)

    call jcup_def_varg( field%varg(a2o_WindStressY_id)%varg_ptr, my_comp%name, 'a2o_WindStressY', OCN_GRID_2D, 1, & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_WindStressY,                                         & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=1 )            ! (in)

    !--
    call jcup_def_varg( field%varg(a2o_LDwRFlx_id)%varg_ptr, my_comp%name, 'a2o_LDwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_LDwRFlx,                                      & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_SDwRFlx_id)%varg_ptr, my_comp%name, 'a2o_SDwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_SDwRFlx,                                      & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_LUwRFlx_id)%varg_ptr, my_comp%name, 'a2o_LUwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_LUwRFlx,                                      & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_SUwRFlx_id)%varg_ptr, my_comp%name, 'a2o_SUwRFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_SUwRFlx,                                      & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)
    
    call jcup_def_varg( field%varg(a2o_LatHFlx_id)%varg_ptr, my_comp%name, 'a2o_LatHFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_LatHFlx,                                      & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_SenHFlx_id)%varg_ptr, my_comp%name, 'a2o_SenHFlx', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_SenHFlx,                                      & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )     ! (in)

    call jcup_def_varg( field%varg(a2o_DSfcHFlxDTs_id)%varg_ptr, my_comp%name, 'a2o_DSfcHFlxDTs', OCN_GRID_2D, 1,  & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_DSfcHFlxDTs,                                          & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=2 )             ! (in)
    !--
    call jcup_def_varg( field%varg(a2o_RainFall_id)%varg_ptr, my_comp%name, 'a2o_RainFall', OCN_GRID_2D, 1,   & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_RainFall,                                        & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=3 )        ! (in)
    
    call jcup_def_varg( field%varg(a2o_SnowFall_id)%varg_ptr, my_comp%name, 'a2o_SnowFall', OCN_GRID_2D, 1,   & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_SnowFall,                                        & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=3 )        ! (in)


    call jcup_def_varg( field%varg(a2o_SfcAirTemp_id)%varg_ptr, my_comp%name, 'a2o_SfcAirTemp', OCN_GRID_2D, 1, & ! (in)
         & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=a2d_SfcAirTemp,                                        & ! (in)
         & RECV_MODE='AVR', INTERVAL=AO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=1, EXCHANGE_TAG=3 )          ! (in)
    
    !- Finish defining variables for jup ---------------------
    
    call jcup_end_var_def()
    
  end subroutine init_jcup_var

  !---------------------------------------------------------------------------------------------------
  
  subroutine init_jcup_interpolate()

    ! モジュール引用; Use statement
    !            
    use jcup_interface, only: &
         & jcup_set_mapping_table, jcup_get_component_name

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
    call init_interpolation(NUM_DCCM_COMP, NUM_OCN_GMAPTAG, my_comp%id)

    
    ! ATM -> OCN grid mapping    **************************
    call jcup_set_mapping_table( my_comp%name,                       & ! (in)
         & atm_comp%name, ATM_GRID_2D, my_comp%name, OCN_GRID_2D,    & ! (in) ATM_GRID_2D -> OCN_GRID_2D
         & GMAPTAG_ATM2D_OCN2D )                                       ! (in)

    call set_operation_index(my_comp%name, atm_comp%name, GMAPTAG_ATM2D_OCN2D)         ! (in)
    
    ! OCN -> ATM grid mapping   *****************************    
    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & my_comp%name, OCN_GRID_2D, atm_comp%name, ATM_GRID_2D,     & ! (in) OCN_GRID_2D -> ATM_GRID_2D
         & GMAPTAG_ATM2D_OCN2D )                                       ! (in)

    call set_A_to_O_coef(GMAPTAG_ATM2D_OCN2D)                    ! (in)
    call set_O_to_A_coef(GMAPTAG_ATM2D_OCN2D)                    ! (in)

  end subroutine init_jcup_interpolate

  !----------------------------------------------------------------------------------
  
  subroutine set_and_put_data( CurrentTimeSec )

    ! モジュール引用; Use statement
    !                

    use mpi
    
    !* Dennou-OGCM / Sea-ice


    use DSIce_Admin_Constants_mod, only: &
         & AlbedoOcean, IceMaskMin

    use DSIce_Admin_Variable_mod, only:  &
         & xya_SIceCon, xya_SIceSfcTemp, &
         & xya_IceThick, xya_SnowThick
    
    use DSIce_Boundary_vars_mod, only: &
         & xy_SfcAlbedoAI

    use DOGCM_Boundary_vars_mod, only: &
         & xy_SeaSfcTemp
    
    use DOGCM_Admin_Variable_mod, only: &
         & xyzaa_TRC,  &
         & TRCID_PTEMP
    
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

    !$omp parallel do private(i)
    do j = JSO, JEO
       do i = ISO, IEO
          if ( xya_SIceCon(i,j,SICE_TLN) > IceMaskMin ) then
             xy_SfcAlbedo4Atm(i,j) = xy_SfcAlbedoAI(i,j)
             xy_SfcTemp4Atm(i,j) = xya_SIceSfcTemp(i,j,SICE_TLN)
             xy_SfcSnow4Atm(i,j) = xya_SnowThick(i,j,SICE_TLN)
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

    
    !* Dennou-OGCM

    use DOGCM_Admin_Constants_mod, only: &
         & LatentHeat
    
    use DOGCM_Boundary_vars_mod, only: &
         & xy_WindStressXAO => xy_WindStressU, &
         & xy_WindStressYAO => xy_WindStressV, &
         & xy_SfcHFlx_ns, xy_SfcHFlx_sr,       &
         & xy_DSfcHFlxDTs,                     &
         & xy_FreshWtFlx, xy_FreshWtFlxS,      &
         & xy_SfcAirTemp

    use DSIce_Admin_Constants_mod, only: &
         & DensFreshWater, LFreeze
    
    use DSIce_Boundary_vars_mod, only: &
         & xy_WindStressXAI => xy_WindStressUAI,            &
         & xy_WindStressYAI => xy_WindStressVAI,            &
         & xy_SDwRFlx, xy_LDwRFlx, xy_LatHFlx, xy_SenHFlx,  &
         & xy_DSfcHFlxAIDTs,                                &
         & xy_RainFall, xy_SnowFall, xy_Evap
    
    ! 宣言文; Declareration statements
    !            
    real(DP), intent(in) :: CurrentTimeSec

    ! 局所変数
    ! Local variables
    !    

    real(DP) :: xy_LUwRFlx(IAO,JAO)
    real(DP) :: xy_SUwRFlx(IAO,JAO)
    
    integer :: i
    integer :: j

    ! 実行文; Executable statement
    !
    
    ! Get oceanic surface temerature send by OGCM.
    !

    call ocn_get_write( a2o_WindStressX_id, "a2o_WindStressX", & ! (in)
         & xy_WindStressXAO )                                    ! (out)

    call ocn_get_write( a2o_WindStressY_id, "a2o_WindStressY", & ! (in)
         & xy_WindStressYAO )                                    ! (out)

    call ocn_get_write( a2o_LDwRFlx_id, "a2o_LDwRFlx",         & ! (in)
         & xy_LDwRFlx )                                          ! (in)

    call ocn_get_write( a2o_SDwRFlx_id, "a2o_SDwRFlx",         & ! (in)
         & xy_SDwRFlx )                                          ! (in)

    call ocn_get_write( a2o_LUwRFlx_id, "a2o_LUwRFlx",         & ! (in)
         & xy_LUwRFlx )                                          ! (in)

    call ocn_get_write( a2o_SUwRFlx_id, "a2o_SUwRFlx",         & ! (in)
         & xy_SUwRFlx )                                          ! (in)

    call ocn_get_write( a2o_LatHFlx_id, "a2o_LatHFlx",         & ! (in)
         & xy_LatHFlx )                                          ! (in)

    call ocn_get_write( a2o_SenHFlx_id, "a2o_SenHFlx",         & ! (in)
         & xy_SenHFlx )                                          ! (in)         

    call ocn_get_write( a2o_DSfcHFlxDTs_id, "a2o_DSfcHFlxDTs", & ! (in)
         & xy_DSfcHFlxDTs )                                      ! (in)         

    call ocn_get_write( a2o_RainFall_id, "a2o_RainFall",       & ! (in)
         & xy_RainFall )                                         ! (in)         

    call ocn_get_write( a2o_SnowFall_id, "a2o_SnowFall",       & ! (in)
         & xy_SnowFall )                                         ! (in)         

    call ocn_get_write( a2o_SfcAirTemp_id, "a2o_SfcAirTemp",   & ! (in)
         & xy_SfcAirTemp )                                       ! (in)         

    !------------------------------------------------------------------
    
    !$omp parallel do private(i)
    do j = JSO, JEO
       do i = ISO, IEO

          ! For ocean model
          xy_SfcHFlx_ns(i,j) =    xy_LUwRFlx(i,j) - xy_LDwRFlx(i,j)  &
               &               +  xy_LatHFlx(i,j) + xy_SenHFlx(i,j)  &
               &               +  LFreeze * xy_SnowFall(i,j)           ! [J/kg * kg.m-2.s-1]

          xy_SfcHFlx_sr(i,j) = xy_SUwRFlx(i,j) - xy_SDwRFlx(i,j)

          xy_Evap(i,j) = xy_LatHFlx(i,j)/(LatentHeat*DensFreshWater)
          
          xy_FreshWtFlxS(i,j) =   (xy_RainFall(i,j) + xy_SnowFall(i,j))/DensFreshWater &
               &                - xy_Evap(i,j)
          xy_FreshWtFlx(i,j) = xy_FreshWtFlxS(i,j)

          ! For sea-ice model
          xy_WindStressXAI(i,j) = xy_WindStressXAO(i,j)
          xy_WindStressYAI(i,j) = xy_WindStressYAO(i,j)          
          xy_DSfcHFlxAIDTs(i,j) = xy_DSfcHFlxDTs(i,j)

       end do
    end do

    
!!$    !
!!$    !
!!$    if( mod(CurrentTimeSec, dble(AO_COUPLING_CYCLE_SEC)) == 0d0 ) then
!!$
!!$    if(my_rank>12) then
!!$       write(*,*) "atm: rank=", my_rank, "SurfTemp=", xy_SurfTemp(0,1:2), "lat=", y_Lat(1:2)/acos(-1d0)*180d0
!!$    end if
!!$       
!!$    end if

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

    use DOGCM_IO_History_mod, only: &
         & DOGCM_IO_History_RegistVar
    
    use mod_common_params
    

    character(TOKEN) :: dims_XYT(3)

    ! 実行文; Executable statement
    !

    dims_XYT = (/ 'lon ', 'lat ', 'time' /)

    call DOGCM_IO_History_RegistVar( 'a2o_WindStressX', 'IJT', 'i axis component of wind stress(ATM->OCN)', 'kg/(m.s2)' )
    call DOGCM_IO_History_RegistVar( 'a2o_WindStressY', 'IJT', 'j axis component of wind stress(ATM->OCN)', 'kg/(m.s2)' )

    call DOGCM_IO_History_RegistVar( 'a2o_LDwRFlx', 'IJT', 'downward long wave radiation', 'W.m-2' )
    call DOGCM_IO_History_RegistVar( 'a2o_SDwRFlx', 'IJT', 'downward short wave radiation', 'W.m-2' )
    call DOGCM_IO_History_RegistVar( 'a2o_LUwRFlx', 'IJT', 'upward long wave radiation', 'W.m-2' )
    call DOGCM_IO_History_RegistVar( 'a2o_SUwRFlx', 'IJT', 'upward short wave radiation', 'W.m-2' )

    call DOGCM_IO_History_RegistVar( 'a2o_LatHFlx', 'IJT', 'Latent heat flux', 'W.m-2')    
    call DOGCM_IO_History_RegistVar( 'a2o_SenHFlx', 'IJT', 'Sensible heat flux', 'W.m-2')
    call DOGCM_IO_History_RegistVar( 'a2o_DSfcHFlxDTs', 'IJT', 'surface temperature dependence of heat flux', 'W.m-2.K-1' )

    call DOGCM_IO_History_RegistVar( 'a2o_RainFall', 'IJT', 'Rain fall', 'm/s')
    call DOGCM_IO_History_RegistVar( 'a2o_SnowFall', 'IJT', 'Snow fall', 'm/s')
    call DOGCM_IO_History_RegistVar( 'a2o_Evap', 'IJT', 'Evarporation', 'm/s')

    call DOGCM_IO_History_RegistVar( 'a2o_SfcAirTemp', 'IJT', 'Surface air temperature', 'K')
    
  end subroutine output_prepare

  subroutine output_var(CurrentTime, varname, data2d)

    ! モジュール引用; Use statement
    !                
    use DOGCM_IO_History_mod, only: &
         & DOGCM_IO_History_HistPut
        
    real(DP), intent(in) :: CurrentTime
    character(*), intent(in) :: varName
    real(DP), intent(in) :: data2d(IAO,JAO)

    ! 実行文; Executable statement
    !
    
    call DOGCM_IO_History_HistPut(varname, data2d(ISO:IEO,JSO:JEO))
    
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

