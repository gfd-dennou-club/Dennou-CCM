!-------------------------------------------------------------
! Copyright (c) 2015-2017 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief A module for surface component
!! 
!! @author Kawai Yuta
!!
!!
module dccm_sfc_mod
   
  ! モジュール引用; Use statements
  !

  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* DCPCM

  use DSFCM_mod, only: &
       & sfcm_main_init => DSFCM_main_Init,                   &
       & sfcm_main_final => DSFCM_main_Final,                 &
       & sfcm_setup => DSFCM_main_setup,                      &
       & sfcm_shutdown => DSFCM_main_shutdown,                &       
       & sfcm_advance_timestep => DSFCM_main_advance_timestep

  use DSFCM_Admin_TInteg_mod, only: &
       & TimeSecB, TimeSecN, TimeSecA, &
       & EndTimeSec => EndTime,        &
       & SFC_TLN => TIMELV_ID_N

  use DSFCM_Admin_Grid_mod, only: &
       & ISS => IS, IES => IE, IMS => IM,   &
       & JSS => JS, JES => JE, JMS => JM,   &
       & IAS => IA, JAS => JA
  
  use component_field, only: &
       & component_field_type

  use dccm_common_params_mod, only: &
       & DEFAULT_DCCM_CONFNAME,      &       
       & NUM_DCCM_COMP,              &
       & NUM_OCN_GMAPTAG,            &
       & NUM_SFC_GMAPTAG,            &
       & COMPNAME_ATM, COMPNAME_OCN, GN25,        &
       & COMPNAME_SFC,                            &
       & ATM_GRID_2D, OCN_GRID_2D, SFC_GRID_2D,   &
       & GMAPTAG_ATM2D_SFC2D,          &
       & GMAPTAG_ATM2D_SFC2D_CONSERVE, &
       & GMAPTAG_SFC2D_OCN2D,          &
       & GMAPTAG_SFC2D_OCN2D_CONSERVE
  
  use mod_common_compdef, only:     &
       & ComponentDef_Init, ComponentDef_Final, &
       & ComponentDef_Share,        &
       & common_read_config,                    &
       & my_comp  => CompDef_sfc,  &
       & ocn_comp  => CompDef_ocn,  &
       & atm_comp  => CompDef_atm,  &
       & sice_comp => CompDef_sice, &
       & GMAPFILENAME_OS,           &
       & GMAPFILENAME_SO,           &
       & GMAPFILENAME_OS_conserve,           &
       & GMAPFILENAME_SO_conserve,           &
       & AS_COUPLING_CYCLE_SEC,     &
       & SO_COUPLING_CYCLE_SEC,     &
       & SFC_GNX, SFC_GNY
       

  
  !* Dennou-OGCM

  implicit none
  private
  
  ! 公開手続き
  ! Public procedure
  !    
  public :: sfc_init
  public :: sfc_run
  public :: sfc_fin

  integer :: tstep

  type(component_field_type) :: field_as
  type(component_field_type) :: field_so
  type(component_field_type) :: field_si

  integer, allocatable :: lnx_field(:)
  integer, allocatable :: lny_field(:)

  real(DP) :: a_Sig1Info(2)
  
  character(*), parameter :: module_name = 'mod_sfc'

  
contains

  subroutine sfc_init()

    ! モジュール引用; Use statements
    !

    use jcup_interface, only: &
         & jcup_init_time, jcup_set_time


    !* Dennou-OGCM

    use OptionParser_mod, only: &
         & OptionParser_Init,   &
         & OptionParser_Final,  &
         & OptionParser_GetInfo
    
    use DSFCM_Admin_TInteg_mod, only: &
         & InitialDate => InitDate, &
         & RestartTime, DelTime

    use mpi
    
    ! 宣言文; Declareration statements
    !    

    character(STRING) :: configNmlFile
    integer :: ierror
    
    ! 実行文; Executable statement
    !
    
    ! Initialze varibale to manage component  ********************
    !

    call common_read_config( DEFAULT_DCCM_CONFNAME )

    call ComponentDef_Init( my_comp,         & ! (inout)
         & COMPNAME_SFC, "Dennou-SFCM" )       ! (in)

    !  Intialize a model associated with  surface component ***
    !
    call MessageNotify( 'M', module_name, &
         & 'Initialize %a .. (MPI_MY_COMM=%d, MY_RANK=%d)', &
         & ca=(/ trim(my_comp%MODELNAME) /), i=(/ my_comp%PRC_comm, my_comp%PRC_rank/))

    call sfcm_main_init()

    call OptionParser_Init()
    call OptionParser_GetInfo( configNmlFile )
    call OptionParser_Final()
    
    call sfcm_setup( configNmlFile, SFC_GNX, SFC_GNY )
    
    !-
    
    my_comp%InitTimeInfo(:) = (/ InitialDate%year, InitialDate%month, InitialDate%day,        &
                        &        InitialDate%hour, InitialDate%min, int(InitialDate%sec)  /)
    my_comp%RestartTimeSec  = RestartTime
    my_comp%DelTime         = DelTime
    my_comp%tstep           = 0    
    call MessageNotify( 'M', module_name, &
         & "Initial Time %4d - %2d - %2d:%2d%2d", i=my_comp%InitTimeInfo )
    call MessageNotify( 'M', module_name, &
         & "RestartTiem = %f [sec]", d=(/ my_comp%RestartTimeSec /) )
    
    my_comp%PRC_NX = 1
    my_comp%GNX    = SFC_GNX
    allocate( lnx_field(my_comp%PRC_NX) )
    lnx_field(:)   = SFC_GNX
    my_comp%LNX    = SFC_GNX

    my_comp%PRC_NY = 1  !size(a_jmax)
    my_comp%GNY    = SFC_GNY    
    allocate( lny_field(my_comp%PRC_NY) )
    lny_field(:)   = SFC_GNY !a_jmax(:)
    my_comp%LNY    = SFC_GNY !a_jmax(my_comp%PRC_rank)

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
    call MessageNotify( 'M', module_name, "Prepare data output..")
!!$    call output_prepare()


    call jcup_init_time(my_comp%InitTimeInfo)    
    
    !--------------------------------------------------

    my_comp%tstep = 0
    my_comp%loop_end_flag = .false.
    
!!$    !
!!$    call sfcm_advance_timestep( &
!!$         & my_comp%tstep,          & ! (in)
!!$         & my_comp%loop_end_flag,  & ! (out)
!!$         & skip_flag = .false. )     ! (in)

    !
!!$    call MessageNotify( 'M', module_name, "Put data ..")
!!$    call set_and_put_data(TimeSecN)

    my_comp%tstep = 1

    call MessageNotify( 'M', module_name, "sfc_init has been finished.")

  end subroutine sfc_init

  
  subroutine sfc_run(loop_flag)

    ! モジュール引用; Use statements
    !    

    !* DCPAM

    use jcup_interface, only: &
         & jcup_set_time, jcup_inc_time, jcup_inc_calendar

    ! 宣言文; Declareration statements
    !        
    logical, intent(inout) :: loop_flag

    ! 局所変数
    ! Local variables
    !
    integer, parameter :: MONITOR_STEPINT =400

    ! 実行文; Executable statement
    !

    my_comp%loop_end_flag = .false.


    do while(.not. my_comp%loop_end_flag)

       my_comp%DelTime = TimeSecA - TimeSecN

       write(*,*) "sfc: set time, TimeSecN=", TimeSecN, "DelTime=", my_comp%DelTime
       
       call jcup_set_time( my_comp%name,                  & ! (in)
            & my_comp%InitTimeInfo, int(my_comp%DelTime) )  ! (in)

       write(*,*) "sfc: get and write date, TimeSecN=", TimeSecN
       call get_and_write_data( TimeSecN )

       !-----------------------------------------------------
       
       if (my_comp%PRC_rank==0 .and. mod(my_comp%tstep, MONITOR_STEPINT) == 0) then
          call MessageNotify( 'M', module_name,            &
               & "TimeSecN=%f (EndTimeSec=%f, tstep=%d)",              &
               & d=(/  TimeSecN, EndTimeSec /), i=(/ my_comp%tstep /)  &
               & )
       end if

       call sfcm_advance_timestep( &
         & my_comp%tstep,          & ! (in)
         & my_comp%loop_end_flag,  & ! (out)
         & skip_flag = .true. )     ! (in)
       
       call set_and_put_data( TimeSecN ) ! (in)
       call jcup_inc_calendar( my_comp%InitTimeInfo, int(my_comp%DelTime) ) ! (in)
       my_comp%tstep = my_comp%tstep + 1       

       if ( EndTimeSec < TimeSecN ) then
          my_comp%loop_end_flag = .true.
       end if

       write(*,*) "<-- sfc: TimeLoopEnd: TimeSecN=", TimeSecN
    end do
    loop_flag = .false.

    
  end subroutine sfc_run

  subroutine sfc_fin()

    ! モジュール引用; Use statement
    !    
    use jcup_interface, only: &
         & jcup_coupling_end

    ! 実行文; Executable statement
    !
    
    call jcup_coupling_end(my_comp%InitTimeInfo, .false.)

    !------------------------
    call MessageNotify( 'M', module_name, "Shutdown surface model..")
    call sfcm_shutdown()

    !-------------------------
    
    call MessageNotify( 'M', module_name, "sfc_fin has been finished. (rank=%d)", &
         & i=(/ my_comp%PRC_rank  /) )
    
  end subroutine sfc_fin

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
    
    use dccm_common_params_mod, only: &
         & SFC_NUM_GRIDTYPE,          &
         & SFC_GRID_2D

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

    call init_field_def(SFC_NUM_GRIDTYPE)

    !-- Define sfc_grid_2d
    call set_field_def( my_comp%name, SFC_GRID_2D, & ! (in)
         & my_comp%GNX, my_comp%GNY, 1, 2,         & ! (in) gnx, gny, gnz, halo
         & my_comp%PRC_NX, my_comp%PRC_NY          & ! (in) 
!         & , lnx_field, lny_field                    & ! (in)
         & )
  
    call get_local_field( component_name=my_comp%name, grid_name=SFC_GRID_2D, & ! (in) 
         & local_is=lis, local_ie=lie, local_js=ljs, local_je=lje             & ! (out)
         & )

    write(*,*) "SFC: rank=", my_comp%PRC_rank, "li=", lis, lie, " ,lj=", ljs, lje

    call init_field(field_as, SFC_GRID_2D, lis, lie, ljs, lje, 1, 1) ! halo=1
    call init_field(field_so, SFC_GRID_2D, lis, lie, ljs, lje, 1, 1) ! halo=1
    call init_field(field_si, SFC_GRID_2D, lis, lie, ljs, lje, 1, 1) ! halo=1
    
    ! Skip to define atm_grid_3d ------------
!!$    call set_field_def(ATM, ATM_GRID_3D, GNXA, GNYA, GNZA, 2, DIV_X, DIV_Y)
!!$    call init_field(field3d, ATM_GRID_3D, lis, lie, ljs, lje, 1, GNZA) ! halo=1    


!!$    call gen_grid_index()
    call cal_grid_index(my_comp%name, SFC_GRID_2D, field_as%grid_index)
    call cal_grid_index(my_comp%name, SFC_GRID_2D, field_so%grid_index)
    call cal_grid_index(my_comp%name, SFC_GRID_2D, field_si%grid_index)

!!$    write(*,*) "SFC: grid_index=", field%grid_index

!    call cal_grid_index(ATM, ATM_GRID_3D, field3d%grid_index)

    call jcup_def_grid(field_as%grid_index, my_comp%name, SFC_GRID_2D) !, GN25)
!    call jcup_def_grid(field3d%grid_index, ATM, ATM_GRID_3D)

    !---- Finish defining grid for jcup ----------------------------------------
    call jcup_end_grid_def()
    
  end subroutine init_jcup_grid
  
  subroutine init_jcup_var()

    ! モジュール引用; Use statement
    !        
    use jcup_interface, only: &
         & jcup_def_varp, jcup_def_varg, jcup_end_var_def

    use component_field, only: &
         & init_field_data

    use mod_common_compdef, only: ComponentDef
    
    use dccm_common_params_mod
    
    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !

    call init_field_data( field_as, num_of_25d=GN25, num_of_varp=NUM_VAR2D_s2a, num_of_varg=NUM_VAR2D_a2s)
    call init_field_data( field_so, num_of_25d=GN25, num_of_varp=NUM_VAR2D_s2o, num_of_varg=NUM_VAR2D_o2s)
    call init_field_data( field_si, num_of_25d=GN25, num_of_varp=0, num_of_varg=NUM_VAR2D_i2s)

    !- Variable put by my own component
    
    call jcup_def_varp( field_as%varp(s2a_WindStressX_id)%varp_ptr, my_comp%name, s2a_WindStressX, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_WindStressY_id)%varp_ptr, my_comp%name, s2a_WindStressY, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_LUwRFlx_id)%varp_ptr, my_comp%name, s2a_LUwRFlx, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_SUwRFlx_id)%varp_ptr, my_comp%name, s2a_SUwRFlx, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_SenHFlx_id)%varp_ptr, my_comp%name, s2a_SenHFlx, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_QVapMFlx_id)%varp_ptr, my_comp%name, s2a_QVapMFlx, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_SfcAlbedo_id)%varp_ptr, my_comp%name, s2a_SfcAlbedo, SFC_GRID_2D )
    call jcup_def_varp( field_as%varp(s2a_SfcRadTemp_id)%varp_ptr, my_comp%name, s2a_SfcRadTemp, SFC_GRID_2D )

    call jcup_def_varp( field_so%varp(s2o_SfcHFlx_ns_id)%varp_ptr, my_comp%name, s2o_SfcHFlx_ns, SFC_GRID_2D )
    call jcup_def_varp( field_so%varp(s2o_SfcHFlx_sr_id)%varp_ptr, my_comp%name, s2o_SfcHFlx_sr, SFC_GRID_2D )    
    call jcup_def_varp( field_so%varp(s2o_SnowFall_id)%varp_ptr, my_comp%name, s2o_SnowFall, SFC_GRID_2D )    
    call jcup_def_varp( field_so%varp(s2o_RainFall_id)%varp_ptr, my_comp%name, s2o_RainFall, SFC_GRID_2D )    
    call jcup_def_varp( field_so%varp(s2o_Evap_id)%varp_ptr, my_comp%name, s2o_Evap, SFC_GRID_2D )    
    call jcup_def_varp( field_so%varp(s2o_WindStressX_id)%varp_ptr, my_comp%name, s2o_WindStressX, SFC_GRID_2D )
    call jcup_def_varp( field_so%varp(s2o_WindStressY_id)%varp_ptr, my_comp%name, s2o_WindStressY, SFC_GRID_2D )

    !- Variable getten from  other components

    call regist_jcup_var_as( a2s_WindU_id, a2s_WindU, 1, GMAPTAG_ATM2D_OCN2D )
    call regist_jcup_var_as( a2s_WindV_id, a2s_WindV, 1, GMAPTAG_ATM2D_OCN2D )
    call regist_jcup_var_as( a2s_SfcAirTemp_id, a2s_SfcAirTemp, 2, GMAPTAG_ATM2D_OCN2D )
    call regist_jcup_var_as( a2s_QVap1_id, a2s_QVap1, 2, GMAPTAG_ATM2D_OCN2D )
    call regist_jcup_var_as( a2s_SfcPress_id, a2s_SfcPress, 2, GMAPTAG_ATM2D_OCN2D )
    call regist_jcup_var_as( a2s_Press1_id, a2s_Press1, 2, GMAPTAG_ATM2D_OCN2D )
    call regist_jcup_var_as( a2s_LDwRFlx_id, a2s_LDwRFlx, 2, GMAPTAG_ATM2D_SFC2D_CONSERVE )
    call regist_jcup_var_as( a2s_SDwRFlx_id, a2s_SDwRFlx, 3, GMAPTAG_ATM2D_SFC2D_CONSERVE )
    call regist_jcup_var_as( a2s_RainFall_id, a2s_RainFall, 3, GMAPTAG_ATM2D_SFC2D_CONSERVE )
    call regist_jcup_var_as( a2s_SnowFall_id, a2s_SnowFall, 3, GMAPTAG_ATM2D_SFC2D_CONSERVE )
    
    call regist_jcup_var_os( o2s_SfcTemp_id, o2s_SfcTemp, 1, GMAPTAG_SFC2D_OCN2D)    
    call regist_jcup_var_os( o2s_SfcAlbedo_id, o2s_SfcAlbedo, 1, GMAPTAG_SFC2D_OCN2D_CONSERVE)

    call regist_jcup_var_is( i2s_SfcTemp_id, i2s_SfcTemp, 2, GMAPTAG_SFC2D_OCN2D)
    call regist_jcup_var_is( i2s_SfcAlbedo_id, i2s_SfcAlbedo, 2, GMAPTAG_SFC2D_OCN2D_CONSERVE)
    call regist_jcup_var_is( i2s_SIceCon_id, i2s_SIceCon, 2, GMAPTAG_SFC2D_OCN2D_CONSERVE)
    
    !- Finish defining variables for jup ---------------------
    
    call jcup_end_var_def()

  contains
    subroutine regist_jcup_var_as(varid, varname, exchange_tag, mapping_tag)
      integer, intent(in) :: varid
      character(*), intent(in) :: varname
      integer, intent(in) :: exchange_tag
      integer, intent(in) :: mapping_tag

      call jcup_def_varg( field_as%varg(varid)%varg_ptr, my_comp%name, varname, SFC_GRID_2D, 1,                      & ! (in)
           & SEND_MODEL_NAME=atm_comp%name, SEND_DATA_NAME=varname, RECV_MODE='SNP',                                 & ! (in)
           & INTERVAL=AS_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=mapping_tag, EXCHANGE_TAG=exchange_tag )   ! (in)
    end subroutine regist_jcup_var_as

    subroutine regist_jcup_var_os(varid, varname, exchange_tag, mapping_tag)
      integer, intent(in) :: varid
      character(*), intent(in) :: varname
      integer, intent(in) :: exchange_tag
      integer, intent(in) :: mapping_tag

      call jcup_def_varg( field_so%varg(varid)%varg_ptr, my_comp%name, varname, SFC_GRID_2D, 1,                      & ! (in)
           & SEND_MODEL_NAME=ocn_comp%name, SEND_DATA_NAME=varname, RECV_MODE='SNP',                                 & ! (in)
           & INTERVAL=SO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=mapping_tag, EXCHANGE_TAG=exchange_tag )   ! (in)
    end subroutine regist_jcup_var_os

    subroutine regist_jcup_var_is(varid, varname, exchange_tag, mapping_tag)
      integer, intent(in) :: varid
      character(*), intent(in) :: varname
      integer, intent(in) :: exchange_tag
      integer, intent(in) :: mapping_tag

      call jcup_def_varg( field_si%varg(varid)%varg_ptr, my_comp%name, varname, SFC_GRID_2D, 1,                      & ! (in)
           & SEND_MODEL_NAME=ocn_comp%name, SEND_DATA_NAME=varname, RECV_MODE='SNP',                                 & ! (in)
           & INTERVAL=SO_COUPLING_CYCLE_SEC, TIME_LAG=-1, MAPPING_TAG=mapping_tag, EXCHANGE_TAG=exchange_tag )   ! (in)
    end subroutine regist_jcup_var_is
    
  end subroutine init_jcup_var
  
  subroutine init_jcup_interpolate()
 
    ! モジュール引用; Use statement
    !            
    use jcup_interface, only: &
         & jcup_set_mapping_table, jcup_get_component_name, &
         & jcup_recv_array

    use jcup_mpi_lib, only: &
         & jml_finalize

    use field_def, only : cal_mn, &
         & set_grid_mapping, set_grid_mapping_3d

    use interpolation_data_latlon_mod, only: &
         & init_interpolation => interpolation_data_latlon_Init, &
         & set_operation_index, set_interpolate_coef

    use grid_mapping_util, only: &
         & set_mappingTable_interpCoef
    
    use grid_mapping_util_jones99, only: &    
         & set_mappingTable_interpCoef_j99 =>  set_mappingTable_interpCoef

    use dccm_common_params_mod
    
    ! 局所変数
    ! Local variables
    !

    integer, allocatable :: send_grid_os(:)
    integer, allocatable :: send_grid_so(:)
    integer, allocatable :: recv_grid_os(:)
    integer, allocatable :: recv_grid_so(:)
    real(DP), allocatable :: coefS_os_global(:)
    real(DP), allocatable :: coefS_so_global(:)
    real(DP), allocatable :: coefS_os_global_conserve(:)
    real(DP), allocatable :: coefS_so_global_conserve(:)
    
    ! 実行文; Executable statement
    !
    
    !
    call init_interpolation(NUM_DCCM_COMP, NUM_SFC_GMAPTAG, my_comp%id)

   
    ! ATM -> SFC grid mapping    **************************
    call jcup_set_mapping_table( my_comp%name,                       & ! (in)
         & atm_comp%name, ATM_GRID_2D, my_comp%name, SFC_GRID_2D,    & ! (in) ATM_GRID_2D -> SFC_GRID_2D
         & GMAPTAG_ATM2D_SFC2D )                                       ! (in)

    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & atm_comp%name, ATM_GRID_2D, my_comp%name, SFC_GRID_2D,     & ! (in) ATM_GRID_2D -> SFC_GRID_2D
         & GMAPTAG_ATM2D_SFC2D_CONSERVE )                               ! (in)

    call set_operation_index(my_comp%name, atm_comp%name, GMAPTAG_ATM2D_SFC2D)            ! (in)
    call set_operation_index(my_comp%name, atm_comp%name, GMAPTAG_ATM2D_SFC2D_CONSERVE)   ! (in)    
    
    !* SFC -> ATM grid mapping   *****************************    
    !
     
    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & my_comp%name, SFC_GRID_2D, atm_comp%name, ATM_GRID_2D,     & ! (in) SFC_GRID_2D -> ATM_GRID_2D
         & GMAPTAG_ATM2D_SFC2D )                                        ! (in)

    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & my_comp%name, SFC_GRID_2D, atm_comp%name, ATM_GRID_2D,     & ! (in) SFC_GRID_2D -> ATM_GRID_2D
         & GMAPTAG_ATM2D_SFC2D_CONSERVE )                               ! (in)


    call MessageNotify('M', module_name, "Set the coefficients for interpolation..")
    
    ! Set the coefficients for ATM -> SFC interpolation     
    call set_interpolate_coef( atm_comp%name, my_comp%name, &
         & atm_comp%name, GMAPTAG_ATM2D_SFC2D )

    call set_interpolate_coef( atm_comp%name, my_comp%name, &
         & atm_comp%name, GMAPTAG_ATM2D_SFC2D_CONSERVE )

    ! Set the coefficients for SFC -> ATM interpolation 
    call set_interpolate_coef( my_comp%name, atm_comp%name, &
         & atm_comp%name, GMAPTAG_ATM2D_SFC2D )

    call set_interpolate_coef( my_comp%name, atm_comp%name, &
         & atm_comp%name, GMAPTAG_ATM2D_SFC2D_CONSERVE )

    
    !* SFC -> OCN grid mapping    **************************
    !

    call MessageNotify( 'M', module_name, "Read the coefficients (SFC -> OCN) of bilinear interpolation.")    

    if(my_comp%PRC_rank==0) then
       call set_mappingTable_interpCoef( &
            & GMAPFILENAME_SO, my_comp%GNX, ocn_comp%GNX,             & ! (in)
            & send_grid_so, recv_grid_so, coefS_so_global             & ! (inout)
            & )
    end if
    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & my_comp%name, SFC_GRID_2D, ocn_comp%name, OCN_GRID_2D,     & ! (in) SFC_GRID_2D -> OCN_GRID_2D
         & GMAPTAG_SFC2D_OCN2D, send_grid_so, recv_grid_so )            ! (in)
    if(my_comp%PRC_rank==0) deallocate( send_grid_so, recv_grid_so )

    call MessageNotify( 'M', module_name, "Read the coefficients (SFC -> OCN) of conservative bilinear interpolation.")    

    if(my_comp%PRC_rank==0) then    
       call set_mappingTable_interpCoef_j99( &
            & GMAPFILENAME_SO_CONSERVE, my_comp%GNX, ocn_comp%GNX,    & ! (in)
            & send_grid_so, recv_grid_so, coefS_so_global_conserve    & ! (inout)
            & )
!       write(*,*) "A2O:", size(send_grid_ao), size(recv_grid_ao), size(coefS_ao_global_conserve)
    end if
    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & my_comp%name, SFC_GRID_2D, ocn_comp%name, OCN_GRID_2D,     & ! (in) SFC_GRID_2D -> OCN_GRID_2D
         & GMAPTAG_SFC2D_OCN2D_CONSERVE, send_grid_so, recv_grid_so )   ! (in)

    !* OCN -> SFC grid mapping   *****************************    
    !

    call MessageNotify( 'M', module_name, "Read the coefficients (OCN -> SFC) of bilinear interpolation.")    

    if(my_comp%PRC_rank==0) then
       call set_mappingTable_interpCoef( &
            & GMAPFILENAME_OS, ocn_comp%GNX, my_comp%GNX,          & ! (in)
            & send_grid_os, recv_grid_os, coefS_os_global          & ! (inout)
            & )
    end if
    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & ocn_comp%name, OCN_GRID_2D, my_comp%name, SFC_GRID_2D,     & ! (in) OCN_GRID_2D -> SFC_GRID_2D
         & GMAPTAG_SFC2D_OCN2D, send_grid_os, recv_grid_os )            ! (in)
    if(my_comp%PRC_rank==0) deallocate( send_grid_os, recv_grid_os )

    call MessageNotify( 'M', module_name, "Read the coefficients (OCN -> SFC) of conservative bilinear interpolation.")    

    if(my_comp%PRC_rank==0) then    
       call set_mappingTable_interpCoef_j99( &
            & GMAPFILENAME_OS_CONSERVE, ocn_comp%GNX, my_comp%GNX, & ! (in)
            & send_grid_os, recv_grid_os, coefS_os_global_conserve & ! (inout)
            & )
    end if    
    call jcup_set_mapping_table( my_comp%name,                        & ! (in)
         & ocn_comp%name, OCN_GRID_2D, my_comp%name, SFC_GRID_2D,     & ! (in) OCN_GRID_2D -> SFC_GRID_2D
         & GMAPTAG_SFC2D_OCN2D_CONSERVE, send_grid_os, recv_grid_os )   ! (in)    
    
    !
    !
    call MessageNotify( 'M', module_name, "Set the operation index ..")    
    
    call set_operation_index(my_comp%name, ocn_comp%name, GMAPTAG_SFC2D_OCN2D)          ! (in)
    call set_operation_index(my_comp%name, ocn_comp%name, GMAPTAG_SFC2D_OCN2D_CONSERVE) ! (in)
    
    
    !***************


    !---

    if(my_comp%PRC_rank==0) then
       call MessageNotify( 'M', module_name, "Set the interpolation coefficients (SFC -> OCN) ..")    

       call set_interpolate_coef( my_comp%name, ocn_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D, coefS_os_global )
    else
       call set_interpolate_coef( my_comp%name, ocn_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D )
    end if
    
    if(my_comp%PRC_rank==0) then
       call set_interpolate_coef( my_comp%name, ocn_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D_CONSERVE, coefS_so_global_conserve )
    else
       call set_interpolate_coef( my_comp%name, ocn_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D_CONSERVE )
    end if
    
    if(my_comp%PRC_rank==0) then
       call MessageNotify( 'M', module_name, "Set the interpolation coefficients (OCN -> SFC) ..")    

       call set_interpolate_coef( ocn_comp%name, my_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D, coefS_os_global )
    else
       call set_interpolate_coef( ocn_comp%name, my_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D )
    end if

    if(my_comp%PRC_rank==0) then
       call set_interpolate_coef( ocn_comp%name, my_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D_CONSERVE, coefS_os_global_conserve )
    else
       call set_interpolate_coef( ocn_comp%name, my_comp%name, &
            & my_comp%name, GMAPTAG_SFC2D_OCN2D_CONSERVE )
    end if

    !
    if(my_comp%PRC_rank==0) then
       call jcup_recv_array( my_comp%name, atm_comp%name, a_Sig1Info )
    end if
    write(*,*) "a_z1Coef:", a_Sig1Info
    
  end subroutine init_jcup_interpolate
  
  !----------------------------------------------------------------------------------
   
  subroutine set_and_put_data( CurrentTimeSec )

    ! モジュール引用; Use statement
    !                
    use mpi

    use jcup_interface, only: &
         & jcup_put_data

    use field_def, only: &
         & set_send_data_2d    

    use dccm_common_params_mod

    use DSFCM_Admin_Variable_mod, only: &
         & xya_WindStressX, xya_WindStressY, &
         & xya_LatHFlx, xya_SenHFlx, xya_QVapMFlx,       &
         & xy_SIceCon, &
         & xya_SfcAlbedo, xya_SfcTemp,                      &
         & xy_RainFall, xy_SnowFall,                        &
         & xya_SUwRFlx, xya_LUwRFlx,                        &
         & xya_SfcHFlx_ns, xya_SfcHFlx_sr, xya_DSfcHFlxDTs, &
         & SFC_PROP_MAX
    
    ! 宣言文; Declareration statements
    !            
    real(DP), intent(in) :: CurrentTimeSec

    ! 局所変数
    ! Local variables
    !

    real(DP) :: xy_SfcHFlx_ns(IAS,JAS)
    real(DP) :: xy_SfcHFlx_sr(IAS,JAS)
    
    integer :: i
    integer :: j


    ! 実行文; Executable statement

    call sfc_set_send_2d( field_as, s2a_WindStressX_id, xya_WindStressX(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_WindStressY_id, xya_WindStressY(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_LUwRFlx_id, xya_LUwRFlx(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_SUwRFlx_id, xya_SUwRFlx(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_SenHFlx_id, xya_SenHFlx(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_QVapMFlx_id, xya_QVapMFlx(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_SfcAlbedo_id, xya_SfcAlbedo(:,:,SFC_PROP_MAX) )
    call sfc_set_send_2d( field_as, s2a_SfcRadTemp_id, xya_SfcTemp(:,:,SFC_PROP_MAX) )

    call sfc_set_send_2d( field_so, s2o_SfcHFlx_ns_id, xya_SfcHFlx_ns(:,:,1) )
    call sfc_set_send_2d( field_so, s2o_SfcHFlx_sr_id, xya_SfcHFlx_sr(:,:,1) )
    call sfc_set_send_2d( field_so, s2o_SnowFall_id, xy_SnowFall(:,:) )
    call sfc_set_send_2d( field_so, s2o_RainFall_id, xy_RainFall(:,:) )
    call sfc_set_send_2d( field_so, s2o_Evap_id, xya_QVapMFlx(:,:,1) )
    call sfc_set_send_2d( field_so, s2o_WindStressX_id, xya_WindStressX(:,:,1) )
    call sfc_set_send_2d( field_so, s2o_WindStressY_id, xya_WindStressY(:,:,1) )
    
  contains
    subroutine sfc_set_send_2d(field, varpID, send_data)
      type(component_field_type), intent(inout) :: field
      integer, intent(in) :: varpID
      real(DP), intent(in) :: send_data(IAS,JAS)
        
      call jcup_put_data( field%varp(varpID)%varp_ptr, &
           & pack(send_data(ISS:IES,JSS:JES), mask=field%mask2d) )
    end subroutine sfc_set_send_2d
    
  end subroutine set_and_put_data

  subroutine get_and_write_data( CurrentTimeSec )

    ! モジュール引用; Use statement
    !                

    !* DCCM
    
    use jcup_interface, only: &
         & jcup_get_data, jcup_recv_array

    use field_def, only: &
         & write_data_2d

    use dccm_common_params_mod

    use DSFCM_Admin_Variable_mod, only: &
         & xya_WindStressX, xya_WindStressY,          &
         & xya_SenHFlx, xya_LatHFlx, xya_QVapMFlx,    &
         & xya_SfcTemp, xya_SfcAlbedo,                &
         & xy_SIceCon, &
         & xya_LUwRFlx, xya_SUwRFlx, xy_LDwRFlx, xy_SDwRFlx,        &
         & xy_RainFall, xy_SnowFall,                  &
         & xya_SfcVelTransCoef, xya_SfcTempTransCoef, &
         & xya_SfcQVapTransCoef,                      &
         & xya_SfcHFlx_sr, xya_SfcHFlx_ns, xya_DSfcHFlxDTs

    use DSFCM_Util_SfcBulkFlux_mod, only: &
         & DSFCM_Util_SfcBulkFlux_Get
    
    ! 宣言文; Declareration statements
    !            
    real(DP), intent(in) :: CurrentTimeSec

    ! 局所変数
    ! Local variables
    !    
    
    integer :: i
    integer :: j

    real(DP) :: xy_WindU(IAS,JAS)
    real(DP) :: xy_WindV(IAS,JAS)
    real(DP) :: xy_SfcAirTemp(IAS,JAS)
    real(DP) :: xy_QVap1(IAS,JAS)
    real(DP) :: xy_Height(IAS,JAS)
    real(DP) :: xy_SfcHeight(IAS,JAS)
    real(DP) :: xy_Press1(IAS,JAS)
    real(DP) :: xy_SfcPress(IAS,JAS)
        
    call sfc_get_write_as(field_as, a2s_WindU_id,a2s_WindU, xy_WindU ) 
    call sfc_get_write_as(field_as, a2s_WindV_id, a2s_WindV, xy_WindV ) 
    call sfc_get_write_as(field_as, a2s_SfcAirTemp_id, a2s_SfcAirTemp, xy_SfcAirTemp )
    call sfc_get_write_as(field_as, a2s_SfcPress_id, a2s_SfcPress, xy_SfcPress ) 
    call sfc_get_write_as(field_as, a2s_QVap1_id, a2s_QVap1, xy_QVap1 )
    call sfc_get_write_as(field_as, a2s_Press1_id, a2s_Press1, xy_Press1 )     
    call sfc_get_write_as(field_as, a2s_LDwRFlx_id, a2s_LDwRFlx, xy_LDwRFlx )
    call sfc_get_write_as(field_as, a2s_SDwRFlx_id, a2s_SDwRFlx, xy_SDwRFlx )    
    call sfc_get_write_as(field_as, a2s_RainFall_id, a2s_RainFall, xy_RainFall )     
    call sfc_get_write_as(field_as, a2s_SnowFall_id, a2s_SnowFall, xy_SnowFall ) 

    call sfc_get_write_so(field_so, o2s_SfcTemp_id, o2s_SfcTemp, xya_SfcTemp(:,:,1) ) 
    call sfc_get_write_so(field_so, o2s_SfcAlbedo_id, o2s_SfcAlbedo, xya_SfcAlbedo(:,:,1) ) 
    call sfc_get_write_so(field_si, i2s_SfcTemp_id, i2s_SfcTemp, xya_SfcTemp(:,:,2) ) 
    call sfc_get_write_so(field_si, i2s_SfcAlbedo_id, i2s_SfcAlbedo, xya_SfcAlbedo(:,:,2) ) 
    call sfc_get_write_so(field_si, i2s_SIceCon_id, i2s_SIceCon, xy_SIceCon ) 

    !---------------------------------------------

    xy_SfcHeight(:,:) = 0d0
    call DSFCM_Util_SfcBulkFlux_Get( &
         & xya_WindStressX, xya_WindStressY,             & ! (out)
         & xya_SenHFlx, xya_QVapMFlx, xya_LatHFlx,       & ! (out)
         & xya_SfcVelTransCoef, xya_SfcTempTransCoef, xya_SfcQVapTransCoef, & ! (out)
         & xya_SUwRFlx, xya_LUwRFlx, &
         & xya_SfcHFlx_ns, xya_SfcHFlx_sr, xya_DSfcHFlxDTs, & ! (out)
         & xy_WindU, xy_WindV, xy_SfcAirTemp, xy_QVap1,  &
         & xy_SDwRFlx, xy_LDwRFlx, &
         & xya_SfcTemp, xya_SfcAlbedo, xy_SIceCon,                &
         & a_Sig1Info, xy_SfcHeight,  xy_SfcPress                 &
         & )
    
  contains
      subroutine sfc_get_write_as(field, vargID, vargName, xy_getdata)
        type(component_field_type), intent(inout) :: field
        integer, intent(in) :: vargID
        character(*), intent(In) :: vargname
        real(DP), intent(out) :: xy_getdata(IAS,JAS)
      
        field%buffer1d(:) = 0d0
        call jcup_get_data(field%varg(vargID)%varg_ptr, field%buffer1d)
      
        if( mod(CurrentTimeSec, dble(AS_COUPLING_CYCLE_SEC)) == 0d0 ) then      
           field%recv_2d(:,:) = unpack(field%buffer1d, field%mask2d, field%recv_2d)
           xy_getdata(ISS:IES,JSS:JES) = field%recv_2d
        end if

      end subroutine sfc_get_write_as

      subroutine sfc_get_write_so(field, vargID, vargName, xy_getdata)
        type(component_field_type), intent(inout) :: field
        integer, intent(in) :: vargID
        character(*), intent(In) :: vargname
        real(DP), intent(out) :: xy_getdata(IAS,JAS)
      
        field%buffer1d(:) = 0d0
        call jcup_get_data(field%varg(vargID)%varg_ptr, field%buffer1d)
      
        if( mod(CurrentTimeSec, dble(SO_COUPLING_CYCLE_SEC)) == 0d0 ) then      
           field%recv_2d(:,:) = unpack(field%buffer1d, field%mask2d, field%recv_2d)
           xy_getdata(ISS:IES,JSS:JES) = field%recv_2d
        end if

      end subroutine sfc_get_write_so
      
  end subroutine get_and_write_data
  
end module dccm_sfc_mod
