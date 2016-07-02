!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module mod_common_compdef

  ! モジュール引用; Use statements
  !

  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use mod_common_params, only: &
       & NUM_DCCM_COMP,                    &
       & JCUP_LOG_LEVEL, JCUP_LOG_STDERROR

  use field_def
  
  use jcup_data, only: &
       & varp_type
  
  use jcup_interface, only: &
       & jcup_initialize,        &
       & jcup_set_new_comp,      &
       & jcup_get_mpi_parameter, &
       & jcup_init_time
  
  
  ! 宣言文; Declareration statements
  !    
  implicit none
  private

  !
  type, public :: ComponentDef
     character(TOKEN) :: NAME
     integer          :: id
     
     character(TOKEN) :: MODELNAME
     
     integer :: GNX
     integer :: GNY
     integer :: GNZ
     integer :: LNX
     integer :: LNY     
     integer :: InitTimeInfo(6) ! Year, Month, Day, Hour, Min, Sec
     real(DP) :: RestartTimeSec

     integer :: PRC_comm
     integer :: PRC_group
     integer :: PRC_rank
     integer :: PRC_size
     integer :: PRC_NX
     integer :: PRC_NY     
     
     real(DP) :: DelTime
     integer  :: tstep
     logical  :: loop_end_flag

  end type ComponentDef
  
  type(ComponentDef), public :: CompDef_atm
  type(ComponentDef), public :: CompDef_ocn
  type(ComponentDef), public :: CompDef_sice
  type(ComponentDef), public :: CompDef_land

  character(TOKEN), public :: GMAPFILENAME_AO
  character(TOKEN), public :: GMAPFILENAME_OA
    
  integer, public :: AO_COUPLING_CYCLE_SEC
  integer, public :: OA_COUPLING_CYCLE_SEC
  integer, public :: SO_COUPLING_CYCLE_SEC
  
  public :: ComponentDef_Init
  public :: ComponentDef_Final
  
  
  !----------------------------------------------------------
  
  character(*), parameter :: module_name = "mod_common_compdef"
  
contains

  subroutine ComponentDef_Init( my_comp, name, modelname )

    type(ComponentDef), intent(inout) :: my_comp
    character(*), intent(in) :: name
    character(*), intent(in) :: modelname

    my_comp%NAME = name
    my_comp%MODELNAME = modelname

    ! Initialize a ATM component for Jcup
    call jcup_set_new_comp( my_comp%name )
    
    call jcup_initialize( my_comp%name,  & ! (in)
         & LOG_LEVEL=JCUP_LOG_LEVEL, LOG_STDERR=JCUP_LOG_STDERROR & ! (in)
         & ) 

    call jcup_get_model_id( my_comp%name, & ! (in)
         & my_comp%id )                     ! (out)
    
    ! Get some information associated with MPI.
    call jcup_get_mpi_parameter( my_comp%name,             & ! (in)
         & my_comp%PRC_comm, my_comp%PRC_group,            & ! (out)
         & my_comp%PRC_size, my_comp%PRC_rank )              ! (out)

    
    !-----------------------------------------------------------------
    
    call MessageNotify( 'M', module_name, &
         & "* Initialize a component ----------------")

    call MessageNotify( 'M', module_name, &
         & "  NAME=%a, id=%d, MODELNAME=%a", &
         & i=(/ my_comp%id /), ca=(/trim(my_comp%name), trim(my_comp%MODELNAME)/) )

    call MessageNotify( 'M', module_name, &
         & "  MPI comm=%d, group=%d, size=%d,  rank=%d", &
         & i=(/ my_comp%PRC_comm, my_comp%PRC_group, my_comp%PRC_size, my_comp%PRC_rank /) &
         & )
    
  end subroutine ComponentDef_Init

  subroutine ComponentDef_Final()
    
  end subroutine ComponentDef_Final

  !----------------------------------------------------------

  subroutine common_read_config( configNmlName )

    use dc_iounit, only: FileOpen
    use dc_types, only: STDOUT
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 局所変数
    ! Local variables
    !    
    NAMELIST /PARAM_DCCM_COMMON/   &
         & AO_COUPLING_CYCLE_SEC,  &
         & GMAPFILENAME_AO,        &
         & GMAPFILENAME_OA


    integer :: unit_nml
    integer :: ierr

    ! 実行文; Executable statement
    !
    
    GMAPFILENAME_OA = 'gmap-OCN_Pl42-ATM_T42.dat'
    GMAPFILENAME_AO = 'gmap-ATM_T42-OCN_Pl42.dat'
    AO_COUPLING_CYCLE_SEC = 7200
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlName, mode = 'r' )  ! (in)

       rewind( unit_nml )
       read( unit_nml, &                  ! (in)
            & nml = PARAM_DCCM_COMMON, &  ! (out)
            & iostat = ierr )             ! (out)       

    end if

    call MessageNotify( 'M', module_name, "AO_COUPLING_CYCLE_SEC=%d [sec]  ", i=(/ AO_COUPLING_CYCLE_SEC /))
    call MessageNotify( 'M', module_name, "GMAPFILENAME_OA =%a             ", ca=(/ GMAPFILENAME_OA /) )
    call MessageNotify( 'M', module_name, "GMAPFILENAME_AO =%a             ", ca=(/ GMAPFILENAME_AO /) )
        
  end subroutine common_read_config
  
end module mod_common_compdef
