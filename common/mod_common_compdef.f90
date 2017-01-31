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
  
  use jcup_data, only: &
       & varp_type
  
  use jcup_interface, only: &
       & jcup_initialize,        &
       & jcup_set_new_comp,      &
       & jcup_get_model_id,      &
       & jcup_get_mpi_parameter, &
       & jcup_init_time

!!$  use field_def
  
  
  ! 宣言文; Declareration statements
  !    
  implicit none
  private

  !
  type, public :: ComponentDef

     character(TOKEN) :: NAME = ''
     integer          :: id   = -1
     
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
  type(ComponentDef), public, pointer :: CompDef_own => null()
  
  character(STRING), public :: GMAPFILENAME_AO
  character(STRING), public :: GMAPFILENAME_AO_CONSERVE
  character(STRING), public :: GMAPFILENAME_OA
  character(STRING), public :: GMAPFILENAME_OA_CONSERVE
    
  integer, public :: AO_COUPLING_CYCLE_SEC
!  integer, public :: OA_COUPLING_CYCLE_SEC
  integer, public :: SO_COUPLING_CYCLE_SEC
  
  public :: ComponentDef_Init
  public :: ComponentDef_Final
  public :: ComponentDef_Share

  public :: common_read_config
  
  !----------------------------------------------------------
  
  character(*), parameter :: module_name = "mod_common_compdef"
  
contains

  subroutine ComponentDef_Init( my_comp, name, modelname )

    
    type(ComponentDef), target, intent(inout) :: my_comp
    character(*), intent(in) :: name
    character(*), intent(in) :: modelname

    CompDef_own => my_comp
    
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

  !-------------------------------------------------------------------
  
  subroutine ComponentDef_Final()
    
  end subroutine ComponentDef_Final

  !-------------------------------------------------------------------
  
  subroutine ComponentDef_Share()

    use mod_common_params
    
!!$    write(*,*) "Share conmponent defininitions.."
!!$    write(*,*) "CompDef_Atm:", CompDef_atm%name
!!$    write(*,*) "CompDef_Ocn:", CompDef_ocn%name
    
    call broadcast_compdef( CompDef_atm, COMPNAME_ATM )
    call broadcast_compdef( CompDef_ocn, COMPNAME_OCN )
!!$    call broadcast_compdef( sice_comp )

!!$    write(*,*) "After CompDef_Atm:", CompDef_atm%name
!!$    write(*,*) "After CompDef_Ocn:", CompDef_ocn%name

  contains
    subroutine broadcast_compdef( comp, comp_name )

      use jcup_comp, only: &
           & get_comp_id_from_name
      
      use jcup_mpi_lib, only: &
           & jml_GetLeaderRank, &
           & jml_BcastGlobal


      type(ComponentDef), intent(inout) :: comp
      character(*), intent(in) :: comp_name

      integer :: src_id
      integer :: src_rank
      integer :: comp_id(1)
      integer :: GRID_NUM(3)
      integer :: ITIME(6)
      
      src_id = get_comp_id_from_name( comp_name )
      src_rank = jml_GetLeaderRank( src_id )

      call jml_BcastGlobal( comp%name, src_rank )
      call jml_BcastGlobal( comp%MODELNAME, src_rank )

      comp_id(:) = comp%id
      call jml_BcastGlobal( comp_id, 1, 1, src_rank )
      comp%id = comp_id(1)
      
      GRID_NUM(1:3) = (/ comp%GNX, comp%GNY, comp%GNZ /)
      call jml_BcastGlobal( GRID_NUM, 1, 3, src_rank )
      comp%GNX = GRID_NUM(1)
      comp%GNY = GRID_NUM(2)
      comp%GNZ = GRID_NUM(3)

      GRID_NUM(1:2) = (/ comp%LNX, comp%LNY /)
      call jml_BcastGlobal( GRID_NUM, 1, 2, src_rank )
      comp%LNX = GRID_NUM(1)
      comp%LNY = GRID_NUM(2)

      call jml_BcastGlobal( comp%InitTimeInfo, 1, 6, src_rank)
      

    end subroutine broadcast_compdef
    
  end subroutine ComponentDef_Share
  
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
         & AO_COUPLING_CYCLE_SEC,    &
         & GMAPFILENAME_AO,          &
         & GMAPFILENAME_OA,          &
         & GMAPFILENAME_AO_CONSERVE, &
         & GMAPFILENAME_OA_CONSERVE


    integer :: unit_nml
    integer :: ierr

    ! 実行文; Executable statement
    !
    
    GMAPFILENAME_OA = 'gmap-OCN_Pl42-ATM_T42.dat'
    GMAPFILENAME_AO = 'gmap-ATM_T42-OCN_Pl42.dat'
    GMAPFILENAME_OA_CONSERVE = 'gmap-OCN_Pl42-ATM_T42_CONSERVE.dat'
    GMAPFILENAME_AO_CONSERVE = 'gmap-ATM_T42-OCN_Pl42_CONSERVE.dat'
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
    call MessageNotify( 'M', module_name, "GMAPFILENAME_OA =%a             ", ca=(/ GMAPFILENAME_OA          /) )
    call MessageNotify( 'M', module_name, "GMAPFILENAME_AO =%a             ", ca=(/ GMAPFILENAME_AO          /) )
    call MessageNotify( 'M', module_name, "GMAPFILENAME_OA_CONSERVE =%a    ", ca=(/ GMAPFILENAME_OA_CONSERVE /) )
    call MessageNotify( 'M', module_name, "GMAPFILENAME_AO_CONSERVE =%a    ", ca=(/ GMAPFILENAME_AO_CONSERVE /) )
        
  end subroutine common_read_config
  
end module mod_common_compdef
