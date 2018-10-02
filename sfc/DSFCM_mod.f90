module DSFCM_mod

  ! モジュール引用; Use statements
  !
 
  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* DSFCM
  
  use DSFCM_Admin_TInteg_mod, only: &
       & DSFCM_Admin_TInteg_Init, DSFCM_Admin_TInteg_Final, &
       & DSFCM_Admin_TInteg_AdvanceLongTStep,               &
       & EndTemporalInteg
  
  use DSFCM_Admin_Variable_mod
  use DSFCM_Admin_Grid_mod

  use DSFCM_Util_SfcBulkFlux_mod, only: &
       & DSFCM_Util_SfcBulkFlux_Init, DSFCM_Util_SfcBulkFlux_Final, &
       & DSFCM_Util_SfcBulkFlux_Get
  
  implicit none
  private  

  public :: DSFCM_main_Init, DSFCM_main_Final
  public :: DSFCM_main_setup
  public :: DSFCM_main_shutdown
  public :: DSFCM_main_advance_timestep

contains

  subroutine DSFCM_main_Init()

    ! 宣言文; Declareration statements
    !            
    
    ! 実行文; Executable statements
    !

    
  end subroutine DSFCM_main_Init

  subroutine DSFCM_main_Final()
    
  end subroutine DSFCM_main_Final

  subroutine DSFCM_main_advance_timestep( &
       & tstep,                           & ! (in)
       & loop_end_flag,                   & ! (inout)
       & skip_flag                        & ! (in)
       & )
           
    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: tstep
    logical, intent(inout) :: loop_end_flag
    logical, intent(in) :: skip_flag
    
    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !

    if( EndTemporalInteg() ) then
       loop_end_flag = .true.; return
    else
       loop_end_flag = .false.
    end if

    if(tstep == 0) then
    else            
       call DSFCM_Admin_TInteg_AdvanceLongTStep()
    end if
    
    !- Output  --------------------------------------------------------------------


    
  end subroutine DSFCM_main_advance_timestep
  
  !----------------------------------
  
  subroutine DSFCM_main_setup( configNmlFileName, GNX, GNY )

    ! モジュール引用; Use statement
    !


    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName
    integer, intent(in) :: GNX
    integer, intent(in) :: GNY
    
    ! 実行文; Executable statement
    !

    call DSFCM_Admin_TInteg_Init(configNmlFileName)
    call DSFCM_Admin_Grid_Init(GNX, GNY)
    call DSFCM_Admin_Variable_Init()
    call DSFCM_Util_SfcBulkFlux_Init()
    
  end subroutine DSFCM_main_setup

  subroutine DSFCM_main_shutdown()

    call DSFCM_Util_SfcBulkFlux_Final()    
    call DSFCM_Admin_Variable_Final()
    call DSFCM_Admin_Grid_Final()
    call DSFCM_Admin_TInteg_Final()
    
  end subroutine DSFCM_main_shutdown

  !----------------------

  
end module DSFCM_mod
