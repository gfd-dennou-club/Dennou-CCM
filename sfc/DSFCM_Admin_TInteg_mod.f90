!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSFCM_Admin_TInteg_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify

  use dc_calendar, only: &
       & DC_CAL, DC_CAL_DATE


  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSFCM_Admin_TInteg_Init, DSFCM_Admin_TInteg_Final
  public :: EndTemporalInteg
  public :: DSFCM_Admin_TInteg_AdvanceLongTStep
  public :: DSFCM_Admin_TInteg_AdvanceShortTStep

  ! 公開変数
  ! Public variable
  !
  integer, parameter, public :: nLongTimeLevel = 3
  integer, save, public :: Bl, Nl, Al
  integer, public :: TIMELV_ID_B = 3
  integer, public :: TIMELV_ID_N = 2
  integer, public :: TIMELV_ID_A = 1

  integer, parameter, public :: nShortTimeLevel = 3
  integer, save, public :: Bs, Ns, As


  real(DP), save, public :: DelTime
  integer, save, public :: SubCycleNum

  real(DP), save, public :: TimeSecA
  real(DP), save, public :: TimeSecN
  real(DP), save, public :: TimeSecB
  
  real(DP), save, public :: CurrentTime

  real(DP), save, public :: RestartTime
  real(DP), save, public :: EndTime

  real(DP), save, public :: IntegTime
  integer, save, public :: CurrentTimeStep
  integer, save, public :: CurrentShortTimeStep
  
  ! Temporal integration scheme

  !
  logical, save, public :: SemiImplicitFlag

  character(TOKEN) :: cal_type  
  type(DC_CAL_DATE), save, public :: InitDate
  type(DC_CAL_DATE), save, public :: RestartDate
  type(DC_CAL_DATE), save, public :: EndDate

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSFCM_Admin_TInteg_mod' !< Module Name

  real(DP), save, public :: ProgMessageInterVal

contains

  !>
  !!
  !!
  subroutine DSFCM_Admin_TInteg_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !


    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 実行文; Executable statements
    !

    ! Set DelTime, SubCycleNum by reading namelist. 
    call read_nmlData( configNmlFileName  )

    DelTime = 1200D0
    CurrentTime = ReStartTime
    TimeSecA = CurrentTime + DelTime
    TimeSecN = CurrentTime
    timeSecB = CurrentTime - DelTime
    CurrentTimeStep = 1
    CurrentShortTimeStep = 1

    Bl = 1; Nl = 2; Al = 3
    Bs = 1; Ns = 2; As = 3

  end subroutine DSFCM_Admin_TInteg_Init

  !>
  !!
  !!
  subroutine DSFCM_Admin_TInteg_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSFCM_Admin_TInteg_Final

  !> @brief 
  !!
  !! @return 
  !!
  function EndTemporalInteg() result(ret)
    
    ! 宣言文; Declaration statement
    !
    logical :: ret

    ! 実行文; Executable statement
    !

    ret = (CurrentTime > EndTime)
    
  end function EndTemporalInteg

  !> @brief 
  !!
  !!
  subroutine DSFCM_Admin_TInteg_AdvanceLongTStep()

    !
    !
    use dc_calendar, only: &
         & DCCalCreate, DCCalDateCreate, &
         & DCCalConvertByUnit, DCCalDateDifference, &
         & DCCalDateInquire, DCCalDateEval
    
    ! 宣言文; Declaration statement
    !
    type(DC_CAL_DATE) :: CurrentDate
    character(TOKEN) :: InitDateStr, RestartDateStr, EndDateStr, CurrentDateStr

    ! 局所変数
    ! Local variables
    !
    integer :: oldBl
    
    ! 実行文; Executable statement
    !
    
    CurrentTime = CurrentTime + DelTime
    CurrentTimeStep = CurrentTimeStep + 1
    CurrentShortTimeStep = 1

    oldBl = Bl
    Bl = Nl
    Nl = Al
    Al = oldBl

    TimeSecB = TimeSecN
    TimeSecN = CurrentTime
    TimeSecA = CurrentTime + DelTime
    
    if( mod(CurrentTime, ProgMessageInterVal) == 0d0 ) then

       call DCCalDateEval(InitDate, CurrentTime, 'sec', date=CurrentDate)
       call DCCalDateInquire(InitDateStr, date=InitDate)
       call DCCalDateInquire(EndDateStr, date=EndDate)
       call DCCalDateInquire(CurrentDateStr, date=CurrentDate)
       call DCCalDateInquire(RestartDateStr, date=RestartDate)

       call MessageNotify( 'M', module_name, "Current Date=%a [%a-%a (Restart Date %a)]", &
            & ca=(/ CurrentDateStr,  InitDateStr, EndDateStr, RestartDateStr /))
    end if

  end subroutine DSFCM_Admin_TInteg_AdvanceLongTStep

  subroutine DSFCM_Admin_TInteg_AdvanceShortTStep()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    integer :: oldBs    
    
    ! 実行文; Executable statement
    !

    oldBs = Bs
    Bs = Ns
    Ns = As
    As = oldBs

    CurrentShortTimeStep = CurrentShortTimeStep + 1

  end subroutine DSFCM_Admin_TInteg_AdvanceShortTStep


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    use dc_calendar, only: &
         & DCCalCreate, DCCalDateCreate, &
         & DCCalConvertByUnit, DCCalDateDifference, &
         & DCCalDateInquire, DCCalDateEval

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
    ! IOSTAT of NAMELIST read

    real(DP) :: DelTimeVal
    character(TOKEN) :: DelTimeUnit
    real(DP) :: IntegTimeVal
    character(TOKEN) :: IntegTimeUnit
    real(DP) :: RestartTimeVal
    character(TOKEN) :: RestartTimeUnit
    real(DP) :: ProgMessageIntVal
    character(TOKEN) :: ProgMessageIntUnit

    integer :: InitYear, InitMonth, InitDay, InitHour, InitMin
    integer :: EndYear, EndMonth, EndDay, EndHour, EndMin
    real(DP) :: InitSec, EndSec
    character(TOKEN) :: InitDateStr, EndDateStr, RestartDateStr
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /sfc_temporalInteg_nml/ &
         & cal_type,                                                 &
         & DelTimeVal, DelTimeUnit, SemiImplicitFlag,                &
         & IntegTimeVal, IntegTimeUnit,                              &
         & InitYear, InitMonth, InitDay, InitHour, InitMin, InitSec, &
         & EndYear, EndMonth, EndDay, EndHour, EndMin, EndSec, &
         & RestartTimeVal, RestartTimeUnit, &
         & SubCycleNum, &
         & ProgMessageIntVal, ProgMessageIntUnit

    
    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    cal_type = 'noleap'
    
    DelTimeVal  = 10d0
    DelTimeUnit = 'min'
    SemiImplicitFlag = .true.
    SubCycleNum = 2
    
    RestartTimeVal  = 0d0
    RestartTimeUnit = 'sec'
    IntegTimeVal  = -1d0
    IntegTimeUnit  = 'sec'

    InitYear=2000; InitMonth=1; InitDay=1; InitHour=0; InitMin=0; InitSec=0d0;
    EndYear=2000; EndMonth=1; EndDay=2; EndHour=0; EndMin=0; EndSec=0d0;

    ProgMessageIntVal = -1d0

    ! Default valuse for Semi-implicit scheme
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                      ! (in)
            & nml = sfc_temporalInteg_nml, &  ! (out)
            & iostat = iostat_nml )           ! (out)
    end if

    ! Determine time to start and finish a temporal integration. 
    !
    call DCCalCreate(cal_type)

    DelTime = DCCalConvertByUnit(DelTimeVal, DelTimeUnit, 'sec')
    RestartTime = DCCalConvertByUnit(RestartTimeVal, RestartTimeUnit, 'sec')

    call DCCalDateCreate( InitYear, InitMonth, InitDay, InitHour, InitMin, InitSec, &
         & InitDate ) !(out)

    call DCCalDateEval( InitDate, RestartTimeVal, RestartTimeUnit, &
         & date=RestartDate) !(out)

    if( IntegTimeVal > 0d0 ) then
       call DCCalDateEval( RestartDate, IntegTimeVal, IntegTimeUnit, &
            & date=EndDate) !(out)
    else
       call DCCalDateCreate(EndYear, EndMonth, EndDay, EndHour, EndMin, EndSec, &
            & EndDate ) !(out)
    end if

    IntegTime = DCCalDateDifference(RestartDate, EndDate)
    EndTime = RestartTime + IntegTime

    ! Determine the time interval to report the progress of temporal integration. 
    !
    if( ProgMessageIntVal < 0.0 ) then
       ProgMessageInterVal = 1d-02*IntegTime
    else
       ProgMessageInterVal = DCCalConvertByUnit(ProgMessageIntVal, ProgMessageIntUnit, 'sec')
    end if

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  cal_type             = %a', ca=(/ cal_type /))
    call MessageNotify( 'M', module_name, '  DelTime              = %f [%c]', d=(/ DelTimeVal /), c1=trim(DelTimeUnit) )
    call MessageNotify( 'M', module_name, '  SemiImplicitFlag     = %b     ', L=(/ SemiImplicitFlag /))
    call MessageNotify( 'M', module_name, '  SubCycleNum          = %d [time]', i=(/ SubCycleNum /))
    call MessageNotify( 'M', module_name, '  RestartTime          = %f [%c]',  d=(/ RestartTimeVal /), c1=trim(RestartTimeUnit) )
    call MessageNotify( 'M', module_name, '  RestartTimeSec       = %f [sec]', d=(/ RestartTime    /) )
    call MessageNotify( 'M', module_name, '  IntegTime            = %f [sec]', d=(/ IntegTime /))

    call DCCalDateInquire(InitDateStr, date=InitDate)
    call DCCalDateInquire(EndDateStr, date=EndDate)
    call DCCalDateInquire(RestartDateStr, date=RestartDate)
    call MessageNotify( 'M', module_name, '  TimeIntegPeriod      = %a - %a', ca=(/InitDateStr, EndDateStr/)) 
    call MessageNotify( 'M', module_name, '  Init/RestartDate     = %c', c1=RestartDateStr )

  end subroutine read_nmlData

end module DSFCM_Admin_TInteg_mod

