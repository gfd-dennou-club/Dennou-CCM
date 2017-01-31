!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief Main program to generate a data file providing the grid-remapping information. 
!! 
!! @author Kawai Yuta
!!
!!
program gmapgen_main

  ! モジュール引用; Use statements
  !  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify
  
  use grid_mapping_util, only: &
       & gen_gridmapfile_lonlat2lonlat

  use grid_mapping_util_jones99, only: &
       & gen_gridmapfile_lonlat2lonlat_j99 => gen_gridmapfile_lonlat2lonlat

  ! 宣言文 ; Declaration statements  
  !  
  implicit none

  character(*), parameter :: DEFAULT_GMAPGEN_CONFIGNML = "gmapgen.conf"
  character(STRING) :: configNmlName
  
  integer :: IMA, JMA, KMA
  integer :: NMA
  character(STRING) :: gmapfile_AO_NAME
  real(DP), allocatable :: x_LonA(:)
  real(DP), allocatable :: y_LatA(:)
  real(DP), allocatable :: x_IntWtLonA(:)
  real(DP), allocatable :: y_IntWtLatA(:)
  
  integer :: IMO, JMO, KMO
  integer :: NMO
  character(STRING) :: gmapfile_OA_NAME
  real(DP), allocatable :: x_LonO(:)
  real(DP), allocatable :: y_LatO(:)
  real(DP), allocatable :: x_IntWtLonO(:)
  real(DP), allocatable :: y_IntWtLatO(:)
  
  character(*), parameter :: PROGRAM_NAME = "gmapgen_main"

  logical :: ConservativeFlag
  
  ! 実行文; Executable statements
  !
  
  ! Read configuration form namelist
  !
  call read_config()
  
  ! Get grid information
  !
  call get_LonLatGrid( x_LonA, y_LatA, x_IntWtLonA, y_IntWtLatA, & ! (out)
       & IMA, JMA, NMA )                 ! (in)
  call get_LonLatGrid( x_LonO, y_LatO, x_IntWtLonO, y_IntWtLatO, & ! (out)
       & IMO, JMO, NMO )                 ! (in)
  
!!$  write(*,*) "=Atm:"
!!$  write(*,*) "*Lon=", x_LonA
!!$  write(*,*) "*Lat=", y_LatA
!!$  write(*,*) "=Ocn:"
!!$  write(*,*) "*Lon=", x_LonO
!!$  write(*,*) "*Lat=", y_LatO

  !
  !
  if (.not. ConservativeFlag) then
     call gen_gridmapfile_lonlat2lonlat( gmapfile_AO_NAME, &
          & x_LonA, y_LatA, x_LonO, y_LatO )

     call gen_gridmapfile_lonlat2lonlat( gmapfile_OA_NAME, &
          & x_LonO, y_LatO, x_LonA, y_LatA ) 
  else
     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_AO_NAME,  &
          & x_LonA, y_LatA, x_LonO, y_LatO,                     &
          & x_IntWtLonA, y_IntWtLatA, x_IntWtLonO, y_IntWtLatO )
  
     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_OA_NAME,  &
          & x_LonO, y_LatO, x_LonA, y_LatA,                     &
          & x_IntWtLonO, y_IntWtLatO, x_IntWtLonA, y_IntWtLatA)
  end if
  
  ! Output some information about grid mapping table
  !
  call check_mappingTable( gmapfile_AO_NAME, IMA, IMO )
  call check_mappingTable( gmapfile_OA_NAME, IMO, IMA )

contains
  subroutine read_config()


    ! モジュール引用; Use statements
    !          
    use dc_iounit, only: FileOpen
    use dc_types, only: STDOUT
    
    use optionparser_mod

    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !    
    NAMELIST /PARAM_DCCM_GRID/   &
         & IMA, JMA, KMA, NMA, &
         & IMO, JMO, KMO, NMO

    NAMELIST /PARAM_GMAPGEN/     &
         & gmapfile_AO_NAME,     &
         & gmapfile_OA_NAME,     &
         & ConservativeFlag


    integer :: unit_nml
    integer :: ierr

    ! 実行文; Executable statement
    !
    
    call OptionParser_Init()
    call OptionParser_GetInfo( configNmlName, DEFAULT_GMAPGEN_CONFIGNML)
    call OptionParser_Final()

    IMA = 64
    JMA = 32
    KMA = 26
    NMA = 21

    IMO = 64
    JMO = 32
    KMO = 26
    NMO = 21

    ConservativeFlag = .false.
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlName) /= '' ) then
       call MessageNotify( 'M', PROGRAM_NAME, "reading namelist '%a'", ca=(/ configNmlName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlName, mode = 'r' )  ! (in)

       rewind( unit_nml )
       read( unit_nml, &                ! (in)
            & nml = PARAM_DCCM_GRID, &  ! (out)
            & iostat = ierr )           ! (out)       

       rewind( unit_nml )
       read( unit_nml, &                ! (in)
            & nml = PARAM_GMAPGEN,   &  ! (out)
            & iostat = ierr )           ! (out)              
    end if

    call MessageNotify( 'M', PROGRAM_NAME, "ATM: (IM, JM, KM)=(%d,%d,%d)", i=(/ IMA, JMA, KMA /) )
    call MessageNotify( 'M', PROGRAM_NAME, "OCN: (IM, JM, KM)=(%d,%d,%d)", i=(/ IMO, JMO, KMO /) )
    call MessageNotify( 'M', PROGRAM_NAME, "gmapfile_OA      =%a        ", ca=(/ gmapfile_OA_NAME /) )
    call MessageNotify( 'M', PROGRAM_NAME, "gmapfile_AO      =%a        ", ca=(/ gmapfile_AO_NAME /) )
    call MessageNotify( 'M', PROGRAM_NAME, "conservativeFlag =%b        ", L=(/ ConservativeFlag /) )
    
  end subroutine read_config
  
  subroutine get_LonLatGrid( &
       & x_Lon, y_Lat, x_LonIntWt, y_LatIntWt, &  ! (inout)
       & iMax, jMax, nMax )                       ! (in)

    ! モジュール引用; Use statements
    !      
    use w_module, only: &
         & w_Initial, w_Finalize,     &
         & xy_Lon, xy_Lat,            &
         & x_Lon_Weight, y_Lat_Weight
    
    use w_zonal_module, only: &
         & w_Initial_zonal => w_Initial, &
         & w_Finalize_zonal => w_Finalize, &       
         & xy_Lon_zonal => xy_Lon, &
         & xy_Lat_zonal => xy_Lat, &
         & x_Lon_Weight_zonal => x_Lon_Weight, &
         & y_Lat_Weight_zonal => y_Lat_Weight

    ! 宣言文; Declaration statement
    !    
    integer, intent(in) :: iMax, jMax, nMax
    real(DP), intent(inout), allocatable :: x_Lon(:)
    real(DP), intent(inout), allocatable :: y_Lat(:)
    real(DP), intent(inout), allocatable :: x_LonIntWt(:)
    real(DP), intent(inout), allocatable :: y_LatIntWt(:)

    ! 実行文; Executable statement
    !
    
    if(iMax==1) then
       call w_Initial_zonal(nMax, iMax, jMax)
    else
       call w_Initial(nMax, iMax, jMax)     
    end if

    allocate(x_Lon(iMax), y_Lat(jMax))
    allocate(x_LonIntWt(iMax), y_LatIntWt(jMax))

    if(iMax==1) then
       x_Lon(:) = xy_Lon_zonal(:,1); y_Lat(:) = xy_Lat_zonal(0,:)
       x_LonIntWt(:) = x_Lon_Weight_zonal
       y_LatIntWt(:) = y_Lat_Weight_zonal
       call w_Finalize_zonal()
    else
       x_Lon(:) = xy_Lon(:,1); y_Lat(:) = xy_Lat(0,:)       
       x_LonIntWt(:) = x_Lon_Weight
       y_LatIntWt(:) = y_Lat_Weight
       call w_Finalize()
    end if
    
  end subroutine get_LonLatGrid
  
  subroutine check_mappingTable(gmapfilename, GNXS, GNXR)

    ! モジュール引用; Use statements
    !          
    use grid_mapping_util, only: set_mappingTable_interpCoef

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: gmapfilename
    integer, intent(in) :: GNXS, GNXR

    integer, allocatable, dimension(:) :: send_index, recv_index
    real(DP), allocatable, dimension(:) :: coef

    ! 実行文; Executable statement
    !
    
    call set_mappingTable_interpCoef( gmapfilename, GNXS, GNXR, &
         & send_index, recv_index, coef )

    write(*,*) "* Set mapping table and coeffecient for interpolation.. file=", trim(gmapfilename)
    write(*,*) "  -- Grid Mapping for 1-4 operation to interpolate ---------"
    write(*,*) "  recv_index=", recv_index(1:4), "send_index=", send_index(1:4)
    write(*,*) "  coef=", coef(1:4)
    
  end subroutine check_mappingTable
  
end program gmapgen_main

