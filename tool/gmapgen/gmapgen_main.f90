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
  real(DP), allocatable :: x_LonA(:)
  real(DP), allocatable :: y_LatA(:)
  real(DP), allocatable :: x_IntWtLonA(:)
  real(DP), allocatable :: y_IntWtLatA(:)

  integer :: IMO, JMO, KMO
  integer :: NMO
  real(DP), allocatable :: x_LonO(:)
  real(DP), allocatable :: y_LatO(:)
  real(DP), allocatable :: x_IntWtLonO(:)
  real(DP), allocatable :: y_IntWtLatO(:)

  integer :: IMS, JMS
  real(DP), save, allocatable :: x_LonS(:)
  real(DP), save, allocatable :: y_LatS(:)
  real(DP), save, allocatable :: x_IntWtLonS(:)
  real(DP), save, allocatable :: y_IntWtLatS(:)
  
  character(STRING) :: gmapfile_AO_NAME
  integer :: interp_order_AO
  character(STRING) :: gmapfile_OA_NAME
  integer :: interp_order_OA
  character(STRING) :: gmapfile_AS_NAME  
  integer :: interp_order_AS
  character(STRING) :: gmapfile_SA_NAME
  integer :: interp_order_SA
  character(STRING) :: gmapfile_OS_NAME
  integer :: interp_order_OS  
  character(STRING) :: gmapfile_SO_NAME
  integer :: interp_order_SO
  
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
  

  call generate_surface_exchage_grid()

  write(*,*) "=Atm:"
  write(*,*) "*Lon=", x_LonA
  write(*,*) "*Lat=", y_LatA

  write(*,*) "=Ocn:"
  write(*,*) "*Lon=", x_LonO
  write(*,*) "*Lat=", y_LatO

  write(*,*) "=Sfc:"
  write(*,*) "*Lon=", x_LonS
  write(*,*) "*Lat=", y_LatS
  
  !
  !
  if (.not. ConservativeFlag) then
     call gen_gridmapfile_lonlat2lonlat( gmapfile_AO_NAME, &
          & x_LonA, y_LatA, x_LonO, y_LatO )

     call gen_gridmapfile_lonlat2lonlat( gmapfile_OA_NAME, &
          & x_LonO, y_LatO, x_LonA, y_LatA ) 

     call gen_gridmapfile_lonlat2lonlat( gmapfile_AS_NAME, &
          & x_LonA, y_LatA, x_LonS, y_LatS )
     
     call gen_gridmapfile_lonlat2lonlat( gmapfile_SA_NAME, &
          & x_LonS, y_LatS, x_LonA, y_LatA )

     call gen_gridmapfile_lonlat2lonlat( gmapfile_OS_NAME, &
          & x_LonO, y_LatO, x_LonS, y_LatS )
     
     call gen_gridmapfile_lonlat2lonlat( gmapfile_SO_NAME, &
          & x_LonS, y_LatS, x_LonO, y_LatO )
     
  else
     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_AO_NAME,  &
          & x_LonA, y_LatA, x_LonO, y_LatO,                     &
          & x_IntWtLonA, y_IntWtLatA, x_IntWtLonO, y_IntWtLatO, &
          & interp_order_AO )
  
     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_OA_NAME,  &
          & x_LonO, y_LatO, x_LonA, y_LatA,                     &
          & x_IntWtLonO, y_IntWtLatO, x_IntWtLonA, y_IntWtLatA, &
          & interp_order_OA )

     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_AS_NAME,  &
          & x_LonA, y_LatA, x_LonS, y_LatS,                     &
          & x_IntWtLonA, y_IntWtLatA, x_IntWtLonS, y_IntWtLatS, &
          & interp_order_AS )
     
     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_SA_NAME,  &
          & x_LonS, y_LatS, x_LonA, y_LatA,                     &
          & x_IntWtLonS, y_IntWtLatS, x_IntWtLonA, y_IntWtLatA, &
          & interp_order_SA )

     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_OS_NAME,  &
          & x_LonO, y_LatO, x_LonS, y_LatS,                     &
          & x_IntWtLonO, y_IntWtLatO, x_IntWtLonS, y_IntWtLatS, &
          & interp_order_OS )
     
     call gen_gridmapfile_lonlat2lonlat_j99( gmapfile_SO_NAME,  &
          & x_LonS, y_LatS, x_LonO, y_LatO,                     &
          & x_IntWtLonS, y_IntWtLatS, x_IntWtLonO, y_IntWtLatO, &
          & interp_order_SO )     
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
         & interp_order_AO,      &
         & gmapfile_OA_NAME,     &
         & interp_order_OA,      &
         & gmapfile_SA_NAME,     &
         & interp_order_SA,      &
         & gmapfile_AS_NAME,     &
         & interp_order_AS,      &         
         & gmapfile_SO_NAME,     &
         & interp_order_SO,      &         
         & gmapfile_OS_NAME,     &
         & interp_order_OS,      &         
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

    interp_order_AO = 2
    interp_order_OA = 2

    interp_order_AS = 2
    interp_order_SA = 1
    interp_order_OS = 1
    interp_order_SO = 1
    
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
    call MessageNotify( 'M', PROGRAM_NAME, "gmapfile_AS      =%a        ", ca=(/ gmapfile_AS_NAME /) )
    call MessageNotify( 'M', PROGRAM_NAME, "gmapfile_SA      =%a        ", ca=(/ gmapfile_SA_NAME /) )
    call MessageNotify( 'M', PROGRAM_NAME, "gmapfile_OS      =%a        ", ca=(/ gmapfile_OS_NAME /) )
    call MessageNotify( 'M', PROGRAM_NAME, "gmapfile_SO      =%a        ", ca=(/ gmapfile_SO_NAME /) )
    call MessageNotify( 'M', PROGRAM_NAME, "conservativeFlag =%b        ", L=(/ ConservativeFlag /) )
    call MessageNotify( 'M', PROGRAM_NAME, "interp_order: (AO,OA,AS,SA,OS,SO)=(%d,%d,%d,%d,%d,%d)", &
         & i=(/ interp_order_AO, interp_order_OA,                                                   &
         & interp_order_AS, interp_order_SA, interp_order_OS, interp_order_SO /) )
    
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

  subroutine generate_surface_exchage_grid()

    integer :: j
    real(DP) :: y_FJA(0:JMA)
    real(DP) :: y_FJO(0:JMO)
    real(DP), allocatable :: y_FJS(:)
    real(DP), allocatable :: y_LatTmp(:)
    real(DP), allocatable :: y_IntWtLatTmp(:)
    integer :: JMS_
    real(DP) :: intWt
    
    real(DP), parameter :: PI = acos(-1d0)

    IMS = IMA
    allocate( x_LonS(IMS), x_IntWtLonS(IMS) )    
    x_LonS(:) = x_LonA;
    x_IntWtLonS(:) = x_IntWtLonA(:)

    if (JMA == JMO) then
       JMS = JMO
       allocate( y_LatS(JMS), y_IntWtLatS(JMS) )
       
       y_LatS(:) = y_LatA
       y_IntWtLatS(:) = y_IntWtLatA
    else
       JMS_ = JMO + JMA - 1
       allocate( y_FJS(0:JMS_), y_LatTmp(JMS_), y_IntWtLatTmp(JMS_) )
       
       y_FJA(0) = - PI/2d0
       do j=1, JMA
          y_FJA(j) = asin( y_IntWtLatA(j) + sin(y_FJA(j-1)) )
       end do
       y_FJA(JMA) = PI/2d0

       y_FJO(0) = - PI/2d0
       do j=1, JMO-1
          y_FJO(j) = asin( y_IntWtLatO(j) + sin(y_FJO(j-1)) )
       end do
       y_FJO(JMO) = PI/2d0
    
       y_FJS(0:JMA-1) = y_FJA(0:JMA-1)
       y_FJS(JMA:JMS_) = y_FJO(1:JMO)
       call sort(y_FJS)
!!$       write(*,*) "FJO:", y_FJO/PI*180d0
!!$       write(*,*) "FJA:", y_FJA/PI*180d0
!!$       write(*,*) "FJS:", y_FJS(0:JMS)/PI*180d0

       JMS = 0
       do j=1, JMS_
          intWt = sin(y_FJS(j)) - sin(y_FJS(j-1))
          if (abs(intWt) > 1d-12) then
             JMS = JMS + 1
             y_LatTmp(JMS) = 0.5d0*(y_FJS(j-1) + y_FJS(j))
             y_IntWtLatTmp(JMS) = intWt
          end if
       end do

       allocate( y_LatS(JMS), y_IntWtLatS(JMS) )
       y_LatS(:) = y_LatTmp(1:JMS)
       y_IntWtLatS(:) = y_IntWtLatTmp(1:JMS)
    end if


    write(*,*) "--Surface exchange grid--"
    write(*,*) "y:", y_LatS(:)*180d0/PI
    write(*,*) x_IntWtLonS(:)
    write(*,*) y_IntWtLatS(:)
    write(*,*) sum(y_IntWtLatS)
    
  end subroutine generate_surface_exchage_grid

  subroutine sort(array) 
    real(DP), intent(inout) :: array(:)
    
    integer :: i
    integer :: j
    integer :: N
    real(DP) :: tmp
    
    N = size(array)
    do i=1, N-1
       do j=i+1,N
          if (array(i) > array(j)) then
             tmp = array(i)
             array(i) = array(j)
             array(j) = tmp
          end if
       end do
    end do
    
  end subroutine sort
  
end program gmapgen_main

