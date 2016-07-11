!-------------------------------------------------------------
! Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module mod_common_params

  implicit none
  public

  !
  !
  character(*), parameter :: DEFAULT_DCCM_CONFNAME = "DCCM.conf"
  integer, parameter :: NUM_DCCM_COMP   = 2

  ! The parameters for jcup
  !
  integer, parameter :: JCUP_LOG_LEVEL    = 0
  logical, parameter :: JCUP_LOG_STDERROR = .false.

  integer, parameter :: GN25              = 25
  
  ! The name of each component
  !
  character(*), parameter :: COMPNAME_ATM   = 'ATM'
  integer, parameter      :: COMPID_ATM     = 1

  character(*), parameter :: COMPNAME_OCN   = 'OCN'
  integer, parameter      :: COMPID_OCN     = 2

  character(*), parameter :: COMPNAME_SICE  = 'SICE'
  integer, parameter      :: COMPID_SICE    = 3
  
  character(*), parameter :: COMPNAME_LAND  = 'LAND'
  integer, parameter      :: COMPID_LAND    = 4


  !* The name of grid with each commponent
  !
  
  character(*), parameter :: ATM_GRID_2D = "atm_grid_2d"
!!$  character(*), parameter :: ATM_GRID_3D = "atm_grid_3d"
  integer, parameter :: NUM_ATM_GMAPTAG  = 1 ! 2

  character(*), parameter :: LAND_GRID   = "land_grid"

  character(*), parameter :: OCN_GRID_2D  = "ocn_grid_2d"
  integer, parameter :: NUM_OCN_GMAPTAG  = 1 ! 2

  character(*), parameter :: SICE_GRID   = "sice_grid"

  character(*), parameter :: CHM_GRID    = "chm_grid"  

  integer, parameter :: GMAPTAG_ATM2D_OCN2D = 1

  integer, parameter :: ATM_NUM_GRIDTYPE  = 1
  integer, parameter :: OCN_NUM_GRIDTYPE  = 1
  
  !* The name of variable exchanged between components
  !
  
  !- Atmosphere -> Ocean ********************************************

  character(*), parameter :: a2d_WindStressX = "a2d_WindStressX"
  integer, parameter :: a2d_WindStressX_id   = 1

  character(*), parameter :: a2d_WindStressY = "a2d_WindStressY"
  integer, parameter :: a2d_WindStressY_id   = 2

  character(*), parameter :: a2d_LDwRFlx     = "a2d_LDwRFlx"
  integer, parameter :: a2d_LDwRFlx_id       = 3

  character(*), parameter :: a2d_SDwRFlx     = "a2d_SDwRFlx"
  integer, parameter :: a2d_SDwRFlx_id       = 4

  character(*), parameter :: a2d_LUwRFlx     = "a2d_LUwRFlx"
  integer, parameter :: a2d_LUwRFlx_id       = 5

  character(*), parameter :: a2d_SUwRFlx     = "a2d_SUwRFlx"
  integer, parameter :: a2d_SUwRFlx_id       = 6
  
  character(*), parameter :: a2d_LatHFlx     = "a2d_LatHFlx"
  integer, parameter :: a2d_LatHFlx_id       = 7
  
  character(*), parameter :: a2d_SenHFlx     = "a2d_SenHFlx"
  integer, parameter :: a2d_SenHFlx_id       = 8
  
  character(*), parameter :: a2d_DSfcHFlxDTs = "a2d_DSfcHFlxDTs"
  integer, parameter :: a2d_DSfcHFlxDTs_id   = 9
  
  character(*), parameter :: a2d_RainFall    = "a2d_RainFall"
  integer, parameter :: a2d_RainFall_id      = 10
  
  character(*), parameter :: a2d_SnowFall    = "a2d_SnowFall"
  integer, parameter :: a2d_SnowFall_id      = 11

  character(*), parameter :: a2d_SfcAirTemp  = "a2d_SfcAirTemp"
  integer, parameter :: a2d_SfcAirTemp_id    = 12

  integer, parameter :: NUM_VAR2D_A2O         = 12  

  ! Ocean  -> Atmosphere ***************************************
  
  character(*), parameter :: o2d_SfcTemp     = "o2d_SfcTemp"
  integer, parameter :: o2d_SfcTemp_id       = 1

  character(*), parameter :: o2d_SfcAlbedo   = "o2d_SfcAlbedo"
  integer, parameter :: o2d_SfcAlbedo_id     = 2

  ! Note: We should send the field of snow depth on sea ice form sea-ice component natuarally. 
  character(*), parameter :: o2d_SfcSnow     = "o2d_SfcSnow" 
  integer, parameter :: o2d_SfcSnow_id       = 3

  integer, parameter :: NUM_VAR2D_O2A         = 3
  
  !* Ocean -> Sea-ice *******************************************
  
  !* Sea-ice -> Atmosphere **************************************
  
  character(*), parameter :: s2d_SfcTemp      = "SfcTemp_S2A"
  character(*), parameter :: s2d_SfcAlbedo    = "SfcAlbedo_S2A"

  integer, parameter :: NUM_VAR2D_S2A         = 0
  
  ! Sea-ice -> Ocean


  !-----------------------------------------------------------------
  !-----------------------------------------------------------------


  !*  Ocean, Sea-ice -> Atmosphere
  
  integer, parameter :: NUM_ATM_PUTVAR2D      = NUM_VAR2D_A2O ! + NUM_VAR2D_A2S
  integer, parameter :: NUM_ATM_GETVAR2D      = NUM_VAR2D_O2A ! + NUM_VAR2D_S2A

  integer :: o2a_SfcTemp_id                   = 1
  integer :: o2a_SfcAlbedo_id                 = 2
  integer :: o2a_SfcSnow_id                   = 3
  
  
  !* Atmosphere -> Ocean(, Sea-ice)
  integer, parameter :: NUM_OCN_PUTVAR2D      = NUM_VAR2D_O2A ! + NUM_VAR2D_O2S
  integer, parameter :: NUM_OCN_GETVAR2D      = NUM_VAR2D_A2O ! + NUM_VAR2D_S2O

  integer, parameter :: a2o_WindStressX_id   = 1
  integer, parameter :: a2o_WindStressY_id   = 2
  integer, parameter :: a2o_LDwRFlx_id       = 3
  integer, parameter :: a2o_SDwRFlx_id       = 4
  integer, parameter :: a2o_LUwRFlx_id       = 5
  integer, parameter :: a2o_SUwRFlx_id       = 6
  integer, parameter :: a2o_LatHFlx_id       = 7
  integer, parameter :: a2o_SenHFlx_id       = 8
  integer, parameter :: a2o_DSfcHFlxDTs_id   = 9
  integer, parameter :: a2o_RainFall_id      = 10
  integer, parameter :: a2o_SnowFall_id      = 11
  integer, parameter :: a2o_SfcAirTemp_id    = 12

end module mod_common_params

