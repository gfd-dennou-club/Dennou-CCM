!-------------------------------------------------------------
! Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DCCM_Common_Params_mod

  implicit none
  public

  !
  !
  character(*), parameter :: DEFAULT_DCCM_CONFNAME = "DCCM.conf"
  integer, parameter :: NUM_DCCM_COMP   = 3

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

  character(*), parameter :: COMPNAME_SFC   = 'SFC'
  integer, parameter      :: COMPID_SFC     = 3
  
  character(*), parameter :: COMPNAME_SICE  = 'SICE'
  integer, parameter      :: COMPID_SICE    = 4
  
  character(*), parameter :: COMPNAME_LAND  = 'LAND'
  integer, parameter      :: COMPID_LAND    = 5

  

  !* The name of grid with each commponent
  !
  
  character(*), parameter :: ATM_GRID_2D = "atm_grid_2d"
!!$  character(*), parameter :: ATM_GRID_3D = "atm_grid_3d"
  integer, parameter :: NUM_ATM_GMAPTAG  = 2

  character(*), parameter :: LAND_GRID   = "land_grid"

  character(*), parameter :: OCN_GRID_2D  = "ocn_grid_2d"
  integer, parameter :: NUM_OCN_GMAPTAG   = 2

  character(*), parameter :: SFC_GRID_2D  = "sfc_grid_2d"
  integer, parameter :: NUM_SFC_GMAPTAG   = 4
  
  character(*), parameter :: SICE_GRID   = "sice_grid"

  character(*), parameter :: CHM_GRID    = "chm_grid"  

  integer, parameter :: GMAPTAG_ATM2D_OCN2D          = 1
  integer, parameter :: GMAPTAG_ATM2D_OCN2D_CONSERVE = 2
  integer, parameter :: GMAPTAG_ATM2D_SFC2D          = 1
  integer, parameter :: GMAPTAG_ATM2D_SFC2D_CONSERVE = 2
  integer, parameter :: GMAPTAG_SFC2D_OCN2D          = 1
  integer, parameter :: GMAPTAG_SFC2D_OCN2D_CONSERVE = 2

  integer, parameter :: ATM_NUM_GRIDTYPE  = 1
  integer, parameter :: OCN_NUM_GRIDTYPE  = 1
  integer, parameter :: SFC_NUM_GRIDTYPE  = 1
  
  !* The name of variable exchanged between components
  !

  !- Atmosphere -> Surface *********************************************

  character(*), parameter :: a2s_WindU       = "a2s_WindU"
  integer, parameter :: a2s_WindU_id         = 1

  character(*), parameter :: a2s_WindV       = "a2s_WindV"
  integer, parameter :: a2s_WindV_id         = 2

  character(*), parameter :: a2s_SfcAirTemp  = "a2s_SfcAirTemp"
  integer, parameter :: a2s_SfcAirTemp_id    = 3

  character(*), parameter :: a2s_SfcPress    = "a2s_SfcPress"
  integer, parameter :: a2s_SfcPress_id      = 4
  
  character(*), parameter :: a2s_QVap1       = "a2s_Qvap1"
  integer, parameter :: a2s_QVap1_id         = 5
  
  character(*), parameter :: a2s_LDwRFlx     = "a2s_LDwRFlx"
  integer, parameter :: a2s_LDwRFlx_id       = 6

  character(*), parameter :: a2s_SDwRFlx     = "a2s_SDwRFlx"
  integer, parameter :: a2s_SDwRFlx_id       = 7

  character(*), parameter :: a2s_RainFall    = "a2s_RainFall"
  integer, parameter :: a2s_RainFall_id      = 8
  
  character(*), parameter :: a2s_SnowFall    = "a2s_SnowFall"
  integer, parameter :: a2s_SnowFall_id      = 9

  character(*), parameter :: a2s_ImplCPLCoef1 = "a2s_ImplCPLCoef1"
  integer, parameter :: a2s_ImplCPLCoef1_id    = 10

  character(*), parameter :: a2s_ImplCPLCoef2 = "a2s_ImplCPLCoef2"
  integer, parameter :: a2s_ImplCPLCoef2_id   = 11

!!$  character(*), parameter :: a2s_Press1      = "a2s_Press1"
!!$  integer, parameter :: a2s_Press1_id        = 12
  
  integer, parameter :: NUM_VAR2D_a2s        = 11
  
  !- Surface -> Atmosphere
 
  character(*), parameter :: s2a_LUwRFlx     = "s2a_LUwRFlx"
  integer, parameter :: s2a_LUwRFlx_id       = 1

  character(*), parameter :: s2a_SUwRFlx     = "s2a_SUwRFlx"
  integer, parameter :: s2a_SUwRFlx_id       = 2
  
  character(*), parameter :: s2a_SenHFlx     = "s2a_SenHFlx"
  integer, parameter :: s2a_SenHFlx_id       = 3

  character(*), parameter :: s2a_QVapMFlx    = "s2a_QVapMFlx"
  integer, parameter :: s2a_QVapMFlx_id      = 4

  character(*), parameter :: s2a_SfcAlbedo   = "s2a_SfcAlbedo"
  integer, parameter :: s2a_SfcAlbedo_id     = 5

  character(*), parameter :: s2a_DelVarImplCPL = "s2a_DelVarImplCPL"
  integer, parameter :: s2a_DelVarImplCPL_id   = 6

!!$  character(*), parameter :: s2a_WindStressX = "s2a_WindStressX"
!!$  integer, parameter :: s2a_WindStressX_id   = 7
!!$
!!$  character(*), parameter :: s2a_WindStressY = "s2a_WindStressY"
!!$  integer, parameter :: s2a_WindStressY_id   = 8  
!!$  character(*), parameter :: s2a_SfcRadTemp   = "s2a_SfcRadTemp"
!!$  integer, parameter :: s2a_SfcRadTemp_id     = 8
  
  integer, parameter :: NUM_VAR2D_s2a          = 6
  
  !- Ocean -> Surface ********************************************

  character(*), parameter :: o2s_SfcTemp      = "o2s_SfcTemp"
  integer, parameter :: o2s_SfcTemp_id        = 1
  character(*), parameter :: o2s_SfcAlbedo    = "o2s_SfcAlbedo"
  integer, parameter :: o2s_SfcAlbedo_id      = 2

  integer, parameter :: NUM_VAR2D_o2s         = 2

  !- Sea-ice -> Surface ********************************************

  character(*), parameter :: i2s_SfcTemp      = "i2s_SfcTemp"
  integer, parameter :: i2s_SfcTemp_id        = 1
  character(*), parameter :: i2s_SfcAlbedo    = "i2s_SfcAlbedo"
  integer, parameter :: i2s_SfcAlbedo_id      = 2
  character(*), parameter :: i2s_SIceCon      = "i2s_SIceCon"
  integer, parameter :: i2s_SIceCon_id        = 3

  integer, parameter :: NUM_VAR2D_i2s         = 3
  
  !- Surface -> Ocean   **************************************

  character(*), parameter :: s2o_SfcHFlx_ns   = "s2o_SfcHFlx_ns"
  integer, parameter :: s2o_SfcHFlx_ns_id     = 1

  character(*), parameter :: s2o_SfcHFlx_sr   = "s2o_SfcHFlx_sr"
  integer, parameter :: s2o_SfcHFlx_sr_id     = 2

  character(*), parameter :: s2o_SnowFall     = "s2o_SnowFall"
  integer, parameter :: s2o_SnowFall_id       = 3

  character(*), parameter :: s2o_RainFall     = "s2o_RainFall"
  integer, parameter :: s2o_RainFall_id       = 4

  character(*), parameter :: s2o_Evap         = "s2o_Evap"
  integer, parameter :: s2o_Evap_id           = 5
  
  character(*), parameter :: s2o_WindStressX  = "s2o_WindStressX"
  integer, parameter :: s2o_WindStressX_id    = 6

  character(*), parameter :: s2o_WindStressY  = "s2o_WindStressY"
  integer, parameter :: s2o_WindStressY_id    = 7

  character(*), parameter :: s2o_DSfcHFlxDTs  = "s2o_DSfcHFlxDTs"
  integer, parameter :: s2o_DSfcHFlxDTs_id    = 8
  
  integer, parameter :: NUM_VAR2D_s2o         = 8

  !- Surface -> Ocean   **************************************

  character(*), parameter :: s2i_SfcHFlx_ns   = "s2i_SfcHFlx_ns"
  integer, parameter :: s2i_SfcHFlx_ns_id     = 1

  character(*), parameter :: s2i_SfcHFlx_sr   = "s2i_SfcHFlx_sr"
  integer, parameter :: s2i_SfcHFlx_sr_id     = 2

  character(*), parameter :: s2i_DSfcHFlxDTs  = "s2i_DSfcHFlxDTs"
  integer, parameter :: s2i_DSfcHFlxDTs_id    = 3

  character(*), parameter :: s2i_Evap         = "s2i_Evap"
  integer, parameter :: s2i_Evap_id           = 4
  
  integer, parameter :: NUM_VAR2D_s2i         = 4
  
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------

end module DCCM_Common_Params_mod
