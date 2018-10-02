module DSFCM_Admin_Variable_mod

  ! モジュール引用; Use statements
  !
 
  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* DSFCM
  
  use DSFCM_Admin_Grid_mod, only: &
       & IA, JA
  
  implicit none
  private


  public :: DSFCM_Admin_Variable_Init
  public :: DSFCM_Admin_Variable_Final
  

  integer, public, parameter :: SFC_PROP_MAX = 3
  
  real(DP), public, allocatable :: xya_WindStressX(:,:,:)
  real(DP), public, allocatable :: xya_WindStressY(:,:,:)
  real(DP), public, allocatable :: xya_LatHFlx(:,:,:)
  real(DP), public, allocatable :: xya_SenHFlx(:,:,:)
  real(DP), public, allocatable :: xya_QvapMFlx(:,:,:)
  real(DP), public, allocatable :: xya_SfcAlbedo(:,:,:)
  real(DP), public, allocatable :: xya_SfcTemp(:,:,:)
  real(DP), public, allocatable :: xy_SDwRFlx(:,:)
  real(DP), public, allocatable :: xy_LDwRFlx(:,:)
  real(DP), public, allocatable :: xya_SUwRFlx(:,:,:)
  real(DP), public, allocatable :: xya_LUwRFlx(:,:,:)
  real(DP), public, allocatable :: xya_DelVarImplCPL(:,:,:)
  
  real(DP), public, allocatable :: xy_SIceCon(:,:)

  real(DP), public, allocatable :: xy_RainFall(:,:)
  real(DP), public, allocatable :: xy_SnowFall(:,:)

  real(DP), public, allocatable :: xya_SfcVelTransCoef(:,:,:)
  real(DP), public, allocatable :: xya_SfcTempTransCoef(:,:,:)
  real(DP), public, allocatable :: xya_SfcQVapTransCoef(:,:,:)

  real(DP), public, allocatable :: xya_SfcHFlx_ns(:,:,:)
  real(DP), public, allocatable :: xya_SfcHFlx_sr(:,:,:)
  real(DP), public, allocatable :: xya_DSfcHFlxDTs(:,:,:)
  
contains

  subroutine DSFCM_Admin_Variable_Init()

    allocate( xya_WindStressX(IA,JA,SFC_PROP_MAX) )
    allocate( xya_WindStressY(IA,JA,SFC_PROP_MAX) )
    allocate( xya_LatHFlx(IA,JA,SFC_PROP_MAX) )
    allocate( xya_SenHFlx(IA,JA,SFC_PROP_MAX) )
    allocate( xya_QVapMFlx(IA,JA,SFC_PROP_MAX) )
    allocate( xya_SfcTemp(IA,JA,SFC_PROP_MAX) )
    allocate( xya_SfcAlbedo(IA,JA,SFC_PROP_MAX) )
    allocate( xy_SIceCon(IA,JA) )

    allocate( xya_SfcVelTransCoef(IA,JA,SFC_PROP_MAX) )
    allocate( xya_SfcTempTransCoef(IA,JA,SFC_PROP_MAX) )
    allocate( xya_SfcQVapTransCoef(IA,JA,SFC_PROP_MAX) )
    
    allocate( xy_RainFall(IA,JA), xy_SnowFall(IA,JA) )
    allocate( xy_SDwRFlx(IA,JA), xy_LDwRFlx(IA,JA) )
    allocate( xya_SUwRFlx(IA,JA,SFC_PROP_MAX), xya_LUwRFlx(IA,JA,SFC_PROP_MAX) )
    
    allocate( xya_SfcHFlx_ns(IA,JA,SFC_PROP_MAX) )
    allocate( xya_SfcHFlx_sr(IA,JA,SFC_PROP_MAX) )
    allocate( xya_DSfcHFlxDTs(IA,JA,SFC_PROP_MAX) )

    allocate( xya_DelVarImplCPL(IA,JA,4) )
    
  end subroutine DSFCM_Admin_Variable_Init

  subroutine DSFCM_Admin_Variable_Final()

    deallocate( xya_WindStressX, xya_WindStressY )
    deallocate( xya_LatHFlx, xya_SenHFlx )
    deallocate( xya_QVapMFlx )
    deallocate( xya_SfcTemp, xya_SfcAlbedo )
    
    deallocate(xy_SIceCon )

    deallocate( xya_SfcVelTransCoef, xya_SfcTempTransCoef )
    deallocate( xya_SfcQVapTransCoef )
    
    deallocate( xy_SDwRFlx, xy_LDwRFlx, xya_SUwRFlx, xya_LUwRFlx )
    deallocate( xy_RainFall, xy_SnowFall )

    deallocate( xya_SfcHFlx_ns, xya_SfcHFlx_sr, xya_DSfcHFlxDTs )

    deallocate( xya_DelVarImplCPL )
    
  end subroutine DSFCM_Admin_Variable_Final
  
end module DSFCM_Admin_Variable_mod
