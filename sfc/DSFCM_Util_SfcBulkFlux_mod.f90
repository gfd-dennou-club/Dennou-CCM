module DSFCM_Util_SfcBulkFlux_mod

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
  
  use DSFCM_Admin_Grid_mod, only: &
       & IA, IS, IE, &
       & JA, JS, JE

  use DSFCM_Admin_Variable_mod, only: &
       & SFC_PROP_MAX

  implicit none
  private

  public :: DSFCM_Util_SfcBulkFlux_Init
  public :: DSFCM_Util_SfcBulkFlux_Final
  public :: DSFCM_Util_SfcBulkFlux_Get  

  real(DP), parameter :: FKarm = 0.4_DP ! Karman constant
  real(DP), parameter, public:: GasRUniv = 8.3144621_DP
  real(DP), parameter, public:: StB = 5.670373e-8_DP
  real(DP), parameter, public:: PI = 3.1415926535897932_DP
  
  real(DP), parameter :: Grav  = 9.8_DP
  real(DP), parameter :: MolWtDry =  1.8D-2
  real(DP), parameter :: MolWtWet =  1.8D-2
  real(DP), parameter :: CpDry    = 1616D0
  real(DP), parameter :: CpWet    = 1616D0
  real(DP), parameter :: LatentHeat = 2425300D0
  real(DP), parameter :: LatentHeatFusion = 334000D0
  real(DP), parameter :: RefPress = 1d5
  
  real(DP), parameter :: GasRDry = GasRUniv / MolWtDry
  real(DP), parameter :: GasRWet = GasRUniv / MolWtWet
  real(DP), parameter :: EpsV    = MolWtWet/MolWtDry

  real(DP), parameter:: Es0 = 611.0_DP
                              ! Saturation water vapor pressure at 0 deg C [Pa]

  real(DP), parameter :: HumidCoef = 1d0
  real(DP), parameter :: RoughLength = 1d-4
  real(DP), parameter :: RoughLenHeatFactor = 1d0
  
  real(DP), save:: VelMinForRi
                            ! Minimum value of velocity for $ R_i $ number
  real(DP), save:: VelMinForVel
                            ! Minimum value of velocity for momentum
  real(DP), save:: VelMinForTemp
                            ! Minimum value of velocity for thermal
  real(DP), save:: VelMinForQVap
                            ! Minimum value of velocity for vapor
  real(DP), save:: VelMaxForVel
                            ! Maximum value of velocity for momentum
  real(DP), save:: VelMaxForTemp
                            ! Maximum value of velocity for thermal
  real(DP), save:: VelMaxForQVap
                            ! Maximum value of velocity for vapor

  ! Bluk coefficients
  !
  logical, save:: FlagConstBulkCoef
                            ! Flag for using constant bulk coefficient
  logical, save:: FlagUseOfBulkCoefInNeutralCond
                            ! Flag for using bulk coefficient in neutral condition
  real(DP), save:: ConstBulkCoef
                            ! Steady value of bulk coefficient
  real(DP), save:: VelBulkCoefMin
                            ! Minimum value of $ u $ bulk coefficient
  real(DP), save:: TempBulkCoefMin
                            ! Minimum value of $ T $ bulk coefficient
  real(DP), save:: QVapBulkCoefMin
                            ! Minimum value of $ q $ bulk coefficient
  real(DP), save:: VelBulkCoefMax
                            ! Maximum value of $ u $ bulk coefficient
  real(DP), save:: TempBulkCoefMax
                            ! Maximum value of $ T $ bulk coefficient
  real(DP), save:: QVapBulkCoefMax
                            ! Maximum value of $ q $ bulk coefficient

  character(*), parameter:: module_name = 'DSFCM_Util_SfcBulkFlux_mod'
  
contains
  subroutine DSFCM_Util_SfcBulkFlux_Init()

    call read_config()
    
  end subroutine DSFCM_Util_SfcBulkFlux_Init

  subroutine DSFCM_Util_SfcBulkFlux_Final()
  end subroutine DSFCM_Util_SfcBulkFlux_Final

  subroutine DSFCM_Util_SfcBulkFlux_Get( &
       & xya_WindStressX, xya_WindStressY, &
       & xya_SenHFlx, xya_QVapMFlx, xya_LatHFlx,  &
       & xya_SfcVelTransCoef, xya_SfcTempTransCoef, xya_SfcQVapTransCoef, & ! (out)
       & xya_DelVarImplCPL,                              &
       & xya_SUwRFlx, xya_LUwRFlx,                      &
       & xya_SfcHFlx_ns, xya_SfcHFlx_sr, xya_DSfcHFlxDTs, &
       & xy_WindU, xy_WindV, xy_SfcAirTemp, xy_QVap1,  &
       & xy_SDwRFlx, xy_LDwRFlx, xya_ImplCplCoef1, xya_ImplCplCoef2,                       &
       & xya_SfcTemp, xya_SfcAlbedo,  xy_SIceCon,      &
       & a_Sig1Info, xy_SfcHeight,                     &
       & xy_SfcPress                                   & 
       & )

    implicit none
    
    real(DP), intent(out) :: xya_WindStressX(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_WindStressY(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_SenHFlx(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_QVapMFlx(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_LatHFlx(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_SfcVelTransCoef(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_SfcTempTransCoef(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_SfcQVapTransCoef(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_SfcHFlx_ns(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_SfcHFlx_sr(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_DSfcHFlxDTs(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_DelVarImplCPL(IA,JA,4)
    real(DP), intent(out) :: xya_SUwRFlx(IA,JA,SFC_PROP_MAX)
    real(DP), intent(out) :: xya_LUwRFlx(IA,JA,SFC_PROP_MAX)    
    real(DP), intent(in) :: xy_WindU(IA,JA)
    real(DP), intent(in) :: xy_WindV(IA,JA)
    real(DP), intent(in) :: xy_SfcAirTemp(IA,JA)
    real(DP), intent(in) :: xy_QVap1(IA,JA)
    real(DP), intent(in) :: xy_SDwRFlx(IA,JA)
    real(DP), intent(in) :: xy_LDwRFlx(IA,JA)
    real(DP), intent(in) :: xya_ImplCplCoef1(IA,JA,4)
    real(DP), intent(in) :: xya_ImplCplCoef2(IA,JA,4)
    real(DP), intent(inout) :: xya_SfcTemp(IA,JA,SFC_PROP_MAX)
    real(DP), intent(inout) :: xya_SfcAlbedo(IA,JA,SFC_PROP_MAX)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)    
    real(DP), intent(in) :: a_Sig1Info(2)
    real(DP), intent(in) :: xy_SfcHeight(IA,JA)
    real(DP), intent(in) :: xy_SfcPress(IA,JA)
    
    integer :: i
    integer :: j
    integer :: n
    
    real(DP) :: xya_SfcMOLength(IA,JA,2)
    
    real(DP) :: xya_SfcRoughLengthMom(IA,JA,2)
    real(DP) :: xya_SfcRoughLengthHeat(IA,JA,2)
    real(DP) :: xya_SfcHumdCoef(IA,JA,SFC_PROP_MAX)
    real(DP) :: xya_SfcBulkRiNum(IA,JA,2)
    real(DP) :: xy_SfcBulkCoefMomInNeutCond(IA,JA)
    real(DP) :: xy_SfcBulkCoefHeatInNeutCond(IA,JA)
    
    real(DP) :: xya_SfcVelBulkCoef (IA,JA,2)
    real(DP) :: xya_SfcTempBulkCoef (IA,JA,2)
    real(DP) :: xya_SfcQVapBulkCoef (IA,JA,2)
    real(DP) :: xy_BetaW(IA,JA)

    real(DP) :: xy_SfcVelAbs(IA,JA)
    real(DP) :: SfcBulkCoefMomInNeutCondTmp

    real(DP) :: xy_VirTemp(IA,JA)
    real(DP) :: xya_SfcVirTemp(IA,JA,2)
    real(DP) :: xya_SfcQVapSat(IA,JA,2)
    real(DP) :: xy_Exner(IA,JA)
    real(DP) :: xy_SfcExner(IA,JA)

    real(DP) :: xya_Frac(IA,JA,SFC_PROP_MAX-1)
    real(DP) :: a_LatentHeatLocal(2)

    real(DP) :: xy_Height(IA,JA)
    real(DP) :: xy_Press1(IA,JA)

    logical :: xy_CalcFlag(IA,JA)
    
    integer :: SPMAX

    real(DP) :: a_f(4)
    real(DP) :: a_Gamma(4)
    real(DP) :: DFsDT1
    
    SPMAX = SFC_PROP_MAX
    
    xya_SfcRoughLengthMom  = RoughLength
    xya_SfcRoughLengthHeat = RoughLenHeatFactor*xya_SfcRoughLengthMom
    xya_SfcHumdCoef        = 1d0

    a_LatentHeatLocal(1) = LatentHeat
    a_LatentHeatLocal(2) = LatentHeat + LatentHeatFusion
    
    !$omp parallel private(n)
    !$omp do collapse(2)
    do j=JS, JE
    do i=IS, IE
       xya_Frac(i,j,1) = 1d0 - xy_SIceCon(i,j)
       xya_Frac(i,j,2) = xy_SIceCon(i,j)
       
       do n=1, SPMAX-1
          xya_SfcQVapSat(i,j,n) = EpsV * Es0 / xy_SfcPress(i,j)                                               &
               & * exp( a_LatentHeatLocal(n)  / GasRWet * (1d0/273d0 - 1d0/xya_SfcTemp(i,j,n)) ) 
          xya_SfcVirTemp(i,j,n) = xya_SfcTemp(i,j,n) * ( 1d0 + ((( 1d0/EpsV ) - 1d0) * xya_SfcQVapSat(i,j,n)) )
       end do
       xy_VirTemp(i,j) = xy_SfcAirTemp(i,j) * ( 1d0 + ((( 1d0/EpsV ) - 1d0) * xy_QVap1(i,j)) )

       xy_Press1(i,j) = xy_SfcPress(i,j) * a_Sig1Info(1)
       xy_Exner(i,j) = (xy_Press1(i,j)/RefPress)**( GasRDry/CpDry )       
       xy_SfcExner(i,j) = (xy_SfcPress(i,j)/RefPress)**( GasRDry/CpDry )

       xy_SfcVelAbs(i,j) = sqrt ( xy_WindU(i,j)**2 + xy_WindV(i,j)**2 )

       xy_Height(i,j) = xy_SfcHeight(i,j) + &
            & GasRDry/Grav * xy_VirTemp(i,j) * (1d0 - a_Sig1Info(1))

       xya_WindStressX(i,j,SPMAX) = 0d0
       xya_WindStressY(i,j,SPMAX) = 0d0
       xya_SenHFlx(i,j,SPMAX) = 0d0
       xya_LatHFlx(i,j,SPMAX) = 0d0
       xya_QvapMFlx(i,j,SPMAX) = 0d0
       xya_SUwRFlx(i,j,SPMAX) = 0d0
       xya_LUwRFlx(i,j,SPMAX) = 0d0
       
       xya_SfcTemp(i,j,SPMAX) = 0d0
       xya_SfcAlbedo(i,j,SPMAX) = 0d0

       xya_SfcVelTransCoef(i,j,SPMAX) = 0d0
       xya_SfcTempTransCoef(i,j,SPMAX) = 0d0
       xya_SfcQVapTransCoef(i,j,SPMAX) = 0d0       
    end do
    end do
    !$omp end parallel

    do n=1, SPMAX-1

       
       !$omp parallel do collapse(2) private(SfcBulkCoefMomInNeutCondTmp)
       do j=JS, JE
       do i=IS, IE
          SfcBulkCoefMomInNeutCondTmp = &
               & ( FKarm                                         &
               & / log (   ( xy_Height(i,j) - xy_SfcHeight(i,j) + xya_SfcRoughLengthMom(i,j,n) ) &
               &         / xya_SfcRoughLengthMom(i,j,n)  ) )
       
          xy_SfcBulkCoefMomInNeutCond(i,j)  = SfcBulkCoefMomInNeutCondTmp**2
          xy_SfcBulkCoefHeatInNeutCond(i,j)  =                   &
               & SfcBulkCoefMomInNeutCondTmp                     &
               & * ( FKarm                                       &
               & / log (   ( xy_Height(i,j) - xy_SfcHeight(i,j) + xya_SfcRoughLengthHeat(i,j,n) ) &
               &         / xya_SfcRoughLengthHeat(i,j,n) ) )

          xya_SfcBulkRiNum(i,j,n) =                                          &
               &   Grav / ( xya_SfcVirTemp(i,j,n) / xy_SfcExner(i,j) )       &
               &   * (   xy_VirTemp(i,j)     / xy_Exner(i,j)                 &
               &       - xya_SfcVirTemp(i,j,n) / xy_SfcExner(i,j) )          &
               &   / max( xy_SfcVelAbs(i,j), VelMinForRi )**2                &
               &   * ( xy_Height(i,j) - xy_SfcHeight(i,j) )

          if (n == 1) then ! Ocean 
             xy_CalcFlag(i,j) = .true.
          else if(n == 2) then ! Sea-ice
             xy_CalcFlag(i,j) = (xya_Frac(i,j,n) > 1d-12)
          end if          

       end do
       end do
  
       call BulkCoefL82( &
            & xya_SfcBulkRiNum(:,:,n),  xya_SfcRoughLengthMom(:,:,n) , xya_SfcRoughLengthHeat(:,:,n),  & ! (in)
            & xy_SfcHeight, xy_Height,                                                                 & ! (in)
            & xy_SfcBulkCoefMomInNeutCond, xy_SfcBulkCoefHeatInNeutCond, xy_CalcFlag,                  & ! (in)
            & xya_SfcVelBulkCoef(:,:,n), xya_SfcTempBulkCoef(:,:,n), xya_SfcQVapBulkCoef(:,:,n),       & ! (out)
            & xy_BetaW, xya_SfcMOLength(:,:,n)                                                         & ! (out)
            & )

       !$omp parallel do collapse(2)
       do j=JS, JE
       do i=IS, IE

          !--
          
          xya_SfcVelTransCoef(i,j,n) = &
               & xya_SfcVelBulkCoef(i,j,n)*xy_SfcPress(i,j)/(GasRDry * xya_SfcVirTemp(i,j,n)) &
               & * min( max( xy_SfcVelAbs(i,j), VelMinForVel), VelMaxForVel)
          
          xya_SfcTempTransCoef(i,j,n) = &
               & xya_SfcTempBulkCoef(i,j,n)*xy_SfcPress(i,j)/(GasRDry * xya_SfcVirTemp(i,j,n)) &
               & * min( max( xy_SfcVelAbs(i,j), VelMinForTemp), VelMaxForTemp)
          
          xya_SfcQVapTransCoef(i,j,n) = &
               & xya_SfcQVapBulkCoef(i,j,n)*xy_SfcPress(i,j)/(GasRDry * xya_SfcVirTemp(i,j,n)) &
               & * min( max( xy_SfcVelAbs(i,j), VelMinForQVap), VelMaxForQVap)

          !--
          if ( xy_CalcFlag(i,j) ) then
             xya_WindStressX(i,j,n) = - xya_SfcVelTransCoef(i,j,n) * xy_WindU(i,j)
             xya_WindStressY(i,j,n) = - xya_SfcVelTransCoef(i,j,n) * xy_WindV(i,j)

             xya_SenHFlx(i,j,n) = - CpDry * xy_SfcExner(i,j) * xya_SfcTempTransCoef(i,j,n) &
                  & * (xy_SfcAirTemp(i,j)/xy_Exner(i,j) - xya_SfcTemp(i,j,n)/xy_SfcExner(i,j))
             
             xya_QVapMFlx(i,j,n) = - xya_SfcHumdCoef(i,j,n) * xya_SfcQVapTransCoef(i,j,n) &
                  & * (xy_QVap1(i,j) - xya_SfcQVapSat(i,j,n))
          
             xya_LatHFlx(i,j,n) = a_LatentHeatLocal(n)*xya_QVapMFlx(i,j,n)

             xya_LUwRFlx(i,j,n) = StB*xya_SfcTemp(i,j,n)**4
             xya_SUwRFlx(i,j,n) = xya_SfcAlbedo(i,j,n)*xy_SDwRFlx(i,j)
                          
             !---          
             xya_SfcTemp(i,j,SPMAX) = xya_SfcTemp(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SfcTemp(i,j,n)**4
             xya_SfcAlbedo(i,j,SPMAX) = xya_SfcAlbedo(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SfcAlbedo(i,j,n)

             xya_WindStressX(i,j,SPMAX) = xya_WindStressX(i,j,SPMAX) + xya_Frac(i,j,n)*xya_WindStressX(i,j,n)
             xya_WindStressY(i,j,SPMAX) = xya_WindStressY(i,j,SPMAX) + xya_Frac(i,j,n)*xya_WindStressY(i,j,n)
             xya_SenHFlx(i,j,SPMAX) = xya_SenHFlx(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SenHFlx(i,j,n)          
             xya_QVapMFlx(i,j,SPMAX) = xya_QVapMFlx(i,j,SPMAX) + xya_Frac(i,j,n)*xya_QVapMFlx(i,j,n)
             xya_LatHFlx(i,j,SPMAX) = xya_LatHFlx(i,j,SPMAX) + xya_Frac(i,j,n)*xya_LatHFlx(i,j,n)
             xya_LUwRFlx(i,j,SPMAX) = xya_LUwRFlx(i,j,SPMAX) + xya_Frac(i,j,n)*xya_LUwRFlx(i,j,n)
             xya_SUwRFlx(i,j,SPMAX) = xya_SUwRFlx(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SUwRFlx(i,j,n)

             xya_SfcVelTransCoef(i,j,SPMAX) = xya_SfcVelTransCoef(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SfcVelTransCoef(i,j,n)
             xya_SfcTempTransCoef(i,j,SPMAX) = xya_SfcTempTransCoef(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SfcTempTransCoef(i,j,n)
             xya_SfcQVapTransCoef(i,j,SPMAX) = xya_SfcQVapTransCoef(i,j,SPMAX) + xya_Frac(i,j,n)*xya_SfcQVapTransCoef(i,j,n)
             
          else
             xya_WindStressX(i,j,n) = 0d0
             xya_WindStressY(i,j,n) = 0d0
             xya_SenHFlx(i,j,n)     = 0d0
             xya_QVapMFlx(i,j,n)    = 0d0
             xya_LatHFlx(i,j,n)     = 0d0
             xya_LUwRFlx(i,j,n)     = 0d0
             xya_SUwRFlx(i,j,n)     = 0d0
          end if
          !------------------------------------------------------------------------
       end do
       end do
    end do

!!$    write(*,*) "WindStressX0:", xya_WindStressX(IS,JS:JE,3)
!!$    write(*,*) "SenHFlx0:", xya_SenHFlx(IS,JS:JE,3)
    
    !$omp parallel private (i, j, n, a_Gamma, a_f, DFsDT1)
    !$omp do collapse(2)
    do j = JS, JE
    do i = IS, IE
       DFsDT1 = - CpDry * xy_SfcExner(i,j) * xya_SfcTempTransCoef(i,j,SPMAX) / xy_Exner(i,j)
       a_Gamma(1) = 1d0 / (xya_ImplCplCoef1(i,j,1) + xya_SfcVelTransCoef(i,j,SPMAX))
       a_Gamma(2) = 1d0 / (xya_ImplCplCoef1(i,j,2) + xya_SfcVelTransCoef(i,j,SPMAX))       
       a_Gamma(3) = 1d0 / (xya_ImplCplCoef1(i,j,3) - DFsDT1)
       a_Gamma(4) = 1d0 / (xya_ImplCplCoef1(i,j,4) + xya_SfcHumdCoef(i,j,SPMAX) * xya_SfcQVapTransCoef(i,j,SPMAX))
 
       a_f(1) = a_Gamma(1)*(xya_WindStressX(i,j,SPMAX) + xya_ImplCplCoef2(i,j,1))
       a_f(2) = a_Gamma(2)*(xya_WindStressY(i,j,SPMAX) + xya_ImplCplCoef2(i,j,2))
       a_f(3) = a_Gamma(3)*(xya_SenHFlx(i,j,SPMAX) + xya_ImplCplCoef2(i,j,3))
       a_f(4) = a_Gamma(4)*(xya_QVapMFlx(i,j,SPMAX) + xya_ImplCplCoef2(i,j,4))

       xya_DelVarImplCPL(i,j,:) = a_f(:)
       
       do n=1, SPMAX
          xya_WindStressX(i,j,n) = xya_WindStressX(i,j,n) - xya_SfcVelTransCoef(i,j,n) * a_f(1)
          xya_WindStressY(i,j,n) = xya_WindStressY(i,j,n) - xya_SfcVelTransCoef(i,j,n) * a_f(2)

          xya_SenHFlx(i,j,n) = xya_SenHFlx(i,j,n) &
               & -  CpDry * xy_SfcExner(i,j) / xy_Exner(i,j) * xya_SfcTempTransCoef(i,j,n) * a_f(3)
          
          xya_QVapMFlx(i,j,n) = xya_QVapMFlx(i,j,n) &
               & -  xya_SfcHumdCoef(i,j,n)*xya_SfcQVapTransCoef(i,j,n) * a_f(4)
          xya_LatHFlx(i,j,n) = a_LatentHeatLocal(n)*xya_QVapMFlx(i,j,n)
       end do       
    end do
    end do

    !$omp do collapse(2)
    do n = 1, SPMAX-1
    do j = JS, JE
    do i = IS, IE

       if (n == 1) then ! Ocean 
          xy_CalcFlag(i,j) = .true.
       else if(n == 2) then ! Sea-ice
          xy_CalcFlag(i,j) = (xya_Frac(i,j,n) > 1d-12)
       end if
       
       if (xy_CalcFlag(i,j)) then
          !--             
          xya_SfcHFlx_ns(i,j,n) = &
               & + xya_LUwRFlx(i,j,n) - xy_LDwRFlx(i,j)      &
               & + xya_LatHFlx(i,j,n) + xya_SenHFlx(i,j,n)
          
          xya_SfcHFlx_sr(i,j,n) = xya_SUwRFlx(i,j,n) - xy_SDwRFlx(i,j)
          
          xya_DSfcHFlxDTs(i,j,n) = &
               & + 4d0*StB*xya_SfcTemp(i,j,n)**3         &
               & + CpDry*xya_SfcTempTransCoef(i,j,n)     &
               & + a_LatentHeatLocal(n)*xya_SfcHumdCoef(i,j,n)*xya_SfcQVapTransCoef(i,j,n)         &
               &     *( a_LatentHeatLocal(n)*xya_SfcQVapSat(i,j,n)/(GasRWet*xya_SfcTemp(i,j,n)**2) )   ! DQsta/DTs
       else
          xya_SfcHFlx_ns(i,j,n)  = 0d0
          xya_SfcHFlx_sr(i,j,n)  = 0d0
          xya_DSfcHFlxDTs(i,j,n) = 0d0
       end if
    end do
    end do
    end do

    !$omp end parallel
    
!!$    write(*,*) "SeaSfcTemp:", xya_SfcTemp(IS,JS:JE,1)
!!$    write(*,*) "LUwRFlx:", xya_LUwRFlx(IS,JS:JE,1)
!!$    write(*,*) "LatHFlx:", xya_LatHFlx(IS,JS:JE,1)
!!$    write(*,*) "SIceSfcTemp:", xya_SfcTemp(IS,JS:JE,2)
!!$    write(*,*) "LUwRFlx:", xya_LUwRFlx(IS,JS:JE,2)
!!$    write(*,*) "LatHFlx:", xya_LatHFlx(IS,JS:JE,2)
!!$    write(*,*) "WindStressX:", xya_WindStressX(IS,JS:JE,2)
!!$    write(*,*) "WindStressX*:", xya_WindStressX(IS,JS:JE,3)
!!$    write(*,*) "SenHFlx*:", xya_SenHFlx(IS,JS:JE,3)
!!$    write(*,*) "DelUImpCPL:", xya_DelVarImplCPL(IS,JS:JE,1)
!!$    write(*,*) "TsO:", xya_SfcTemp(IS,JS:JE,1)
!!$    write(*,*) "T1:", xy_SfcAirTemp(IS,JS:JE)
!!$    write(*,*) "LD:", sum(xy_LDWRFlx(IS:IE,JS:JE),1)/64d0
!!$    write(*,*) "Lat:", sum(xya_LatHFlx(IS:IE,JS:JE,1),1)/64d0
!!$    write(*,*) "Sen:", sum(xya_SenHFlx(IS:IE,JS:JE,1),1)/64d0
!!$    write(*,*) "LU:", sum(xya_LUWRFlx(IS:IE,JS:JE,1),1)/64d0
!!$    write(*,*) "SD:", sum(xy_SDWRFlx(IS:IE,JS:JE),1)/64d0
!!$    write(*,*) "AO:", sum(xya_SfcHFlx_ns(IS:IE,JS:JE,1)+xya_SfcHFlx_sr(IS:IE,JS:JE,1),1)/64d0
!!$    stop
    
  end subroutine DSFCM_Util_SfcBulkFlux_Get

  !-----------------------------------------

  subroutine BulkCoefL82( &
    & xy_SfcBulkRiNum,  xy_SfcRoughLengthMom , xy_SfcRoughLengthHeat,          & ! (in)
    & xy_SfcHeight, xy_Height,                                                 & ! (in)
    & xy_SfcBulkCoefMomInNeutCond, xy_SfcBulkCoefHeatInNeutCond, xy_CalcFlag,  & ! (in)
    & xy_SfcVelBulkCoef, xy_SfcTempBulkCoef, xy_SfcQVapBulkCoef,               & ! (out)
    & xy_BetaW, xy_SfcMOLength                                                 & ! (out)
    & )

    ! 宣言文 ; Declaration statements
    !
    real(DP), intent(in):: xy_SfcBulkRiNum (IA,JA)
    real(DP), intent(in):: xy_SfcRoughLengthMom (IA,JA)
    real(DP), intent(in):: xy_SfcRoughLengthHeat(IA,JA)
    real(DP), intent(in):: xy_SfcHeight(IA,JA)
    real(DP), intent(in):: xy_Height (IA,JA)
    real(DP), intent(in):: xy_SfcBulkCoefMomInNeutCond (IA,JA)
    real(DP), intent(in):: xy_SfcBulkCoefHeatInNeutCond(IA,JA)
    logical, intent(in) :: xy_CalcFlag(IA,JA)
    real(DP), intent(out):: xy_SfcVelBulkCoef (IA,JA)
    real(DP), intent(out):: xy_SfcTempBulkCoef (IA,JA)
    real(DP), intent(out):: xy_SfcQVapBulkCoef (IA,JA)
    real(DP), intent(out):: xy_BetaW(IA,JA)
    real(DP), intent(out):: xy_SfcMOLength(IA,JA)

    ! 作業変数
    ! Work variables
    !
    real(DP) :: SfcBulkRiNum
    real(DP) :: xy_MOLength(IA,JA)
    integer :: i
    integer :: j
    
    ! 実行文 ; Executable statement
    !

    xy_BetaW(:,:) = 0d0

    ! Calculate bulk coefficients in non-neutral condition
    !
    ! Parameterization by Louis et al. (1982)
    !
    
    !$omp parallel do collapse(2) private(SfcBulkRinum)
    do j = JS, JE
    do i = IS, IE

       if (xy_CalcFlag(i,j)) then
          if ( xy_SfcBulkRiNum(i,j) > 0.0_DP ) then 

             xy_SfcVelBulkCoef(i,j) =                                       &
                  &   xy_SfcBulkCoefMomInNeutCond(i,j)                         &
                  &   / (   1.0_DP                                              &
                  &       + 10.0_DP * xy_SfcBulkRiNum(i,j)                     &
                  &         / sqrt( 1.0_DP + 5.0_DP * xy_SfcBulkRiNum(i,j) )   &
                  &     )

             xy_SfcTempBulkCoef(i,j) =                                      &
                  &   xy_SfcBulkCoefHeatInNeutCond(i,j)                        &
                  &   / (   1.0_DP                                              &
                  &       + 15.0_DP * xy_SfcBulkRiNum(i,j)                     &
                  &         * sqrt( 1.0_DP + 5.0_DP * xy_SfcBulkRiNum(i,j) )   &
                  &     )

             xy_SfcQVapBulkCoef(i,j) = xy_SfcTempBulkCoef(i,j)
             
          else
          
             xy_SfcVelBulkCoef(i,j) =                                              &
                  &   xy_SfcBulkCoefMomInNeutCond(i,j)                                &
                  &   * (   1.0_DP                                                     &
                  &       - 10.0_DP * xy_SfcBulkRiNum(i,j)                            &
                  &         / (   1.0_DP                                               &
                  &             + 75.0_DP * xy_SfcBulkCoefMomInNeutCond(i,j)          &
                  &               * sqrt( - ( xy_Height(i,j) - xy_SfcHeight(i,j) + xy_SfcRoughLengthMom(i,j) ) &
                  &                         / xy_SfcRoughLengthMom(i,j)               &
                  &                         * xy_SfcBulkRiNum(i,j)                    &
                  &                     )                                              &
                  &           )                                                        &
                  &     )

             xy_SfcTempBulkCoef(i,j) =                                             &
                  &   xy_SfcBulkCoefHeatInNeutCond(i,j)                               &
                  &   * (   1.0_DP                                                     &
                  &       - 15.0_DP * xy_SfcBulkRiNum(i,j)                            &
                  &         / (   1.0_DP                                               &
                  &             + 75.0_DP * xy_SfcBulkCoefHeatInNeutCond(i,j)         &
                  &               * sqrt( - ( xy_Height(i,j) - xy_SfcHeight(i,j) + xy_SfcRoughLengthHeat(i,j) ) &
                  &                         / xy_SfcRoughLengthHeat(i,j)              &
                  &                         * xy_SfcBulkRiNum(i,j)                    &
                  &                     )                                              &
                  &           )                                                        &
                  &     )

             xy_SfcQVapBulkCoef(i,j) = xy_SfcTempBulkCoef(i,j)
          end if
       else
          xy_SfcVelBulkCoef(i,j) = 0d0
          xy_SfcTempBulkCoef(i,j) = 0d0
          xy_SfcQVapBulkCoef(i,j) = 0d0          
       end if

       ! Calculation of Monin-Obukhov length
       SfcBulkRiNum = xy_SfcBulkRiNum(i,j)
       if ( SfcBulkRiNum == 0.0_DP ) SfcBulkRiNum = 1.0e-10_DP
       xy_MOLength(i,j) =                               &
            &   ( xy_Height(i,j) - xy_SfcHeight(i,j) ) &
            & / ( FKarm * SfcBulkRiNum )                  &
            & * xy_SfcVelBulkCoef(i,j)**1.5_DP / xy_SfcTempBulkCoef(i,j)

       
       ! Measure maximum/minimum
       !       
       xy_SfcVelBulkCoef(i,j)  = &
            & max( min( xy_SfcVelBulkCoef(i,j), VelBulkCoefMax ), &
            &      VelBulkCoefMin )

       xy_SfcTempBulkCoef(i,j) = &
            & max( min( xy_SfcTempBulkCoef(i,j), TempBulkCoefMax ), &
            &      TempBulkCoefMin )

       xy_SfcQVapBulkCoef(i,j) = &
            & max( min( xy_SfcQVapBulkCoef(i,j), QVapBulkCoefMax ), &
            &      QVapBulkCoefMin )
       
    end do
    end do

    xy_SfcMOLength = xy_MOLength
    
  end subroutine BulkCoefL82
  
  subroutine read_config()
    
    VelMinForRi   = 0.01_DP
    VelMinForVel  = 0.01_DP
    VelMinForTemp = 0.01_DP
    VelMinForQVap = 0.01_DP
    VelMaxForVel  = 1000.0_DP
    VelMaxForTemp = 1000.0_DP
    VelMaxForQVap = 1000.0_DP

    VelBulkCoefMin  =  0.0_DP
    TempBulkCoefMin =  0.0_DP
    QVapBulkCoefMin =  0.0_DP
    VelBulkCoefMax  =  1.0_DP
    TempBulkCoefMax =  1.0_DP
    QVapBulkCoefMax =  1.0_DP
    
  end subroutine read_config
  
end module DSFCM_Util_SfcBulkFlux_mod
