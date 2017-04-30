!= dcpam 主プログラム
!
!= dcpam main program
!
! Authors::   Yasuhiro Morikawa, Satoshi Noda, Yoshiyuki O. Takahashi
! Version::   $Id: dcpam_main.f90,v 1.67 2015/03/11 04:54:54 yot Exp $ 
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008-2010. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module dcpam_main_mod
  
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  ! モデルの使い方については {チュートリアル}[link:../../../doc/tutorial/rakuraku/] を
  ! 参照してください.
  !
  ! See {Tutorial}[link:../../../doc/tutorial/rakuraku/index.htm.en] for usage of the 
  ! model. 
  !

  ! モジュール引用 ; USE statements
  !

  ! 力学過程 (スペクトル法, Arakawa and Suarez (1983))
  ! Dynamical process (Spectral method, Arakawa and Suarez (1983))
  !
  use dynamics_hspl_vas83, only: DynamicsHsplVAS83, VerticalFilterAdjust, SurfPresChangeWithWtVap

  ! 物理過程のみの計算のための力学過程
  ! A dynamics for calculation with physical processes only
  !
  use dynamics_physicsonly, only: DynamicsPhysicsOnly

  !
  ! Dynamical process for TWP-ICE experiment
  !
  use dynamics_twpice_scm_exp, only : DynamicsTWPICESCMExp

  ! Held and Suarez (1994) による強制と散逸
  ! Forcing and dissipation suggested by Held and Suarez (1994)
  !
  use held_suarez_1994, only: HS94Forcing

  ! 簡単金星計算のための強制
  ! forcing for simple Venus calculation
  !
  use yt2003_forcing, only: YT2003Forcing

  ! Schneider and Liu (2009) による鉛直混合課程
  ! Vertical diffusion by Schneider and Liu (2009)
  !
  use sl09_diffusion, only : SL09Diffusion

  ! 放射フラックス (GFD 電脳倶楽部開発の放射モデル)
  ! Radiation flux (radiation model developed by GFD Dennou Club)
  !
  use rad_DennouAGCM, only: RadDennouAGCMFlux

  ! 放射関連ルーチン
  ! Routines for radiation calculation
  !
  use rad_utils, only : RadDTempDt, RadFluxOutput

  ! 地球大気向け放射モデル Ver. 2
  ! radiation model for the Earth's atmosphere Ver. 2
  !
  use rad_Earth_V2, only: RadEarthV2Flux

  ! wrapper of RRTMG
  ! wrapper of RRTMG
  !
  use rad_rrtmg_wrapper, only: RadRRTMGWrapperFlux

  ! 火星大気向け放射モデル Ver. 1
  ! radiation model for the Mars' atmosphere Ver. 1
  !
  use rad_Mars_V1, only: RadMarsV1Flux

!!$  ! (火星大気向け) Non-LTE 放射モデル
!!$  ! Non-NLTE radiation model (for the Mars' atmosphere)
!!$  !
!!$  use rad_15m_NLTE, only: rad15mNLTEMergeHR

  ! 火星計算用近赤外加熱率計算
  ! Calculation of near infrared heating rate in the case of Mars
  !
  use rad_Mars_NIR, only : RadMarsNIRINOUT

  ! Schneider and Liu (2009) の放射モデル
  ! Radiation model by Schneider and Liu (2009)
  !
  use rad_SL09, only : RadSL09Flux

  ! 簡単放射モデル
  ! Simple radiation model
  !
  use rad_simple, only : RadSimpleFlux

  ! 何もしない放射モデル
  ! radiation model with no absorption and no scattering
  !
  use rad_none, only : RadNoneFlux

  ! 鉛直拡散フラックス
  ! Vertical diffusion flux
  !
  use vdiffusion_my, only: VDiffusionMY25, VDiffusion, VDiffusionOutput
  use vdiffusion_my, only: VDiffusionMY251DWrapper3D
  use vdiffusion_my, only: VDiffusionMY25GBT94

  ! JMA 乱流混合モジュール
  ! JMA turbulent mixing module
  !
  use vdiffusion_jma_my_wrapper, only : VDiffusionJMAMYWrapper3D

  ! 地下における熱の鉛直拡散
  ! Vertical diffusion of heat under the ground
  !
  use subsurface_diffusion_heat, only: &
!!$    & SubsurfaceDiffusionFlagLand, &
    & SubsurfaceDiffusion

  ! Gravity wave drag by McFarlane (1987)
  ! Gravity wave drag by McFarlane (1987)
  !
  use gwd_m1987, only : GWDM1987

  ! 積雲パラメタリゼーション (対流調節)
  ! Cumulus parameterization (convection adjust)
  !
  use moist_conv_adjust, only: MoistConvAdjust, MoistConvAdjustI98

  ! Relaxed Arakawa-Schubert scheme
  ! Relaxed Arakawa-Schubert scheme
  !
  use relaxed_arakawa_schubert, only : RAS
  use relaxed_arakawa_schubert, only : RAS1DWrapper3D
  use relaxed_arakawa_schubert, only : RASWithIce1DWrapper3DWrapper

  ! 大規模凝結 (非対流性凝結)
  ! Large scale condensation
  !
  use lscond, only: LScaleCond
  use lscond, only: LScaleCond1D3DWrapper
  use lscond, only: LScaleCondLL911D3DWrapper

  ! 大規模凝結 (非対流性凝結) (Le Treut and Li, 1991)
  ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
  !
!!$  use lscond_LL91, only : LScaleCondLL91

  !
  != Saturation adjustment
  !
  use saturation_adjust, only : SaturationAdjust

  ! 地表面フラックス
  ! Surface flux
  use surface_flux_bulk, only: SurfaceFlux, SurfaceFluxOutput

  !
  ! set dust flux
  !
!!$  use set_dust_flux, only : SetDustFlux

  ! 下部境界フラックス
  ! Lower boundary flux
  !
  use lb_flux_simple, only : LBFluxSimple

  ! 乾燥対流調節
  ! Dry convective adjustment
  !
  use dry_conv_adjust, only: DryConvAdjust

  ! 質量の補正
  ! Mass fixer
  !
  use mass_fixer, only: MassFixerColumn

  ! 温度の半整数σレベルの補間, 気圧と高度の算出
  ! Interpolate temperature on half sigma level, 
  ! and calculate pressure and height
  !
  use auxiliary, only: AuxVars

  ! 陰解法による時間積分のためのルーチン
  ! Routines for time integration with implicit scheme
  !
  use phy_implicit_utils, only : PhyImplEvalRadLFluxA

  ! 陰解法のための行列処理 (一部の物理過程用)
  ! Matrices handling for implicit scheme (for a part of physical processes)
  !
  use phy_implicit, only: PhyImplTendency

  ! 陰解法のための行列処理 (一部の物理過程用)
  ! Matrices handling for implicit scheme (for a part of physical processes)
  !
!!$  use phy_implicit_sdh, only:              &
!!$    & PhyImplSDHSetMethodFromMatthews,     &
!!$    & PhyImplSDHTendency,                  &
!!$    & PhyImplSDHCorSOTempBySnowMelt

!!$  ! 陰解法のための行列処理 (一部の物理過程用)
!!$  ! Matrices handling for implicit scheme (for a part of physical processes)
!!$  !
!!$  use phy_implicit_sdh_V2, only:       &
!!$    & PhyImplSDHV2SetMethodMatthews,   &
!!$    & PhyImplSDHV2Tendency,            &
!!$    & PhyImplSDHV2CorSOTempBySnowMelt

  ! 陰解法のための行列処理 (一部の物理過程用)
  ! Matrices handling for implicit scheme (for a part of physical processes)
  !
  use phy_implicit_sdh_V3, only:       &
    & PhyImplSDHV3SetMethodMatthews,   &
    & PhyImplSDHV3Tendency,            &
    & PhyImplSDHV3CorSOTempBySnowMelt

  ! 陰解法による時間積分 (大気のみ / 惑星表面温度・土壌温度計算なし)
  ! Time integration by using implicit scheme in case without calculation of surface and soil temperature
  !
  use phy_implicit_atmonly, only : PhyImplAtmOnlyTendency

  ! 地面温度の時間積分・地表面放射補正
  ! Time integration of surface temperature, correction of flux on surface
  !
  use intg_surftemp, only: IntegralSurfTemp

  ! バケツモデル
  ! Bucket model
  !
  use Bucket_Model, only :            &
    & BucketSetFlagOceanFromMatthews, &
    & BucketModEvapAndLatentHeatFlux, &
    & BucketIntegration,              &
    & BucketPRCPAdjust

  ! タイムフィルター (Asselin, 1972)
  ! Time filter (Asselin, 1972)
  !
  use timefilter_asselin1972, only: TimeFilter, TimeFilterSurfVars

  ! 時間フィルター (Williams, 2009)
  ! Time filter (Williams, 2009)
  !
  use timefilter_williams2009, only: &
    & TimeFilterWilliams2009, TimeFilterWilliams2009SurfVars

  ! 時刻管理
  ! Time control
  !
  use timeset, only: TimesetProgress, &
    & TimeB, &                ! ステップ $ t - \Delta t $ の時刻. 
                              ! Time of step $ t - \Delta t $. 
    & TimeN, &                ! ステップ $ t $ の時刻. 
                              ! Time of step $ t $. 
    & TimeA, &                ! ステップ $ t + \Delta t $ の時刻. 
                              ! Time of step $ t + \Delta t $. 
    & EndTime, &              ! 計算終了時刻. 
                              ! End time of calculation
    & DelTime                 ! $ \Delta t $ [s]

  ! リスタートデータ入出力
  ! Restart data input/output
  !
  use restart_file_io, only: RestartFileOutPut

  ! 地表面温度リスタートデータ入出力
  ! Restart data of surface temperature input/output
  !
  use restart_surftemp_io, only: RestartSurfTempOutPut

  ! 惑星表面特性の設定
  ! Setting of surface properties
  !
  use surface_properties, only: SetSurfaceProperties

  ! 雪, 氷の割合
  ! snow/ice fraction
  !
  use snowice_frac, only : CalcSnowFrac

  ! ヒストリデータ出力
  ! History data output
  !
  use gtool_historyauto, only: HistoryAutoPut, HistoryAutoAllVarFix

  ! 組成に関わる配列の設定
  ! Settings of array for atmospheric composition
  !
  use composition, only: &
    &                    ncmax, &
                              ! 成分の数
                              ! Number of composition
    &                    a_QMixName, &
                              ! 成分の変数名
                              ! Name of variables for composition
    &                    a_QMixLongName, &
                              ! 成分の長い変数名
                              ! Long name of variables for composition
    &                    IndexH2OVap, &
    &                    IndexH2OLiq, &
    &                    IndexH2OSol, &
    &                    IndexTKE,    &
    &                    CompositionInqIndex


  ! 格子点設定
  ! Grid points settings
  !
  use gridset, only: imax, &  ! 経度格子点数. 
                              ! Number of grid points in longitude
    &                jmax, &  ! 緯度格子点数. 
                              ! Number of grid points in latitude
    &                kmax     ! 鉛直層数. 
                              ! Number of vertical level

  ! 物理定数設定
  ! Physical constants settings
  !
  use constants, only:  &
    & LatentHeat      , &
    & LatentHeatFusion, &
    & Grav
                              ! $ g $ [m s-2]. 
                              ! 重力加速度. 
                              ! Gravitational acceleration

  !
  !
  use dc_message, only: MessageNotify

  ! 種別型パラメタ
  ! Kind type parameter
  !
  use dc_types, only: DP, &      ! 倍精度実数型. Double precision. 
    &                 STRING, &  ! 文字列.       Strings. 
    &                 TOKEN      ! キーワード.   Keywords. 

  ! 雲なしモデル
  ! No cloud model
  !
  use cloud_none, only : &
    & CloudNone,         &
    & CloudNoneWithIce

  ! 簡単雲モデル
  ! Simple cloud
  !
  use cloud_simple, only:               &
    & CloudSimple,                      &
    & CloudSimpleWithIce,               &
    & CloudSimpleCalcCloudCover,        &
    & CloudSimpleDivideWatAndIce,       &
    & CloudSimpleCalcPRCPKeyLLTemp3D

  ! Tiedtke (1993) に基づく雲モデル
  ! Cloud model based on Tiedtke (1993)
  !
  use cloud_T1993base, only : CloudT1993baseWithIce

  ! 火星 H2O 雲モデル
  ! Mars H2O cloud model
  !
  use cloud_mars_h2o, only : CloudMarsH2O

  ! 重力沈降過程
  ! Gravitational sedimentation process
  !
  use grav_sed, only : GravSed

  ! 主成分相変化
  ! Phase change of atmospheric major component
  !
  use major_comp_phase_change, only : &
    & MajorCompPhaseChangeInAtm,      &
    & MajorCompPhaseChangeOnGround

  ! 予報変数の値の確認
  ! Check values of prognostic variables
  !
  use check_prog_vars, only: CheckProgVars

  ! Output frequently used variables
  ! Output frequently used variables
  !
  use output_freq_used_vars, only : OutputFreqUsedVars

  use ProfUtil_mod

  ! 宣言文 ; Declaration statements
  !
  implicit none
  
  !
  public :: dcpam_main_Init, dcpam_main_Final
  public :: MainInit, MainTerminate
  public :: dcpam_advance_timestep
  public :: dcpam_UpdateSurfaceProperties
  
  character(*), parameter:: prog_name = 'dcpam_main'
                            ! 主プログラム名. 
                            ! Main program name

  ! 予報変数 (ステップ $ t-\Delta t $ , $ t $ , $ t+\Delta t $ )
  ! Prediction variables  (Step $ t-\Delta t $ , $ t $ , $ t+\Delta t $ )
  !
  real(DP), allocatable:: xyz_UB (:,:,:)
                              ! $ u (t-\Delta t) $ .   東西風速. Eastward wind (m s-1)
  real(DP), allocatable:: xyz_VB (:,:,:)
                              ! $ v (t-\Delta t) $ .   南北風速. Northward wind (m s-1)
  real(DP), allocatable:: xyz_TempB (:,:,:)
                              ! $ T (t-\Delta t) $ .   温度. Temperature (K)
  real(DP), allocatable:: xyzf_QMixB(:,:,:,:)
                              ! $ q (t-\Delta t) $ .   混合比. Mass mixing ratio of constituents (1)
  real(DP), allocatable:: xy_PsB (:,:)
                              ! $ p_s (t-\Delta t) $ . 地表面気圧. Surface pressure (Pa)
  real(DP), allocatable:: xyz_UN (:,:,:)
                              ! $ u (t) $ .     東西風速. Eastward wind (m s-1)
  real(DP), allocatable:: xyz_VN (:,:,:)
                              ! $ v (t) $ .     南北風速. Northward wind (m s-1)
  real(DP), allocatable:: xyz_TempN (:,:,:)
                              ! $ T (t) $ .     温度. Temperature (K)
  real(DP), allocatable:: xyzf_QMixN(:,:,:,:)
                              ! $ q (t) $ .     混合比. Mass mixing ratio of constituents (1)
  real(DP), allocatable:: xy_PsN (:,:)
                              ! $ p_s (t) $ .   地表面気圧. Surface pressure (Pa)
  real(DP), allocatable:: xyz_UA (:,:,:)
                              ! $ u (t+\Delta t) $ .   東西風速. Eastward wind (m s-1)
  real(DP), allocatable:: xyz_VA (:,:,:)
                              ! $ v (t+\Delta t) $ .   南北風速. Northward wind (m s-1)
  real(DP), allocatable:: xyz_TempA (:,:,:)
                              ! $ T (t+\Delta t) $ .   温度. Temperature (K)
  real(DP), allocatable:: xyzf_QMixA(:,:,:,:)
                              ! $ q (t+\Delta t) $ .   混合比. Mass mixing ratio of constituents (1)
  real(DP), allocatable:: xy_PsA (:,:)
                              ! $ p_s (t+\Delta t) $ . 地表面気圧. Surface pressure (Pa)

  real(DP), allocatable:: xy_SurfMajCompIceB (:,:)
                              ! $ M_mcs (t-\Delta t) $ . (kg m-2)
                              ! Surface major component ice amount (kg m-2)
  real(DP), allocatable:: xy_SoilMoistB (:,:)
                              ! $ M_ws (t-\Delta t) $ . 土壌水分 (kg m-2)
                              ! Soil moisture (kg m-2)
  real(DP), allocatable:: xy_SurfSnowB (:,:)
                              ! $ M_ss (t-\Delta t) $ . 積雪量 (kg m-2)
                              ! Surface snow amount (kg m-2)
  real(DP), allocatable:: xy_SurfMajCompIceN (:,:)
                              ! $ M_mcs (t) $ . (kg m-2)
                              ! Surface major component ice amount (kg m-2)
  real(DP), allocatable:: xy_SoilMoistN (:,:)
                              ! $ M_ws (t) $          . 土壌水分 (kg m-2)
                              ! Soil moisture (kg m-2)
  real(DP), allocatable:: xy_SurfSnowN (:,:)
                              ! $ M_ss (t) $ . 積雪量 (kg m-2)
                              ! Surface snow amount (kg m-2)
  real(DP), allocatable:: xy_SurfMajCompIceA (:,:)
                              ! $ M_mcs (t+\Delta t) $ . (kg m-2)
                              ! Surface major component ice amount (kg m-2)
  real(DP), allocatable:: xy_SoilMoistA (:,:)
                              ! $ M_ws (t+\Delta t) $ . 土壌水分 (kg m-2)
                              ! Soil moisture (kg m-2)
  real(DP), allocatable:: xy_SurfSnowA (:,:)
                              ! $ M_ss (t+\Delta t) $ . 積雪量 (kg m-2)
                              ! Surface snow amount (kg m-2)


  ! 診断変数, 他
  ! Diagnostic variables, etc.
  !
  real(DP), allocatable:: xyz_DUDt (:,:,:)
                              ! $ \DP{u}{t} $ . 東西風速変化 (m s-2)
                              ! Eastward wind tendency (m s-2)
  real(DP), allocatable:: xyz_DVDt (:,:,:)
                              ! $ \DP{v}{t} $ . 南北風速変化 (m s-2)
                              ! Northward wind tendency (m s-2)
  real(DP), allocatable:: xyz_DTempDt (:,:,:)
                              ! $ \DP{T}{t} $ . 温度変化 (K s-1)
                              ! Temperature tendency (K s-1)
  real(DP), allocatable:: xyzf_DQMixDt(:,:,:,:)
                              ! $ \DP{q}{t} $ . 混合比変化 (s-1)
                              ! Mass mixing ratio tendency (s-1)

  real(DP), allocatable:: xyz_DTurKinEneDt(:,:,:)
                              ! 
                              ! Turbulent kinetic energy tendency (m2 s-3)

  real(DP), allocatable:: xy_SurfHeight (:,:)
                              ! $ z_s $ . 地表面高度 (m)
                              ! Surface height (m)
  real(DP), allocatable:: xy_SurfHeightStd (:,:)
                              ! 
                              ! Surface height standard deviation(m)

  real(DP), allocatable:: xy_SurfTemp (:,:)
                              ! 地表面温度 (K)
                              ! Surface temperature (K)
  real(DP), allocatable:: xyz_SoilTemp(:,:,:)
                              ! 土壌温度 (K)
                              ! Soil temperature (K)


  real(DP), allocatable:: xy_SurfAlbedo (:,:)
                              ! 地表アルベド (1)
                              ! Surface albedo (1)
  real(DP), allocatable:: xy_SurfHumidCoef (:,:)
                              ! 地表湿潤度 (1)
                              ! Surface humidity coefficient (1)
  real(DP), allocatable:: xy_SurfRoughLenMom (:,:)
                              ! 地表粗度長 (m)
                              ! Surface rough length for momentum (m)
  real(DP), allocatable:: xy_SurfRoughLenHeat(:,:)
                              ! 地表粗度長 (m)
                              ! Surface rough length for heat (m)
  real(DP), allocatable:: xy_SurfHeatCapacity (:,:)
                              ! 地表熱容量 (J m-2 K-1)
                              ! Surface heat capacity (J m-2 K-1)
  real(DP), allocatable:: xy_SeaIceConc(:,:)
                              ! 海氷密度 (0 <= xy_SeaIceConc <= 1) (1)
                              ! Sea ice concentration (0 <= xy_SeaIceConc <= 1) (1)
  integer , allocatable:: xy_SurfCond (:,:)
                              ! 惑星表面状態 (0: 固定, 1: 可変) (1)
                              ! Surface condition (0: fixed, 1: variable) (1)
  integer , allocatable:: xy_SurfType (:,:)
                              ! 惑星表面タイプ (土地利用, Matthews 分布) (1)
                              ! Surface type (land use type classified by Matthews) (1)
  real(DP), allocatable:: xy_DeepSubSurfHeatFlux (:,:)
                              ! 地中熱フラックス (W m-2)
                              ! "Deep subsurface heat flux" (W m-2)
                              ! Heat flux at the bottom of surface/soil layer.
  real(DP), allocatable:: xy_SoilHeatCap (:,:)
                              ! 土壌熱容量 (J K-1 kg-1)
                              ! Specific heat of soil (J K-1 kg-1)
  real(DP), allocatable:: xy_SoilHeatDiffCoef (:,:)
                              ! 土壌熱伝導係数 (J m-3 K-1)
                              ! Heat conduction coefficient of soil (J m-3 K-1)

  real(DP), allocatable:: xy_SnowFrac(:,:)
                              ! 
                              ! Snow fraction (1)

!!$  logical , allocatable:: xy_FlagMatthewsLand(:,:)
!!$                              !
!!$                              ! Flag for land grid point based on Matthews
  integer , allocatable:: xy_PhyImplSDHIndexCalcMethod(:,:)
                              !
                              ! Index for calculation method used in PhyImplSDHTendency
  logical , allocatable:: xy_BucketFlagOceanGrid(:,:)
                              !
                              ! Flag for ocean grid point used in bucket model

  real(DP), allocatable:: xyr_Temp (:,:,:)
                              ! $ \hat{T} $ . 温度 (半整数レベル) (K)
                              ! Temperature (half level) (K)
  real(DP), allocatable:: xyz_VirTemp (:,:,:)
                              ! $ T_v $ . 仮温度 (K)
                              ! Virtual temperature (K)
  real(DP), allocatable:: xyr_VirTemp (:,:,:)
                              ! $ \hat{T}_v $ . 仮温度 (半整数レベル) (K)
                              ! Virtual temperature (half level) (K)
  real(DP), allocatable:: xy_SurfVirTemp (:,:)
                              ! $ \hat{T}_{v,s} $ . 仮温度 (惑星表面) (K)
                              ! Virtual temperature (surface) (K)
  real(DP), allocatable:: xyz_Press (:,:,:)
                              ! $ p $ . 気圧 (整数レベル) (Pa)
                              ! Air pressure (full level) (Pa)
  real(DP), allocatable:: xyr_Press (:,:,:)
                              ! $ \hat{p} $ . 気圧 (半整数レベル) (Pa)
                              ! Air pressure (half level) (Pa)
  real(DP), allocatable:: xyz_Height (:,:,:)
                              ! 高度 (整数レベル) (m)
                              ! Height (full level) (m)
  real(DP), allocatable:: xyr_Height (:,:,:)
                              ! 高度 (半整数レベル) (m)
                              ! Height (half level) (m)
  real(DP), allocatable:: xyz_Exner (:,:,:)
                              ! Exner 関数 (整数レベル) (1)
                              ! Exner function (full level) (1)
  real(DP), allocatable:: xyr_Exner (:,:,:)
                              ! Exner 関数 (半整数レベル) (1)
                              ! Exner function (half level) (1)

  real(DP), allocatable:: xyr_RadLUwFlux (:,:,:)
                              ! 長波フラックス (W m-2)
                              ! Upward longwave flux (W m-2)
  real(DP), allocatable:: xyr_RadLDwFlux (:,:,:)
                              ! 長波フラックス (W m-2)
                              ! Downward longwave flux (W m-2)
  real(DP), allocatable:: xyr_RadLFlux  (:,:,:)
  real(DP), allocatable:: xyr_RadLFluxA (:,:,:)
                              ! 陰解法で解いた地表面熱収支と整合的な $ t+\Delta t $ に
                              ! おける長波フラックスの計算
                              !   * ここで計算された値が直接次のステップ $ t $ における
                              !     長波フラックスとして用いられるわけではない
                              !   * 現在の時間ステップにおける長波放射加熱率の計算に使
                              !     われる
                              !
                              ! Evaluate longwave flux at $ t+\Delta t $ consistent 
                              ! with surface energy balance solved with implicit method
                              !   * The evaluated value is not used directly as Longwave
                              !     flux at next step $ t $.
                              !   * The evaluated value is used to calculate long wave 
                              !     radiative heating rate in the current time step.
  real(DP), allocatable:: xyr_RadSFlux (:,:,:)
  real(DP), allocatable:: xyr_RadSUwFlux (:,:,:)
                              ! 短波 (日射) フラックス (W m-2)
                              ! Upward shortwave flux (W m-2)
  real(DP), allocatable:: xyr_RadSDwFlux (:,:,:)
                              ! 短波 (日射) フラックス (W m-2)
                              ! Downward shortwave flux (W m-2)
  real(DP), allocatable:: xyra_DelRadLFlux   (:,:,:,:)
                              ! 長波地表温度変化 (W m-2)
  real(DP), allocatable:: xyra_DelRadLUwFlux (:,:,:,:)
                              ! 長波地表温度変化 (W m-2)
                              ! 
  real(DP), allocatable:: xyra_DelRadLDwFlux (:,:,:,:)
                              ! 長波地表温度変化 (W m-2)
                              ! 

  real(DP), allocatable:: xyr_MomFluxX (:,:,:)
                              ! 東西方向運動量フラックス
                              ! Eastward momentum flux
  real(DP), allocatable:: xyr_MomFluxY (:,:,:)
                              ! 南北方向運動量フラックス. 
                              ! Northward momentum flux
  real(DP), allocatable:: xyr_HeatFlux (:,:,:)
                              ! 熱フラックス. 
                              ! Heat flux
  real(DP), allocatable:: xyrf_QMixFlux(:,:,:,:)
                              ! 成分質量フラックス. 
                              ! Mass flux of compositions

  real(DP), allocatable:: xy_SurfMomFluxX (:,:)
                              ! 惑星表面東西方向運動量フラックス
                              ! Eastward momentum flux at surface
  real(DP), allocatable:: xy_SurfMomFluxY (:,:)
                              ! 惑星表面南北方向運動量フラックス. 
                              ! Northward momentum flux at surface
  real(DP), allocatable:: xy_SurfHeatFlux (:,:)
                              ! 惑星表面熱フラックス. 
                              ! Heat flux at surface
  real(DP), allocatable:: xyf_SurfQMixFlux(:,:,:)
                              ! 惑星表面成分質量フラックス. 
                              ! Mass flux of compositions at surface

  real(DP), allocatable:: xy_SurfH2OVapFluxA(:,:)
                              ! 惑星表面水蒸気フラックス.
                              ! Water vapor flux at the surface
  real(DP), allocatable:: xy_SurfLatentHeatFluxA(:,:)
                              ! 惑星表面潜熱フラックス.
                              ! Latent heat flux at the surface
        ! NOTE:
        ! Only if the evaporation of liquid water is considered, a variable, 
        ! xy_SurfLatentHeatFlux is not required, since a latent heat flux
        ! at the surface is equal to water mass flux times latent heat. 
        ! But, if the evaporation of snow is considered, that is not the case 
        ! and a variable for the latent heat flux is required in addition to 
        ! that for the water mass flux.
        !

  real(DP), allocatable:: xyr_SoilHeatFlux (:,:,:)
                              ! 土壌中の熱フラックス (W m-2)
                              ! Heat flux in sub-surface soil (W m-2)

  real(DP), allocatable:: xyr_VelDiffCoef (:,:,:)
                              ! 拡散係数：運動量. 
                              ! Diffusion coefficient: velocity
  real(DP), allocatable:: xyr_TempDiffCoef (:,:,:)
                              ! 拡散係数：温度. 
                              ! Transfer coefficient: temperature
  real(DP), allocatable:: xyr_QMixDiffCoef (:,:,:)
                              ! 拡散係数：比湿. 
                              ! Diffusion coefficient: specific humidity


  real(DP), allocatable:: xy_SurfVelTransCoef (:,:)
                              ! 輸送係数：運動量. 
                              ! Diffusion coefficient: velocity
  real(DP), allocatable:: xy_SurfTempTransCoef (:,:)
                              ! 輸送係数：温度. 
                              ! Transfer coefficient: temperature
  real(DP), allocatable:: xy_SurfQVapTransCoef (:,:)
                              ! 輸送係数：水蒸気
                              ! Transfer coefficient: water vapor
  real(DP), allocatable:: xyr_SoilTempTransCoef (:,:,:)
                              ! 輸送係数：土壌温度.
                              ! Transfer coefficient: soil temperature

  real(DP), allocatable:: xy_SurfMOLength(:,:)
                              ! 
                              ! Monin-Obukov length

  real(DP), allocatable:: xy_DSurfTempDt (:,:)
                              ! 地表面温度変化率. 
                              ! Surface temperature tendency
  real(DP), allocatable:: xyz_DSoilTempDt (:,:,:)
                              ! $ \DP{Tg}{t} $ . 土壌温度変化 (K s-1)
                              ! Temperature tendency (K s-1)

  real(DP), allocatable:: xy_DPsDt (:,:)
                              !  (Pa s-1)
                              ! Surface pressure tendency (Pa s-1)

  real(DP), allocatable:: xy_DSurfMajCompIceDt (:,:)
                              !  (kg m-2 s-1)
                              ! Surface major component ice tendency (kg m-2 s-1)
  real(DP), allocatable:: xy_DSoilMoistDt (:,:)
                              ! 土壌水分時間変化率 (kg m-2 s-1)
                              ! Soil temperature tendency (kg m-2 s-1)
  real(DP), allocatable:: xy_DSurfSnowDt (:,:)
                              ! 積雪時間変化率 (kg m-2 s-1)
                              ! Surface snow amount tendency (kg m-2 s-1)

  real(DP), allocatable:: xyz_DTempDtVDiff(:,:,:)
                              ! 鉛直拡散による加熱率 (K s-1)
                              ! Temperature tendency due to vertical diffusion (K s-1)

  real(DP), allocatable:: xyz_DTempDtRadL (:,:,:)
                              ! 長波加熱率. 
                              ! Temperature tendency with longwave
  real(DP), allocatable:: xyz_DTempDtRadS (:,:,:)
                              ! 短波加熱率. 
                              ! Temperature tendency with shortwave

  real(DP), allocatable:: xyz_DUDtGWD (:,:,:)
                              ! $ \DP{u}{t} $ . 東西風速変化 (m s-2)
                              ! Eastward wind tendency by gravity wave drag (m s-2)
  real(DP), allocatable:: xyz_DVDtGWD (:,:,:)
                              ! $ \DP{v}{t} $ . 南北風速変化 (m s-2)
                              ! Northward wind tendency by gravity wave drag (m s-2)

  real(DP), allocatable:: xyz_OMG (:,:,:)
                              ! Vertical velocity in pressure coordinate

  real(DP), allocatable:: xy_Rain (:,:)
                              ! 降水量. 
                              ! Precipitation
  real(DP), allocatable:: xy_RainCumulus  (:,:)
                              ! 
                              ! Rain due to moist convection
  real(DP), allocatable:: xy_RainLsc      (:,:)
                              ! 
                              ! Rain due to non-convective condensation

  real(DP), allocatable:: xy_Snow (:,:)
                              ! 
                              ! Snow fall
  real(DP), allocatable:: xy_SnowCumulus  (:,:)
                              ! 
                              ! Snow fall due to moist convection
  real(DP), allocatable:: xy_SnowLsc      (:,:)
                              ! 
                              ! Snow fall due to non-convective condensation

  real(DP), allocatable:: xyz_DTempDtCum(:,:,:)
  real(DP), allocatable:: xyz_DQVapDtCum(:,:,:)
  real(DP), allocatable:: xyz_DQH2OLiqDtCum(:,:,:)
                              ! Production rate of cloud water in the layer 
                              ! due to condensation in cumulus convection 
                              ! parameterization (kg kg-1)
  real(DP), allocatable:: xyz_DQH2OSolDtCum(:,:,:)
                              ! Production rate of cloud ice in the layer 
                              ! due to condensation in cumulus convection 
                              ! parameterization (kg kg-1)
  real(DP), allocatable:: xyz_DUDtCum(:,:,:)
  real(DP), allocatable:: xyz_DVDtCum(:,:,:)

  real(DP), allocatable:: xyz_MoistConvDetTend       (:,:,:)
  real(DP), allocatable:: xyz_MoistConvSubsidMassFlux(:,:,:)

  real(DP), allocatable:: xyz_DQH2OLiqDtLSC(:,:,:)
                              ! Production rate of cloud water in the layer 
                              ! due to condensation in large scale condensation 
                              ! (kg kg-1)
  real(DP), allocatable:: xyz_DQH2OSolDtLSC(:,:,:)
                              ! Production rate of cloud ice in the layer 
                              ! due to condensation in large scale condensation 
                              ! (kg kg-1)

  real(DP), allocatable:: xyz_QH2OLiqforRad(:,:,:)
                              ! Array for liquid water for radiation calculation
                              ! (kg kg-1)
  real(DP), allocatable:: xyz_QH2OSolforRad(:,:,:)
                              ! Array for solid water (ice) for radiation calculation
                              ! (kg kg-1)
  real(DP), allocatable:: xyz_CloudCoverforRad(:,:,:)
                              ! Cloud cover (cloud fraction)

  real(DP), allocatable:: xy_SurfDustGravSedFlux(:,:)


  !* For coupler
  real(DP), allocatable, dimension(:,:) :: &
       & xy_TauXAtm, xy_TauYAtm, xy_SensAtm, xy_LatentAtm,           &
       & xy_LDWRFlxAtm, xy_LUWRFlxAtm, xy_SDWRFlxAtm, xy_SUWRFlxAtm, &
       & xy_SurfAirTemp, xy_DSurfHFlxDTs, xy_DSurfLatentFlxDTs, &
       & xy_RainAtm, xy_SnowAtm
  
  ! 作業変数
  ! Work variables
  !
  integer           :: IDDynMode                 ! 使用する力学過程
                                                 ! Dynamics used for an experiment
  !
  integer, parameter:: IDDynModeHSPLVAS83       = 0
  integer, parameter:: IDDynModeNoHorAdv        = 1
  integer, parameter:: IDDynModeTWPICE          = 2

  integer           :: IDPhysMode                 ! 使用する物理過程
                                                 ! Physics used for an experiment
  !
  integer, parameter:: IDPhysModeNoPhysics       = 0
  integer, parameter:: IDPhysModeFullPhysics     = 1
  integer, parameter:: IDPhysModeHS94            = 2
  integer, parameter:: IDPhysModeVenusSimple     = 3
  integer, parameter:: IDPhysModeJupiterSimple   = 4
  integer, parameter:: IDPhysModeJupiterSimpleV2 = 5

  integer           :: IDPhyTendMethod           ! 物理過程による変化率の計算方法
                                                 ! Method calculating physics tendency
  !
  integer, parameter:: IDPhyTendMethodImp1LayModel = 10
  integer, parameter:: IDPhyTendMethodImpSoilModel = 11
  integer, parameter:: IDPhyTendMethodImpAtmOnly   = 12

  integer           :: IDRadMethod               ! 放射過程の計算方法
                                                 ! Method for radiation
  !
  integer, parameter:: IDRadMethodDennouAGCM  = 20
  integer, parameter:: IDRadMethodEarthV2     = 21
  integer, parameter:: IDRadMethodMarsV1      = 22
  integer, parameter:: IDRadMethodSL09        = 23
!!$  integer, parameter:: IDRadMethodVenusSimple = 24
  integer, parameter:: IDRadMethodSimple      = 25
  integer, parameter:: IDRadMethodRRTMG       = 26
  integer, parameter:: IDRadMethodNone        = 27


  integer           :: IDSfcFluxMethod          ! 
                                                ! Method for surface flux evaluation
  !
  integer, parameter:: IDSfcFluxMethodL82     = 90
  integer, parameter:: IDSfcFluxMethodBH91B94 = 91


  integer           :: IDVDiffMethod            ! 
                                                ! Method for vertical diffusion
  !
  integer, parameter:: IDVDiffMethodMY2   = 80
  integer, parameter:: IDVDiffMethodMY25  = 81
  integer, parameter:: IDVDiffMethodJMA   = 82

  integer           :: IDMCMethod                ! 湿潤対流の計算方法
                                                 ! Method for moist convection

  integer, parameter:: IDMCMethodNone        = 30
  integer, parameter:: IDMCMethodMCA         = 31
  integer, parameter:: IDMCMethodRAS         = 32
  integer, parameter:: IDMCMethodRASWithIce  = 33
  integer, parameter:: IDMCMethodMCAI98      = 34

  integer           :: IDLSCMethod               ! 大規模凝結 (非対流性凝結) の計算方法
                                                 ! Method for large scale condensation
                                                 ! (non-convective condensation)
  !
  integer, parameter:: IDLSCMethodNone        = 40
  integer, parameter:: IDLSCMethodM65         = 41
  integer, parameter:: IDLSCMethodLL91        = 42
  integer, parameter:: IDLSCMethodSatAdjM65   = 43
  integer, parameter:: IDLSCMethodM65WithIce  = 44
  integer, parameter:: IDLSCMethodLL91WithIce = 45

  integer           :: IDCloudMethod             ! 雲の計算方法
                                                 ! Method for cloud
  !
  integer, parameter:: IDCloudMethodNone          = 50
  integer, parameter:: IDCloudMethodSimple        = 51
  integer, parameter:: IDCloudMethodSimpleWithIce = 52
  integer, parameter:: IDCloudMethodT1993WithIce  = 53
  integer, parameter:: IDCloudMethodMarsH2OCloud  = 54

  integer           :: IDSfcMoistMethod          ! 惑星表面水分の計算方法
                                                 ! Method for surface moisture calculation
  !
  integer, parameter:: IDSfcMoistMethodNone   = 60
  integer, parameter:: IDSfcMoistMethodBucket = 61

  integer           :: IDDCMethod                ! 乾燥対流の計算方法
                                                 ! Method for dry convection
  !
  integer, parameter:: IDDCMethodNone = 70
  integer, parameter:: IDDCMethodDCA  = 71

  integer           :: IDGWDMethod               ! 
                                                 ! Method for gravity wave drag
  !
  integer, parameter:: IDGWDMethodNone  = 80
  integer, parameter:: IDGWDMethodM1987 = 81


  logical:: firstloop = .true.
                              ! 初回のループであることを示すフラグ. 
                              ! Flag implying first loop

  logical:: flag_initial
                              ! 内部サブルーチン MainInit で設定されます. 
                              ! リスタートデータを読み込む場合には, 
                              ! .false. が, 初期値データを読み込む場合には
                              ! .true. が設定されます. 
                              ! 
                              ! This variable is set in an internal 
                              ! subroutine "MainInit". 
                              ! If restart data is loaded, .false. is set. 
                              ! On the other hand, if initial data is loaded, 
                              ! .true. is set.


  integer :: MPI_MY_COMM

  real(DP) :: WtMassA, WtMassB
  
contains

  !------------------------------------------

  subroutine dcpam_main_Init(MY_COMM)
    integer, intent(in), optional :: MY_COMM

    MPI_MY_COMM = -1
    if(present(MY_COMM)) MPI_MY_COMM = MY_COMM
    
  end subroutine dcpam_main_Init

  subroutine dcpam_main_Final()
  end subroutine dcpam_main_Final
  
  !------------------------------------------------------------------

  ! Interface for coupler to update variables of surface properties in DCPAM.
  subroutine dcpam_UpdateSurfaceProperties( &
       & xy_SurfTempRecv, xy_SurfAlbedoRecv, xy_SeaIceConcRecv, xy_SurfSnowRecv,  &
       & xy_SfcEngyFlxModRecv                                                     &
       & )
    use constants, only: Grav, CpDry
    
    real(DP), intent(in), optional :: xy_SurfTempRecv(0:iMax-1,jMax)
    real(DP), intent(in), optional :: xy_SurfAlbedoRecv(0:iMax-1,jMax)
    real(DP), intent(in), optional :: xy_SeaIceConcRecv(0:iMax-1,jMax)
    real(DP), intent(in), optional :: xy_SurfSnowRecv(0:iMax-1,jMax)
    real(DP), intent(in), optional :: xy_SfcEngyFlxModRecv(0:iMax-1,jMax)

    if(present(xy_SurfTempRecv))  xy_SurfTemp(:,:) = xy_SurfTempRecv
    if(present(xy_SurfAlbedoRecv)) xy_SurfAlbedo(:,:) = xy_SurfAlbedoRecv
    if(present(xy_SeaIceConcRecv)) xy_SeaIceConc(:,:) = xy_SeaIceConcRecv
    if(present(xy_SurfSnowRecv)) xy_SurfSnowB(:,:) = xy_SurfSnowRecv

    if (present(xy_SfcEngyFlxModRecv)) then
       call AuxVars( &
            & xy_PsN, xyz_TempN, xyzf_QMixN(:,:,:,IndexH2OVap), & ! (in )
            & xyr_Press = xyr_Press                             & ! (out) optional
            & )
       
       xyz_TempN(:,:,1) = xyz_TempN(:,:,1) + &
            & (xy_SfcEngyFlxModRecv(:,:) - 0d0) / (xyr_Press(:,:,0) - xyr_Press(:,:,1)) &
            & * Grav / CpDry
    end if
    
  end subroutine dcpam_UpdateSurfaceProperties

  subroutine dcpam_StoreAtmSurfFlxInfo()

    use constants, only: LatentHeat, CpDry

    use composition, only: IndexH2OVap
    
    use saturate, only: &
      & xy_CalcQVapSatOnLiq,       &
      & xy_CalcQVapSatOnSol,       &
      & xy_CalcDQVapSatDTempOnLiq, &
      & xy_CalcDQVapSatDTempOnSol

    use intavr_operate
    use mpi_wrapper
    
    real(DP), dimension(0:iMax-1,jMax) :: &
         & xy_SurfQVapSatOnLiq, xy_SurfQVapSatOnSol, xy_SurfQVapSat, &
         & xy_SurfDQVapSatDTempOnLiq, xy_SurfDQVapSatDTempOnSol, xy_SurfDQVapSatDTemp

    real(DP) :: delta_t
    real(DP) :: xy_Tmp(0:iMax-1,jMax)
    real(DP) :: avr_QMixFlx, avr_QMixFlx2, avr_Tmp
    real(DP) :: avr_LUW, avr_LDW, avr_SUW, avr_SDW, avr_Sens, avr_Lat
    
    xy_SurfQVapSatOnLiq  = &
      & xy_CalcQVapSatOnLiq( xy_SurfTemp, xyr_Press(:,:,0) )
    xy_SurfQVapSatOnSol  = &
      & xy_CalcQVapSatOnSol( xy_SurfTemp, xyr_Press(:,:,0) )
    xy_SurfQVapSat       = &
      &   ( 1.0_DP - xy_SnowFrac ) * xy_SurfQVapSatOnLiq &
      & + xy_SnowFrac              * xy_SurfQVapSatOnSol
    xy_SurfDQVapSatDTempOnLiq = &
      & xy_CalcDQVapSatDTempOnLiq( xy_SurfTemp, xy_SurfQVapSatOnLiq )
    xy_SurfDQVapSatDTempOnSol = &
      & xy_CalcDQVapSatDTempOnSol( xy_SurfTemp, xy_SurfQVapSatOnSol )
    xy_SurfDQVapSatDTemp = &
      &   ( 1.0_DP - xy_SnowFrac ) * xy_SurfDQVapSatDTempOnLiq &
      & + xy_SnowFrac              * xy_SurfDQVapSatDTempOnSol

    !
    delta_t = DelTime

    
    !$omp parallel 
    !$omp workshare
    xy_TauXAtm(:,:) = xy_SurfMomFluxX - xy_SurfVelTransCoef*xyz_DUDt(:,:,1)*2d0*delta_t
    xy_TauYAtm(:,:) = xy_SurfMomFluxY - xy_SurfVelTransCoef*xyz_DVDt(:,:,1)*2d0*delta_t

    xy_SensAtm(:,:) = xyr_HeatFlux(:,:,0) - &
            CpDry*xyr_Exner(:,:,0)*xy_SurfTempTransCoef &
          * (xyz_DTempDtVDiff(:,:,1)/xyz_Exner(:,:,1) - xy_DSurfTempDt/xyr_Exner(:,:,0)) &
          * 2d0*delta_t
    xy_LatentAtm(:,:) = LatentHeat*( &
         & xyrf_QMixFlux(:,:,0,IndexH2OVap)  &
         &   - xy_SurfHumidCoef*xy_SurfQVapTransCoef*( &
         &     xyzf_DQMixDt(:,:,1,IndexH2OVap) - xy_SurfDQVapSatDTemp*xy_DSurfTempDt &
         &   )* 2d0*delta_t &
         & )
!!$    xy_LatentAtm(:,:) = xy_SurfLatentHeatFluxA
    
    xy_LDWRFlxAtm(:,:) = xyr_RadLDwFlux(:,:,0) + 2d0*delta_t*( &
         &    xy_DSurfTempDt * xyra_DelRadLDwFlux(:,:,0,0)            &
         & +  xyz_DTempDtVDiff(:,:,1) * xyra_DelRadLDwFlux(:,:,0,1)   &
         & )
    xy_LUWRFlxAtm(:,:) = xyr_RadLUwFlux(:,:,0) + 2d0*delta_t*( &
         &    xy_DSurfTempDt * xyra_DelRadLUwFlux(:,:,0,0)            &
         & +  xyz_DTempDtVDiff(:,:,1) * xyra_DelRadLUwFlux(:,:,0,1)   &
         & )

    xy_SDWRFlxAtm(:,:) = xyr_RadSDwFlux(:,:,0)
    xy_SUWRFlxAtm(:,:) = xyr_RadSUwFlux(:,:,0)

    
    !
    xy_SurfAirTemp(:,:) = xyr_Exner(:,:,0)/xyz_Exner(:,:,1)*xyz_TempN(:,:,1)
    xy_DSurfLatentFlxDTs(:,:) = LatentHeat*xy_SurfHumidCoef*xy_SurfQVapTransCoef*xy_SurfDQVapSatDTemp
    xy_DSurfHFlxDTs(:,:) = &
         &   CpDry*xy_SurfTempTransCoef      &
         & + xy_DSurfLatentFlxDTs            &
         & + 0d0*xyra_DelRadLUwFlux(:,:,0,0) &
         & - xyra_DelRadLDwFlux(:,:,0,0)
    !$omp end workshare
    !$omp end parallel


!!$    xy_Tmp = 1d0
!!$    avr_Tmp = IntLonLat_xy(xy_Rain + xy_Snow)
!!$    avr_QMixFlx = IntLonLat_xy(xyrf_QMixFlux(:,:,0,IndexH2OVap))
!!$    avr_QMixFlx2 = IntLonLat_xy(xy_LatentAtm)/LatentHeat
!!$    if (myrank == 0) then
!!$       write(*,*) "QMix net flux=",  avr_QMixFlx / (4d0*acos(-1d0)), avr_QMixFlx2 / (4d0*acos(-1d0))
!!$    end if
!!$    
!!$    avr_QMixFlx = IntLonLat_xy(xyrf_QMixFlux(:,:,0,IndexH2OVap) - xy_SurfHumidCoef*xy_SurfQVapTransCoef*( &
!!$         &     xyzf_DQMixDt(:,:,1,IndexH2OVap) )*2d0*delta_t)
!!$    if (myrank == 0) then
!!$       write(*,*) "QMix BudgetA=", avr_QMixFlx * LatentHeat / (4d0*acos(-1d0)), &
!!$            avr_QMixFlx2 / (4d0*acos(-1d0))            
!!$    end if
!!$
!!$    xy_Tmp = &
!!$         - ( xy_LDWRFlxAtm - xy_LUWRFlxAtm) &
!!$         - ( xy_SDWRFlxAtm - xy_SUWRFlxAtm) &
!!$         + xy_SensAtm + xy_LatentAtm
!!$    avr_Tmp = IntLonLat_xy(xy_Tmp)
!!$    if (myrank == 0) then
!!$       write(*,*) "SurfHFlxO=", avr_Tmp/(4d0*acos(-1d0))
!!$    end if
    
  end subroutine dcpam_StoreAtmSurfFlxInfo

  !------------------------------------------------------------------
  
  subroutine dcpam_advance_timestep(tstep, end_loop_flag, skip_flag)

    integer, intent(in) :: tstep
    logical, intent(inout) :: end_loop_flag
    logical, intent(in) :: skip_flag

    integer:: n                 ! 組成方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in dimension of constituents

    integer:: i
    integer:: j
    integer:: k

    real(DP) :: xy_Tmp1(0:imax-1,jmax)
    
    ! 時間積分
    ! Time integration
    !
    !loop_time : do while ( TimeB < EndTime )
    end_loop_flag = (TimeB >= EndTime)
    if(end_loop_flag) return

    call ProfUtil_RapStart('TimeLoop', 0)
    
!    write(*,*) "ATM :advance_timestep..", tstep, TimeN
if (.not. skip_flag) then    
   
    ! 地表面高度の設定
    ! Set surface height
    !
    call SetSurfaceProperties( &
      & xy_SurfHeight = xy_SurfHeight  & ! (inout) optional
      & )


    select case ( IDPhysMode )
    case ( IDPhysModeNoPhysics )

      xyz_DUDt     = 0.0_DP
      xyz_DVDt     = 0.0_DP
      xyz_DTempDt  = 0.0_DP
      xy_DPsDt     = 0.0_DP
      xyzf_DQMixDt = 0.0_DP

    case ( IDPhysModeHS94 )

      ! Held and Suarez (1994) による強制と散逸
      ! Forcing and dissipation suggested by Held and Suarez (1994)
      !
      call HS94Forcing( &
        & xyz_UB, xyz_VB, xyz_TempB, xy_PsB, & ! (in)
        & xyz_DUDt, xyz_DVDt, xyz_DTempDt )    ! (out)
      xy_DPsDt     = 0.0_DP
      xyzf_DQMixDt = 0.0_DP

    case ( IDPhysModeVenusSimple )

      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level,
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsB, xyz_TempB, xyzf_QMixB(:,:,:,IndexH2OVap),    & ! (in )
        & xyr_Temp, xyz_VirTemp, xyr_VirTemp,                  & ! (out) optional
        & xyr_Press     = xyr_Press,                           & ! (out) optional
        & xyz_Press     = xyz_Press,                           & ! (out) optional
        & xy_SurfHeight = xy_SurfHeight,                       & ! (in ) optional
        & xyz_Height    = xyz_Height,                          & ! (out) optional
        & xyr_Height    = xyr_Height,                          & ! (out) optional
        & xyz_Exner     = xyz_Exner,  xyr_Exner  = xyr_Exner   & ! (out) optional
        & )

      ! 簡単金星計算のための強制
      ! forcing for simple Venus calculation
      !
      call YT2003Forcing(                                      &
        & xy_SurfHeight,                                       & ! (in )
        & xyz_UB, xyz_VB, xyz_TempB, xyz_VirTemp, xyr_VirTemp, & ! (in )
        & xy_PsB, xyz_Press, xyr_Press, xyr_Temp,              & ! (in )
        & xyz_Height, xyr_Height, xyz_Exner, xyr_Exner,        & ! (in )
        & xyz_DUDt, xyz_DVDt, xyz_DTempDt                      & ! (out)
        & )

      xy_DPsDt     = 0.0_DP
      xyzf_DQMixDt = 0.0_DP


    case ( IDPhysModeJupiterSimple )

      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level,
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsB, xyz_TempB, xyzf_QMixB(:,:,:,IndexH2OVap),  & ! (in )
        & xyr_VirTemp   = xyr_VirTemp,                       & ! (out) optional
        & xyr_Press     = xyr_Press,                         & ! (out) optional
        & xyz_Press     = xyz_Press,                         & ! (out) optional
        & xyz_Exner     = xyz_Exner,                         & ! (out) optional
        & xy_SurfHeight = xy_SurfHeight,                     & ! (out) optional
        & xyz_Height    = xyz_Height                         & ! (out) optional
        & )


      ! Schneider and Liu (2009) による鉛直混合課程
      ! Vertical diffusion by Schneider and Liu (2009)
      !
      call SL09Diffusion(                                     &
        & xy_SurfHeight, xyz_Height,                          &
        & xyz_UB, xyz_VB, xyzf_QMixB, xyr_Press, xyr_VirTemp, & ! (in)
        & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt  & ! (out)
        & )

      ! Schneider and Liu (2009) の放射モデル
      ! Radiation model by Schneider and Liu (2009)
      !
      call RadSL09Flux(                           &
        & xyr_Press, xyz_Press, xyz_TempB,        & ! (in)
        & xyr_RadSUwFlux, xyr_RadSDwFlux,         & ! (out)
        & xyr_RadLUwFlux, xyr_RadLDwFlux,         & ! (out)
        & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux  & ! (out)
        & )

      ! Net flux is calculated.
      !
      xyr_RadSFlux     = xyr_RadSUwFlux     - xyr_RadSDwFlux
      xyr_RadLFlux     = xyr_RadLUwFlux     - xyr_RadLDwFlux
      xyra_DelRadLFlux = xyra_DelRadLUwFlux - xyra_DelRadLDwFlux


      ! 放射による温度変化率
      ! Temperature tendency with radiation
      !
      call RadDTempDt(                           &
        & xyr_RadLFlux, xyr_RadSFlux, xyr_Press, &   ! (in)
        & xyz_DTempDtRadL, xyz_DTempDtRadS )          ! (out)
      xyr_RadLDwFlux     = 0.0_DP
      xyra_DelRadLDwFlux = 0.0_DP

      ! 非断熱加熱率の総和の計算
      ! Sum all diabatic heating rates
      !
      xyz_DTempDt = xyz_DTempDtVDiff + xyz_DTempDtRadL + xyz_DTempDtRadS

      xy_DPsDt     = 0.0_DP


    case ( IDPhysModeJupiterSimpleV2 )

      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level,
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsB, xyz_TempB, xyzf_QMixB(:,:,:,IndexH2OVap),  & ! (in )
        & xyr_VirTemp   = xyr_VirTemp,                       & ! (out) optional
        & xyr_Press     = xyr_Press,                         & ! (out) optional
        & xyz_Press     = xyz_Press,                         & ! (out) optional
        & xy_SurfHeight = xy_SurfHeight,                     & ! (out) optional
        & xyz_Height    = xyz_Height,                        & ! (out) optional
        & xyz_Exner     = xyz_Exner,  xyr_Exner  = xyr_Exner & ! (out) optional
        & )


      xyr_MomFluxX (:,:,1:kmax)   = 0.0_DP
      xyr_MomFluxY (:,:,1:kmax)   = 0.0_DP
      xyr_HeatFlux (:,:,1:kmax)   = 0.0_DP
      xyrf_QMixFlux(:,:,1:kmax,:) = 0.0_DP

      xyr_VelDiffCoef  = 0.0_DP
      xyr_TempDiffCoef = 0.0_DP
      xyr_QMixDiffCoef = 0.0_DP

      ! 下部境界フラックス
      ! Lower boundary flux
      !
      call LBFluxSimple( &
        & xyz_UB, xyz_VB, xyz_TempB, xyr_VirTemp, xyzf_QMixB,         & ! (in)
        & xyr_Press, xy_SurfHeight, xyz_Height, xyz_Exner, xyr_Exner, & ! (in)
        & xyr_MomFluxX(:,:,0:0), xyr_MomFluxY(:,:,0:0), xyr_HeatFlux(:,:,0:0), xyrf_QMixFlux(:,:,0:0,:),    & ! (out)
        & xy_SurfVelTransCoef, xy_SurfTempTransCoef,                  & ! (out)
        & xy_SurfQVapTransCoef                                        & ! (out)
        & )

      call PhyImplAtmOnlyTendency(                                 &
        & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
        & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
        & xyr_VirTemp, xyz_Height,                                 & ! (in)
        & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
        & xy_SurfVelTransCoef, xy_SurfTempTransCoef,               & ! (in)
        & xy_SurfQVapTransCoef,                                    & ! (in)
        & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt       & ! (out)
        & )

      xy_DPsDt             = 0.0_DP

      ! Schneider and Liu (2009) の放射モデル
      ! Radiation model by Schneider and Liu (2009)
      !
      call RadSL09Flux(                           &
        & xyr_Press, xyz_Press, xyz_TempB,        & ! (in)
        & xyr_RadSUwFlux, xyr_RadSDwFlux,         & ! (out)
        & xyr_RadLUwFlux, xyr_RadLDwFlux,         & ! (out)
        & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux  & ! (out)
        & )

      ! Net flux is calculated.
      !
      xyr_RadSFlux     = xyr_RadSUwFlux     - xyr_RadSDwFlux
      xyr_RadLFlux     = xyr_RadLUwFlux     - xyr_RadLDwFlux
      xyra_DelRadLFlux = xyra_DelRadLUwFlux - xyra_DelRadLDwFlux

      call ProfUtil_RapEnd('Rad', 1)

      ! 放射による温度変化率
      ! Temperature tendency with radiation
      !
      call RadDTempDt(                           &
        & xyr_RadLFlux, xyr_RadSFlux, xyr_Press, &   ! (in)
        & xyz_DTempDtRadL, xyz_DTempDtRadS )          ! (out)
      xyr_RadLDwFlux     = 0.0_DP
      xyra_DelRadLDwFlux = 0.0_DP

      ! 非断熱加熱率の総和の計算
      ! Sum all diabatic heating rates
      !
      xyz_DTempDt = xyz_DTempDtVDiff + xyz_DTempDtRadL + xyz_DTempDtRadS


    case ( IDPhysModeFullPhysics )

      ! 地表面条件の設定
      ! Configure surface conditions
      !
      call SetSurfaceProperties(                           &  
        & xy_SurfMajCompIceB, xy_SoilMoistB,            & ! (in)    optional
        & xy_SurfSnowB,                      & ! (in)    optional
!!$        & xy_SurfTemp            = xy_SurfTemp,            & ! (inout) optional
!!$        & xy_SurfAlbedo          = xy_SurfAlbedo,          & ! (inout) optional
        & xy_SurfHumidCoef       = xy_SurfHumidCoef,       & ! (inout) optional
        & xy_SurfRoughLenMom     = xy_SurfRoughLenMom,     & ! (inout) optional
        & xy_SurfRoughLenHeat    = xy_SurfRoughLenHeat,    & ! (inout) optional
        & xy_SurfHeatCapacity    = xy_SurfHeatCapacity,    & ! (inout) optional
        & xy_DeepSubSurfHeatFlux = xy_DeepSubSurfHeatFlux, & ! (inout) optional
        & xy_SurfCond            = xy_SurfCond,            & ! (inout) optional
        & xy_SurfType            = xy_SurfType,            & ! (inout) optional
        & xy_SeaIceConc          = xy_SeaIceConc,          & ! (inout) optional
        & xy_SoilHeatCap         = xy_SoilHeatCap,         & ! (inout) optional
        & xy_SoilHeatDiffCoef    = xy_SoilHeatDiffCoef,    & ! (inout) optional
        & xy_SurfHeightStd       = xy_SurfHeightStd        & ! (inout) optional
        & )

      ! 雪, 氷の割合
      ! snow/ice fraction
      !
      call CalcSnowFrac(       &
!!$        & xy_FlagLand, xy_SurfSnow,  & ! (in )
        & xy_SurfSnowB,              & ! (in )
        & xy_SnowFrac                & ! (out)
        & )

      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level, 
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsB, xyz_TempB, xyzf_QMixB(:,:,:,IndexH2OVap),    & ! (in )
        & xyr_Temp, xyz_VirTemp, xyr_VirTemp, xy_SurfVirTemp,  & ! (out) optional
        & xyr_Press     = xyr_Press,                           & ! (out) optional
        & xyz_Press     = xyz_Press,                           & ! (out) optional
        & xy_SurfHeight = xy_SurfHeight, xy_SurfTemp = xy_SurfTemp, & ! (in ) optional
        & xyz_Height    = xyz_Height, xyr_Height = xyr_Height, & ! (out) optional
        & xyz_Exner     = xyz_Exner,  xyr_Exner  = xyr_Exner   & ! (out) optional
        & )


      call ProfUtil_RapStart('Rad', 1)      
      select case ( IDRadMethod )
      case ( IDRadMethodDennouAGCM )

        ! 放射フラックス (GFD 電脳倶楽部開発の放射モデル)
        ! Radiation flux (radiation model developed by GFD Dennou Club)
        !
        call RadDennouAGCMFlux(                                  &
          & xyz_TempB, xyzf_QMixB(:,:,:,IndexH2OVap), xyr_Press, &   ! (in)
          & xy_SurfTemp, xy_SurfAlbedo,      &   ! (in)
          & xyr_RadSUwFlux, xyr_RadSDwFlux,            & ! (out)
          & xyr_RadLUwFlux, xyr_RadLDwFlux,            & ! (out)
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux     & ! (out)
          & )


      case ( IDRadMethodEarthV2 )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        ! 
        ! Cloud model
        !
        select case ( IDCloudMethod )
        case ( IDCloudMethodNone )
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = 0.0_DP
        case ( IDCloudMethodSimple )
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = 0.0_DP
        case ( IDCloudMethodSimpleWithIce )
          if ( IndexH2OSol <= 0 ) then
            call MessageNotify( 'E', prog_name, &
              & trim(a_QMixName(IndexH2OSol))//' is not found.' )
          end if
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = xyzf_QMixB(:,:,:,IndexH2OSol)
        case ( IDCloudMethodT1993WithIce )
          if ( IndexH2OSol <= 0 ) then
            call MessageNotify( 'E', prog_name, &
              & trim(a_QMixName(IndexH2OSol))//' is not found.' )
          end if
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = xyzf_QMixB(:,:,:,IndexH2OSol)
        end select

!!$        call CloudSimpleDivideWatAndIce(                           &
!!$          & xyz_TempB,                                             & ! (in )
!!$          & xyzf_QMixB(:,:,:,IndexH2OLiq),                         & ! (in )
!!$          & xyz_QH2OLiqforRad, xyz_QH2OSolforRad                   & ! (out)
!!$          & )

        ! Cloud cover is calculated.
        !
        select case ( IDCloudMethod )
        case ( IDCloudMethodNone )
          xyz_CloudCoverforRad = 1.0_DP
        case ( IDCloudMethodSimple )
          call CloudSimpleCalcCloudCover(                            &
            & xyz_Press, xyz_TempB,                                  & ! (in )
!!$            & xyzf_QMixB(:,:,:,IndexH2OVap)                       &
!!$            & + xyzf_QMixB(:,:,:,IndexH2OLiq),                    & ! (in )
            & xyzf_QMixB(:,:,:,IndexH2OVap),                         & ! (in )
            & xyz_CloudCoverforRad                                   & ! (out)
            & )
        case ( IDCloudMethodSimpleWithIce )
          call CloudSimpleCalcCloudCover(                            &
            & xyz_Press, xyz_TempB,                                  & ! (in )
!!$            & xyzf_QMixB(:,:,:,IndexH2OVap)      &
!!$            & + xyzf_QMixB(:,:,:,IndexH2OLiq)    & ! (in )
!!$            & + xyzf_QMixB(:,:,:,IndexH2OSol),   & ! (in )
            & xyzf_QMixB(:,:,:,IndexH2OVap),                         & ! (in )
            & xyz_CloudCoverforRad                                   & ! (out)
            & )
        case ( IDCloudMethodT1993WithIce )
          xyz_CloudCoverforRad = xyzf_QMixB(:,:,:,CompositionInqIndex( 'CloudCover' ))
        end select

        call RadEarthV2Flux(                              &
          & xy_SurfAlbedo,                                &
          & xyz_Press, xyr_Press, xyz_TempB,              &
          & xyzf_QMixB(:,:,:,IndexH2OVap),                &
          & xyz_QH2OLiqforRad, xyz_QH2OSolforRad,         &
          & xyz_CloudCoverforRad,                         &
          & xy_SurfTemp,                                  &
          & xyr_RadSUwFlux, xyr_RadSDwFlux,               &
          & xyr_RadLUwFlux, xyr_RadLDwFlux,               & ! (out)
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux        & ! (out)
          & )

      case ( IDRadMethodRRTMG )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        ! 
        ! Cloud model
        !
        select case ( IDCloudMethod )
        case ( IDCloudMethodNone )
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = 0.0_DP
        case ( IDCloudMethodSimple )
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = 0.0_DP
        case ( IDCloudMethodSimpleWithIce )
          if ( IndexH2OSol <= 0 ) then
            call MessageNotify( 'E', prog_name, &
              & trim(a_QMixName(IndexH2OSol))//' is not found.' )
          end if
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = xyzf_QMixB(:,:,:,IndexH2OSol)
        case ( IDCloudMethodT1993WithIce )
          if ( IndexH2OSol <= 0 ) then
            call MessageNotify( 'E', prog_name, &
              & trim(a_QMixName(IndexH2OSol))//' is not found.' )
          end if
          xyz_QH2OLiqforRad = xyzf_QMixB(:,:,:,IndexH2OLiq)
          xyz_QH2OSolforRad = xyzf_QMixB(:,:,:,IndexH2OSol)
        end select

        ! Cloud cover is calculated.
        !
        select case ( IDCloudMethod )
        case ( IDCloudMethodNone )
          xyz_CloudCoverforRad = 1.0_DP
        case ( IDCloudMethodSimple )
          call CloudSimpleCalcCloudCover(                            &
            & xyz_Press, xyz_TempB,                                  & ! (in )
!!$            & xyzf_QMixB(:,:,:,IndexH2OVap)                          &
!!$            & + xyzf_QMixB(:,:,:,IndexH2OLiq),   & ! (in )
            & xyzf_QMixB(:,:,:,IndexH2OVap),                         & ! (in )
            & xyz_CloudCoverforRad                                   & ! (out)
            & )
        case ( IDCloudMethodSimpleWithIce )
          call CloudSimpleCalcCloudCover(                            &
            & xyz_Press, xyz_TempB,                                  & ! (in )
!!$            & xyzf_QMixB(:,:,:,IndexH2OVap)                          &
!!$            & + xyzf_QMixB(:,:,:,IndexH2OLiq)    & ! (in )
!!$            & + xyzf_QMixB(:,:,:,IndexH2OSol),   & ! (in )
            & xyzf_QMixB(:,:,:,IndexH2OVap),                         & ! (in )
            & xyz_CloudCoverforRad                                   & ! (out)
            & )
        case ( IDCloudMethodT1993WithIce )
          xyz_CloudCoverforRad = xyzf_QMixB(:,:,:,CompositionInqIndex( 'CloudCover' ))
        end select


        ! wrapper of RRTMG
        ! wrapper of RRTMG
        !
        call RadRRTMGWrapperFlux(               &
          & xy_SurfAlbedo,                                &
          & xyz_Press, xyr_Press, xyz_TempB, xyr_Temp,    &
          & xyzf_QMixB(:,:,:,IndexH2OVap),                &
          & xyz_QH2OLiqforRad, xyz_QH2OSolforRad,         &
          & xyz_CloudCoverforRad,                         &
          & xy_SurfTemp,                                  &
          & xyr_RadSUwFlux, xyr_RadSDwFlux,               &
          & xyr_RadLUwFlux, xyr_RadLDwFlux,               & ! (out)
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux        & ! (out)
          & )

      case ( IDRadMethodMarsV1 )

        if ( CompositionInqIndex( 'QDust' ) <= 0 ) then
          call MessageNotify( 'E', prog_name, 'QDust is not found.' )
        end if

        call RadMarsV1Flux(                                          &
          & xy_SurfType, xy_SurfMajCompIceB,                         &
          & xy_SurfAlbedo,                                           &
          & xyz_Press, xyr_Press, xyz_TempB, xyr_Temp, xy_SurfTemp,  &
          & xyzf_QMixB(:,:,:,CompositionInqIndex( 'QDust' )),        & ! (in )
          & xyr_RadSUwFlux, xyr_RadSDwFlux,                          &
          & xyr_RadLUwFlux, xyr_RadLDwFlux,                          &
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux                   &
          & )

      case ( IDRadMethodSL09 )

        ! Schneider and Liu (2009) の放射モデル
        ! Radiation model by Schneider and Liu (2009)
        !
        call RadSL09Flux(                           &
          & xyr_Press, xyz_Press, xyz_TempB,        & ! (in)
          & xyr_RadSUwFlux, xyr_RadSDwFlux,         & ! (out)
          & xyr_RadLUwFlux, xyr_RadLDwFlux,         & ! (out)
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux  & ! (out)
          & )

      case ( IDRadMethodSimple )

        ! 簡単放射モデル
        ! Simple radiation model
        !
        call RadSimpleFlux(                                              &
          & xy_SurfAlbedo, xy_SurfTemp, xyr_Press, xyz_Press, xyz_TempB, &
          & xyzf_QMixB(:,:,:,IndexH2OVap),                               &
          & xyr_RadSUwFlux, xyr_RadSDwFlux,                              &
          & xyr_RadLUwFlux, xyr_RadLDwFlux,                              &
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux                       &
          & )

      case ( IDRadMethodNone )

        ! 何もしない放射モデル
        ! radiation model with no absorption and no scattering
        !
        call RadNoneFlux(                                      &
          & xy_SurfAlbedo,                                           & ! (in)
          & xy_SurfTemp,                                             & ! (in)
          & xyr_RadSUwFlux, xyr_RadSDwFlux,                          & ! (out)
          & xyr_RadLUwFlux, xyr_RadLDwFlux,                          & ! (out)
          & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux                   & ! (out)
          & )

      end select

      ! Net flux is calculated.
      !
      xyr_RadSFlux     = xyr_RadSUwFlux     - xyr_RadSDwFlux
      xyr_RadLFlux     = xyr_RadLUwFlux     - xyr_RadLDwFlux
      xyra_DelRadLFlux = xyra_DelRadLUwFlux - xyra_DelRadLDwFlux

      call ProfUtil_RapEnd('Rad', 1)

      call ProfUtil_RapStart('SfcFlx', 1)      
      ! 地表面フラックス
      ! Surface flux
      !
      select case ( IDSfcFluxMethod )
      case ( IDSfcFluxMethodL82 )

        call SurfaceFlux(                                                  &
          & 'L82',                                                         & ! (in)
          & xyz_UB, xyz_VB,                                                & ! (in)
          & xyz_TempB, xyr_Temp, xyz_VirTemp, xyr_VirTemp, xy_SurfVirTemp, & ! (in)
          & xyzf_QMixB,                                                    & ! (in)
          & xyr_Press, xy_SurfHeight, xyz_Height, xyz_Exner, xyr_Exner,   & ! (in)
          & xy_SurfTemp, xy_SurfHumidCoef,                                & ! (in)
          & xy_SurfRoughLenMom, xy_SurfRoughLenHeat,                      & ! (in)
          & xy_SnowFrac,                                                  & ! (in)
          & xy_SurfMomFluxX, xy_SurfMomFluxY, xy_SurfHeatFlux, xyf_SurfQMixFlux, & ! (out)
          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,                    & ! (out)
          & xy_SurfQVapTransCoef,                                         & ! (out)
          & xy_SurfMOLength                                               & ! (out)
          & )

      case ( IDSfcFluxMethodBH91B94 )

        call SurfaceFlux(                                                 &
          & 'BH91B94',                                                    & ! (in)
          & xyz_UB, xyz_VB,                                               & ! (in)
          & xyz_TempB, xyr_Temp, xyz_VirTemp, xyr_VirTemp, xy_SurfVirTemp, & ! (in)
          & xyzf_QMixB,                                                    & ! (in)
          & xyr_Press, xy_SurfHeight, xyz_Height, xyz_Exner, xyr_Exner,   & ! (in)
          & xy_SurfTemp, xy_SurfHumidCoef,                                & ! (in)
          & xy_SurfRoughLenMom, xy_SurfRoughLenHeat,                      & ! (in)
          & xy_SnowFrac,                                                  & ! (in)
          & xy_SurfMomFluxX, xy_SurfMomFluxY, xy_SurfHeatFlux, xyf_SurfQMixFlux, & ! (out)
          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,                    & ! (out)
          & xy_SurfQVapTransCoef,                                         & ! (out)
          & xy_SurfMOLength                                               & ! (out)
          & )

      end select
      !
      ! set dust flux
      !   This is ad hoc treatment now (yot, 2013/09/28)
      !
!!$      if ( CompositionInqIndex('QDust') > 0 ) then
!!$        call SetDustFlux(                                       &
!!$          & xyf_SurfQMixFlux(:,:,CompositionInqIndex('QDust'))  & ! (out)
!!$          & )
!!$      end if

      call ProfUtil_RapEnd('SfcFlx', 1)      

      call ProfUtil_RapStart('VDiffMY', 1)      
      ! 鉛直拡散フラックス
      ! Vertical diffusion flux
      !
      select case ( IDVDiffMethod )
      case ( IDVDiffMethodMY2 )

        call VDiffusion(                                              &
          & xyz_UB,     xyz_VB,     xyzf_QMixB,                       & ! (in)
          & xyz_TempB, xyr_Temp, xyz_VirTemp, xyr_VirTemp, xyr_Press, & ! (in)
          & xy_SurfHeight,                                            & ! (in)
          & xyz_Height, xyr_Height, xyz_Exner,    xyr_Exner,          & ! (in)
          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux,  & ! (out)
          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef       & ! (out)
          & )

      case ( IDVDiffMethodMY25 )

        if ( IndexTKE <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexTKE))//' is not found.' )
        end if

        call VDiffusionMY25(                                          &
!        call VDiffusionMY251DWrapper3D(                               &
!        call VDiffusionMY25GBT94(                                     &
          & xyz_UB,     xyz_VB,     xyzf_QMixB,                       & ! (in)
          & xyz_TempB, xyr_Temp, xyz_VirTemp, xyr_VirTemp, xyr_Press, & ! (in)
          & xy_SurfHeight,                                            & ! (in)
          & xyz_Height, xyr_Height, xyz_Exner, xyr_Exner,             & ! (in)
          & xyzf_QMixB(:,:,:,IndexTKE),                               & ! (in)
          & xy_SurfMomFluxX, xy_SurfMomFluxY,                         & ! (in)
          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux,  & ! (out)
          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,      & ! (out)
          & xyz_DTurKinEneDt                                          & ! (out)
          & )

      case ( IDVDiffMethodJMA )

        if ( IndexTKE <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexTKE))//' is not found.' )
        end if

        ! JMA 乱流混合モジュール
        ! JMA turbulent mixing module
        !
        call VDiffusionJMAMYWrapper3D(                          &
          & xyz_UB,     xyz_VB,     xyzf_QMixB,                       & ! (in)
          & xyz_TempB, xyr_Temp, xyz_VirTemp, xyr_VirTemp,            & ! (in)
          & xyz_Press, xyr_Press,                                     & ! (in)
          & xy_SurfHeight,                                            & ! (in)
          & xyz_Height, xyr_Height, xyz_Exner, xyr_Exner,             & ! (in)
          & xy_SurfMOLength,                                          & ! (in)
          & xyzf_QMixB(:,:,:,IndexTKE),                               & ! (in)
          & xy_SurfMomFluxX, xy_SurfMomFluxY,                         & ! (in)
          & xy_SurfHeatFlux, xyf_SurfQMixFlux,                        & ! (in)
          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux,  & ! (out)
          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,      & ! (out)
          & xyz_DTurKinEneDt                                          & ! (out)
          & )

      end select

      xyr_MomFluxX (:,:,0)   = xy_SurfMomFluxX
      xyr_MomFluxY (:,:,0)   = xy_SurfMomFluxY
      xyr_HeatFlux (:,:,0)   = xy_SurfHeatFlux
      xyrf_QMixFlux(:,:,0,:) = xyf_SurfQMixFlux
      call ProfUtil_RapEnd('VDiffMY', 1)


      call ProfUtil_RapStart('PhysImplSDH', 1)      
      ! 一部の物理過程の時間変化率の計算 (陰解法)
      ! Calculate tendency by a part of physical processes (implicit)
      !
      select case ( IDPhyTendMethod )
      case ( IDPhyTendMethodImp1LayModel )

        call PhyImplTendency(                                        &
          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
          & xyr_RadSFlux, xyr_RadLFlux,                              & ! (in)
          & xy_DeepSubSurfHeatFlux,                                  & ! (in)
          & xy_SurfTemp, xy_SurfHumidCoef, xy_SurfCond,              & ! (in)
          & xy_SurfHeatCapacity,                                     & ! (in)
          & xyra_DelRadLFlux,                                        & ! (in)
          & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
          & xyr_VirTemp, xyz_Height,                                 & ! (in)
          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,               & ! (in)
          & xy_SurfQVapTransCoef,                                    & ! (in)
          & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt,      & ! (out)
          & xy_DSurfTempDt                                           & ! (out)
          & )

        xy_SurfH2OVapFluxA     = 0.0_DP
        xy_SurfLatentHeatFluxA = 0.0_DP

        xyz_DSoilTempDt      = 0.0_DP
        xy_DPsDt             = 0.0_DP
        xy_DSurfMajCompIceDt = 0.0_DP
        xy_DSoilMoistDt      = 0.0_DP
        xy_DSurfSnowDt       = 0.0_DP

      case ( IDPhyTendMethodImpSoilModel )


        ! 地下における熱の鉛直拡散
        ! Vertical diffusion of heat under the ground
        !
!!$        call SubsurfaceDiffusionFlagLandMatthews(  &
!!$          & xy_SurfType,                           & ! (in)
!!$          & xy_FlagMatthewsLand                    & ! (out)
!!$          & )
        call SubsurfaceDiffusion(                    &
          & xy_DeepSubSurfHeatFlux,                  &          ! (in )
          & xy_SoilHeatCap, xy_SoilHeatDiffCoef,     &          ! (in )
          & xy_SurfTemp, xyz_SoilTemp, xy_SurfSnowB, &          ! (in )
          & xyr_SoilTempTransCoef, xyr_SoilHeatFlux  &          ! (out)
          & )


        ! This value is not correct, if snow evaporates. 
        ! If a bucket model is used, this surface latent heat flux is 
        ! corrected in BucketModEvapAndLatentHeatFlux by the use of surface 
        ! moisture and surface snow amount. 
        !
!!$        xy_SurfLatentHeatFluxA = LatentHeat * xyrf_QMixFlux(:,:,0,IndexH2OVap)


!!$        select case ( IDSfcMoistMethod )
!!$        case ( IDSfcMoistMethodBucket )
!!$          ! バケツモデルのための地表面フラックス修正
!!$          ! Modification of surface flux for bucket model
!!$          !
!!$          call BucketSetFlagOceanFromMatthews( &
!!$            & xy_SurfType,                     & ! (in)
!!$            & xy_BucketFlagOceanGrid           & ! (out)
!!$            & )
!!$          call BucketModEvapAndLatentHeatFlux(                         &
!!$            & xy_BucketFlagOceanGrid, xy_SoilMoistB, xy_SurfSnowB,     & ! (in   )
!!$            & xyrf_QMixFlux(:,:,0,IndexH2OVap), xy_SurfLatentHeatFluxA & ! (inout)
!!$            & )
!!$        end select



!!$        call PhyImplSDHSetMethodFromMatthews(  &
!!$          & xy_SurfType, xy_SeaIceConc,        & ! (in)
!!$          & xy_PhyImplSDHIndexCalcMethod       & ! (out)
!!$          & )
!!$        call PhyImplSDHTendency(                                      &
!!$          & xy_PhyImplSDHIndexCalcMethod,                             & ! (in)
!!$          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux,  & ! (in)
!!$          & xy_SurfLatentHeatFluxA,                                   & ! (in)
!!$          & xyr_SoilHeatFlux,                                         & ! (in)
!!$          & xyr_RadSFlux, xyr_RadLFlux,                               & ! (in)
!!$          & xy_DeepSubSurfHeatFlux,                                   & ! (in)
!!$          & xy_SurfTemp, xyz_SoilTemp,                                & ! (in)
!!$          & xy_SurfHumidCoef,                                         & ! (in)
!!$          & xy_SurfHeatCapacity,                                      & ! (in)
!!$          & xy_SoilHeatCap, xy_SoilHeatDiffCoef,                      & ! (in)
!!$          & xyra_DelRadLFlux,                                         & ! (in)
!!$          & xyr_Press, xyz_Exner, xyr_Exner,                          & ! (in)
!!$          & xyr_VelTransCoef, xyr_TempTransCoef,                      & ! (in)
!!$          & xyr_QMixTransCoef,                                        & ! (in)
!!$          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,                & ! (in)
!!$          & xy_SurfQVapTransCoef,                                     & ! (in)
!!$          & xyr_SoilTempTransCoef,                                    & ! (in)
!!$          & xy_SurfMajCompIceB,                                       & ! (in)
!!$          & xy_SurfSnowB,                                             & ! (in)
!!$          & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt,       & ! (out)
!!$          & xy_DSurfTempDt,                                           & ! (out)
!!$          & xyz_DSoilTempDt,                                          & ! (out)
!!$          & xy_DPsDt, xy_DSurfMajCompIceDt,                           & ! (out)
!!$          & xy_DSoilMoistDt,                                          & ! (out)
!!$          & xy_DSurfSnowDt                                            & ! (out)
!!$          & )

!!$        ! This is temporal treatment.
!!$        xy_SurfH2OVapFluxA = xyrf_QMixFlux(:,:,0,IndexH2OVap)

!!$        call PhyImplSDHV2SetMethodMatthews(  &
!!$          & xy_SurfType, xy_SeaIceConc,      & ! (in)
!!$          & xy_PhyImplSDHIndexCalcMethod     & ! (out)
!!$          & )
!!$        select case ( IDSfcMoistMethod )
!!$        case ( IDSfcMoistMethodBucket )
!!$          call BucketSetFlagOceanFromMatthews( &
!!$            & xy_SurfType,                     & ! (in)
!!$            & xy_BucketFlagOceanGrid           & ! (out)
!!$            & )
!!$        case default
!!$          xy_BucketFlagOceanGrid = .true.
!!$        end select
!!$        call PhyImplSDHV2Tendency(                                   &
!!$          & xy_PhyImplSDHIndexCalcMethod, xy_BucketFlagOceanGrid,    & ! (in)
!!$          & xy_SnowFrac,                                             & ! (in)
!!$          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
!!$          & xy_SurfH2OVapFluxA, xy_SurfLatentHeatFluxA,              & ! (out)
!!$          & xyr_SoilHeatFlux,                                        & ! (in)
!!$          & xyr_RadSFlux, xyr_RadLFlux,                              & ! (in)
!!$          & xy_DeepSubSurfHeatFlux,                                  & ! (in)
!!$          & xyz_TempB, xy_SurfTemp, xyz_SoilTemp,                    & ! (in)
!!$          & xyzf_QMixB,                                              & ! (in)
!!$          & xy_SurfHumidCoef,                                        & ! (in)
!!$          & xy_SurfHeatCapacity,                                     & ! (in)
!!$          & xy_SoilHeatCap, xy_SoilHeatDiffCoef,                     & ! (in)
!!$          & xyra_DelRadLFlux,                                        & ! (in)
!!$          & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
!!$          & xyr_VirTemp, xyz_Height,                                 & ! (in)
!!$          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
!!$          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,               & ! (in)
!!$          & xy_SurfQVapTransCoef,                                    & ! (in)
!!$          & xyr_SoilTempTransCoef,                                   & ! (in)
!!$          & xy_SurfMajCompIceB,                                      & ! (in)
!!$          & xy_SoilMoistB, xy_SurfSnowB,                             & ! (in)
!!$          & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt,      & ! (out)
!!$          & xy_DSurfTempDt,                                          & ! (out)
!!$          & xyz_DSoilTempDt,                                         & ! (out)
!!$          & xy_DPsDt, xy_DSurfMajCompIceDt,                          & ! (out)
!!$          & xy_DSoilMoistDt,                                         & ! (out)
!!$          & xy_DSurfSnowDt                                           & ! (out)
!!$          & )

        call PhyImplSDHV3SetMethodMatthews(          &
          & xy_SurfCond, xy_SurfType, xy_SeaIceConc, & ! (in)
          & xy_PhyImplSDHIndexCalcMethod             & ! (out)
          & )
        select case ( IDSfcMoistMethod )
        case ( IDSfcMoistMethodBucket )
          call BucketSetFlagOceanFromMatthews( &
            & xy_SurfType,                     & ! (in)
            & xy_BucketFlagOceanGrid           & ! (out)
            & )
        case default
          xy_BucketFlagOceanGrid = .true.
        end select
        call PhyImplSDHV3Tendency(                                   &
          & xy_PhyImplSDHIndexCalcMethod, xy_BucketFlagOceanGrid,    & ! (in)
          & xy_SnowFrac,                                             & ! (in)
          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
          & xy_SurfH2OVapFluxA, xy_SurfLatentHeatFluxA,              & ! (out)
          & xyr_SoilHeatFlux,                                        & ! (in)
          & xyr_RadSFlux, xyr_RadLFlux,                              & ! (in)
          & xy_DeepSubSurfHeatFlux,                                  & ! (in)
          & xyz_TempB, xy_SurfTemp, xyz_SoilTemp,                    & ! (in)
          & xyzf_QMixB,                                              & ! (in)
          & xy_SurfHumidCoef,                                        & ! (in)
!!$          & xy_SurfHeatCapacity,                                     & ! (in)
          & xy_SoilHeatCap, xy_SoilHeatDiffCoef,                     & ! (in)
          & xyra_DelRadLFlux,                                        & ! (in)
          & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
          & xyr_VirTemp, xyz_Height,                                 & ! (in)
          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,               & ! (in)
          & xy_SurfQVapTransCoef,                                    & ! (in)
          & xyr_SoilTempTransCoef,                                   & ! (in)
          & xy_SurfMajCompIceB,                                      & ! (in)
          & xy_SoilMoistB, xy_SurfSnowB,                             & ! (in)
          & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt,      & ! (out)
          & xy_DSurfTempDt,                                          & ! (out)
          & xyz_DSoilTempDt,                                         & ! (out)
          & xy_DPsDt, xy_DSurfMajCompIceDt,                          & ! (out)
          & xy_DSoilMoistDt,                                         & ! (out)
          & xy_DSurfSnowDt                                           & ! (out)
          & )

      case ( IDPhyTendMethodImpAtmOnly )

        call PhyImplAtmOnlyTendency(                           &
          & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
          & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
          & xyr_VirTemp, xyz_Height,                                 & ! (in)
          & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
          & xy_SurfVelTransCoef, xy_SurfTempTransCoef,               & ! (in)
          & xy_SurfQVapTransCoef,                                    & ! (in)
          & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt       & ! (out)
          & )

        xy_DSurfTempDt       = 0.0_DP
        xyz_DSoilTempDt      = 0.0_DP
        xy_DPsDt             = 0.0_DP
        xy_DSurfMajCompIceDt = 0.0_DP
        xy_DSoilMoistDt      = 0.0_DP
        xy_DSurfSnowDt       = 0.0_DP

        xy_SurfH2OVapFluxA     = 0.0_DP
        xy_SurfLatentHeatFluxA = 0.0_DP

      end select


      ! Overwrite tendency of turbulent kinetic energy
      !
      select case ( IDVDiffMethod )
      case ( IDVDiffMethodMY25, IDVDiffMethodJMA )
        if ( IndexTKE <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexTKE))//' is not found.' )
        end if
        xyzf_DQMixDt(:,:,:,IndexTKE) = &
          & xyz_DTurKinEneDt
      end select
      call ProfUtil_RapEnd('PhysImplSDH', 1)


      ! Gravity wave drag
      ! Gravity wave drag
      !
      select case ( IDGWDMethod )
      case ( IDGWDMethodNone )

        xyz_DUDtGWD = 0.0_DP
        xyz_DVDtGWD = 0.0_DP

      case ( IDGWDMethodM1987 )

        ! Gravity wave drag by McFarlane (1987)
        ! Gravity wave drag by McFarlane (1987)
        !
        call GWDM1987(                                       &
          & xyz_UB, xyz_VB, xyz_TempB,                       & ! (in)
          & xyz_Press, xyr_Press, xyz_Exner, xyz_Height,     & ! (in)
          & xy_SurfHeight, xy_SurfHeightStd,                 & ! (in)
          & xyz_DUDtGWD, xyz_DVDtGWD                         & ! (out)
          & )
        xyz_DUDt = xyz_DUDt + xyz_DUDtGWD
        xyz_DVDt = xyz_DVDt + xyz_DVDtGWD

      end select


      ! 陰解法で解いた地表面熱収支と整合的な $ t+\Delta t $ における長波
      ! フラックスの計算
      !   * ここで計算された値が直接次のステップ $ t $ における長波フラックス
      !     として
      !   * 用いられるわけではない
      !   * 現在の時間ステップにおける長波放射加熱率の計算に使われる
      !
      ! Evaluate longwave flux at $ t+\Delta t $ consistent with surface 
      ! energy balance solved with implicit method
      !   * The evaluated value is not used directly as Longwave flux at next 
      !     step $t$.
      !   * The evaluated value is used to calculate long wave radiative 
      !     heating rate 
      !     in the current time step.
      !
      call PhyImplEvalRadLFluxA( &
        & xyr_RadLFlux, &                                  ! (in)
        & xyz_DTempDtVDiff, xy_DSurfTempDt, xyra_DelRadLFlux, & ! (in)
        & xyr_RadLFluxA )                                  ! (out)

      ! 放射による温度変化率
      ! Temperature tendency with radiation
      !
      call RadDTempDt(                            &
        & xyr_RadLFluxA, xyr_RadSFlux, xyr_Press, &   ! (in)
        & xyz_DTempDtRadL, xyz_DTempDtRadS )          ! (out)

!!$      select case ( IDRadMethod )
!!$      case ( IDRadMethodMarsV1 )
!!$        ! (火星大気向け) Non-LTE 放射モデル
!!$        ! Non-NLTE radiation model (for the Mars' atmosphere)
!!$        !
!!$        call rad15mNLTEMergeHR(                &
!!$          & xyz_Press, xyz_TempB, xyz_VirTemp, &
!!$          & xyz_DTempDtRadL                    &
!!$          & )
!!$      end select

      ! 非断熱加熱率の総和の計算
      ! Sum all diabatic heating rates
      !
      xyz_DTempDt = xyz_DTempDtVDiff + xyz_DTempDtRadL + xyz_DTempDtRadS


      select case ( IDRadMethod )
      case ( IDRadMethodMarsV1 )
        ! 火星計算用近赤外加熱率計算
        ! Calculation of near infrared heating rate in the case of Mars
        !
        call RadMarsNIRINOUT(     &
          & xyz_Press,            &  ! (in)
          & xyz_DTempDt           &  ! (inout)
          & )
      end select


      ! 鉛直拡散フラックスの出力 
      !   * 出力のみのサブルーチンであり, 計算には影響しない
      ! 
      ! Output Vertical diffusion fluxes
      !   * This subroutine works for output only, 
      !     so it does not influence a calculation.
      !
      call VDiffusionOutput(                                       &
        & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
        & xyz_DUDt,  xyz_DVDt,  xyz_DTempDtVDiff,  xyzf_DQMixDt,   & ! (in)
        & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
        & xyr_VirTemp, xyz_Height,                                 & ! (in)
        & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef      & ! (in)
        & )

      ! 地表面フラックスの出力 
      !   * 出力のみのサブルーチンであり, 計算には影響しない
      ! 
      ! Output surface fluxes
      !   * This subroutine works for output only, 
      !     so it does not influence a calculation.
      !
      call SurfaceFluxOutput(                                       &
        & xy_SnowFrac,                                              & ! (in)
        & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux,  & ! (in)
        & xy_SurfH2OVapFluxA, xy_SurfLatentHeatFluxA,               & ! (in)
        & xyz_DUDt, xyz_DVDt, xyz_DTempDtVDiff, xyzf_DQMixDt,       & ! (in)
        & xy_SurfTemp, xy_DSurfTempDt,                              & ! (in)
        & xyr_Press, xyz_Exner, xyr_Exner, xy_SurfHumidCoef,        & ! (in)
        & xy_SurfVelTransCoef, xy_SurfTempTransCoef,                & ! (in)
        & xy_SurfQVapTransCoef                                      & ! (in)
        & )

      ! 放射フラックスの出力 
      !   * 出力のみのサブルーチンであり, 計算には影響しない
      ! 
      ! Output radiation fluxes
      !   * This subroutine works for output only, 
      !     so it does not influence a calculation.
      !
      call RadFluxOutput(                                &
        & xyr_RadSUwFlux, xyr_RadSDwFlux,                & ! (in)
        & xyr_RadLUwFlux, xyr_RadLDwFlux,                & ! (in)
        & xyra_DelRadLUwFlux, xyra_DelRadLDwFlux,        & ! (in)
        & xy_DSurfTempDt, xyz_DTempDtVDiff               & ! (in)
        & )

    end select

    !* For coupler ******
    call dcpam_StoreAtmSurfFlxInfo()
    

    call ProfUtil_RapStart('Dyn', 1)    
    ! 力学過程
    ! Dynamical core
    !
    select case ( IDDynMode )
    case ( IDDynModeHSPLVAS83 )
      call DynamicsHSplVAS83( &
        & xyz_UB,   xyz_VB,   xyz_TempB,   xyzf_QMixB,   xy_PsB, & ! (in)
        & xyz_UN,   xyz_VN,   xyz_TempN,   xyzf_QMixN,   xy_PsN, & ! (in)
        & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt,         & ! (in)
        & xy_SurfHeight,                                         & ! (in)
        & xyz_UA,   xyz_VA,   xyz_TempA,   xyzf_QMixA,   xy_PsA, & ! (out)
        & xyz_ArgOMG = xyz_OMG                                   & ! (out) optional
        & )
    case ( IDDynModeNoHorAdv )
      call DynamicsPhysicsOnly(                          &
        & xyz_Exner, xy_SurfHeight, xyz_Height,          & ! (in)
        & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt, & ! (in)
        & xy_PsB, xyz_UB, xyz_VB, xyz_TempB, xyzf_QMixB, & ! (in)
        & xy_PsN, xyz_UN, xyz_VN, xyz_TempN, xyzf_QMixN, & ! (in)
        & xy_PsA, xyz_UA, xyz_VA, xyz_TempA, xyzf_QMixA  & ! (out)
        & )
      xyz_OMG = 0.0_DP
    case ( IDDynModeTWPICE )
      !
      ! Dynamical process for TWP-ICE experiment
      !
      call DynamicsTWPICESCMExp( &
        & TimeN,                                         &
        & xyz_Press, xyz_Exner, xyz_Height,              &
        & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt, &
        & xy_PsB, xyz_UB, xyz_VB, xyz_TempB, xyzf_QMixB, &
        & xy_PsA, xyz_UA, xyz_VA, xyz_TempA, xyzf_QMixA  &
        & )
      xyz_OMG = 0.0_DP
    end select
    call ProfUtil_RapEnd('Dyn', 1)


    ! 地面温度・土壌温度・土壌水分・積雪量の積分
    ! Time integration of surface temperature, soil temperature, soil 
    ! moisture, and surface snow amount
    !
    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )

      ! 地面温度・土壌温度の時間積分
      ! Time integration of surface temperature and soil temperature
      !
      call IntegralSurfTemp( &
        & xy_DSurfTempDt, xyz_DSoilTempDt, &     ! (in)
        & xy_SurfTemp   , xyz_SoilTemp     &     ! (inout)
        & )

      select case ( IDSfcMoistMethod )
      case ( IDSfcMoistMethodNone )
        xy_SoilMoistA = xy_SoilMoistB
        xy_SurfSnowA  = xy_SurfSnowB
      case ( IDSfcMoistMethodBucket )
        ! 土壌水分・地面積雪量の時間積分
        ! Time integration of soil moisture and snow amount
        !
        call BucketSetFlagOceanFromMatthews( &
          & xy_SurfType,                     & ! (in)
          & xy_BucketFlagOceanGrid           & ! (out)
          & )
        call BucketIntegration(              &
          & xy_BucketFlagOceanGrid,          & ! (in )
          & xy_DSoilMoistDt, xy_DSurfSnowDt, & ! (in )
          & xy_SoilMoistB, xy_SurfSnowB,     & ! (in )
          & xy_SoilMoistA, xy_SurfSnowA      & ! (out)
          & )
      end select

      xy_SurfMajCompIceA = xy_SurfMajCompIceB
!!$      call MajorCompPhaseChangeOnGround( &
!!$        & xy_DPsDt, xy_DSurfMajCompIceDt,         & ! (in)
!!$        & xy_PsA, xyzf_QMixA, xy_SurfMajCompIceA  & ! (inout)
!!$        & )
      call MajorCompPhaseChangeOnGround( &
        & xy_DPsDt, xy_DSurfMajCompIceDt,         & ! (in)
        & xy_PsA, xyzf_QMixA, xyz_UA, xyz_VA,     & ! (inout)
        & xy_SurfMajCompIceA                      & ! (inout)
        & )

    case default

!!$      xy_SurfMajCompIceA = xy_SurfMajCompIceB + xy_DSurfMajCompIceDt * ( 2.0_DP * DelTime )
!!$      xy_SurfMajCompIceA = xy_SurfMajCompIceA + xy_DSurfMajCompIceDt * ( 2.0_DP * DelTime )
      !
      ! Surface pressure is adjusted, here.
      !
      xy_PsA = xy_PsA + xy_DPsDt * ( 2.0_DP * DelTime )

    end select



    select case ( IDPhysMode )
    case ( IDPhysModeJupiterSimple )
      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level, 
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsA, xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),  & ! (in)
        & xyz_Press = xyz_Press,                             & ! (out) optional
        & xyr_Press = xyr_Press                              & ! (out) optional
        & )

      select case ( IDDCMethod )
      case ( IDDCMethodDCA )
        ! 乾燥対流調節
        ! Dry convective adjustment
        !
        call DryConvAdjust(                        &
          & xyz_TempA, xyz_UA, xyz_VA, xyzf_QMixA, &  ! (inout)
          & xyz_Press, xyr_Press                   &  ! (in)
          & )
      end select

    case ( IDPhysModeJupiterSimpleV2 )
      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level, 
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsA, xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),  & ! (in)
        & xyz_Press = xyz_Press,                             & ! (out) optional
        & xyr_Press = xyr_Press                              & ! (out) optional
        & )

      ! 湿潤対流調節
      ! Moist convective adjustment
      !
      call MoistConvAdjust( &
        & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),      & ! (inout)
        & xyz_Press, xyr_Press,                          & ! (in)
        & xyz_DQH2OLiqDtCum                              & ! (out)
        & )

      ! 大規模凝結 (非対流性凝結) (Manabe, 1965)
      ! Large scale condensation (non-convective condensation) 
      ! (Manabe, 1965)
      !
      call LScaleCond( &
        & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),        & ! (inout)
        & xyz_Press, xyr_Press,                            & ! (in)
        & xyz_DQH2OLiqDtLSC                                & ! (out)
        & )

      call CloudSimpleCalcPRCPKeyLLTemp3D(           &
        & xyr_Press, xyz_Press, xyz_DQH2OLiqDtCum,   &  ! (in )
        & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),  &  ! (in )
        & xy_RainCumulus, xy_SnowCumulus             &  ! (out)
        & )
      call CloudSimpleCalcPRCPKeyLLTemp3D(           &
        & xyr_Press, xyz_Press, xyz_DQH2OLiqDtLsc,   &  ! (in )
        & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),  &  ! (in )
        & xy_RainLsc, xy_SnowLsc                     &  ! (out)
        & )

      xy_Rain = xy_RainCumulus + xy_RainLsc
      xy_Snow = xy_SnowCumulus + xy_SnowLsc

      select case ( IDDCMethod )
      case ( IDDCMethodDCA )
        ! 乾燥対流調節
        ! Dry convective adjustment
        !
        call DryConvAdjust(                        &
          & xyz_TempA, xyz_UA, xyz_VA, xyzf_QMixA, &  ! (inout)
          & xyz_Press, xyr_Press                   &  ! (in)
          & )
      end select

      ! 成分の質量の補正
      ! Fix masses of constituents
      !
      call MassFixerColumn( &
        & xyr_Press,  & ! (in)
        & xyzf_QMixA  & ! (inout)
        & )

    case ( IDPhysModeFullPhysics )

      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level, 
      ! and calculate pressure and height
      !
      call AuxVars( &
        & xy_PsA, xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap), & ! (in )
        & xyr_Press     = xyr_Press,                        & ! (out) optional
        & xyz_Press     = xyz_Press,                        & ! (out) optional
        & xyz_Exner     = xyz_Exner, xyr_Exner = xyr_Exner, & ! (out) optional
        & xy_SurfHeight = xy_SurfHeight,                    & ! (in ) optional
        & xyr_Height    = xyr_Height,                       & ! (out) optional
        & xyz_Height    = xyz_Height                        & ! (out) optional
        & )

      call ProfUtil_RapStart('Cloud_Cumulus', 1)
      ! 積雲パラメタリゼーション
      ! Cumulus parameterization
      !
      select case ( IDMCMethod )
      case ( IDMCMethodNone )

        xyz_DQH2OLiqDtCum = 0.0_DP
        xyz_DQH2OSolDtCum = 0.0_DP

        xyz_MoistConvDetTend        = 0.0_DP
        xyz_MoistConvSubsidMassFlux = 0.0_DP

        xy_RainCumulus = 0.0_DP
        xy_SnowCumulus = 0.0_DP

      case ( IDMCMethodMCA, IDMCMethodMCAI98 )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        ! 湿潤対流調節
        ! Moist convective adjustment
        !
        if ( IDMCMethod == IDMCMethodMCA ) then
           call MoistConvAdjust( &
                & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),      & ! (inout)
                & xyz_Press, xyr_Press,                          & ! (in)
                & xyz_DQH2OLiqDtCum                              & ! (out)
                & )
        else
           call MoistConvAdjustI98( &
                & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),      & ! (inout)
                & xyz_Press, xyr_Press,                          & ! (in)
                & xyz_DQH2OLiqDtCum                              & ! (out)
                & )
        end if
        
        ! It would be better that lines below are included in 
        ! MoistConvAdjust subroutine.
        xyzf_QMixA(:,:,:,IndexH2OLiq) = &
          &   xyzf_QMixA(:,:,:,IndexH2OLiq) &
          & + xyz_DQH2OLiqDtCum * 2.0_DP * DelTime
        xyz_DQH2OLiqDtCum = 0.0_DP
        xyz_DQH2OSolDtCum = 0.0_DP

        xyz_MoistConvDetTend        = 0.0_DP
        xyz_MoistConvSubsidMassFlux = 0.0_DP

        xy_RainCumulus = 0.0_DP
        xy_SnowCumulus = 0.0_DP

      case ( IDMCMethodRAS )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        ! Relaxed Arakawa-Schubert scheme
        ! Relaxed Arakawa-Schubert scheme
        !
        call RAS(                                        &
          & xy_SurfTemp,                                 &  ! (in)
          & xyz_Press, xyr_Press, xyz_Exner, xyr_Exner,  &  ! (in)
          & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),    &  ! (inout)
          & xyz_DQH2OLiqDtCum,                           &  ! (out) optional
          & xyz_MoistConvDetTend,                        &  ! (out) optional
          & xyz_MoistConvSubsidMassFlux                  &  ! (out) optional
          & )

        ! It would be better that lines below are included in 
        ! RAS subroutine.
        xyzf_QMixA(:,:,:,IndexH2OLiq) = &
          &   xyzf_QMixA(:,:,:,IndexH2OLiq) &
          & + xyz_DQH2OLiqDtCum * 2.0_DP * DelTime
        xyz_DQH2OLiqDtCum = 0.0_DP
        xyz_DQH2OSolDtCum = 0.0_DP

        xy_RainCumulus = 0.0_DP
        xy_SnowCumulus = 0.0_DP

      case ( IDMCMethodRASWithIce )
        ! Relaxed Arakawa-Schubert scheme
        ! Relaxed Arakawa-Schubert scheme
        !
        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        if ( IndexH2OSol <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OSol))//' is not found.' )
        end if

        ! It should be noted that H2OLiq and H2OSol have updated. 
        ! But, DQH2OLiqDtCum and DQH2OSolDtCum are not zero, so those
        ! will be set zero below.
        call RASWithIce1DWrapper3DWrapper( &
          & xy_SurfTemp,                                 & ! (in)
          & xyz_Press, xyr_Press, xyz_Exner, xyr_Exner,  & ! (in)
          & xyz_TempA,                                   & ! (in)
          & xyzf_QMixA(:,:,:,IndexH2OVap),               & ! (in)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),               & ! (in)
          & xyzf_QMixA(:,:,:,IndexH2OSol),               & ! (in)
          & xyz_UA, xyz_VA,                                       & ! (in)
          & xyz_DTempDtCum,                                       & ! (out)
          & xyz_DQVapDtCum, xyz_DQH2OLiqDtCum, xyz_DQH2OSolDtCum, & ! (out)
          & xyz_DUDtCum, xyz_DVDtCum,                             & ! (out)
          & xy_RainCumulus, xy_SnowCumulus,                       & ! (out)
          & xyz_MoistConvDetTend,                        & ! (out) optional
          & xyz_MoistConvSubsidMassFlux                  & ! (out) optional
          & )

        xyz_TempA =     &
          &   xyz_TempA &
          & + xyz_DTempDtCum * ( 2.0_DP * DelTime )
        xyzf_QMixA(:,:,:,IndexH2OVap) =     &
          &   xyzf_QMixA(:,:,:,IndexH2OVap) &
          & + xyz_DQVapDtCum * ( 2.0_DP * DelTime )
        xyzf_QMixA(:,:,:,IndexH2OLiq) =     &
          &   xyzf_QMixA(:,:,:,IndexH2OLiq) &
          & + xyz_DQH2OLiqDtCum * ( 2.0_DP * DelTime )
        xyzf_QMixA(:,:,:,IndexH2OSol) =     &
          &   xyzf_QMixA(:,:,:,IndexH2OSol) &
          & + xyz_DQH2OSolDtCum * ( 2.0_DP * DelTime )
        xyz_UA =     &
          &   xyz_UA &
          & + xyz_DUDtCum * ( 2.0_DP * DelTime )
        xyz_VA =     &
          &   xyz_VA &
          & + xyz_DVDtCum * ( 2.0_DP * DelTime )

        ! 成分の質量の補正
        ! Fix masses of constituents
        !
        call MassFixerColumn( &
          & xyr_Press,  & ! (in)
          & xyzf_QMixA  & ! (inout)
          & )

        xyz_DTempDtCum    = 0.0_DP
        xyz_DQVapDtCum    = 0.0_DP
        xyz_DQH2OLiqDtCum = 0.0_DP
        xyz_DQH2OSolDtCum = 0.0_DP
        xyz_DUDtCum       = 0.0_DP
        xyz_DVDtCum       = 0.0_DP


      end select
      call ProfUtil_RapEnd('Cloud_Cumulus', 1)


      call ProfUtil_RapStart('Cloud_Lscond', 1)      
      ! 大規模凝結 (非対流性凝結)
      ! Large scale condensation
      !
      select case ( IDLSCMethod )
      case ( IDLSCMethodNone )

        xyz_DQH2OLiqDtLSC = 0.0_DP
        xyz_DQH2OSolDtLSC = 0.0_DP

      case ( IDLSCMethodM65 )
        ! 大規模凝結 (非対流性凝結) (Manabe, 1965)
        ! Large scale condensation (non-convective condensation) 
        ! (Manabe, 1965)
        !
        call LScaleCond( &
!!$        call LScaleCond1D3DWrapper( &
          & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),        & ! (inout)
          & xyz_Press, xyr_Press,                            & ! (in)
          & xyz_DQH2OLiqDtLSC                                & ! (out)
          & )

        ! It would be better that lines below are included in LScaleCond
        ! subroutine.
        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        xyzf_QMixA(:,:,:,IndexH2OLiq) = &
          &   xyzf_QMixA(:,:,:,IndexH2OLiq) &
          & + xyz_DQH2OLiqDtLSC * 2.0_DP * DelTime
        xyz_DQH2OLiqDtLSC = 0.0_DP
        xyz_DQH2OSolDtLSC = 0.0_DP

!!$      case ( IDLSCMethodLL91 )
!!$          ! 大規模凝結 (非対流性凝結) (Le Treut and Li, 1991)
!!$          ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
!!$          !
!!$          call LScaleCondLL91(                                 &
!!$            & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),        &  ! (inout)
!!$            & xyz_DTempDtCond, xyz_DQVapDtCond,                &  ! (inout)
!!$            & xyz_Press, xyr_Press,                            &  ! (in)
!!$            & xyz_DQH2OLiqDtLSC                                &  ! (out)
!!$            & )
      case ( IDLSCMethodSatAdjM65 )
        !
        != Saturation adjustment
        !
        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        call SaturationAdjust(                                 &
          & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),          & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),                     & ! (in)
          & xyz_Press, xyr_Press,                              & ! (in)
          & xyz_DQH2OLiqDtLSC                                  & ! (out)
          & )

      case ( IDLSCMethodM65WithIce )
        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        if ( IndexH2OSol <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OSol))//' is not found.' )
        end if

        ! 大規模凝結 (非対流性凝結) (Manabe, 1965)
        ! Large scale condensation (non-convective condensation) 
        ! (Manabe, 1965)
        !
        ! It should be noted that H2OLiq and H2OSol have updated in above 
        ! subroutine.

        call LScaleCond1D3DWrapper(                            &
          & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),          & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),                     & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OSol),                     & ! (inout)
          & xyz_Press, xyr_Press,                              & ! (in)
          & xyz_DQH2OLiqDtLSC, xyz_DQH2OSolDtLSC               & ! (out)
          & )

        ! Temporal treatment
        xyz_DQH2OLiqDtLSC = 0.0_DP
        xyz_DQH2OSolDtLSC = 0.0_DP

      case ( IDLSCMethodLL91WithIce )
        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        if ( IndexH2OSol <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OSol))//' is not found.' )
        end if

        ! 大規模凝結 (非対流性凝結) (Le Treut and Li, 1991)
        ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
        !
        ! It should be noted that H2OLiq and H2OSol have updated in above 
        ! subroutine.

        call LScaleCondLL911D3DWrapper(                        &
          & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),          & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),                     & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OSol),                     & ! (inout)
          & xyz_Press, xyr_Press,                              & ! (in)
          & xyz_DQH2OLiqDtLSC, xyz_DQH2OSolDtLSC               & ! (out)
!!$          & xyz_DQH2OLiqDtLSC, xyz_DQH2OSolDtLSC,              & ! (out)
!!$          & xyz_CloudCover = xyzf_QMixA(:,:,:,CompositionInqIndex( 'CloudCover' )) & ! (inout)
          & )

        ! Temporal treatment
        xyz_DQH2OLiqDtLSC = 0.0_DP
        xyz_DQH2OSolDtLSC = 0.0_DP

      end select
      call ProfUtil_RapEnd('Cloud_Lscond', 1)


      call ProfUtil_RapStart('Cloud_clmodel', 1)            
      ! 
      ! Cloud model
      !
      select case ( IDCloudMethod )
      case ( IDCloudMethodNone )

        ! 雲なしモデル
        ! No cloud model
        !
        if ( IndexH2OSol <= 0 ) then
          if ( IndexH2OLiq <= 0 ) then
            xy_RainLsc = 0.0_DP
            xy_SnowLsc = 0.0_DP
          else
            call CloudNone(                      &
              & xyr_Press, xyz_Press,                             & ! (in   )
              & xyzf_QMixA(:,:,:,IndexH2OLiq)                     &
              &   / ( 2.0_DP * DelTime ),                         & ! (in )
              & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),         & ! (in )
              & xy_RainLsc                                        & ! (out)
              & )
            xyzf_QMixA(:,:,:,IndexH2OLiq) = 0.0_DP
            xy_SnowLsc = 0.0_DP
          end if
        else
          call CloudNoneWithIce(                                &
            & xyr_Press, xyz_Press,                             & ! (in   )
            & xyzf_QMixA(:,:,:,IndexH2OLiq)                     &
            &   / ( 2.0_DP * DelTime ),                         & ! (in )
            & xyzf_QMixA(:,:,:,IndexH2OSol)                     &
            &   / ( 2.0_DP * DelTime ),                         & ! (in )
            & xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),         & ! (in )
            & xy_RainLsc, xy_SnowLsc                            & ! (out)
            & )
          xyzf_QMixA(:,:,:,IndexH2OLiq) = 0.0_DP
          xyzf_QMixA(:,:,:,IndexH2OSol) = 0.0_DP
        end if

      case ( IDCloudMethodSimple )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        ! Update cloud water
        !
        call CloudSimple(                                    &
          & xyr_Press, xyz_Press,                            & ! (in)
          & xyz_TempA,                                       & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OVap),                   & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),                   & ! (inout)
          & xy_RainLsc, xy_SnowLsc                           & ! (out)
          & )

      case ( IDCloudMethodSimpleWithIce )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        if ( IndexH2OSol <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OSol))//' is not found.' )
        end if

        ! Update cloud water
        !
        call CloudSimpleWithIce(                             &
          & xyr_Press, xyz_Press,                            & ! (in)
          & xyz_TempA,                                       & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OVap),                   & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),                   & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OSol),                   & ! (inout)
          & xy_RainLsc, xy_SnowLsc                           & ! (out)
          & )

!!$        call PhyImplSDHV2SetMethodMatthews(  &
!!$          & xy_SurfType, xy_SeaIceConc,      & ! (in)
!!$          & xy_PhyImplSDHIndexCalcMethod     & ! (out)
!!$          & )
!!$        call PhyImplSDHV2CorSOTempBySnowMelt(  &
!!$          & xy_PhyImplSDHIndexCalcMethod,    & ! (in   )
!!$          & xy_Snow,                         & ! (in   )
!!$          & xy_SurfTemp                      & ! (inout)
!!$          & )

      case ( IDCloudMethodT1993WithIce )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        if ( IndexH2OSol <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OSol))//' is not found.' )
        end if
        if ( CompositionInqIndex( 'CloudCover' ) <= 0 ) then
          call MessageNotify( 'E', prog_name, 'CloudCover is not found.' )
        end if

        ! 温度の半整数σレベルの補間, 気圧と高度の算出
        ! Interpolate temperature on half sigma level, 
        ! and calculate pressure and height
        !
        call AuxVars( &
          & xy_PsA, xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap), & ! (in )
          & xyr_Press     = xyr_Press,                        & ! (out) optional
          & xyz_Press     = xyz_Press,                        & ! (out) optional
          & xyz_VirTemp   = xyz_VirTemp                       & ! (out) optional
          & )

        ! Tiedtke (1993) に基づく雲モデル
        ! Cloud model based on Tiedtke (1993)
        !
        call CloudT1993baseWithIce(                              &
          & xyz_Press, xyr_Press, xyz_VirTemp,                   & ! (in)
          & xyz_DQH2OLiqDtCum, xyz_DQH2OSolDtCum,                & ! (in)
          & xyz_MoistConvDetTend,                                & ! (in)
          & xyz_OMG, xyz_MoistConvSubsidMassFlux, xyz_DTempDt,   & ! (in)
          & xyz_TempA,                                           & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OVap),                       & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OLiq),                       & ! (inout)
          & xyzf_QMixA(:,:,:,IndexH2OSol),                       & ! (inout)
          & xyzf_QMixA(:,:,:,CompositionInqIndex('CloudCover')), & ! (inout)
          & xy_RainLsc, xy_SnowLsc                               & ! (out)
          & )

      case ( IDCloudMethodMarsH2OCloud )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        ! 温度の半整数σレベルの補間, 気圧と高度の算出
        ! Interpolate temperature on half sigma level, 
        ! and calculate pressure and height
        !
        call AuxVars( &
          & xy_PsA, xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),& ! (in )
          & xyz_VirTemp = xyz_VirTemp,                       & ! (out) optional
          & xyr_Press     = xyr_Press,                       & ! (out) optional
          & xy_SurfHeight = xy_SurfHeight,                   & ! (in ) optional
          & xyr_Height = xyr_Height                          & ! (out) optional
          & )

        ! 火星 H2O 雲モデル
        ! Mars H2O cloud model
        !
        call CloudMarsH2O(                                   &
          & xyr_Press, xyz_VirTemp, xyr_Height,              & ! (in)
          & xyz_DQH2OSolDtCum, xyz_DQH2OSolDtLSC,            & ! (in)
          & xyzf_QMixA(:,:,:,IndexH2OSol),                   & ! (inout)
          & xy_RainLsc, xy_SnowLsc                           & ! (out)
          & )

      end select


      ! Sum of cumulus and non-convective (large scale) condensation
      xy_Rain = xy_RainCumulus + xy_RainLsc
      xy_Snow = xy_SnowCumulus + xy_SnowLsc
      call ProfUtil_RapEnd('Cloud_clmodel', 1)      


      select case ( IDSfcMoistMethod )
      case ( IDSfcMoistMethodBucket )
        ! バケツモデル, 降水に伴う地表水量変化の計算
        ! bucket model, calculation of change of surface water due to 
        ! precipitation
        !
        call BucketSetFlagOceanFromMatthews( &
          & xy_SurfType,                     & ! (in)
          & xy_BucketFlagOceanGrid           & ! (out)
          & )
        call BucketPRCPAdjust(                  &
          & xy_BucketFlagOceanGrid, xy_Rain, xy_Snow, &  ! (in )
          & xy_SoilMoistA, xy_SurfSnowA         &  ! (inout)
          & )
      end select


      select case ( IDDCMethod )
      case ( IDDCMethodDCA )
        ! 乾燥対流調節
        ! Dry convective adjustment
        !
        call DryConvAdjust(                        &
          & xyz_TempA, xyz_UA, xyz_VA, xyzf_QMixA, &  ! (inout)
          & xyz_Press, xyr_Press                   &  ! (in)
          & )
      end select



      ! 重力沈降過程
      ! Gravitational sedimentation process
      !
      if ( CompositionInqIndex('QDust') > 0 ) then
        call GravSed(                                       &
          & 'MarsDust',                                     & ! (in )
          & xyr_Press, xyz_VirTemp, xyr_Height,             & ! (in )
          & xyzf_QMixA(:,:,:,CompositionInqIndex('QDust')), & ! (out)
          & xy_SurfDustGravSedFlux                          & ! (out) optional
          & )
      end if

!!$      if ( IDPhysMode == IDPhysModeFullPhysics ) then

        ! 主成分相変化
        ! Phase change of atmospheric major component
        !
!!$        call MajorCompPhaseChangeInAtm(                        &
!!$          & xyr_Press, xyz_Press,                              &  ! (in)
!!$          & xy_PsA, xyz_TempA, xyzf_QMixA, xy_SurfMajCompIceA  &  ! (inout)
!!$          & )
        call MajorCompPhaseChangeInAtm(                        &
          & xyr_Press, xyz_Press,                              &  ! (in)
          & xy_PsA, xyz_TempA, xyzf_QMixA, xyz_UA, xyz_VA,     &  ! (inout)
          & xy_SurfMajCompIceA                                 &  ! (inout)
          & )

!!$      end if

      ! 温度の半整数σレベルの補間, 気圧と高度の算出
      ! Interpolate temperature on half sigma level, 
      ! and calculate pressure and height
      !
       call AuxVars( &
        & xy_PsA, xyz_TempA, xyzf_QMixA(:,:,:,IndexH2OVap),& ! (in )
        & xyr_Press = xyr_Press                            & ! (out) optional
        & )
      ! 成分の質量の補正
      ! Fix masses of constituents
      !
      call MassFixerColumn( &
        & xyr_Press,  & ! (in)
        & xyzf_QMixA  & ! (inout)
        & )

    end select

#ifdef INTH98_MODIFY
    call VerticalFilterAdjust( &
         & xyz_UA, xyz_VA, xyz_TempA, &
         & xy_PsA )

    call SurfPresChangeWithWtVap( &
         & xy_PsA, xyzf_QMixA(:,:,:,IndexH2Ovap),  &
         & xy_SurfH2OVapFluxA, xy_Rain + xy_Snow,  &
         & xy_PsB, xyzf_QMixB(:,:,:,IndexH2Ovap) )
#endif
    
    ! 時間フィルター (Asselin, 1972)
    ! Time filter (Asselin, 1972)
    !
!!$    if ( .not. flag_initial .or. .not. firstloop ) then
!!$      call TimeFilter( &
!!$        & xyz_UB, xyz_VB, xyz_TempB, xyzf_QMixB, xy_PsB, &   ! (in)
!!$        & xyz_UN, xyz_VN, xyz_TempN, xyzf_QMixN, xy_PsN, &   ! (inout)
!!$        & xyz_UA, xyz_VA, xyz_TempA, xyzf_QMixA, xy_PsA  )   ! (in)
!!$
!!$      select case ( IDPhysMode )
!!$      case ( IDPhysModeFullPhysics )
!!$        call TimeFilterSurfVars(                             &
!!$          & xy_SurfMajCompIceB, xy_SoilMoistB, xy_SurfSnowB, &   ! (in)
!!$          & xy_SurfMajCompIceN, xy_SoilMoistN, xy_SurfSnowN, &   ! (inout)
!!$          & xy_SurfMajCompIceA, xy_SoilMoistA, xy_SurfSnowA  &   ! (in)
!!$          & )
!!$      end select
!!$    end if

    call ProfUtil_RapStart('TimeAdvance', 1)    
    ! 時間フィルター (Williams, 2009)
    ! Time filter (Williams, 2009)
    !
    if ( .not. flag_initial .or. .not. firstloop ) then

       
      call TimeFilterWilliams2009(                 &
        & xyz_UB, xyz_VB, xyz_TempB, xyzf_QMixB, xy_PsB, &   ! (in)
        & xyz_UN, xyz_VN, xyz_TempN, xyzf_QMixN, xy_PsN, &   ! (inout)
        & xyz_UA, xyz_VA, xyz_TempA, xyzf_QMixA, xy_PsA  &   ! (inout)
        & )

      select case ( IDPhysMode )
      case ( IDPhysModeFullPhysics )
        call TimeFilterWilliams2009SurfVars(           &
          & xy_SurfMajCompIceB, xy_SoilMoistB, xy_SurfSnowB, & ! (in)
          & xy_SurfMajCompIceN, xy_SoilMoistN, xy_SurfSnowN, & ! (inout)
          & xy_SurfMajCompIceA, xy_SoilMoistA, xy_SurfSnowA, & ! (inout)
          & xy_PsA                                           & ! (inout)
          & )
      end select
    end if


    ! 予報変数の値の確認
    ! Check values of prognostic variables
    !
    call CheckProgVars( &
      & xy_PsA, xyz_UA, xyz_VA, xyz_TempA, xyzf_QMixA   & ! (in)
      & )

    ! ヒストリデータ出力
    ! History data output
    !
    call HistoryAutoPut( TimeA, 'U',    xyz_UA )
    call HistoryAutoPut( TimeA, 'V',    xyz_VA )
    call HistoryAutoPut( TimeA, 'Temp', xyz_TempA )
    do n = 1, ncmax
      call HistoryAutoPut( TimeA, a_QMixName(n), xyzf_QMixA(:,:,:,n) )
    end do
    ! Backward compatibility
    n = IndexH2OVap
    if ( a_QMixName(n) /= 'QVap' ) then
      call HistoryAutoPut( TimeA, 'QVap', xyzf_QMixA(:,:,:,n) )
    end if
    call HistoryAutoPut( TimeA, 'Ps',   xy_PsA )

    ! Output frequently used variables
    ! Output frequently used variables
    !
    call OutputFreqUsedVars(           &
      & xy_PsA, xyz_TempA, xyzf_QMixA, & ! (in)
      & xy_SurfHeight                  & ! (in)
      & )

    select case ( IDPhysMode )
    case ( IDPhysModeJupiterSimpleV2 )
      call HistoryAutoPut( TimeN, 'Rain', xy_Rain )
      call HistoryAutoPut( TimeN, 'Snow', xy_Snow )
      call HistoryAutoPut( TimeN, 'PRCP', xy_Rain + xy_Snow )

    case ( IDPhysModeFullPhysics )
      call HistoryAutoPut( TimeN, 'SurfTemp', xy_SurfTemp )
      if ( size( xyz_SoilTemp ) /= 0 ) then
        call HistoryAutoPut( TimeN, 'SoilTemp', xyz_SoilTemp )
      end if

      call HistoryAutoPut( TimeA, 'SurfMajCompIce', xy_SurfMajCompIceA )
      call HistoryAutoPut( TimeA, 'SoilMoist'     , xy_SoilMoistA      )
      call HistoryAutoPut( TimeA, 'SurfSnow'      , xy_SurfSnowA       )

!!$      call HistoryAutoPut( TimeN, 'Rain', xy_Rain * LatentHeat )
!!$      call HistoryAutoPut( TimeN, 'Rain', ( xy_Rain + xy_Snow ) * LatentHeat )
      call HistoryAutoPut( TimeN, 'Rain', xy_Rain )
      call HistoryAutoPut( TimeN, 'Snow', xy_Snow )
      call HistoryAutoPut( TimeN, 'PRCP', xy_Rain + xy_Snow )
      call HistoryAutoPut( TimeN, 'RainCum', xy_RainCumulus )
      call HistoryAutoPut( TimeN, 'SnowCum', xy_SnowCumulus )
      call HistoryAutoPut( TimeN, 'PRCPCum', xy_RainCumulus + xy_SnowCumulus )
      call HistoryAutoPut( TimeN, 'RainLsc', xy_RainLsc )
      call HistoryAutoPut( TimeN, 'SnowLsc', xy_SnowLsc )
      call HistoryAutoPut( TimeN, 'PRCPLsc', xy_RainLsc + xy_SnowLsc )

      call HistoryAutoPut( TimeN, 'SeaIceConc'      , xy_SeaIceConc       )
      call HistoryAutoPut( TimeN, 'SurfAlbedo'      , xy_SurfAlbedo       )
      call HistoryAutoPut( TimeN, 'SurfRoughLenMom' , xy_SurfRoughLenMom  )
      call HistoryAutoPut( TimeN, 'SurfRoughLenHeat', xy_SurfRoughLenHeat )

      if ( CompositionInqIndex('QDust') > 0 ) then
        call HistoryAutoPut( TimeN, 'SurfDustGravSedFlux', xy_SurfDustGravSedFlux )
      end if
    end select

    ! 次の時間ステップに向けて予報変数の入れ替え
    ! Exchange prediction variables for the next time step
    !
    !$omp parallel
    !$omp workshare
    xyz_UB     = xyz_UN     ; xyz_UN     = xyz_UA     ; xyz_UA     = 0.
    xyz_VB     = xyz_VN     ; xyz_VN     = xyz_VA     ; xyz_VA     = 0.
    xyz_TempB  = xyz_TempN  ; xyz_TempN  = xyz_TempA  ; xyz_TempA  = 0.
    !$omp end workshare
    !$omp workshare
    xyzf_QMixB = xyzf_QMixN ; xyzf_QMixN = xyzf_QMixA ; xyzf_QMixA = 0.
    !$omp end workshare
    !$omp workshare
    xy_PsB     = xy_PsN     ; xy_PsN     = xy_PsA     ; xy_PsA     = 0.
    !$omp end workshare
    !$omp end parallel

    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )
      xy_SurfMajCompIceB = xy_SurfMajCompIceN
      xy_SurfMajCompIceN = xy_SurfMajCompIceA
      xy_SurfMajCompIceA = 0.0_DP

      xy_SoilMoistB = xy_SoilMoistN
      xy_SoilMoistN = xy_SoilMoistA
      xy_SoilMoistA = 0.0_DP

      xy_SurfSnowB  = xy_SurfSnowN
      xy_SurfSnowN  = xy_SurfSnowA
      xy_SurfSnowA  = 0.0_DP
    end select
    call ProfUtil_RapEnd('TimeAdvance', 1)    

endif ! end if for skip_flag
 
    ! 時刻の進行
    ! Progress time
    !
    call TimesetProgress

    ! NAMELIST から読み込んだ変数名に無効なものが存在したかどうかをチェック
    ! HistoryAutoAddVariable で登録した変数名を印字
    !
    ! Check that invalid variable names are loaded from NAMELIST or not
    ! Print registered variable names by "HistoryAutoAddVariable"
    !
    !!! if ( firstloop ) call HistoryAutoAllVarFix

    ! リスタートデータ出力
    ! Restart data output
    !
    call RestartFileOutput(                            &
      & xyz_UB, xyz_VB, xyz_TempB, xyzf_QMixB, xy_PsB, &  ! (in)
      & xyz_UN, xyz_VN, xyz_TempN, xyzf_QMixN, xy_PsN  &  ! (in)
      & )

    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )
      ! 地表面温度リスタートデータ出力
      ! Restart data of surface temperature output
      !
      call RestartSurfTempOutput(                             &
        & xy_SurfTemp, xyz_SoilTemp,                          & ! (in)
        & xy_SurfMajCompIceB, xy_SoilMoistB, xy_SurfSnowB,    & ! (in) optional
        & xy_SurfMajCompIceN, xy_SoilMoistN, xy_SurfSnowN     & ! (in) optional
        & )
    end select

    firstloop = .false.

  ! 時間積分終了
  ! Time integration is finished
  !
  !end do loop_time

    call ProfUtil_RapEnd('TimeLoop', 0)
    
  end subroutine dcpam_advance_timestep

  !-------------------------------------------------------------------

  subroutine MainInit
    !
    ! 主プログラムの初期化手続き. 
    !
    ! Initialization procedure for the main program. 
    !

    ! MPI
    !
    use mpi_wrapper, only : MPIWrapperInit

    use dc_message, only: MessageNotify

    ! コマンドライン引数処理
    ! Command line option parser
    !
    use option_parser, only: OptParseInit

    ! NAMELIST ファイル入力に関するユーティリティ
    ! Utilities for NAMELIST file input
    !
    use namelist_util, only: NmlutilInit, NmlutilMsg

    ! 時刻管理
    ! Time control
    !
    use timeset, only: TimesetInit, TimesetDelTimeHalf, &
      & TimeN                 ! ステップ $ t $ の時刻. Time of step $ t $. 

    ! 出力ファイルの基本情報管理
    ! Management basic information for output files
    ! 
    use fileset, only: FilesetInit

    ! 格子点設定
    ! Grid points settings
    !
    use gridset, only: GridsetInit, &
      &                  imax, & ! 経度格子点数. 
                                 ! Number of grid points in longitude
      &                  jmax, & ! 緯度格子点数. 
                                 ! Number of grid points in latitude
      &                  kmax, & ! 鉛直層数. 
                                 ! Number of vertical level
      &                kslmax    ! 地下の鉛直層数. 
                                 ! Number of subsurface vertical level

    ! 組成に関わる配列の設定
    ! Settings of array for atmospheric composition
    !
    use composition, only: CompositionInit

    ! 物理定数設定
    ! Physical constants settings
    !
    use constants, only: ConstantsInit

    ! 雪と海氷の定数の設定
    ! Setting constants of snow and sea ice
    !
    use constants_snowseaice, only: ConstantsSnowSeaIceInit

    ! 座標データ設定
    ! Axes data settings
    !
    use axesset, only: AxessetInit

    ! リスタートデータ入出力
    ! Restart data input/output
    !
    use restart_file_io, only: RestartFileInit, RestartFileOpen, RestartFileGet

    ! 地表面温度リスタートデータ入出力
    ! Restart data of surface temperature input/output
    !
    use restart_surftemp_io, only: RestartSurfTempInit, RestartSurfTempOpen, RestartSurfTempGet

    ! ヒストリデータ出力
    ! History data output
    !
    use history_file_io, only: HistoryFileOpen
    use gtool_historyauto, only: HistoryAutoAddVariable, HistoryAutoPut

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: &
      & STDOUT               ! 標準出力の装置番号. Unit number of standard output

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 時系列データの読み込み
    ! Reading time series
    !
    use read_time_series, only : ReadTimeSeriesInit

    ! 地表面データ提供
    ! Prepare surface data
    !
    use surface_data, only : SurfDataInit

    ! ファイルから 1 次元プロファイルを読んで設定する. 
    ! read 1-D profile from a file and set it 
    !
    use set_1d_profile, only : Set1DProfileInit


    ! 補助的な変数を計算するサブルーチン・関数群
    ! Subroutines and functions for calculating auxiliary variables
    !
    use auxiliary, only : AuxVarsInit, AuxVars

    ! Held and Suarez (1994) による強制と散逸
    ! Forcing and dissipation suggested by Held and Suarez (1994)
    !
    use held_suarez_1994, only : HS94Init

    ! Yamamoto and Takahashi (2003) に従った簡単金星計算のための強制
    ! forcing for simple Venus calculation following Yamamoto and Takahashi (2003)
    !
    use yt2003_forcing, only : YT2003ForcingInit

    ! 惑星表面データの設定
    ! Setting planetary surface properties
    !
    use surface_properties, only : SurfacePropertiesInit

    ! 雪, 氷の割合
    ! snow/ice fraction
    !
    use snowice_frac, only : SnowIceFracInit

    ! 地下における熱の鉛直拡散
    ! Vertical diffusion of heat under the ground
    !
    use subsurface_diffusion_heat, only : SubsurfaceDiffusionInit

    ! 陰解法による時間積分
    ! Time integration with implicit scheme
    !
    use phy_implicit, only : PhyImplInit

    ! 地下熱伝導モデルを用いた場合の陰解法による時間積分
    !
    ! Time integration by using implicit scheme in case using subsurface thermal diffusion model
    use phy_implicit_sdh, only : PhyImplSDHInit

!!$    ! 地下熱伝導モデルを用いた場合の陰解法による時間積分
!!$    !
!!$    ! Time integration by using implicit scheme in case using subsurface thermal diffusion model
!!$    use phy_implicit_sdh_V2, only : PhyImplSDHV2Init
    ! 地下熱伝導モデルを用いた場合の陰解法による時間積分
    !
    ! Time integration by using implicit scheme in case using subsurface thermal diffusion model
    use phy_implicit_sdh_V3, only : PhyImplSDHV3Init

    ! 陰解法による時間積分 (大気のみ / 惑星表面温度・土壌温度計算なし)
    ! Time integration by using implicit scheme in case without calculation of surface and soil temperature
    !
    use phy_implicit_atmonly, only : PhyImplAtmOnlyInit

    ! 陰解法による時間積分のためのルーチン
    ! Routines for time integration with implicit scheme
    !
    use phy_implicit_utils, only : PhyImplUtilsInit

    ! バケツモデル
    ! Bucket model
    !
    use Bucket_Model, only : BucketModelInit

    ! 放射フラックス (GFD 電脳倶楽部開発の放射モデル)
    ! Radiation flux (radiation model developed by GFD Dennou Club)
    !
    use rad_DennouAGCM, only : RadDennouAGCMInit

    ! 地球大気向け放射モデル Ver. 2
    ! radiation model for the Earth's atmosphere Ver. 2
    !
    use rad_Earth_V2, only : RadEarthV2Init

    ! wrapper of RRTMG
    ! wrapper of RRTMG
    !
    use rad_rrtmg_wrapper, only : RadRRTMGWrapperInit

    ! 火星大気向け放射モデル Ver. 1
    ! radiation model for the Mars' atmosphere Ver. 1
    !
    use rad_Mars_V1, only : RadMarsV1Init

!!$    ! (火星大気向け) Non-LTE 放射モデル
!!$    ! Non-NLTE radiation model (for the Mars' atmosphere)
!!$    !
!!$    use rad_15m_NLTE, only: Rad15mNLTEInit

    ! 火星計算用近赤外加熱率計算
    ! Calculation of near infrared heating rate in the case of Mars
    !
    use rad_Mars_NIR, only : RadMarsNIRInit

    ! Schneider and Liu (2009) の放射モデル
    ! Radiation model by Schneider and Liu (2009)
    !
    use rad_SL09, only : RadSL09Init

    ! 簡単放射モデル
    ! Simple radiation model
    !
    use rad_simple, only : RadSimpleInit

    ! 何もしない放射モデル
    ! radiation model with no absorption and no scattering
    !
    use rad_none, only : RadNoneInit

    ! 放射関連ルーチン
    ! Routines for radiation calculation
    !
    use rad_utils, only : RadUtilsInit

    ! 鉛直拡散フラックス
    ! Vertical diffusion flux
    !
    use vdiffusion_my, only : VDiffusionInit

    ! JMA 乱流混合モジュール
    ! JMA turbulent mixing module
    !
    use vdiffusion_jma_my_wrapper, only : VDiffusionJMAInit

    ! Schneider and Liu (2009) による鉛直混合課程
    ! Vertical diffusion by Schneider and Liu (2009)
    !
    use sl09_diffusion, only : SL09DiffusionInit

    ! Gravity wave drag by McFarlane (1987)
    ! Gravity wave drag by McFarlane (1987)
    !
    use gwd_m1987, only : GWDM1987Init

    ! 雲なしモデル
    ! No cloud model
    !
    use cloud_none, only : CloudNoneInit

    ! 簡単雲モデル
    ! Simple cloud
    !
    use cloud_simple, only: CloudSimpleInit

    ! Tiedtke (1993) に基づく雲モデル
    ! Cloud model based on Tiedtke (1993)
    !
    use cloud_T1993base, only : CloudT1993baseInit

    ! 火星 H2O 雲モデル
    ! Mars H2O cloud model
    !
    use cloud_mars_h2o, only : CloudMarsH2OInit

    ! 地表面フラックス (バルク法)
    ! Surface flux (Bulk method)
    !
    use surface_flux_bulk, only : SurfaceFluxInit

    !
    ! set dust flux
    !
!!$    use set_dust_flux, only : SetDustFluxInit

    ! 下部境界フラックス
    ! Lower boundary flux
    !
    use lb_flux_simple, only : LBFluxSimpleInit

    ! 地面温度, 土壌温度の時間積分
    ! Time integration of surface temperature and soil temperature
    !
    use intg_surftemp, only : IntgSurfTempInit

    ! 力学過程 (スペクトル法, Arakawa and Suarez (1983))
    ! Dynamical process (Spectral method, Arakawa and Suarez (1983))
    !
    use dynamics_hspl_vas83, only : DynamicsHSplVAS83Init

    ! 物理過程のみの計算のための力学過程
    ! A dynamics for calculation with physical processes only
    !
    use dynamics_physicsonly, only : DynamicsPhysicsOnlyInit

    !
    ! Dynamical process for TWP-ICE experiment
    !
    use dynamics_twpice_scm_exp, only : DynamicsTWPICESCMExpInit

    ! 湿潤対流調節
    ! Moist convective adjustment
    !
    use moist_conv_adjust, only : MoistConvAdjustInit

    ! Relaxed Arakawa-Schubert scheme
    ! Relaxed Arakawa-Schubert scheme
    !
    use relaxed_arakawa_schubert, only : RASInit

    ! 大規模凝結 (非対流性凝結)
    ! Large scale condensation (non-convective condensation)
    !
    use lscond, only : LScaleCondInit

    ! 大規模凝結 (非対流性凝結) (Le Treut and Li, 1991)
    ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
    !
!!$    use lscond_LL91, only : LScaleCondLL91Init

    !
    != Saturation adjustment
    !
    use saturation_adjust, only : SaturationAdjustInit

    ! 乾燥対流調節
    ! Dry convective adjustment
    !
    use dry_conv_adjust, only : DryConvAdjustInit

    ! 重力沈降過程
    ! Gravitational sedimentation process
    !
    use grav_sed, only : GravSedInit

    ! 主成分相変化
    ! Phase change of atmospheric major component
    !
    use major_comp_phase_change, only : MajorCompPhaseChangeInit

    ! タイムフィルター (Asselin, 1972)
    ! Time filter (Asselin, 1972)
    !
    use timefilter_asselin1972, only : TimeFiltInit

    ! 時間フィルター (Williams, 2009)
    ! Time filter (Williams, 2009)
    !
    use timefilter_williams2009, only: TimeFilterWilliams2009Init

    ! 予報変数の値の確認
    ! Check values of prognostic variables
    !
    use check_prog_vars, only : CheckProgVarsInit

    ! 質量の補正
    ! Mass fixer
    !
    use mass_fixer, only : MassFixerInit

    ! Output frequently used variables
    ! Output frequently used variables
    !
    use output_freq_used_vars, only : OutputFreqUsedVarsInit


    use omp_wrapper, only: OMPWrapperInit
    
    ! 宣言文 ; Declaration statements
    !
    implicit none

    character(*), parameter:: version = &
      & '$Name:  $' // &
      & '$Id: dcpam_main.f90,v 1.67 2015/03/11 04:54:54 yot Exp $'
                              ! 主プログラムのバージョン
                              ! Main program version

    character(STRING)      :: namelist_filename
                              ! NAMELIST ファイルの名称. 
                              ! NAMELIST file name

    character(STRING)      :: DynMode
                                ! Dynamics used in calculation
    character(STRING)      :: PhysMode
                                ! Physics used in calculation
    character(STRING)      :: RadModel
                                ! Radiation model used in calculation
    character(STRING)      :: SfcFluxMethod
                                ! Method for surface flux evaluation used in calculation
    character(STRING)      :: VDiffMethod
                                ! Method for vertical diffusion evaluation used in calculation
    character(STRING)      :: PhysImpMode
                                ! Mode for implicit method used in calculation
    character(STRING)      :: MCMethod
                                ! Moist convection parameterization
    character(STRING)      :: LSCMethod
                                ! Large scale condensation parameterization
    character(STRING)      :: CloudMethod
                                ! Cloud model
    character(STRING)      :: SfcMoistMethod
                                ! Surface moist model
    character(STRING)      :: GWDMethod
                                ! Gravity wave drag
    character(STRING)      :: DCMethod
                                ! Dry convection

    logical                :: FlagPhysImpSoilModelSO
                                ! flag for use of slab ocean
    logical                :: FlagSnow
                                ! flag for treating snow
    logical                :: FlagMajCompPhaseChange
                                ! flag for use of major component phase change

    character(STRING):: CondMajCompName
                                ! name of condensable major component

    character(STRING):: briefexpldyn
                              ! 実行ファイルの簡潔な説明 (力学過程)
                              ! Brief account of executable file (dynamics)
    character(STRING):: briefexplphy
                              ! 実行ファイルの簡潔な説明 (物理過程)
                              ! Brief account of executable file (physics)
    character(STRING):: briefexplrad
                              ! 実行ファイルの簡潔な説明 (放射過程)
                              ! Brief account of executable file (radiation)
    character(STRING):: briefexplsfcflux
                              ! 実行ファイルの簡潔な説明 (惑星表面フラックス)
                              ! Brief account of executable file (surface flux)
    character(STRING):: briefexplvdiff
                              ! 実行ファイルの簡潔な説明 (鉛直拡散過程)
                              ! Brief account of executable file (vertical diffusion)
    character(STRING):: briefexplimp
                              ! 実行ファイルの簡潔な説明 (陰解法方程式構築方法)
                              ! Brief account of executable file (implicit method)

    logical:: FlagBucketModel

    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
                              ! Unit number for NAMELIST file open
    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
                              ! IOSTAT of NAMELIST read

    integer:: n               ! 組成方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in dimension of constituents

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /dcpam_main_nml/                                                 &
      & DynMode, PhysMode, RadModel, SfcFluxMethod, VDiffMethod, PhysImpMode, &
      & MCMethod, LSCMethod, CloudMethod, SfcMoistMethod, GWDMethod, DCMethod,&
      & FlagSnow, FlagMajCompPhaseChange, CondMajCompName
          !
          ! デフォルト値については初期化手続 "main/dcpam_main.F90#MainInit" 
          ! のソースコードを参照のこと. 
          !
          ! Refer to source codes in the initialization procedure
          ! "main/dcpam_main.F90#MainInit" for the default values. 
          !

    ! 実行文 ; Executable statement
    !

    ! Initialize MPI
    !
    if(MPI_MY_COMM /= -1) then
       call MPIWrapperInit(MPI_MY_COMM)
    else
       call MPIWrapperInit()
    end if

    ! コマンドライン引数処理
    ! Command line option parser
    !
    call OptParseInit(       &
      & namelist_filename,    & ! (out)
      & prog_name            & ! (in )
      & )

    ! NAMELIST ファイル名入力
    ! Input NAMELIST file name
    !
    call NmlutilInit( &
      & namelist_filename  & ! (in)
      & )

    ! デフォルト値の設定
    ! Default values settings
    !
    DynMode                 = 'HSPLVAS83'
!!$    DynMode                 = 'NoHorAdv'

    PhysMode                = 'FullPhysics'
!!$    PhysMode                = 'HS94'
!!$    PhysMode                = 'VenusSimple'
!!$    PhysMode                = 'JupiterSimple'
!!$    PhysMode                = 'NoPhysics'

    RadModel                = 'DennouAGCM'
!!$    RadModel                = 'Earth'
!!$    RadModel                = 'Mars'
!!$    RadModel                = 'SL09'

    SfcFluxMethod           = 'L82'
!!$    SfcFluxMethod           = 'BH91B94'

    VDiffMethod             = 'MY2'
!!$    VDiffMethod             = 'MY2.5'

    PhysImpMode             = '1LayModel'
!!$    PhysImpMode             = 'SoilModel'
!!$    PhysImpMode             = 'SoilModelSO'
!!$    PhysImpMode             = 'AtmOnly'

!!$    MCMethod                = 'None'
    MCMethod                = 'MCA'
!!$    MCMethod                = 'RAS'

!!$    LSCMethod               = 'None'
    LSCMethod               = 'M65'
!!$    LSCMethod               = 'M65WithIce'
!!$    LSCMethod               = 'LL91'
!!$    LSCMethod               = 'LL91WithIce'

    CloudMethod             = 'None'
!!$    CloudMethod             = 'Simple'
!!$    CloudMethod             = 'T1993WithIce'

    SfcMoistMethod          = 'None'
!!$    SfcMoistMethod          = 'Bucket'

    GWDMethod               = 'None'
!!$    GWDMethod               = 'M1987'

!!$    DCMethod                = 'None'
    DCMethod                = 'DCA'


    FlagSnow                = .false.

    FlagMajCompPhaseChange  = .false.

    CondMajCompName         = ''
!!$    CondMajCompName         = 'CO2'


    ! 計算モードの設定
    ! Configure calculation mode
    !
    if ( trim(namelist_filename) /= '' ) then
      call FileOpen( unit_nml, &          ! (out)
        & namelist_filename, mode = 'r' ) ! (in)

      rewind( unit_nml )
      read( unit_nml, &            ! (in)
        & nml = dcpam_main_nml, &  ! (out)
        & iostat = iostat_nml )    ! (out)
      close( unit_nml )

      call NmlutilMsg( iostat_nml, prog_name ) ! (in)
      if ( iostat_nml == 0 ) write( STDOUT, nml = dcpam_main_nml )
    end if


    ! Identification of calculation method for dynamics
    !
    call MessageNotify( 'M', prog_name, &
      & 'DynMode=<%c>.', &
      & c1 = trim(DynMode) )
    !
    select case ( DynMode )
    case ( 'HSPLVAS83' )
      IDDynMode = IDDynModeHSPLVAS83
    case ( 'NoHorAdv' )
      IDDynMode = IDDynModeNoHorAdv
    case ( 'TWPICE' )
      IDDynMode = IDDynModeTWPICE
    case default
      call MessageNotify( 'E', prog_name, &
        & 'DynMode=<%c> is not supported.', &
        & c1 = trim(DynMode) )
    end select


    ! Identification of calculation method for physics
    !
    call MessageNotify( 'M', prog_name, &
      & 'PhysMode=<%c>.', &
      & c1 = trim(PhysMode) )
    !
    select case ( PhysMode )
    case ( 'FullPhysics' )
      IDPhysMode = IDPhysModeFullPhysics
    case ( 'HS94' )
      IDPhysMode = IDPhysModeHS94
    case ( 'VenusSimple' )
      IDPhysMode = IDPhysModeVenusSimple
    case ( 'JupiterSimple' )
      IDPhysMode = IDPhysModeJupiterSimple
    case ( 'JupiterSimpleV2' )
      IDPhysMode = IDPhysModeJupiterSimpleV2
    case ( 'NoPhysics' )
      IDPhysMode = IDPhysModeNoPhysics
    case default
      call MessageNotify( 'E', prog_name, &
        & 'PhysMode=<%c> is not supported.', &
        & c1 = trim(PhysMode) )
    end select


    ! Identification of calculation method for radiation
    !
    call MessageNotify( 'M', prog_name, &
      & 'RadModel=<%c>.', &
      & c1 = trim(RadModel) )
    !
    select case ( RadModel )
    case ( 'DennouAGCM' )
      IDRadMethod = IDRadMethodDennouAGCM
    case ( 'Earth' )
      IDRadMethod = IDRadMethodEarthV2
    case ( 'Mars' )
      IDRadMethod = IDRadMethodMarsV1
    case ( 'SL09' )
      IDRadMethod = IDRadMethodSL09
!!$    case ( 'VenusSimple' )
!!$      IDRadMethod = IDRadMethodVenusSimple
    case ( 'Simple' )
      IDRadMethod = IDRadMethodSimple
    case ( 'RRTMG' )
      IDRadMethod = IDRadMethodRRTMG
    case ( 'none' )
      IDRadMethod = IDRadMethodNone
    case default
      call MessageNotify( 'E', prog_name, &
        & 'RadModel=<%c> is not supported.', &
        & c1 = trim(RadModel) )
    end select

    ! Identification of calculation method for surface flux
    !
    call MessageNotify( 'M', prog_name, &
      & 'SfcFluxMethod=<%c>.', &
      & c1 = trim(SfcFluxMethod) )
    !
    select case ( SfcFluxMethod )
    case ( 'L82' )
      IDSfcFluxMethod = IDSfcFluxMethodL82
    case ( 'BH91B94' )
      IDSfcFluxMethod = IDSfcFluxMethodBH91B94
    case default
      call MessageNotify( 'E', prog_name, &
        & 'SfcFluxMethod=<%c> is not supported.', &
        & c1 = trim(SfcFluxMethod) )
    end select

    ! Identification of calculation method for vertical diffusion
    !
    call MessageNotify( 'M', prog_name, &
      & 'VDiffMethod=<%c>.', &
      & c1 = trim(VDiffMethod) )
    !
    select case ( VDiffMethod )
    case ( 'MY2' )
      IDVDiffMethod = IDVDiffMethodMY2
    case ( 'MY2.5' )
      IDVDiffMethod = IDVDiffMethodMY25
    case ( 'JMA' )
      IDVDiffMethod = IDVDiffMethodJMA
    case default
      call MessageNotify( 'E', prog_name, &
        & 'VDiffMethod=<%c> is not supported.', &
        & c1 = trim(VDiffMethod) )
    end select

    ! Identification of calculation method for solving simultaneous linear equations 
    ! of physics
    !
    call MessageNotify( 'M', prog_name, &
      & 'PhysImpMode=<%c>.', &
      & c1 = trim(PhysImpMode) )
    !
    select case ( PhysImpMode )
    case ( '1LayModel' )
      IDPhyTendMethod = IDPhyTendMethodImp1LayModel
      FlagPhysImpSoilModelSO = .false.
    case ( 'SoilModel' )
      IDPhyTendMethod = IDPhyTendMethodImpSoilModel
      FlagPhysImpSoilModelSO = .false.
    case ( 'SoilModelSO' )
      IDPhyTendMethod = IDPhyTendMethodImpSoilModel
      FlagPhysImpSoilModelSO = .true.
    case ( 'AtmOnly' )
      IDPhyTendMethod = IDPhyTendMethodImpAtmOnly
      FlagPhysImpSoilModelSO = .false.
    case default
      call MessageNotify( 'E', prog_name, &
        & 'PhysImpMode=<%c> is not supported.', &
        & c1 = trim(PhysImpMode) )
    end select
    !
    !   Check for value of FlagFullPhysics
    !
    if ( ( IDPhysMode /= IDPhysModeFullPhysics ) .and. &
      &  ( IDPhyTendMethod == IDPhyTendMethodImpSoilModel ) ) then
      call MessageNotify( 'E', prog_name, &
        & 'PhysMode has to be "FullPhysics" true, when PhyImpMode is "SoilModel" or "SoilModelSO".' )
    end if
    if ( ( IDPhysMode /= IDPhysModeFullPhysics ) .and. &
      &  ( IDPhyTendMethod == IDPhyTendMethodImpAtmOnly ) ) then
      call MessageNotify( 'E', prog_name, &
        & 'PhysMode has to be "FullPhysics" true, when PhyImpMode is "AtmOnly".' )
    end if


    ! Identification of calculation method for moist convection
    !
    call MessageNotify( 'M', prog_name, 'MCMethod=<%c>.', c1 = trim(MCMethod) )
    !
    select case ( MCMethod )
    case ( 'None' )
      IDMCMethod = IDMCMethodNone
    case ( 'MCA' )
      IDMCMethod = IDMCMethodMCA
    case ( 'MCAI98' )
      IDMCMethod = IDMCMethodMCAI98
    case ( 'RAS' )
      IDMCMethod = IDMCMethodRAS
    case ( 'RASWithIce' )
      IDMCMethod = IDMCMethodRASWithIce
    case default
      call MessageNotify( 'E', prog_name, 'MCMethod=<%c> is not supported.', &
        & c1 = trim(MCMethod) )
    end select


    ! Identification of calculation method for large scale condensation
    !
    call MessageNotify( 'M', prog_name, 'LSCMethod=<%c>.', &
      & c1 = trim(LSCMethod) )
    !
    select case ( LSCMethod )
    case ( 'None' )
      IDLSCMethod = IDLSCMethodNone
    case ( 'M65' )
      IDLSCMethod = IDLSCMethodM65
!!$    case ( 'LL91' )
!!$      IDLSCMethod = IDLSCMethodLL91
    case ( 'SatAdjM65' )
      IDLSCMethod = IDLSCMethodSatAdjM65
    case ( 'M65WithIce' )
      IDLSCMethod = IDLSCMethodM65WithIce
    case ( 'LL91WithIce' )
      IDLSCMethod = IDLSCMethodLL91WithIce
    case default
      call MessageNotify( 'E', prog_name, 'LSCMethod=<%c> is not supported.', &
        & c1 = trim(LSCMethod) )
    end select


    ! Identification of calculation method for cloud
    !
    call MessageNotify( 'M', prog_name, 'CloudMethod=<%c>.', &
      & c1 = trim(CloudMethod) )
    !
    select case ( CloudMethod )
    case ( 'None' )
      IDCloudMethod = IDCloudMethodNone
    case ( 'Simple' )
      IDCloudMethod = IDCloudMethodSimple
    case ( 'SimpleWithIce' )
      IDCloudMethod = IDCloudMethodSimpleWithIce
    case ( 'T1993WithIce' )
      IDCloudMethod = IDCloudMethodT1993WithIce
    case ( 'MarsH2OCloud' )
      IDCloudMethod = IDCloudMethodMarsH2OCloud
    case default
      call MessageNotify( 'E', prog_name, 'LSCloudMethod=<%c> is not supported.', &
        & c1 = trim(CloudMethod) )
    end select
    ! check a non-convective condensation
    if ( IDCloudMethod == IDCloudMethodT1993WithIce ) then
      if ( IDLSCMethod /= IDLSCMethodNone ) then
        call MessageNotify( 'E', prog_name, 'If CloudMethod="T1993WithIce", LscMethod has to be "None".' )
      end if
    end if

    ! Identification of calculation method for surface moisture
    !
    call MessageNotify( 'M', prog_name, 'SfcMoistMethod=<%c>.', &
      & c1 = trim(SfcMoistMethod) )
    !
    select case ( SfcMoistMethod )
    case ( 'None' )
      IDSfcMoistMethod = IDSfcMoistMethodNone
      FlagBucketModel  = .false.
    case ( 'Bucket' )
      IDSfcMoistMethod = IDSfcMoistMethodBucket
      FlagBucketModel  = .true.
    case default
      call MessageNotify( 'E', prog_name, 'SfcMoistMethod=<%c> is not supported.', &
        & c1 = trim(SfcMoistMethod) )
    end select


    ! Identification of calculation method for gravity wave drag
    !
    call MessageNotify( 'M', prog_name, 'GWDMethod=<%c>.', c1 = trim(GWDMethod) )
    !
    select case ( GWDMethod )
    case ( 'None' )
      IDGWDMethod = IDGWDMethodNone
    case ( 'M1987' )
      IDGWDMethod = IDGWDMethodM1987
    case default
      call MessageNotify( 'E', prog_name, 'GWDMethod=<%c> is not supported.', &
        & c1 = trim(GWDMethod) )
    end select


    ! Identification of calculation method for dry convection
    !
    call MessageNotify( 'M', prog_name, 'DCMethod=<%c>.', &
      & c1 = trim(DCMethod) )
    !
    select case ( DCMethod )
    case ( 'None' )
      IDDCMethod = IDDCMethodNone
    case ( 'DCA' )
      IDDCMethod = IDDCMethodDCA
    case default
      call MessageNotify( 'E', prog_name, 'DCMethod=<%c> is not supported.', &
        & c1 = trim(DCMethod) )
    end select



    ! 計算モードの表示
    ! Display calculation mode
    !
    select case ( IDDynMode )
    case ( IDDynModeHSPLVAS83 )
      briefexpldyn = 'used'
    case ( IDDynModeNoHorAdv )
      briefexpldyn = 'not used'
    case ( IDDynModeTWPICE )
      briefexpldyn = 'TWP-ICE'
    end select

    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )
      briefexplphy = 'parameterization suite is used'
    case ( IDPhysModeHS94 )
      briefexplphy = 'forcing for Held and Suarez (1994) dynamical core test'
    case ( IDPhysModeVenusSimple )
      briefexplphy = 'simple forcing for a Venus-like planet'
    case ( IDPhysModeJupiterSimple )
      briefexplphy = 'simple forcing for a Jupiter-like planet'
    case ( IDPhysModeJupiterSimpleV2 )
      briefexplphy = 'simple forcing for a Jupiter-like planet Ver. 2'
    case ( IDPhysModeNoPhysics )
      briefexplphy = 'not used'
    end select

    if ( IDPhysMode == IDPhysModeFullPhysics ) then
      select case ( IDRadMethod )
      case ( IDRadMethodDennouAGCM )
        briefexplrad = 'dennou AGCM5 default'
      case ( IDRadMethodEarthV2 )
        briefexplrad = 'Earth'
      case ( IDRadMethodMarsV1 )
        briefexplrad = 'Mars'
      case ( IDRadMethodSL09 )
        briefexplrad = 'Schneider and Liu (2009) (Jupiter-like planet)'
!!$      case ( IDRadMethodVenusSimple )
!!$        briefexplrad = 'Venus simple radiation model (tentative)'
      case ( IDRadMethodSimple )
        briefexplrad = 'Simple radiation model (tentative)'
      case ( IDRadMethodRRTMG )
        briefexplrad = 'RRTMG'
      case ( IDRadMethodNone )
        briefexplrad = 'none'
      case default
        call MessageNotify( 'E', 'Unexpected error in setting briefexplrad', '' )
      end select

      select case ( IDSfcFluxMethod )
      case ( IDSfcFluxMethodL82 )
        briefexplsfcflux = 'surface flux by the method of Louis et al. (1982)'
      case ( IDSfcFluxMethodBH91B94 )
        briefexplsfcflux = 'surface flux by the method of Beljaars and Holtslag (1991), Beljaars (1994)'
      case default
        call MessageNotify( 'E', 'Unexpected error in setting briefexplsfcflux', '' )
      end select

      select case ( IDVDiffMethod )
      case ( IDVDiffMethodMY2 )
        briefexplvdiff = 'vertical diffusion with the method of Mellor and Yamada level 2'
      case ( IDVDiffMethodMY25 )
        briefexplvdiff = 'vertical diffusion with the method of Mellor and Yamada level 2.5'
      case ( IDVDiffMethodJMA )
        briefexplvdiff = 'vertical diffusion with the method of JMA physical library of Mellor and Yamada'! level 2.5'
      case default
        call MessageNotify( 'E', 'Unexpected error in setting briefexplvdiff', '' )
      end select

      select case ( IDPhyTendMethod )
      case ( IDPhyTendMethodImp1LayModel )
        briefexplimp = 'system with surface 1 layer model'
      case ( IDPhyTendMethodImpSoilModel )
        briefexplimp = 'system with thermal diffusion soil model'
      case ( IDPhyTendMethodImpAtmOnly )
        briefexplimp = 'system only with atmosphere'
      case default
        call MessageNotify( 'E', 'Unexpected error in setting briefexplimp', '' )
      end select
    else
      briefexplrad = ''
    end if

    call MessageNotify( 'M', prog_name, '' )
    call MessageNotify( 'M', prog_name,   '+-------------------------------------' )
    call MessageNotify( 'M', prog_name,   '|  Dynamics: %c', c1 = trim(briefexpldyn) )
    call MessageNotify( 'M', prog_name,   '|  Physics : %c', c1 = trim(briefexplphy) )
    if ( IDPhysMode == IDPhysModeFullPhysics ) then
      call MessageNotify( 'M', prog_name, '|    Radiation model : %c', c1 = trim(briefexplrad) )
      call MessageNotify( 'M', prog_name, '|    Surface flux method : %c', c1 = trim(briefexplsfcflux) )
      call MessageNotify( 'M', prog_name, '|    Vertical diffusion method : %c', c1 = trim(briefexplvdiff) )
      call MessageNotify( 'M', prog_name, '|    Implicit method : %c', c1 = trim(briefexplimp) )
      call MessageNotify( 'M', prog_name, '|    Major component phase change : %b', l = (/ FlagMajCompPhaseChange /) )
    end if
    call MessageNotify( 'M', prog_name,   '|  -- version = %c', c1 = trim(version) )
    call MessageNotify( 'M', prog_name,   '+-------------------------------------' )
    call MessageNotify( 'M', prog_name, '' )


    ! Initialization of modules used in this module
    !

    call ProfUtil_Init( namelist_filename )    
    call ProfUtil_RapStart('Setup', 0)

    ! 時刻管理
    ! Time control
    !
    call TimesetInit

    ! 出力ファイルの基本情報管理
    ! Management basic information for output files
    ! 
    call FilesetInit

    ! 格子点設定
    ! Grid points settings
    !
    call GridsetInit

    ! 組成に関わる配列の設定
    ! Settings of array for atmospheric composition
    !
    call CompositionInit

    ! 物理定数設定
    ! Physical constants settings
    !
    call ConstantsInit

    ! 雪と海氷の定数の設定
    ! Setting constants of snow and sea ice
    !
    call ConstantsSnowSeaIceInit

    ! 座標データ設定
    ! Axes data settings
    !
    call AxessetInit



    ! 時系列データの読み込み
    ! Reading time series
    !
    call ReadTimeSeriesInit

    ! 地表面データ提供
    ! Prepare surface data
    !
    call SurfDataInit

    ! ファイルから 1 次元プロファイルを読んで設定する. 
    ! read 1-D profile from a file and set it 
    !
    call Set1DProfileInit

    !
    !
    call OMPWrapperInit



    ! 予報変数の割付
    ! Allocation of prediction variables
    !
    allocate( xyz_UB    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_VB    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_TempB (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyzf_QMixB(0:imax-1, 1:jmax, 1:kmax, 1:ncmax) )
    allocate( xy_PsB    (0:imax-1, 1:jmax) )

    allocate( xyz_UN    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_VN    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_TempN (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyzf_QMixN(0:imax-1, 1:jmax, 1:kmax, 1:ncmax) )
    allocate( xy_PsN    (0:imax-1, 1:jmax) )

    allocate( xyz_UA    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_VA    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_TempA (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyzf_QMixA(0:imax-1, 1:jmax, 1:kmax, 1:ncmax) )
    allocate( xy_PsA    (0:imax-1, 1:jmax) )


    ! リスタートデータ入力
    ! Restart data input
    !
    call RestartFileInit
    call RestartFileGet( &
      & xyz_UB, xyz_VB, xyz_TempB, xyzf_QMixB, xy_PsB, & ! (out)
      & xyz_UN, xyz_VN, xyz_TempN, xyzf_QMixN, xy_PsN, & ! (out)
      & flag_initial )                                   ! (out) optional


    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )

      ! 地表面温度, 土壌温度の割付
      ! Allocation of surface temperature and soil temperature
      !
      allocate( xy_SurfTemp (0:imax-1, 1:jmax) )
      allocate( xyz_SoilTemp(0:imax-1, 1:jmax, 1:kslmax) )

      allocate( xy_SurfMajCompIceB(0:imax-1, 1:jmax) )
      allocate( xy_SoilMoistB     (0:imax-1, 1:jmax) )
      allocate( xy_SurfSnowB      (0:imax-1, 1:jmax) )

      allocate( xy_SurfMajCompIceN(0:imax-1, 1:jmax) )
      allocate( xy_SoilMoistN     (0:imax-1, 1:jmax) )
      allocate( xy_SurfSnowN      (0:imax-1, 1:jmax) )

      allocate( xy_SurfMajCompIceA(0:imax-1, 1:jmax) )
      allocate( xy_SoilMoistA     (0:imax-1, 1:jmax) )
      allocate( xy_SurfSnowA      (0:imax-1, 1:jmax) )


      ! 地表面温度リスタートデータ入力
      ! Restart data of surface temperature input
      !
!!$      call RestartSurfTempGet( &
!!$        & xy_SurfTemp  )          ! (out)
      call RestartSurfTempInit
      call RestartSurfTempGet(                                      &
        & xy_SurfTemp,                                              & ! (out)
        & xyz_SoilTemp,                                             & ! (out)
        & xy_SurfMajCompIceB, xy_SoilMoistB, xy_SurfSnowB,          & ! (out)
        & xy_SurfMajCompIceN, xy_SoilMoistN, xy_SurfSnowN           & ! (out)
        & )

!!$      xy_SurfSnowN = 1.0d1
!!$      xy_SurfSnowB = xy_SurfSnowN

      ! THIS IS A TEMPORARY LINE.
!!$      ! 土壌温度, ..., の初期値設定, いずれリスタートファイルから読むようにする
!!$      ! Setting of initial values of soil temperature, ..., these values are input from restart file in near future
!!$      !
!!$      do k = 1, kslmax
!!$        xyz_SoilTemp(:,:,k) = xy_SurfTemp
!!$      end do
!!$
!!$      xy_SoilMoistN = 0.0_DP
!!$      xy_SurfSnowN  = 0.0_DP
!!$
!!$      xy_SoilMoistB = xy_SoilMoistN
!!$      xy_SurfSnowB  = xy_SurfSnowN
!!$
!!$      xy_SoilMoistA = 0.0_DP
!!$      xy_SurfSnowA  = 0.0_DP

!!$      xy_SurfMajCompIceN = 0.0_DP
!!$      xy_SurfMajCompIceB = xy_SurfMajCompIceN

    end select

    ! リスタートデータファイルの初期化
    ! Initialization of restart data file
    !
    call RestartFileOpen

    ! ヒストリデータファイルの初期化
    ! Initialization of history data files
    !
    call HistoryFileOpen

    ! ヒストリデータ出力のためのへの変数登録
    ! Register of variables for history data output
    !
    call HistoryAutoAddVariable( 'U' , &         ! (in)
      & (/ 'lon ', 'lat ', 'sig ', 'time' /), &  ! (in)
      & 'eastward wind', 'm s-1' )               ! (in)

    call HistoryAutoAddVariable( 'V' , &         ! (in)
      & (/ 'lon ', 'lat ', 'sig ', 'time' /), &  ! (in)
      & 'northward wind', 'm s-1' )              ! (in)

    call HistoryAutoAddVariable( 'Temp' , &      ! (in)
      & (/ 'lon ', 'lat ', 'sig ', 'time' /), &  ! (in)
      & 'temperature', 'K' )                     ! (in)

    do n = 1, ncmax
      call HistoryAutoAddVariable( a_QMixName(n) , & ! (in)
        & (/ 'lon ', 'lat ', 'sig ', 'time' /),    & ! (in)
        & a_QMixLongName(n), 'kg kg-1' )             ! (in)
    end do
    n = IndexH2OVap
    if ( a_QMixName(n) /= 'QVap' ) then
      call HistoryAutoAddVariable( 'QVap' ,        & ! (in)
        & (/ 'lon ', 'lat ', 'sig ', 'time' /),    & ! (in)
        & a_QMixLongName(n), 'kg kg-1' )             ! (in)
    end if

    call HistoryAutoAddVariable( 'Ps' , &        ! (in)
      & (/ 'lon ', 'lat ', 'time' /), &          ! (in)
      & 'surface pressure', 'Pa' )               ! (in)


    ! ヒストリデータ出力 (スタート時刻)
    ! History data output (Start time)
    !
    call HistoryAutoPut( TimeN, 'U', xyz_UN )
    call HistoryAutoPut( TimeN, 'V', xyz_VN )
    call HistoryAutoPut( TimeN, 'Temp', xyz_TempN )
    do n = 1, ncmax
      call HistoryAutoPut( TimeN, a_QMixName(n), xyzf_QMixN(:,:,:,n) )
    end do
    call HistoryAutoPut( TimeN, 'Ps', xy_PsN )
    ! Backward compatibility
    n = IndexH2OVap
    if ( a_QMixName(n) /= 'QVap' ) then
      call HistoryAutoPut( TimeN, 'QVap', xyzf_QMixN(:,:,:,n) )
    end if

    ! Output frequently used variables
    ! Output frequently used variables
    !
    call OutputFreqUsedVarsInit


    select case ( IDPhysMode )
    case ( IDPhysModeJupiterSimpleV2 )

      ! ヒストリデータ出力のためのへの変数登録
      ! Register of variables for history data output
      !
      call HistoryAutoAddVariable( 'Rain', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'precipitation (rain)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'Snow', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'precipitation (snow)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'PRCP', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'precipitation', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'RainCum', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'convective precipitation (rain)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'SnowCum', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'convective precipitation (snow)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'PRCPCum', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'convective precipitation', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'RainLsc', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'non-convective precipitation (rain)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'SnowLsc', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'non-convective precipitation (snow)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'PRCPLsc', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'non-convective precipitation', 'kg m-2 s-1' )

    case ( IDPhysModeFullPhysics )

      ! 地表面温度リスタートデータファイルの初期化
      ! Initialization of restart data file of surface temperature
      !
      call RestartSurfTempOpen

      ! ヒストリデータ出力のためのへの変数登録
      ! Register of variables for history data output
      !
      call HistoryAutoAddVariable( 'SurfTemp' , &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'surface temperature', 'K' )

      call HistoryAutoAddVariable( 'SoilTemp', &
        & (/ 'lon ', 'lat ', 'ssz ', 'time' /), &
        & 'soil temperature', 'K' )

      call HistoryAutoAddVariable( 'SurfMajCompIce' , & ! (in)
        & (/ 'lon ', 'lat ', 'time' /),               & ! (in)
        & 'surface major component ice', 'kg m-2' )     ! (in)

      call HistoryAutoAddVariable( 'SoilMoist' , & ! (in)
        & (/ 'lon ', 'lat ', 'time' /), &          ! (in)
        & 'soil moisture', 'kg m-2' )              ! (in)

      call HistoryAutoAddVariable( 'SurfSnow'  , & ! (in)
        & (/ 'lon ', 'lat ', 'time' /), &          ! (in)
        & 'surface snow amount', 'kg m-2' )        ! (in)

      call HistoryAutoAddVariable( 'Rain', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'precipitation (rain)', 'kg m-2 s-1' )

      call HistoryAutoAddVariable( 'Snow', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'precipitation (snow)', 'kg m-2 s-1' )

      call HistoryAutoAddVariable( 'PRCP', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'precipitation', 'kg m-2 s-1' )

      call HistoryAutoAddVariable( 'RainCum', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'convective precipitation (rain)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'SnowCum', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'convective precipitation (snow)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'PRCPCum', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'convective precipitation', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'RainLsc', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'non-convective precipitation (rain)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'SnowLsc', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'non-convective precipitation (snow)', 'kg m-2 s-1' )
      call HistoryAutoAddVariable( 'PRCPLsc', &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'non-convective precipitation', 'kg m-2 s-1' )

      ! ヒストリデータ出力 (スタート時刻)
      ! History data output (Start time)
      !
      call HistoryAutoPut( TimeN, 'SurfMajCompIce', xy_SurfMajCompIceN )
      call HistoryAutoPut( TimeN, 'SoilMoist'     , xy_SoilMoistN      )
      call HistoryAutoPut( TimeN, 'SurfSnow'      , xy_SurfSnowN       )
    end select


    ! 診断変数の割付
    ! Allocation of diagnostic variables
    !
    allocate( xyz_DUDt    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_DVDt    (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_DTempDt (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xy_DPsDt    (0:imax-1, 1:jmax) )
    allocate( xyzf_DQMixDt(0:imax-1, 1:jmax, 1:kmax, 1:ncmax) )

    allocate( xyz_DTurKinEneDt(0:imax-1, 1:jmax, 1:kmax) )

    allocate( xyz_OMG(0:imax-1, 1:jmax, 1:kmax) )

    allocate( xy_SurfHeight(0:imax-1, 1:jmax) )
    allocate( xyz_Height   (0:imax-1, 1:jmax, 1:kmax) )
    allocate( xyz_Exner    (0:imax-1, 1:jmax, 1:kmax) )

    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )
      allocate( xy_SurfAlbedo         (0:imax-1, 1:jmax) )
      allocate( xy_SurfHumidCoef      (0:imax-1, 1:jmax) )
      allocate( xy_SurfRoughLenMom    (0:imax-1, 1:jmax) )
      allocate( xy_SurfRoughLenHeat   (0:imax-1, 1:jmax) )
      allocate( xy_SurfHeatCapacity   (0:imax-1, 1:jmax) )
      allocate( xy_SeaIceConc         (0:imax-1, 1:jmax) )
      allocate( xy_SurfCond           (0:imax-1, 1:jmax) )
      allocate( xy_SurfType           (0:imax-1, 1:jmax) )
      allocate( xy_DeepSubSurfHeatFlux(0:imax-1, 1:jmax) )
      allocate( xy_SoilHeatCap        (0:imax-1, 1:jmax) )
      allocate( xy_SoilHeatDiffCoef   (0:imax-1, 1:jmax) )
      allocate( xy_SurfHeightStd      (0:imax-1, 1:jmax) )
      allocate( xy_SnowFrac           (0:imax-1, 1:jmax) )

!!$      allocate( xy_FlagMatthewsLand         (0:imax-1, 1:jmax) )
      allocate( xy_PhyImplSDHIndexCalcMethod(0:imax-1, 1:jmax) )
      allocate( xy_BucketFlagOceanGrid      (0:imax-1, 1:jmax) )

      allocate( xyr_Temp   (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyz_VirTemp   (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_VirTemp   (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xy_SurfVirTemp(0:imax-1, 1:jmax) )
      allocate( xyz_Press  (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_Press  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_Height (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_Exner  (0:imax-1, 1:jmax, 0:kmax) )

      allocate( xyr_RadLFlux      (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLFluxA     (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLUwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLDwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSFlux      (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSUwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSDwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyra_DelRadLFlux  (0:imax-1, 1:jmax, 0:kmax,  0:1) )
      allocate( xyra_DelRadLUwFlux(0:imax-1, 1:jmax, 0:kmax,  0:1) )
      allocate( xyra_DelRadLDwFlux(0:imax-1, 1:jmax, 0:kmax,  0:1) )

      allocate( xyr_MomFluxX  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_MomFluxY  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_HeatFlux  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyrf_QMixFlux (0:imax-1, 1:jmax, 0:kmax, 1:ncmax) )

      allocate( xy_SurfMomFluxX  (0:imax-1, 1:jmax) )
      allocate( xy_SurfMomFluxY  (0:imax-1, 1:jmax) )
      allocate( xy_SurfHeatFlux  (0:imax-1, 1:jmax) )
      allocate( xyf_SurfQMixFlux (0:imax-1, 1:jmax, 1:ncmax) )

      allocate( xy_SurfH2OVapFluxA    (0:imax-1, 1:jmax) )
      allocate( xy_SurfLatentHeatFluxA(0:imax-1, 1:jmax) )

      allocate( xyr_SoilHeatFlux(0:imax-1, 1:jmax, 0:kslmax) )

      allocate( xyr_VelDiffCoef (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_TempDiffCoef(0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_QMixDiffCoef(0:imax-1, 1:jmax, 0:kmax) )

      allocate( xy_SurfVelTransCoef (0:imax-1, 1:jmax) )
      allocate( xy_SurfTempTransCoef(0:imax-1, 1:jmax) )
      allocate( xy_SurfQVapTransCoef(0:imax-1, 1:jmax) )

      allocate( xy_SurfMOLength     (0:imax-1, 1:jmax) )

      allocate( xyr_SoilTempTransCoef (0:imax-1, 1:jmax, 0:kslmax) )

      allocate( xy_DSurfTempDt (0:imax-1, 1:jmax) )
      allocate( xyz_DSoilTempDt(0:imax-1, 1:jmax, 1:kslmax) )

      allocate( xy_DSurfMajCompIceDt(0:imax-1, 1:jmax) )
      allocate( xy_DSoilMoistDt     (0:imax-1, 1:jmax) )
      allocate( xy_DSurfSnowDt      (0:imax-1, 1:jmax) )

      allocate( xyz_DTempDtVDiff(0:imax-1, 1:jmax, 1:kmax) )

      allocate( xyz_DUDtGWD(0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyz_DVDtGWD(0:imax-1, 1:jmax, 1:kmax) )

      allocate( xyz_DTempDtRadL (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyz_DTempDtRadS (0:imax-1, 1:jmax, 1:kmax) )


      allocate( xy_Rain         (0:imax-1, 1:jmax) )
      allocate( xy_RainCumulus  (0:imax-1, 1:jmax) )
      allocate( xy_RainLsc      (0:imax-1, 1:jmax) )
      allocate( xy_Snow         (0:imax-1, 1:jmax) )
      allocate( xy_SnowCumulus  (0:imax-1, 1:jmax) )
      allocate( xy_SnowLsc      (0:imax-1, 1:jmax) )

      allocate( xyz_DTempDtCum    (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQVapDtCum    (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQH2OLiqDtCum (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQH2OSolDtCum (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DUDtCum       (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DVDtCum       (0:imax-1,1:jmax,1:kmax) )

      allocate( xyz_DQH2OLiqDtLSC (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQH2OSolDtLSC (0:imax-1,1:jmax,1:kmax) )

      allocate( xyz_MoistConvDetTend       (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyz_MoistConvSubsidMassFlux(0:imax-1, 1:jmax, 1:kmax) )

      allocate( xyz_QH2OLiqforRad   (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_QH2OSolforRad   (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_CloudCoverforRad(0:imax-1,1:jmax,1:kmax) )

      allocate( xy_SurfDustGravSedFlux(0:imax-1,1:jmax) )
      
      ! ヒストリデータ出力のためのへの変数登録
      ! Register of variables for history data output
      !
      call HistoryAutoAddVariable( 'SeaIceConc' , &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'sea ice concentration', '1' )

      call HistoryAutoAddVariable( 'SurfAlbedo' , &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'surface albedo', '1' )

      call HistoryAutoAddVariable( 'SurfRoughLenMom' , &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'surface roughness length for momentum', 'm' )

      call HistoryAutoAddVariable( 'SurfRoughLenHeat' , &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'surface roughness length for heat', 'm' )

      call HistoryAutoAddVariable( 'SurfDustGravSedFlux' , &
        & (/ 'lon ', 'lat ', 'time' /), &
        & 'surface roughness length', 'm' )

    case ( IDPhysModeVenusSimple )

      allocate( xyr_Temp   (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyz_VirTemp   (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_VirTemp   (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyz_Press  (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_Press  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_Height (0:imax-1, 1:jmax, 0:kmax) )
!!$      allocate( xyz_Exner  (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_Exner  (0:imax-1, 1:jmax, 0:kmax) )

    case ( IDPhysModeJupiterSimple )

      allocate( xyz_Press  (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_Press  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_VirTemp(0:imax-1, 1:jmax, 0:kmax) )

      allocate( xyz_DTempDtVDiff(0:imax-1, 1:jmax, 1:kmax) )

      allocate( xyr_HeatFlux      (0:imax-1, 1:jmax, 0:kmax) )

      allocate( xyr_RadLFlux      (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLUwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLDwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSFlux      (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSUwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSDwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyz_DTempDtRadL   (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyz_DTempDtRadS   (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyra_DelRadLFlux  (0:imax-1, 1:jmax, 0:kmax,  0:1) )
      allocate( xyra_DelRadLUwFlux(0:imax-1, 1:jmax, 0:kmax,  0:1) )
      allocate( xyra_DelRadLDwFlux(0:imax-1, 1:jmax, 0:kmax,  0:1) )


    case ( IDPhysModeJupiterSimpleV2 )

      allocate( xyz_Press  (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_Press  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_VirTemp(0:imax-1, 1:jmax, 0:kmax) )
!!$      allocate( xyz_Exner  (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyr_Exner  (0:imax-1, 1:jmax, 0:kmax) )

      allocate( xyz_DTempDtVDiff(0:imax-1, 1:jmax, 1:kmax) )

      allocate( xyr_MomFluxX  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_MomFluxY  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_HeatFlux  (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyrf_QMixFlux (0:imax-1, 1:jmax, 0:kmax, 1:ncmax) )

      allocate( xyr_VelDiffCoef (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_TempDiffCoef(0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_QMixDiffCoef(0:imax-1, 1:jmax, 0:kmax) )

      allocate( xy_SurfVelTransCoef (0:imax-1, 1:jmax) )
      allocate( xy_SurfTempTransCoef(0:imax-1, 1:jmax) )
      allocate( xy_SurfQVapTransCoef(0:imax-1, 1:jmax) )

      allocate( xyr_RadLFlux      (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLUwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadLDwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSFlux      (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSUwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyr_RadSDwFlux    (0:imax-1, 1:jmax, 0:kmax) )
      allocate( xyz_DTempDtRadL   (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyz_DTempDtRadS   (0:imax-1, 1:jmax, 1:kmax) )
      allocate( xyra_DelRadLFlux  (0:imax-1, 1:jmax, 0:kmax,  0:1) )
      allocate( xyra_DelRadLUwFlux(0:imax-1, 1:jmax, 0:kmax,  0:1) )
      allocate( xyra_DelRadLDwFlux(0:imax-1, 1:jmax, 0:kmax,  0:1) )

      allocate( xy_Rain         (0:imax-1, 1:jmax) )
      allocate( xy_RainCumulus  (0:imax-1, 1:jmax) )
      allocate( xy_RainLsc      (0:imax-1, 1:jmax) )
      allocate( xy_Snow         (0:imax-1, 1:jmax) )
      allocate( xy_SnowCumulus  (0:imax-1, 1:jmax) )
      allocate( xy_SnowLsc      (0:imax-1, 1:jmax) )

      allocate( xyz_DTempDtCum    (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQVapDtCum    (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQH2OLiqDtCum (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQH2OSolDtCum (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DUDtCum       (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DVDtCum       (0:imax-1,1:jmax,1:kmax) )

      allocate( xyz_DQH2OLiqDtLSC (0:imax-1,1:jmax,1:kmax) )
      allocate( xyz_DQH2OSolDtLSC (0:imax-1,1:jmax,1:kmax) )

    end select

    !* For coupler ----------------------------------------------------------------
    allocate( xy_TauXAtm(0:imax-1,jmax), xy_TauYAtm(0:imax-1,jmax) )
    allocate( xy_SensAtm(0:imax-1,jmax), xy_LatentAtm(0:imax-1,jmax) )
    allocate( xy_LDWRFlxAtm(0:imax-1,jmax), xy_LUWRFlxAtm(0:imax-1,jmax) )
    allocate( xy_SDWRFlxAtm(0:imax-1,jmax), xy_SUWRFlxAtm(0:imax-1,jmax) )
    allocate( xy_SurfAirTemp(0:imax-1,jmax) )
    allocate( xy_DSurfHFlxDTs(0:imax-1,jmax), xy_DSurfLatentFlxDTs(0:imax-1,jmax) )
    allocate( xy_RainAtm(0:imax-1,jmax), xy_SnowAtm(0:imax-1,jmax) )
    !-------------------------------------------------------------------------------


    ! 惑星表面データの設定
    ! Setting planetary surface properties
    !
    select case ( IDSfcMoistMethod )
    case ( IDSfcMoistMethodNone )
      call SurfacePropertiesInit(                    &
        & FlagPhysImpSoilModelSO, .false., FlagSnow  & ! (in)
        & )
    case ( IDSfcMoistMethodBucket )
      call SurfacePropertiesInit(                    &
        & FlagPhysImpSoilModelSO, .true. , FlagSnow  & ! (in)
        & )
    end select


    select case ( IDPhysMode )
    case ( IDPhysModeNoPhysics )

      ! Nothing to do

    case ( IDPhysModeHS94 )

      ! Held and Suarez (1994) による強制と散逸
      ! Forcing and dissipation suggested by Held and Suarez (1994)
      !
      call HS94Init

    case ( IDPhysModeVenusSimple )

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      ! Yamamoto and Takahashi (2003) に従った簡単金星計算のための強制
      ! forcing for simple Venus calculation following Yamamoto and Takahashi (2003)
      !
      call YT2003ForcingInit

    case ( IDPhysModeJupiterSimple )

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      ! Schneider and Liu (2009) による鉛直混合課程
      ! Vertical diffusion by Schneider and Liu (2009)
      !
      call SL09DiffusionInit

      ! Schneider and Liu (2009) の放射モデル
      ! Radiation model by Schneider and Liu (2009)
      !
      call RadSL09Init

      ! 放射関連ルーチン
      ! Routines for radiation calculation
      !
      call RadUtilsInit

    case ( IDPhysModeJupiterSimpleV2 )

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      ! 下部境界フラックス
      ! Lower boundary flux
      !
      call LBFluxSimpleInit

      ! 陰解法による時間積分 (大気のみ / 惑星表面温度・土壌温度計算なし)
      ! Time integration by using implicit scheme in case without calculation of surface and soil temperature
      !
      call PhyImplAtmOnlyInit

      ! Schneider and Liu (2009) の放射モデル
      ! Radiation model by Schneider and Liu (2009)
      !
      call RadSL09Init

      ! 放射関連ルーチン
      ! Routines for radiation calculation
      !
      call RadUtilsInit

    case ( IDPhysModeFullPhysics )

      ! 惑星表面データの設定
      ! Setting planetary surface properties
      !
      select case ( IDSfcMoistMethod )
      case ( IDSfcMoistMethodNone )
        call SurfacePropertiesInit(                    &
          & FlagPhysImpSoilModelSO, .false., FlagSnow  & ! (in)
          & )
      case ( IDSfcMoistMethodBucket )
        call SurfacePropertiesInit(                    &
          & FlagPhysImpSoilModelSO, .true. , FlagSnow  & ! (in)
          & )
      end select

      ! 雪, 氷の割合
      ! snow/ice fraction
      !
      call SnowIceFracInit

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      select case ( IDRadMethod )
      case ( IDRadMethodDennouAGCM )

        ! 放射フラックス (GFD 電脳倶楽部開発の放射モデル)
        ! Radiation flux (radiation model developed by GFD Dennou Club)
        !
        call RadDennouAGCMInit( flag_rst = .not. flag_initial )

      case ( IDRadMethodEarthV2 )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        select case ( IDCloudMethod )
        case ( IDCloudMethodNone )
          ! 雲なしモデル
          ! No cloud model
          !
          call CloudNoneInit( &
            & FlagSnow             &
            & )
        case ( IDCloudMethodSimple )
          ! 簡単雲モデル
          ! Simple cloud
          !
          call CloudSimpleInit( &
            & FlagSnow          &
            & )
        case ( IDCloudMethodSimpleWithIce )
          ! 簡単雲モデル
          ! Simple cloud
          !
          call CloudSimpleInit( &
            & FlagSnow          &
            & )
        case ( IDCloudMethodT1993WithIce )
          ! Tiedtke (1993) に基づく雲モデル
          ! Cloud model based on Tiedtke (1993)
          !
          call CloudT1993baseInit( &
            & FlagSnow                &
            & )
        end select

        ! 地球大気向け放射モデル Ver. 2
        ! radiation model for the Earth's atmosphere Ver. 2
        !
        call RadEarthV2Init( &
          & FlagSnow         &
          & )

      case ( IDRadMethodRRTMG )

        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if

        select case ( IDCloudMethod )
        case ( IDCloudMethodNone )
          ! 雲なしモデル
          ! No cloud model
          !
          call CloudNoneInit( &
            & FlagSnow             &
            & )
        case ( IDCloudMethodSimple )
          ! 簡単雲モデル
          ! Simple cloud
          !
          call CloudSimpleInit( &
            & FlagSnow          &
            & )
        case ( IDCloudMethodSimpleWithIce )
          ! 簡単雲モデル
          ! Simple cloud
          !
          call CloudSimpleInit( &
            & FlagSnow          &
            & )
        case ( IDCloudMethodT1993WithIce )
          ! Tiedtke (1993) に基づく雲モデル
          ! Cloud model based on Tiedtke (1993)
          !
          call CloudT1993baseInit( &
            & FlagSnow                &
            & )
        end select

        ! wrapper of RRTMG
        ! wrapper of RRTMG
        !
        call RadRRTMGWrapperInit(    &
          & FlagSnow                  &
          & )

      case ( IDRadMethodMarsV1 )

        ! 火星大気向け放射モデル Ver. 1
        ! radiation model for the Mars' atmosphere Ver. 1
        !
        call RadMarsV1Init

      case ( IDRadMethodSL09 )

        ! Schneider and Liu (2009) の放射モデル
        ! Radiation model by Schneider and Liu (2009)
        !
        call RadSL09Init

      case ( IDRadMethodSimple )

        ! 簡単放射モデル
        ! Simple radiation model
        !
        call RadSimpleInit

      case ( IDRadMethodNone )

        ! 何もしない放射モデル
        ! radiation model with no absorption and no scattering
        !
        call RadNoneInit

      end select


      ! 鉛直拡散フラックス
      ! Vertical diffusion flux
      !
      select case ( IDVDiffMethod )
      case ( IDVDiffMethodMY2, IDVDiffMethodMY25 )

        ! 鉛直拡散フラックス (Mellor and Yamada, 1974)
        ! Vertical diffusion flux (Mellor and Yamada, 1974)
        !
        call VDiffusionInit

      case ( IDVDiffMethodJMA )
        ! JMA 乱流混合モジュール
        ! JMA turbulent mixing module
        !
        call VDiffusionJMAInit

      end select

      ! 地表面フラックス (バルク法)
      ! Surface flux (Bulk method)
      !
      call SurfaceFluxInit

      !
      ! set dust flux
      !
!!$      call SetDustFluxInit

      ! 一部の物理過程の時間変化率の計算 (陰解法)
      ! Calculate tendency by a part of physical processes (implicit)
      !
      select case ( IDPhyTendMethod )
      case ( IDPhyTendMethodImp1LayModel )

        ! 陰解法による時間積分
        ! Time integration with implicit scheme
        !
        call PhyImplInit(             &
          & FlagBucketModel, FlagSnow &
          & )

      case ( IDPhyTendMethodImpSoilModel )

        select case ( IDSfcMoistMethod )
        case ( IDSfcMoistMethodBucket )
          ! バケツモデル
          ! Bucket model
          !
          call BucketModelInit( &
            & FlagSnow          & ! (in)
            & )

        end select

        ! 地下における熱の鉛直拡散
        ! Vertical diffusion of heat under the ground
        !
        call SubsurfaceDiffusionInit

        ! 地下熱伝導モデルを用いた場合の陰解法による時間積分
        !
        ! Time integration by using implicit scheme in case using subsurface thermal diffusion model
!!$        call PhyImplSDHInit(                                                 &
!!$          & FlagBucketModel, FlagSnow,                                       &
!!$          & FlagPhysImpSoilModelSO, FlagMajCompPhaseChange, CondMajCompName  &
!!$          & )

!!$        call PhyImplSDHV2Init(                                                &
!!$          & FlagBucketModel, FlagSnow,                                       &
!!$          & FlagPhysImpSoilModelSO, FlagMajCompPhaseChange, CondMajCompName  &
!!$          & )
        call PhyImplSDHV3Init(                                                &
          & FlagBucketModel, FlagSnow,                                       &
          & FlagPhysImpSoilModelSO, FlagMajCompPhaseChange, CondMajCompName  &
          & )

      case ( IDPhyTendMethodImpAtmOnly )

        ! 陰解法による時間積分 (大気のみ / 惑星表面温度・土壌温度計算なし)
        ! Time integration by using implicit scheme in case without calculation of surface and soil temperature
        !
        call PhyImplAtmOnlyInit

      end select


      ! Gravity wave drag by McFarlane (1987)
      ! Gravity wave drag by McFarlane (1987)
      !
      call GWDM1987Init

      ! 陰解法による時間積分のためのルーチン
      ! Routines for time integration with implicit scheme
      !
      call PhyImplUtilsInit

      ! 放射関連ルーチン
      ! Routines for radiation calculation
      !
      call RadUtilsInit


!!$      select case ( IDRadMethod )
!!$      case ( IDRadMethodMarsV1 )
!!$        ! (火星大気向け) Non-LTE 放射モデル
!!$        ! Non-NLTE radiation model (for the Mars' atmosphere)
!!$        !
!!$        call Rad15mNLTEInit
!!$      end select

      select case ( IDRadMethod )
      case ( IDRadMethodMarsV1 )
        ! 火星計算用近赤外加熱率計算
        ! Calculation of near infrared heating rate in the case of Mars
        !
        call RadMarsNIRInit
      end select


      ! 鉛直拡散フラックス (Mellor and Yamada, 1974)
      ! Vertical diffusion flux (Mellor and Yamada, 1974)
      !
      call VDiffusionInit

      ! 地表面フラックス (バルク法)
      ! Surface flux (Bulk method)
      !
      call SurfaceFluxInit

      ! 放射関連ルーチン
      ! Routines for radiation calculation
      !
      call RadUtilsInit

    end select



    ! 地面温度・土壌温度・土壌水分・積雪量の積分
    ! Time integration of surface temperature, soil temperature, soil moisture, 
    ! and surface snow amount
    !
    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )

      ! 地面温度, 土壌温度の時間積分
      ! Time integration of surface temperature and soil temperature
      !
      call IntgSurfTempInit

      select case ( IDSfcMoistMethod )
      case ( IDSfcMoistMethodBucket )
        ! バケツモデル
        ! Bucket model
        !
        call BucketModelInit( &
          & FlagSnow          & ! (in)
          & )
      end select

    end select




    ! 力学過程
    ! Dynamical core
    !
    select case ( IDDynMode )
    case ( IDDynModeHSPLVAS83 )
      ! 力学過程 (スペクトル法, Arakawa and Suarez (1983))
      ! Dynamical process (Spectral method, Arakawa and Suarez (1983))
      !
      call DynamicsHSplVAS83Init
    case ( IDDynModeNoHorAdv )
      ! 物理過程のみの計算のための力学過程
      ! A dynamics for calculation with physical processes only
      !
      call DynamicsPhysicsOnlyInit
    case ( IDDynModeTWPICE )
      !
      ! Dynamical process for TWP-ICE experiment
      !
      call DynamicsTWPICESCMExpInit
    end select


    select case ( IDPhysMode )
    case ( IDPhysModeJupiterSimple )

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      select case ( IDDCMethod )
      case ( IDDCMethodDCA )
        ! 乾燥対流調節
        ! Dry convective adjustment
        !
        call DryConvAdjustInit
      end select

    case ( IDPhysModeJupiterSimpleV2 )

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      ! 湿潤対流調節
      ! Moist convective adjustment
      !
      call MoistConvAdjustInit

      ! 大規模凝結 (非対流性凝結) (Manabe, 1965)
      ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
      !
      call LScaleCondInit( &
        & FlagSnow &
        & )

      ! 簡単雲モデル
      ! Simple cloud
      !
      call CloudSimpleInit( &
        & FlagSnow          &
        & )

      select case ( IDDCMethod )
      case ( IDDCMethodDCA )
        ! 乾燥対流調節
        ! Dry convective adjustment
        !
        call DryConvAdjustInit
      end select

      ! 質量の補正
      ! Mass fixer
      !
      call MassFixerInit

    case ( IDPhysModeFullPhysics )

      ! 補助的な変数を計算するサブルーチン・関数群
      ! Subroutines and functions for calculating auxiliary variables
      !
      call AuxVarsInit

      select case ( IDMCMethod )
      case ( IDMCMethodMCA, IDMCMethodMCAI98 )
        ! 湿潤対流調節
        ! Moist convective adjustment
        !
        call MoistConvAdjustInit
      case ( IDMCMethodRAS )
        ! Relaxed Arakawa-Schubert scheme
        ! Relaxed Arakawa-Schubert scheme
        !
        call RASInit
      case ( IDMCMethodRASWithIce )
        ! Relaxed Arakawa-Schubert scheme
        ! Relaxed Arakawa-Schubert scheme
        !
        call RASInit
      end select

      select case ( IDLSCMethod )
      case ( IDLSCMethodM65 )
        ! 大規模凝結 (非対流性凝結) (Manabe, 1965)
        ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
        !
        call LScaleCondInit( &
          & FlagSnow &
          & )
!!$      case ( IDLSCMethodLL91 )
!!$        ! 大規模凝結 (非対流性凝結) (Le Treut and Li, 1991)
!!$        ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
!!$        !
!!$        call LScaleCondLL91Init
      case ( IDLSCMethodSatAdjM65 )
        !
        != Saturation adjustment
        !
        if ( IndexH2OLiq <= 0 ) then
          call MessageNotify( 'E', prog_name, &
            & trim(a_QMixName(IndexH2OLiq))//' is not found.' )
        end if
        call SaturationAdjustInit( FlagSnow )
      case ( IDLSCMethodM65WithIce )
        ! 大規模凝結 (非対流性凝結) (Manabe, 1965)
        ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
        !
        call LScaleCondInit( &
          & FlagSnow &
          & )
      case ( IDLSCMethodLL91WithIce )
        ! 大規模凝結 (非対流性凝結) (Le Treut and Li, 1991)
        ! Large scale condensation (non-convective condensation) (Le Treut and Li, 1991)
        !
        call LScaleCondInit( &
          & FlagSnow &
          & )
      end select


      ! 
      ! Cloud model
      !
      select case ( IDCloudMethod )
      case ( IDCloudMethodNone )
        ! 雲なしモデル
        ! No cloud model
        !
        call CloudNoneInit( &
          & FlagSnow             &
          & )
      case ( IDCloudMethodSimple )
        ! 簡単雲モデル
        ! Simple cloud
        !
        call CloudSimpleInit( &
          & FlagSnow          &
          & )
      case ( IDCloudMethodSimpleWithIce )
        ! 簡単雲モデル
        ! Simple cloud
        !
        call CloudSimpleInit( &
          & FlagSnow          &
          & )
      case ( IDCloudMethodT1993WithIce )
        ! Tiedtke (1993) に基づく雲モデル
        ! Cloud model based on Tiedtke (1993)
        !
        call CloudT1993baseInit( &
          & FlagSnow                &
          & )
      case ( IDCloudMethodMarsH2OCloud )
        ! 火星 H2O 雲モデル
        ! Mars H2O cloud model
        !
        call CloudMarsH2OInit
      end select


      select case ( IDSfcMoistMethod )
      case ( IDSfcMoistMethodBucket )
        ! バケツモデル
        ! Bucket model
        !
        call BucketModelInit( &
          & FlagSnow          & ! (in)
          & )
      end select

      select case ( IDDCMethod )
      case ( IDDCMethodDCA )
        ! 乾燥対流調節
        ! Dry convective adjustment
        !
        call DryConvAdjustInit
      end select

      ! 重力沈降過程
      ! Gravitational sedimentation process
      !
      call GravSedInit

      ! 主成分相変化
      ! Phase change of atmospheric major component
      !
      call MajorCompPhaseChangeInit(                &
        & FlagMajCompPhaseChange, CondMajCompName   & ! (in)
        & )

      ! 質量の補正
      ! Mass fixer
      !
      call MassFixerInit

    end select


    ! タイムフィルター (Asselin, 1972)
    ! Time filter (Asselin, 1972)
    !
!!$    call TimeFiltInit

    ! 時間フィルター (Williams, 2009)
    ! Time filter (Williams, 2009)
    !
    call TimeFilterWilliams2009Init

    ! 予報変数の値の確認
    ! Check values of prognostic variables
    !
    call CheckProgVarsInit

    ! 補助的な変数を計算するサブルーチン・関数群
    ! Subroutines and functions for calculating auxiliary variables
    !
    call AuxVarsInit

    !
    ! End of initialization
    !
    call ProfUtil_RapEnd('Setup', 0)


    ! 初回だけはオイラー法を用いるため, Δt を半分に
    ! Delta t is reduced to half in order to use Euler method at initial step
    !
    if ( flag_initial ) then
      call TimesetDelTimeHalf
    end if

  end subroutine MainInit

  !-------------------------------------------------------------------

  subroutine MainTerminate
    !
    ! 主プログラムの終了処理手続き. 
    !
    ! Termination procedure for the main program. 
    !

    ! モジュール引用 ; USE statements
    !

    ! MPI
    !
    use mpi_wrapper, only : MPIWrapperFinalize

    ! 力学過程 (スペクトル法, Arakawa and Suarez (1983))
    ! Dynamical process (Spectral method, Arakawa and Suarez (1983))
    !
    use dynamics_hspl_vas83, only : DynamicsHSplVAS83Finalize

    !
    ! Dynamical process for TWP-ICE experiment
    !
    use dynamics_twpice_scm_exp, only : DynamicsTWPICESCMExpFinalize

    ! Held and Suarez (1994) による強制と散逸
    ! Forcing and dissipation suggested by Held and Suarez (1994)
    !
    use held_suarez_1994, only: Hs94Finalize

    ! 放射フラックス (バンドモデル)
    ! Radiation flux (band model)
    !
    use rad_DennouAGCM, only: RadDennouAGCMFinalize

    ! 座標データ設定
    ! Axes data settings
    !
    use axesset, only: AxessetFinalize

    ! 温度の半整数σレベルの補間, 気圧と高度の算出
    ! Interpolate temperature on half sigma level, 
    ! and calculate pressure and height
    !
    use auxiliary, only: AuxVarsFinalize

    ! 時刻管理
    ! Time control
    !
    use timeset, only: TimesetClose

    ! リスタートデータ入出力
    ! Restart data input/output
    !
    use restart_file_io, only: RestartFileClose

    ! 地表面温度リスタートデータ入出力
    ! Restart data of surface temperature input/output
    !
    use restart_surftemp_io, only: RestartSurfTempClose

    ! ヒストリデータ出力
    ! History data output
    !
    use history_file_io, only: HistoryFileClose

    use dc_string, only: CPrintf
    use mpi_wrapper, only: myrank
    
    ! 宣言文 ; Declaration statements
    !
    implicit none

    character(STRING) :: profileName
    
    ! 実行文 ; Executable statement
    !

    call ProfUtil_RapStart('Shutdown', 0)    
    ! リスタートデータファイルクローズ
    ! Close restart data file
    !
    call RestartFileClose

    if ( IDPhysMode == IDPhysModeFullPhysics ) then
      ! 地表面温度リスタートデータファイルクローズ
      ! Close restart data file of surface temperature
      !
      call RestartSurfTempClose
    end if

    ! ヒストリデータファイルクローズ
    ! Close history data files
    !
    call HistoryFileClose

    ! 予報変数の割付解除
    ! Deallocation of prediction variables
    !
    deallocate( xyz_UB     )
    deallocate( xyz_VB     )
    deallocate( xyz_TempB  )
    deallocate( xyzf_QMixB )
    deallocate( xy_PsB     )

    deallocate( xyz_UN     )
    deallocate( xyz_VN     )
    deallocate( xyz_TempN  )
    deallocate( xyzf_QMixN )
    deallocate( xy_PsN     )

    deallocate( xyz_UA     )
    deallocate( xyz_VA     )
    deallocate( xyz_TempA  )
    deallocate( xyzf_QMixA )
    deallocate( xy_PsA     )

    ! 診断変数の割付解除
    ! Dellocation of diagnostic variables
    !
    deallocate( xyz_DUDt     )
    deallocate( xyz_DVDt     )
    deallocate( xyz_DTempDt  )
    deallocate( xyzf_DQMixDt )

    deallocate( xyz_DTurKinEneDt )

    deallocate( xyz_OMG )

    deallocate( xy_SurfHeight )
    deallocate( xyz_Height )

    if ( IDPhysMode == IDPhysModeFullPhysics ) then
      deallocate( xy_SurfAlbedo          )
      deallocate( xy_SurfHumidCoef       )
      deallocate( xy_SurfRoughLenMom     )
      deallocate( xy_SurfRoughLenHeat    )
      deallocate( xy_SurfHeatCapacity    )
      deallocate( xy_SeaIceConc          )
      deallocate( xy_SurfCond            )
      deallocate( xy_SurfType            )
      deallocate( xy_DeepSubSurfHeatFlux )
      deallocate( xy_SoilHeatCap         )
      deallocate( xy_SoilHeatDiffCoef    )
      deallocate( xy_SurfHeightStd       )
      deallocate( xy_SnowFrac            )

!!$      deallocate( xy_FlagMatthewsLand          )
      deallocate( xy_PhyImplSDHIndexCalcMethod )
      deallocate( xy_BucketFlagOceanGrid       )

      deallocate( xyr_Temp   )
      deallocate( xyz_VirTemp )
      deallocate( xyr_VirTemp )
      deallocate( xy_SurfVirTemp )
      deallocate( xyz_Press  )
      deallocate( xyr_Press  )
      deallocate( xyr_Height )
      deallocate( xyz_Exner  )
      deallocate( xyr_Exner  )

      deallocate( xyr_RadLFlux     )
      deallocate( xyr_RadLFluxA    )
      deallocate( xyr_RadLUwFlux   )
      deallocate( xyr_RadLDwFlux   )
      deallocate( xyr_RadSFlux     )
      deallocate( xyr_RadSUwFlux     )
      deallocate( xyr_RadSDwFlux     )
      deallocate( xyra_DelRadLUwFlux )
      deallocate( xyra_DelRadLDwFlux )

      deallocate( xyr_MomFluxX  )
      deallocate( xyr_MomFluxY  )
      deallocate( xyr_HeatFlux  )
      deallocate( xyrf_QMixFlux )

      deallocate( xy_SurfMomFluxX  )
      deallocate( xy_SurfMomFluxY  )
      deallocate( xy_SurfHeatFlux  )
      deallocate( xyf_SurfQMixFlux )

      deallocate( xy_SurfH2OVapFluxA     )
      deallocate( xy_SurfLatentHeatFluxA )

      deallocate( xyr_SoilHeatFlux )

      deallocate( xyr_VelDiffCoef  )
      deallocate( xyr_TempDiffCoef )
      deallocate( xyr_QMixDiffCoef )

      deallocate( xy_SurfVelTransCoef  )
      deallocate( xy_SurfTempTransCoef )
      deallocate( xy_SurfQVapTransCoef )

      deallocate( xy_SurfMOLength )

      deallocate( xy_DSurfTempDt )
      deallocate( xyz_DSoilTempDt )

      deallocate( xy_DPsDt )

      deallocate( xy_DSurfMajCompIceDt )
      deallocate( xy_DSoilMoistDt      )
      deallocate( xy_DSurfSnowDt       )

      deallocate( xyz_DTempDtVDiff )

      deallocate( xyz_DUDtGWD )
      deallocate( xyz_DVDtGWD )

      deallocate( xyz_DTempDtRadL )
      deallocate( xyz_DTempDtRadS )


      deallocate( xy_Rain          )
      deallocate( xy_RainCumulus   )
      deallocate( xy_RainLsc       )
      deallocate( xy_SnowCumulus   )
      deallocate( xy_SnowLsc       )

      deallocate( xyz_DTempDtCum     )
      deallocate( xyz_DQVapDtCum     )
      deallocate( xyz_DQH2OLiqDtCum  )
      deallocate( xyz_DQH2OSolDtCum  )
      deallocate( xyz_DUDtCum        )
      deallocate( xyz_DVDtCum        )

      deallocate( xyz_DQH2OLiqDtLSC  )
      deallocate( xyz_DQH2OSolDtLSC  )

      deallocate( xyz_MoistConvDetTend        )
      deallocate( xyz_MoistConvSubsidMassFlux )

      deallocate( xyz_QH2OLiqforRad    )
      deallocate( xyz_QH2OSolforRad    )
      deallocate( xyz_CloudCoverforRad )
    end if ! FlagFullPhysics

    !* For coupler ----------------------------------------
    deallocate( xy_TauXAtm, xy_TauYAtm )
    deallocate( xy_SensAtm, xy_LatentAtm )
    deallocate( xy_LDWRFlxAtm, xy_LUWRFlxAtm )
    deallocate( xy_SDWRFlxAtm, xy_SUWRFlxAtm )
    deallocate( xy_SurfAirTemp )    
    deallocate( xy_DSurfHFlxDTs, xy_DSurfLatentFlxDTs )
    deallocate( xy_RainAtm, xy_SnowAtm )
    !-------------------------------------------------------
    
    select case ( IDDynMode )
    case ( IDDynModeHSPLVAS83 )
      ! 各モジュール内の変数の割付解除
      ! Dellocation of variables in modules
      !
      call DynamicsHSplVAS83Finalize
    case ( IDDynModeNoHorAdv )
    case ( IDDynModeTWPICE )
      !
      ! Dynamical process for TWP-ICE experiment
      !
      call DynamicsTWPICESCMExpFinalize
    end select


    call AuxVarsFinalize

    select case ( IDPhysMode )
    case ( IDPhysModeHS94 )
      call Hs94Finalize
    end select

    select case ( IDPhysMode )
    case ( IDPhysModeFullPhysics )
      ! 割付解除とリスタートファイルの終了処理
      ! Dellocation and close a restart file
      !
      call RadDennouAGCMFinalize
    end select

    call AxessetFinalize

    ! 時刻管理終了処理
    ! Termination of time control
    !
    call TimesetClose

    call ProfUtil_RapEnd('Shutdown', 0)

    ProfileName = CPrintf("dcpam_rank%06d.prof", i=(/myrank/))
    call ProfUtil_RapReport(ProfileName)
    call ProfUtil_Final()
    
    ! Finalize MPI
    !
    call MPIWrapperFinalize


  end subroutine MainTerminate

end module dcpam_main_mod
