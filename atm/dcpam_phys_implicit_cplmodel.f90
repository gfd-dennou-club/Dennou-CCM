!= 陰解法による時間積分 (大気のみ / 惑星表面温度・土壌温度計算なし)
!
!= Time integration by using implicit scheme in case without calculation of surface and soil temperature
!
! Authors::   Yoshiyuki O. Takahashi
! Version::   $Id: phy_implicit_cplmodel.f90,v 1.5 2015/01/29 12:05:01 yot Exp $
! Tag Name::  $Name:  $
! Copyright:: Copyright (C) GFD Dennou Club, 2008-2009. All rights reserved.
! License::   See COPYRIGHT[link:../../../COPYRIGHT]
!

module phy_implicit_cplmodel
  !
  != 陰解法による時間積分 (大気のみ / 惑星表面温度・土壌温度計算なし)
  !
  != Time integration by using implicit scheme in case without calculation of surface and soil temperature
  !
  ! <b>Note that Japanese and English are described in parallel.</b>
  !
  !== Procedures List
  !
  ! PhyImplTendency      :: 時間変化率の計算
  ! ------------         :: ------------
  ! PhyImplTendency      :: Calculate tendency
  !
  !--
  !== NAMELIST
  !
  ! NAMELIST#phy_implicit_nml
  !++

  ! モジュール引用 ; USE statements
  !

  ! 格子点設定
  ! Grid points settings
  !
  use gridset, only:   imax, & ! 経度格子点数. 
                               ! Number of grid points in longitude
    &                  jmax, & ! 緯度格子点数. 
                               ! Number of grid points in latitude
    &                  kmax, & ! 鉛直層数. 
                               ! Number of vertical level
    &                  kslmax  ! 地下の鉛直層数. 
                               ! Number of subsurface vertical level

  ! 組成に関わる配列の設定
  ! Settings of array for atmospheric composition
  !
  use composition, only: ncmax, IndexH2OVap

  ! 種別型パラメタ
  ! Kind type parameter
  !
  use dc_types, only: DP, &      ! 倍精度実数型. Double precision. 
    &                 STRING     ! 文字列.       Strings. 

  ! メッセージ出力
  ! Message output
  !
  use dc_message, only: MessageNotify

  ! 宣言文 ; Declaration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: PhyImplCplmodelTendency_Forward
  public :: PhyImplCplmodelTendency_Backward
  public :: PhyImplCplmodelInit


  ! 公開変数
  ! Public variables
  !

  ! 非公開変数
  ! Private variables
  !
  logical, save :: FlagPresSurfTemp
  logical, save :: FlagPresSurfQMix


  logical, save :: phy_implicit_cplmodel_inited = .false.
                              ! 初期設定フラグ. 
                              ! Initialization flag

  character(*), parameter:: module_name = 'phy_implicit_cplmodel'
                              ! モジュールの名称. 
                              ! Module name
  character(*), parameter:: version = &
    & '$Name:  $' // &
    & '$Id: phy_implicit_cplmodel.f90,v 1.5 2015/01/29 12:05:01 yot Exp $'
                              ! モジュールのバージョン
                              ! Module version

  real(DP), allocatable :: xyza_UVLUMtx    (:,:,:,:)
  real(DP), allocatable :: xyza_TempLUMtx    (:,:,:,:)
  real(DP), allocatable :: xyza_QMixLUMtx    (:,:,:,:)
 
contains

  subroutine PhyImplCplmodelTendency_Backward(                 &
    & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt            & ! (out)
    & )

    ! 時刻管理
    ! Time control
    !
    use timeset, only: &
         & DelTime
    
    real(DP), intent(inout):: xyz_DUDt (0:imax-1, 1:jmax, 1:kmax)
                              ! $ \DP{u}{t} $ . 東西風速変化. 
                              ! Eastward wind tendency
    real(DP), intent(inout):: xyz_DVDt (0:imax-1, 1:jmax, 1:kmax)
                              ! $ \DP{v}{t} $ . 南北風速変化. 
                              ! Northward wind tendency
    real(DP), intent(inout):: xyz_DTempDt (0:imax-1, 1:jmax, 1:kmax)
                              ! $ \DP{T}{t} $ . 温度変化. 
                              ! Temperature tendency
    real(DP), intent(inout):: xyzf_DQMixDt(0:imax-1, 1:jmax, 1:kmax, 1:ncmax)
                              ! $ \DP{q}{t} $ . 質量混合比変化. 
                              ! Mass mixing ratio tendency

    ! モジュール引用 ; USE statements
    !

    ! 作業変数
    ! Work variables
    !    
    integer:: k               ! 鉛直方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in vertical direction
!!$    integer:: l               ! 行列用 DO ループ用作業変数
!!$                              ! Work variables for DO loop of matrices
    integer:: n               ! 組成方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in dimension of constituents

    ! 実行文 ; Executable statement
    !

    call PhyImplLUSolve3_Backward( &
         & xyz_DUDt, &            ! (inout)
         & xyza_UVLUMtx, &        ! (in)
         & 1, imax * jmax, kmax ) ! (in)
    call kreverse_and_tendency(xyz_DUDt)

    call PhyImplLUSolve3_Backward( &
         & xyz_DVDt, &            ! (inout)
         & xyza_UVLUMtx, &        ! (in)
         & 1, imax * jmax, kmax ) ! (in)
    call kreverse_and_tendency(xyz_DVDt)

    call PhyImplLUSolve3_Backward(           &
         & xyz_DTempDt,                & ! (inout)
         & xyza_TempLUMtx,             & ! (in)
         & 1, imax * jmax , kmax       & ! (in)
         & )
    call kreverse_and_tendency(xyz_DTempDt)

    do n = 1, ncmax
      call PhyImplLUSolve3_Backward(      &
        & xyzf_DQMixDt(:,:,:,n),          & ! (inout)
        & xyza_QMixLUMtx,                 & ! (in)
        & 1, imax * jmax , kmax           & ! (in)
        & )
      call kreverse_and_tendency(xyzf_DQMixDt)
    end do

  contains
    subroutine kreverse_and_tendency(xyz)
      real(DP), intent(inout) :: xyz(0:imax-1,jmax,1:kmax)

      real(DP) :: xyz_(0:imax-1,jmax,1:kmax)
      integer :: k
      
      xyz_(:,:,:) = xyz
      do k=1, kmax
         xyz(:,:,kmax-k+1) = xyz_(:,:,k)/( 2.0_DP * DelTime)
      end do
    end subroutine kreverse_and_tendency
    
  end subroutine PhyImplCplmodelTendency_Backward
  
  subroutine PhyImplCplmodelTendency_Forward(                           &
    & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
    & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
    & xyr_VirTemp, xyz_Height,                                 & ! (in)
    & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
    & xy_SurfVelTransCoef, xy_SurfTempTransCoef,               & ! (in)
    & xy_SurfQVapTransCoef,                                    & ! (in)
    & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt            & ! (out)
    & )
    !
    ! 時間変化率の計算を行います. 
    !
    ! Calculate tendencies. 
    !

    ! モジュール引用 ; USE statements
    !

    ! 物理定数設定
    ! Physical constants settings
    !
    use constants, only: &
      & Grav, &               ! $ g $ [m s-2]. 
                              ! 重力加速度. 
                              ! Gravitational acceleration
      & CpDry, &
                              ! $ C_p $ [J kg-1 K-1]. 
                              ! 乾燥大気の定圧比熱. 
                              ! Specific heat of air at constant pressure
      & LatentHeat, &
                              ! $ L $ [J kg-1] . 
                              ! 凝結の潜熱. 
                              ! Latent heat of condensation
      & GasRDry
                              ! $ R $ [J kg-1 K-1]. 
                              ! 乾燥大気の気体定数. 
                              ! Gas constant of air

    ! 時刻管理
    ! Time control
    !
    use timeset, only: &
      & DelTime, &            ! $ \Delta t $ [s]
      & TimeN, &              ! ステップ $ t $ の時刻. Time of step $ t $. 
      & TimesetClockStart, TimesetClockStop

    ! 陰解法による時間積分のためのルーチン
    ! Routines for time integration with implicit scheme
    !
    use phy_implicit_utils, only : PhyImplLUDecomp3, PhyImplLUSolve3

    ! 宣言文 ; Declaration statements
    !
    implicit none

    real(DP), intent(in):: xyr_MomFluxX (0:imax-1, 1:jmax, 0:kmax)
                              ! 東西方向運動量フラックス. 
                              ! Eastward momentum flux
    real(DP), intent(in):: xyr_MomFluxY (0:imax-1, 1:jmax, 0:kmax)
                              ! 南北方向運動量フラックス. 
                              ! Northward momentum flux
    real(DP), intent(in):: xyr_HeatFlux (0:imax-1, 1:jmax, 0:kmax)
                              ! 熱フラックス. 
                              ! Heat flux
    real(DP), intent(in):: xyrf_QMixFlux(0:imax-1, 1:jmax, 0:kmax, 1:ncmax)
                              ! 比湿フラックス. 
                              ! Specific humidity flux

    real(DP), intent(in):: xyr_Press (0:imax-1, 1:jmax, 0:kmax)
                              ! $ \hat{p} $ . 気圧 (半整数レベル). 
                              ! Air pressure (half level)
    real(DP), intent(in):: xyz_Exner (0:imax-1, 1:jmax, 1:kmax)
                              ! Exner 関数 (整数レベル). 
                              ! Exner function (full level)
    real(DP), intent(in):: xyr_Exner (0:imax-1, 1:jmax, 0:kmax)
                              ! Exner 関数 (半整数レベル). 
                              ! Exner function (half level)

    real(DP), intent(in):: xyr_VirTemp (0:imax-1, 1:jmax, 0:kmax)
                              ! $ \hat{T}_v $ . 仮温度 (半整数レベル). 
                              ! Virtual temperature (half level)
    real(DP), intent(in):: xyz_Height (0:imax-1, 1:jmax, 1:kmax)
                              ! 高度 (整数レベル). 
                              ! Height (full level)

    real(DP), intent(in):: xyr_VelDiffCoef (0:imax-1, 1:jmax, 0:kmax)
                              ! 拡散係数：運動量. 
                              ! Diffusion coefficient: velocity
    real(DP), intent(in):: xyr_TempDiffCoef (0:imax-1, 1:jmax, 0:kmax)
                              ! 拡散係数：温度. 
                              ! Transfer coefficient: temperature
    real(DP), intent(in):: xyr_QMixDiffCoef (0:imax-1, 1:jmax, 0:kmax)
                              ! 拡散係数：比湿. 
                              ! Diffusion coefficient: specific humidity

    real(DP), intent(in):: xy_SurfVelTransCoef (0:imax-1, 1:jmax)
                              ! 輸送係数：運動量. 
                              ! Diffusion coefficient: velocity
    real(DP), intent(in):: xy_SurfTempTransCoef (0:imax-1, 1:jmax)
                              ! 輸送係数：温度. 
                              ! Transfer coefficient: temperature
    real(DP), intent(in):: xy_SurfQVapTransCoef (0:imax-1, 1:jmax)
                              ! 輸送係数：比湿. 
                              ! Transfer coefficient: specific humidity

    real(DP), intent(out):: xyz_DUDt (0:imax-1, 1:jmax, 1:kmax)
                              ! $ \DP{u}{t} $ . 東西風速変化. 
                              ! Eastward wind tendency
    real(DP), intent(out):: xyz_DVDt (0:imax-1, 1:jmax, 1:kmax)
                              ! $ \DP{v}{t} $ . 南北風速変化. 
                              ! Northward wind tendency
    real(DP), intent(out):: xyz_DTempDt (0:imax-1, 1:jmax, 1:kmax)
                              ! $ \DP{T}{t} $ . 温度変化. 
                              ! Temperature tendency
    real(DP), intent(out):: xyzf_DQMixDt(0:imax-1, 1:jmax, 1:kmax, 1:ncmax)
                              ! $ \DP{q}{t} $ . 質量混合比変化. 
                              ! Mass mixing ratio tendency

    ! 作業変数
    ! Work variables
    !

    real(DP) :: xyr_VelTransCoef (0:imax-1, 1:jmax, 0:kmax)
                              ! 輸送係数：運動量. 
                              ! Transfer coefficient: velocity
    real(DP) :: xyr_TempTransCoef (0:imax-1, 1:jmax, 0:kmax)
                              ! 輸送係数：温度. 
                              ! Transfer coefficient: temperature
    real(DP) :: xyr_QMixTransCoef(0:imax-1, 1:jmax, 0:kmax)
                              ! 輸送係数：質量. 
                              ! Transfer coefficient: mass of constituents

    real(DP):: xyza_UVMtx  (0:imax-1, 1:jmax, 1:kmax, -1:1)
                              ! 速度陰解行列. 
                              ! Implicit matrix about velocity 
    real(DP):: xyza_TempMtx(0:imax-1, 1:jmax, 1:kmax, -1:1)
                              ! 温度陰解行列. 
                              ! Implicit matrix about temperature
    real(DP):: xyza_QMixMtx(0:imax-1, 1:jmax, 1:kmax, -1:1)
                              ! 質量混合比陰解行列. 
                              ! Implicit matrix about mass mixing ratio

    real(DP):: xyz_DelQMixLUVec(0:imax-1, 1:jmax, 1:kmax)
                              ! $ q $ の時間変化.
                              ! Tendency of $ q $

                              ! LU 行列.
                              ! LU matrix
    real(DP):: xya_DelTempLUVec(0:imax-1, 1:jmax, 1:kmax)
                              ! $ T, Tg $ の時間変化.
                              ! Tendency of $ T $ and $ Tg |


!!$    integer:: i               ! 経度方向に回る DO ループ用作業変数
!!$                              ! Work variables for DO loop in longitude
!!$    integer:: j               ! 緯度方向に回る DO ループ用作業変数
!!$                              ! Work variables for DO loop in latitude
    integer:: k               ! 鉛直方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in vertical direction
!!$    integer:: l               ! 行列用 DO ループ用作業変数
!!$                              ! Work variables for DO loop of matrices
    integer:: n               ! 組成方向に回る DO ループ用作業変数
                              ! Work variables for DO loop in dimension of constituents

    ! 実行文 ; Executable statement
    !

    ! 初期化確認
    ! Initialization check
    !
    if ( .not. phy_implicit_cplmodel_inited ) then
      call MessageNotify( 'E', module_name, 'This module has not been initialized.' )
    end if


    ! 計算時間計測開始
    ! Start measurement of computation time
    !
    call TimesetClockStart( module_name )


    ! 輸送係数の計算
    ! Calculate transfer coefficient
    !
    xyr_VelTransCoef (:,:,0)    = 0.0_DP
    xyr_VelTransCoef (:,:,kmax) = 0.0_DP
    xyr_TempTransCoef(:,:,0)    = 0.0_DP
    xyr_TempTransCoef(:,:,kmax) = 0.0_DP
    xyr_QMixTransCoef(:,:,0)    = 0.0_DP
    xyr_QMixTransCoef(:,:,kmax) = 0.0_DP

    do k = 1, kmax-1
      xyr_VelTransCoef(:,:,k) =                                     &
        &   xyr_VelDiffCoef(:,:,k)                                  &
        &     * xyr_Press(:,:,k) / ( GasRDry * xyr_VirTemp(:,:,k) ) &
        &     / ( xyz_Height(:,:,k+1) - xyz_Height(:,:,k) )

      xyr_TempTransCoef(:,:,k) =                                    &
        &   xyr_TempDiffCoef(:,:,k)                                 &
        &     * xyr_Press(:,:,k) / ( GasRDry * xyr_VirTemp(:,:,k) ) &
        &     / ( xyz_Height(:,:,k+1) - xyz_Height(:,:,k) )

      xyr_QMixTransCoef(:,:,k) =                                    &
        &   xyr_QMixDiffCoef(:,:,k)                                 &
        &     * xyr_Press(:,:,k) / ( GasRDry * xyr_VirTemp(:,:,k) ) &
        &     / ( xyz_Height(:,:,k+1) - xyz_Height(:,:,k) )
    end do


    ! 陰解法のための行列作成
    ! Create matrices for implicit scheme
    !

    ! 鉛直拡散スキームの輸送係数から陰解行列の計算 (速度)
    ! Calculate implicit matrices from transfer coefficient of vertical diffusion scheme (velocity)
    !
    k = 1
    xyza_UVMtx  (:,:,k,-1) = 0.0_DP
    xyza_UVMtx  (:,:,k, 0) =                                                  &
      & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
      & + xy_SurfVelTransCoef(:,:)                                            &
      & + xyr_VelTransCoef(:,:,k  )
    xyza_UVMtx  (:,:,k, 1) = &
      & - xyr_VelTransCoef(:,:,k)

    do k = 2, kmax-1
      xyza_UVMtx  (:,:,k,-1) = &
        & - xyr_VelTransCoef(:,:,k-1)
      xyza_UVMtx  (:,:,k, 0) =                                                  &
        & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + xyr_VelTransCoef(:,:,k-1)                                           &
        & + xyr_VelTransCoef(:,:,k  )
      xyza_UVMtx  (:,:,k, 1) = &
        & - xyr_VelTransCoef(:,:,k)
    end do

    k = kmax
    xyza_UVMtx  (:,:,k,-1) = &
      & - xyr_VelTransCoef(:,:,k-1)
    xyza_UVMtx  (:,:,k, 0) =                                                  &
      & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
      & + xyr_VelTransCoef(:,:,k-1)
    xyza_UVMtx  (:,:,k, 1) = 0.0_DP


    ! 鉛直拡散スキームの輸送係数から陰解行列の計算 (温度)
    ! Calculate implicit matrices from transfer coefficient of vertical diffusion scheme (temperature)
    !
    k = 1
    xyza_TempMtx(:,:,k,-1) = 0.0_DP
    if ( FlagPresSurfTemp ) then
      ! Prescribe surface temperature
      xyza_TempMtx(:,:,k, 0) =                                                          &
        & - CpDry * ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + CpDry * xyr_Exner(:,:,k-1) / xyz_Exner(:,:,k  ) * xy_SurfTempTransCoef(:,:) &
        & + CpDry * xyr_Exner(:,:,k  ) / xyz_Exner(:,:,k  ) * xyr_TempTransCoef(:,:,k  )
    else
      ! Prescribe surface flux
      xyza_TempMtx(:,:,k, 0) =                                                          &
        & - CpDry * ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + CpDry * xyr_Exner(:,:,k  ) / xyz_Exner(:,:,k  ) * xyr_TempTransCoef(:,:,k  )
    end if
    xyza_TempMtx(:,:,k, 1) = &
      & - CpDry * xyr_Exner(:,:,k  ) / xyz_Exner(:,:,k+1) * xyr_TempTransCoef(:,:,k  )

    do k = 2, kmax-1
      xyza_TempMtx(:,:,k,-1) = &
        & - CpDry * xyr_Exner(:,:,k-1) / xyz_Exner(:,:,k-1) * xyr_TempTransCoef(:,:,k-1)
      xyza_TempMtx(:,:,k, 0) =                                                          &
        & - CpDry * ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + CpDry * xyr_Exner(:,:,k-1) / xyz_Exner(:,:,k  ) * xyr_TempTransCoef(:,:,k-1)&
        & + CpDry * xyr_Exner(:,:,k  ) / xyz_Exner(:,:,k  ) * xyr_TempTransCoef(:,:,k  )
      xyza_TempMtx(:,:,k, 1) = &
        & - CpDry * xyr_Exner(:,:,k  ) / xyz_Exner(:,:,k+1) * xyr_TempTransCoef(:,:,k  )
    end do

    k = kmax
    xyza_TempMtx(:,:,k,-1) = &
      & - CpDry * xyr_Exner(:,:,k-1) / xyz_Exner(:,:,k-1) * xyr_TempTransCoef(:,:,k-1)
    xyza_TempMtx(:,:,k, 0) =                                                          &
      & - CpDry * ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
      & + CpDry * xyr_Exner(:,:,k-1) / xyz_Exner(:,:,k  ) * xyr_TempTransCoef(:,:,k-1)
    xyza_TempMtx(:,:,k, 1) = 0.0_DP



    ! 鉛直拡散スキームの輸送係数から陰解行列の計算 (比湿)
    ! Calculate implicit matrices from transfer coefficient of vertical diffusion scheme (specific humidity)
    !

    k = 1
    xyza_QMixMtx(:,:,k,-1) =                                                  &
      & 0.0_DP
    if ( FlagPresSurfQMix ) then
      ! Prescribe surface mixing ratio
      xyza_QMixMtx(:,:,k, 0) =                                                  &
        & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + xy_SurfQVapTransCoef(:,:)                                           &
        & + xyr_QMixTransCoef(:,:,k  )
    else
      ! Prescribe surface flux
      xyza_QMixMtx(:,:,k, 0) =                                                  &
        & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + xyr_QMixTransCoef(:,:,k  )
    end if
    xyza_QMixMtx(:,:,k, 1) =                                                  &
      & - xyr_QMixTransCoef(:,:,k  )

    do k = 2, kmax-1
      xyza_QMixMtx(:,:,k,-1) =                                                  &
        & - xyr_QMixTransCoef(:,:,k-1)
      xyza_QMixMtx(:,:,k, 0) =                                                  &
        & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
        & + xyr_QMixTransCoef(:,:,k-1)                                          &
        & + xyr_QMixTransCoef(:,:,k  )
      xyza_QMixMtx(:,:,k, 1) =                                                  &
        & - xyr_QMixTransCoef(:,:,k  )
    end do

    k = kmax
    xyza_QMixMtx(:,:,k,-1) =                                                  &
      & - xyr_QMixTransCoef(:,:,k-1)
    xyza_QMixMtx(:,:,k, 0) =                                                  &
      & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
      & + xyr_QMixTransCoef(:,:,k-1)
    xyza_QMixMtx(:,:,k, 1) = 0.0_DP




    ! 東西風速, 南北風速の計算
    ! Calculate eastward and northward wind
    !
    do k=1, kmax
       xyza_UVLUMtx(:,:,kmax-k+1,:) = xyza_UVMtx(:,:,k,:)
    end do
    call PhyImplLUDecomp3( &
      & xyza_UVLUMtx, &     ! (inout)
      & imax * jmax, kmax ) ! (in)

    do k = 1, kmax
      xyz_DUDt(:,:,kmax-k+1) = - ( xyr_MomFluxX(:,:,k) - xyr_MomFluxX(:,:,k-1) )
      xyz_DVDt(:,:,kmax-k+1) = - ( xyr_MomFluxY(:,:,k) - xyr_MomFluxY(:,:,k-1) )
    end do

    call PhyImplLUSolve3_Forward( &
      & xyz_DUDt, &            ! (inout)
      & xyza_UVLUMtx, &        ! (in)
      & 1, imax * jmax, kmax ) ! (in)

    call PhyImplLUSolve3_Forward( &
      & xyz_DVDt, &            ! (inout)
      & xyza_UVLUMtx, &        ! (in)
      & 1, imax * jmax, kmax ) ! (in)

!!$    do k = 1, kmax
!!$      xyz_DUDt(:,:,k) = xyz_DUDt(:,:,k) / ( 2.0_DP * DelTime )
!!$      xyz_DVDt(:,:,k) = xyz_DVDt(:,:,k) / ( 2.0_DP * DelTime )
!!$    end do


    ! 温度の計算
    ! Calculate temperature
    !
    xyza_TempLUMtx = xyza_TempMtx

    call PhyImplLUDecomp3( &
      & xyza_TempLUMtx,    & ! (inout)
      & imax * jmax, kmax  & ! (in)
      )

    do k = 1, kmax
      xyz_DTempDt(:,:,kmax-k+1) = &
        - ( xyr_HeatFlux(:,:,k) - xyr_HeatFlux(:,:,k-1) )
    end do

    call PhyImplLUSolve3_Forward(           &
         & xyz_DTempDt,                & ! (inout)
         & xyza_TempLUMtx,             & ! (in)
         & 1, imax * jmax , kmax       & ! (in)
         & )

!!$    xyz_DTempDt = xya_DelTempLUVec / ( 2.0_DP * DelTime )


    ! 比湿の計算
    ! Calculate specific humidity
    !
    xyza_QMixLUMtx = xyza_QMixMtx

    call PhyImplLUDecomp3( &
      & xyza_QMixLUMtx,    &   ! (inout)
      & imax * jmax, kmax  &   ! (in)
      & )

    do n = 1, ncmax
      do k = 1, kmax
         xyzf_DQMixDt(:,:,kmax-k+1,n) = &
              & - ( xyrf_QMixFlux(:,:,k,n) - xyrf_QMixFlux(:,:,k-1,n) )
      end do

      call PhyImplLUSolve3_Forward(      &
        & xyzf_DQMixDt(:,:,:,n),         & ! (inout)
        & xyza_QMixLUMtx,        & ! (in)
        & 1, imax * jmax , kmax  & ! (in)
        & )

!!$      xyzf_DQMixDt(:,:,:,n) = xyz_DelQMixLUVec(:,:,:) / ( 2.0_DP * DelTime )
    end do


    ! 計算時間計測一時停止
    ! Pause measurement of computation time
    !
    call TimesetClockStop( module_name )

  end subroutine PhyImplCplmodelTendency_Forward

  !-------------------------------------------------------------------

  subroutine PhyImplCplmodelInit
    !
    ! phy_implicit_cplmodel モジュールの初期化を行います. 
    ! NAMELIST#phy_implicit_cplmodel_nml の読み込みはこの手続きで行われます. 
    !
    ! "phy_implicit_cplmodel" module is initialized. 
    ! "NAMELIST#phy_implicit_cplmodel_nml" is loaded in this procedure. 
    !

    ! モジュール引用 ; USE statements
    !

    ! NAMELIST ファイル入力に関するユーティリティ
    ! Utilities for NAMELIST file input
    !
    use namelist_util, only: namelist_filename, NmlutilMsg, NmlutilAryValid

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    ! 文字列操作
    ! Character handling
    !
    use dc_string, only: StoA

    ! 宣言文 ; Declaration statements
    !
    implicit none

    ! 作業変数
    ! Work variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
                              ! Unit number for NAMELIST file open
    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
                              ! IOSTAT of NAMELIST read

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /phy_implicit_cplmodel_nml/ &
      & FlagPresSurfTemp, FlagPresSurfQMix
          !
          ! デフォルト値については初期化手続 "phy_implicit_cplmodel#PhyImplInit" 
          ! のソースコードを参照のこと. 
          !
          ! Refer to source codes in the initialization procedure
          ! "phy_implicit_cplmodel#PhyImplInit" for the default values. 
          !

    ! 実行文 ; Executable statement
    !

    if ( phy_implicit_cplmodel_inited ) return

    ! デフォルト値の設定
    ! Default values settings
    !
    FlagPresSurfTemp = .false.
    FlagPresSurfQMix = .false.


    ! NAMELIST の読み込み
    ! NAMELIST is input
    !
    if ( trim(namelist_filename) /= '' ) then
      call FileOpen( unit_nml, &          ! (out)
        & namelist_filename, mode = 'r' ) ! (in)

      rewind( unit_nml )
      read( unit_nml,                     &  ! (in)
        & nml = phy_implicit_cplmodel_nml, &  ! (out)
        & iostat = iostat_nml )              ! (out)
      close( unit_nml )

      call NmlutilMsg( iostat_nml, module_name ) ! (in)
    end if

    allocate( xyza_UVLUMtx(0:imax-1,jmax,kmax,-1:1) )
    allocate( xyza_TempLUMtx(0:imax-1,jmax,kmax,-1:1) )
    allocate( xyza_QMixLUMtx(0:imax-1,jmax,kmax,-1:1) )
    
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  FlagPresSurfTemp = %b', l = (/ FlagPresSurfTemp /) )
    call MessageNotify( 'M', module_name, '  FlagPresSurfQMix = %b', l = (/ FlagPresSurfQMix /) )
    call MessageNotify( 'M', module_name, '-- version = %c', c1 = trim(version) )

    phy_implicit_cplmodel_inited = .true.

  end subroutine PhyImplCplmodelInit

  !-------------------------------------------------------------------

  subroutine PhyImplLUSolve3_Forward( &
    & ijn_Vector, &       ! (inout)
    & jna_LUMtx, &        ! (in)
    & IDim, JDim, NDim &  ! (in)
    & )
    !
    ! LU 分解による解の計算 (3重対角行列用) を行います.
    !
    ! Solve with LU decomposition (For triple diagonal matrix). 
    !

    ! 宣言文 ; Declaration statements
    !
    implicit none
    integer, intent(in):: IDim
    integer, intent(in):: JDim
    integer, intent(in):: NDim
    real(DP), intent(in):: jna_LUMtx(JDim, NDim, -1:1)
                              ! LU 行列. 
                              ! LU matrix
    real(DP), intent(inout):: ijn_Vector(IDim, JDim, NDim)
                              ! 右辺ベクトル / 解. 
                              ! Right-hand side vector / solution

    ! 作業変数
    ! Work variables
    ! 
    integer:: i, j, n         ! DO ループ用作業変数
                              ! Work variables for DO loop

    ! 実行文 ; Executable statement
    !

    ! 前進代入
    ! Forward substitution
    !
    do i = 1, IDim
      do j = 1, JDim
        ijn_Vector(i,j,1) = ijn_Vector(i,j,1) / jna_LUMtx(j,1,0)
      end do
    end do

    do n = 2, NDim
      do i = 1, IDim
        do j = 1, JDim
          ijn_Vector(i,j,n) = (   ijn_Vector(i,j,n) &
            &                   - ijn_Vector(i,j,n-1) * jna_LUMtx(j,n,-1) &
            &                  ) / jna_LUMtx(j,n,0)
        end do
      end do
    end do

  end subroutine PhyImplLUSolve3_Forward

  subroutine PhyImplLUSolve3_Backward( &
    & ijn_Vector, &       ! (inout)
    & jna_LUMtx, &        ! (in)
    & IDim, JDim, NDim &  ! (in)
    & )
    !
    ! LU 分解による解の計算 (3重対角行列用) を行います.
    !
    ! Solve with LU decomposition (For triple diagonal matrix). 
    !

    ! 宣言文 ; Declaration statements
    !
    implicit none
    integer, intent(in):: IDim
    integer, intent(in):: JDim
    integer, intent(in):: NDim
    real(DP), intent(in):: jna_LUMtx(JDim, NDim, -1:1)
                              ! LU 行列. 
                              ! LU matrix
    real(DP), intent(inout):: ijn_Vector(IDim, JDim, NDim)
                              ! 右辺ベクトル / 解. 
                              ! Right-hand side vector / solution

    ! 作業変数
    ! Work variables
    ! 
    integer:: i, j, n         ! DO ループ用作業変数
                              ! Work variables for DO loop
    
    ! 実行文 ; Executable statement
    !

    ! 後退代入
    ! Backward substitution
    !
    do n = NDim-1, 1, -1
      do i = 1, IDim
        do j = 1, JDim
          ijn_Vector(i,j,n) =   ijn_Vector(i,j,n) &
            &                 - ijn_Vector(i,j,n+1) * jna_LUMtx(j,n,1)
        end do
     end do
    end do
    
  end subroutine PhyImplLUSolve3_Backward
  
end module phy_implicit_cplmodel
