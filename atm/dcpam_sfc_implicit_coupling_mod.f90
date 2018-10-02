module dcpam_sfc_implicit_coupling_mod

  use gridset, only: &
       & imax, jmax, kmax
  use composition, only: ncmax, IndexH2Ovap
  use dc_types, only: DP, STRING
  use dc_message, only: MessageNotify

  implicit none
  private

  public :: dcpam_sfc_implicit_coupling_Init
  public :: SfcImplicitCoupling_VDiffBackward
  public :: SfcImplicitCoupling_VDiffForward
  
  real(DP), save, allocatable :: xyza_UVMtx(:,:,:,:)
  real(DP), save, allocatable :: xyza_TempMtx(:,:,:,:)
  real(DP), save, allocatable :: xyza_QMixMtx(:,:,:,:)

  character(*), parameter :: module_name = 'dcpam_sfc_impl_coupling_mod'
  
contains


  subroutine SfcImplicitCoupling_VDiffBackward( &
    & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt            & ! (out)
    & )

    use timeset, only: &
      & DelTime    
    use mpi_wrapper
    
    implicit none
    
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

    integer :: k
    integer :: n
    
    call Solve_TriDiagSystem_Backward( xyza_UVMtx, xyz_DUDt )
    call Solve_TriDiagSystem_Backward( xyza_UVMtx, xyz_DVDt )
    call Solve_TriDiagSystem_Backward( xyza_TempMtx, xyz_DTempDt )
    do k=1, kmax
       xyz_DUDt(:,:,k) = xyz_DUDt(:,:,k) / (2d0 * DelTime)
       xyz_DVDt(:,:,k) = xyz_DVDt(:,:,k) / (2d0 * DelTime)
       xyz_DTempDt(:,:,k) = xyz_DTempDt(:,:,k) / (2d0 * DelTime)       
    end do

    do n=1, ncmax
       call Solve_TriDiagSystem_Backward( xyza_QMixMtx, xyzf_DQMixDt(:,:,:,n) )
       xyzf_DQMixDt(:,:,:,n) = xyzf_DQMixDt(:,:,:,n) / (2d0 * DelTime)
    end do

!!$    write(*,*) myrank, ":U", sum(xyz_DUDt(0,1,:)*2d0*DelTime)
!!$    write(*,*) myrank, ":Temp", sum(xyz_DTempDt(0,1,:)*2d0*DelTime)
!!$    write(*,*) myrank, ":Vap", sum(xyzf_DQMixDt(0,1,:,IndexH2OVap)*2d0*DelTime)
!!$    stop
    
  end subroutine SfcImplicitCoupling_VDiffBackward
  
  subroutine SfcImplicitCoupling_VDiffForward( &
    & xyr_MomFluxX, xyr_MomFluxY, xyr_HeatFlux, xyrf_QMixFlux, & ! (in)
    & xyr_Press, xyz_Exner, xyr_Exner,                         & ! (in)
    & xyr_VirTemp, xyz_Height,                                 & ! (in)
    & xyr_VelDiffCoef, xyr_TempDiffCoef, xyr_QMixDiffCoef,     & ! (in)
    & xyz_DUDt, xyz_DVDt, xyz_DTempDt, xyzf_DQMixDt,           & ! (out)
    & xya_ImplCplCoef1, xya_ImplCplCoef2                                         & ! (out)
    & )
   
    
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

    use mpi_wrapper, only: myrank
    
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
    real(DP), intent(out):: xya_ImplCplCoef1(0:imax-1,1:jmax,4)
    real(DP), intent(out):: xya_ImplCplCoef2(0:imax-1,1:jmax,4)
    
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    
    real(DP) :: xyr_VelTransCoef(0:imax-1,jmax,0:kmax)
    real(DP) :: xyr_TempTransCoef(0:imax-1,jmax,0:kmax)
    real(DP) :: xyr_QMixTransCoef(0:imax-1,jmax,0:kmax)
    real(DP) :: xy_Tmp(0:imax-1,jmax)

    real(DP) :: DFADUV1    
    real(DP) :: DFADUV2
    real(DP) :: DFADT1    
    real(DP) :: DFADT2
    real(DP) :: DFADQ1    
    real(DP) :: DFADQ2
    real(DP) :: xyza_QMixMtxSave(0:imax-1,jmax,1:kmax,-1:1)
    real(DP) :: xyza_UVMtxSave(0:imax-1,jmax,1:kmax,-1:1)
    
!!$    call MessageNotify('M', module_name, "Set Mtx ..")        
    xyr_VelTransCoef (:,:,0)    = 0.0_DP
    xyr_VelTransCoef (:,:,kmax) = 0.0_DP
    xyr_TempTransCoef(:,:,0)    = 0.0_DP
    xyr_TempTransCoef(:,:,kmax) = 0.0_DP
    xyr_QMixTransCoef(:,:,0)    = 0.0_DP
    xyr_QMixTransCoef(:,:,kmax) = 0.0_DP

    do k = 1, kmax-1
       xy_Tmp(:,:) = xyr_Press(:,:,k) / (GasRDry * xyr_VirTemp(:,:,k) ) &
            & / ( xyz_Height(:,:,k+1) - xyz_Height(:,:,k) )

       xyr_VelTransCoef(:,:,k) = xyr_VelDiffCoef(:,:,k)*xy_Tmp(:,:)
       xyr_TempTransCoef(:,:,k) = xyr_TempDiffCoef(:,:,k)*xy_Tmp(:,:)
       xyr_QMixTransCoef(:,:,k) = xyr_QMixDiffCoef(:,:,k)*xy_Tmp(:,:)
    end do

    !------------
    
    k = 1
    xyza_UVMtx  (:,:,k,-1) = 0.0_DP
    xyza_UVMtx  (:,:,k, 0) =                                                  &
      & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
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

    !-----------

    k = 1
    xyza_TempMtx(:,:,k,-1) = 0.0_DP
    ! Prescribe surface flux
    xyza_TempMtx(:,:,k, 0) =                                                          &
         & - CpDry * ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
         & + CpDry * xyr_Exner(:,:,k  ) / xyz_Exner(:,:,k  ) * xyr_TempTransCoef(:,:,k  )
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

    !-----------

    k = 1
    xyza_QMixMtx(:,:,k,-1) =                                                  &
      & 0.0_DP
    ! Prescribe surface flux
    xyza_QMixMtx(:,:,k, 0) =                                                  &
         & - ( xyr_Press(:,:,k) - xyr_Press(:,:,k-1) ) / Grav / ( 2.0_DP * DelTime ) &
         & + xyr_QMixTransCoef(:,:,k  )
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

    !--------

    !$omp parallel
    !$omp do
    do k = 1, kmax
       xyz_DUDt(:,:,k) = - (xyr_MomFluxX(:,:,k) - xyr_MomFluxX(:,:,k-1))
       xyz_DVDt(:,:,k) = - (xyr_MomFluxY(:,:,k) - xyr_MomFluxY(:,:,k-1))
       xyz_DTempDt(:,:,k) = - (xyr_HeatFlux(:,:,k) - xyr_HeatFlux(:,:,k-1))
    end do
    !$omp do collapse(2)
    do n = 1, ncmax
    do k = 1, kmax
       xyzf_DQMixDt(:,:,k,n) = &
            & - ( xyrf_QMixFlux(:,:,k,n) - xyrf_QMixFlux(:,:,k-1,n) )
    end do
    end do
    !$omp end parallel

    xya_ImplCplCoef2(:,:,1) = xyz_DUDt(:,:,1)
    xya_ImplCplCoef2(:,:,2) = xyz_DVDt(:,:,1)
    xya_ImplCplCoef2(:,:,3) = xyz_DTempDt(:,:,1)
    xya_ImplCplCoef2(:,:,4) = xyzf_DQMixDt(:,:,1,IndexH2OVap)

!!$    call MessageNotify('M', module_name, "TriDiagsystem forward..")
!!$    if (myrank == 0)then
!!$       do k=1, kmax
!!$          write(*,*) k, ":", xyza_UVMtx(1,1,k,-1:1), ":", xyz_DUDt(1,1,k)
!!$       end do
!!$    end if

    xyza_UVMtxSave = xyza_UVMtx
    call Solve_TriDiagSystem_Forward( xyza_UVMtx, xyz_DUDt )
    call Solve_TriDiagSystem_Forward( xyza_UVMtxSave, xyz_DVDt )
    call Solve_TriDiagSystem_Forward( xyza_TempMtx, xyz_DTempDt )

    xyza_QMixMtxSave = xyza_QMixMtx
    do n=1, ncmax
       xyza_QMixMtx = xyza_QMixMtxSave
       call Solve_TriDiagSystem_Forward( xyza_QMixMtx, xyzf_DQMixDt(:,:,:,n) )
    end do

!!$    if (myrank == 0)then
!!$       do k=1, kmax
!!$          write(*,*) k, ":", xyza_UVMtx(1,1,k,-1:1), ":", xyz_DUDt(1,1,k)
!!$       end do
!!$       stop
!!$    end if
!!$    call MessageNotify('M', module_name, "Set ImplCplCoef..")    

    do j=1, jmax
    do i=0, imax-1
       xy_Tmp(i,j) = - ( xyr_Press(i,j,1) - xyr_Press(i,j,0) ) / Grav / ( 2.0_DP * DelTime )

       !--
       
       DFADUV1 =   xyr_VelTransCoef(i,j,1)
       DFADUV2 = - xyr_VelTransCoef(i,j,1)
       
       xya_ImplCplCoef1(i,j,1) = xy_Tmp(i,j) + DFADUV1 - DFADUV2/xyza_UVMtx(i,j,2,0)
       xya_ImplCplCoef2(i,j,1) = xya_ImplCplCoef2(i,j,1) - DFADUV2*xyz_DUDt(i,j,2)/xyza_UVMtx(i,j,2,0)

       xya_ImplCplCoef1(i,j,2) = xy_Tmp(i,j) + DFADUV1 - DFADUV2/xyza_UVMtx(i,j,2,0)
       xya_ImplCplCoef2(i,j,2) = xya_ImplCplCoef2(i,j,2) - DFADUV2*xyz_DVDt(i,j,2)/xyza_UVMtx(i,j,2,0)
 
       !--
       
       DFADT1 =    CpDry*xyr_Exner(i,j,1)*xyr_TempTransCoef(i,j,1)/xyz_Exner(i,j,1)
       DFADT2 =  - CpDry*xyr_Exner(i,j,1)*xyr_TempTransCoef(i,j,1)/xyz_Exner(i,j,2)

       xya_ImplCplCoef1(i,j,3) = CpDry*xy_Tmp(i,j) + DFADT1 - DFADT2/xyza_TempMtx(i,j,2,0)
       xya_ImplCplCoef2(i,j,3) = xya_ImplCplCoef2(i,j,3) - DFADT2*xyz_DTempDt(i,j,2)/xyza_TempMtx(i,j,2,0)

       !--

       DFADQ1 =   xyr_QMixTransCoef(i,j,1)
       DFADQ2 = - xyr_QMixTransCoef(i,j,1)

       xya_ImplCplCoef1(i,j,4) = xy_Tmp(i,j) + DFADQ1 - DFADQ2/xyza_QMixMtx(i,j,2,0)
       xya_ImplCplCoef2(i,j,4) = xya_ImplCplCoef2(i,j,4) - DFADQ2*xyzf_DQMixDt(i,j,2,IndexH2OVap)/xyza_QMixMtx(i,j,2,0)

    end do
    end do
    
  end subroutine SfcImplicitCoupling_VDiffForward

  subroutine Solve_TriDiagSystem_Forward( &
       & xyza_Mtx, xyz_RHS )

    real(DP), intent(inout) :: xyza_Mtx(0:imax-1,jmax,1:kmax,-1:1)
    real(DP), intent(inout) :: xyz_RHS(0:imax-1,jmax,1:kmax)

    integer :: k

    k = kmax
    xyza_Mtx(:,:,k,0) = xyza_Mtx(:,:,k,0)/xyza_Mtx(:,:,k,-1)
    xyz_RHS(:,:,k) = xyz_RHS(:,:,k)/xyza_Mtx(:,:,k,-1)
    xyza_Mtx(:,:,k,-1) = 1d0

    do k=kmax-1, 1, -1
       xyza_Mtx(:,:,k,0) = (xyza_Mtx(:,:,k,0)*xyza_Mtx(:,:,k+1,0) - xyza_Mtx(:,:,k,1)) &
            &              / (xyza_Mtx(:,:,k,-1)*xyza_Mtx(:,:,k+1,0))
       xyz_RHS(:,:,k) = (xyz_RHS(:,:,k)*xyza_Mtx(:,:,k+1,0) - xyza_Mtx(:,:,k,1)*xyz_RHS(:,:,k+1)) &
            &              / (xyza_Mtx(:,:,k,-1)*xyza_Mtx(:,:,k+1,0))
       xyza_Mtx(:,:,k,-1) = 1d0
       xyza_Mtx(:,:,k, 1) = 0d0       
    end do
    
  end subroutine Solve_TriDiagSystem_Forward
  
  subroutine Solve_TriDiagSystem_Backward( &
       & xyza_Mtx, xyz_RHS )

    real(DP), intent(inout) :: xyza_Mtx(0:imax-1,jmax,1:kmax,-1:1)
    real(DP), intent(inout) :: xyz_RHS(0:imax-1,jmax,1:kmax)
    
    integer :: k

!!$    k = 1
!!$    xyz_RHS(:,:,k) = xyz_RHS(:,:,k)/xyza_Mtx(:,:,k,0)
    do k = 2, kmax
       xyz_RHS(:,:,k) = (xyz_RHS(:,:,k) - xyza_Mtx(:,:,k,-1)*xyz_RHS(:,:,k-1))/xyza_Mtx(:,:,k,0)
    end do
    
  end subroutine Solve_TriDiagSystem_Backward

  subroutine dcpam_sfc_implicit_coupling_Init()
    
    allocate( xyza_UVMtx(0:imax-1,jmax,1:kmax,-1:1) )
    allocate( xyza_TempMtx(0:imax-1,jmax,1:kmax,-1:1) )
    allocate( xyza_QMixMtx(0:imax-1,jmax,1:kmax,-1:1) )    

  end subroutine dcpam_sfc_implicit_coupling_Init
  
end module dcpam_sfc_implicit_coupling_mod
