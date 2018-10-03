program EBM_main
  use dc_types
  use l_module
  implicit none

  integer, parameter :: jm = 1024
  integer, parameter :: nm = jm
  integer, parameter :: lm = nm + 1
  real(DP), parameter :: PI = acos(-1d0)

!!$  real(DP), parameter :: IR_A = 298.1d0
!!$  real(DP), parameter :: IR_B = 2.89d0
!!$  real(DP), parameter :: COALBO_P0 = 1d0
!!$  real(DP), parameter :: COALBO_P2 = 0d0
!!$  real(DP), parameter :: COALBI    = 0.5d0
!!$  real(DP), parameter :: D = 0.136d0

!!$  real(DP), parameter :: IR_A = 212d0
!!$  real(DP), parameter :: IR_B = 1.7d0
!!$  real(DP), parameter :: COALBO_P0 = 0.7d0
!!$  real(DP), parameter :: COALBO_P2 = -0.078d0
!!$  real(DP), parameter :: COALBI    = 0.38d0
!!$  real(DP), parameter :: D = 2.2d13/(2d0*PI*6.373d6**2)
  
!!$  real(DP), parameter :: IR_A = 211.2d0
!!$  real(DP), parameter :: IR_B = 1.55d0*1.5d0
!!$  real(DP), parameter :: COALBO_P0 = 0.697d0
!!$  real(DP), parameter :: COALBO_P2 = -0.0779d0
!!$  real(DP), parameter :: COALBI    = 0.38d0
!!$  real(DP), parameter :: D = 0.382d0*0.41d0
!!$  real(DP), parameter :: SR_P0 = 1d0
!!$  real(DP), parameter :: SR_P2 = -0.482d0

  real(DP), parameter :: IR_A = 300d0
  real(DP), parameter :: IR_B = 3.5d0
  real(DP), parameter :: COALBO_P0 = 1d0
  real(DP), parameter :: COALBO_P2 = 0d0
  real(DP), parameter :: COALBI    = 0.5d0
  real(DP), parameter :: D = 1.5d13/(2d0*PI*6.373d6**2)
  
  integer, parameter :: nParamExp = 200
  real(DP) :: IceLineLatParams(nParamExp)
  real(DP) :: IceLineTempParams(nParamExp)
  integer :: j
  integer :: n
  integer :: m

  call initialize()

!!$  IceLineLatParams(:) = asin(0.95d0) !PI/2d0
  do m=1, nParamExp
     IceLineLatParams(m) = 0d0 + (PI/2d0/dble(nParamExp))*m
  end do
  IceLineTempParams(:) = - 10d0
  do m=1, nParamExp-1
     call solve_EBM( IceLineLatParams(m), IceLineTempParams(m) )
  end do
  
  call shutdown()
  
contains

  subroutine initialize()
    call l_initial(nm, jm)
  end subroutine initialize

  subroutine shutdown()
    call l_finalize()
  end subroutine shutdown

  subroutine solve_EBM(iceline_lat, iceline_temp)
    real(DP), intent(in) :: iceline_lat
    real(DP), intent(in) :: iceline_temp
    
    real(DP) :: Is
    real(DP) :: l_H(lm)
    real(DP) :: tmp
    real(DP) :: l_I(lm)
    
    integer :: j
    integer :: l    
    real(DP) :: x(jm)
    real(DP) :: P2(jm)
    integer :: nm(2)
    real(DP) :: n
    real(DP) :: Pn_xs(lm)
    real(DP) :: Q
    
    x(:) = sin(y_Lat)
    P2(:) = (3d0*x**2 - 1d0)*0.5d0
    call calc_Hn_AlbMod( & 
!!$    call calc_Hn( & 
!!$    call calc_Hn_exact( &
         & iceline_lat, &         
         & l_H )
    
    Pn_xs(1) = 1d0
    Pn_xs(2) = sin(iceline_lat)
    do l=3, lm
!       nm = nm_l(l); n = dble(nm(1))
       n = dble(l-1)
       Pn_xs(l) = ( (2d0*n - 1d0)*Pn_xs(2)*Pn_xs(l-1) - (n-1)*Pn_xs(l-2) )/dble(n)
    end do
!    write(*,*) Pn_xs

    Is = IR_A + IR_B*iceline_temp
    tmp = 0d0    
    do l=1, lm
       n = dble(l-1)
       tmp = tmp + l_H(l) * Pn_xs(l)/( n*(n+1d0)*D + 1d0 )
!       write(*,*) "H", l_H(l)
    end do
!!$    write(*,*) "check", sum(y_Lat_weight(:)*y_H(:))*0.5d0
!!$    write(*,*) "check", sum(y_Lat_weight(:)*y_H(:)*P2)*0.5d0*5d0
!!$    write(*,*) "check*P2", sum(y_Lat_weight(:)*P2)*0.5d0
!!$    write(*,*) "check*P2*P2", sum(y_Lat_weight(:)*P2*P2)
!!$    write(*,*) "check*P2*P2*P2", sum(y_Lat_weight(:)*P2*P2*P2)*0.5d0

    Q = Is/tmp
    
    do l=1, lm
       !       nm = nm_l(l); n = dble(nm(1))
       n = dble(l-1)
       l_I(l) = Q*l_H(l)/( n*(n+1d0)*D + 1d0 )
!!$       write(*,*) "I", l_I(l)
       
       l_I(l) = l_I(l)/sqrt(2d0*n + 1d0)
    end do

    write(*,*) Q*4d0, iceline_lat*180d0/PI
!!$    write(*,*) "T:", (y_l(l_I) - IR_A)/IR_B
!!$    write(*,*) "ASR:", Q*y_S*y_CoAlbedo

    do j=1, jm
!       write(*,*) "lat", y_Lat(j)*180d0/PI, "T:", 
    end do
  end subroutine solve_EBM

  subroutine calc_Hn(iceline_lat, l_H)
    use l_module, only: y_Lat_Weight
    
    real(DP), intent(in) :: iceline_lat
    real(DP), intent(out) :: l_H(lm)

    real(DP) :: y_CoAlbedo(jm)
    real(DP) :: y_S(jm)
    real(DP) :: y_H(jm)
    
    real(DP) :: x(jm)
    real(DP) :: P2(jm)
    integer :: j
    integer :: l
    real(DP) :: n
    
    x(:) = sin(y_Lat)
    P2(:) = (3d0*x**2 - 1d0)*0.5d0
    
    do j=1, jm
       if ( abs(y_Lat(j)) >= iceline_lat ) then
          y_CoAlbedo(j) = COALBI
       else
          y_CoAlbedo(j) = COALBO_P0 + COALBO_P2*P2(j)
       end if
       y_S(j) = SR_P0 + SR_P2*P2(j)

       y_H(j) = y_S(j)*y_CoAlbedo(j)
!       write(*,*) "lat", y_Lat(j), "CoAlb", y_CoAlbedo(j), "S", y_S(j)
    end do

    l_H = l_y(y_H)
    do l=1, lm
       n = dble(l-1)
       l_H(l) = sqrt(2d0*n + 1d0)*l_H(l)
    end do
    
  end subroutine calc_Hn

  subroutine calc_Hn_AlbMod(iceline_lat, l_H)
    use l_module, only: y_Lat_Weight
    
    real(DP), intent(in) :: iceline_lat
    real(DP), intent(out) :: l_H(lm)

    real(DP) :: y_CoAlbedo(jm)
    real(DP) :: y_S(jm)
    real(DP) :: y_H(jm)
    
    real(DP) :: y_x(jm)
    real(DP) :: v_x(0:jm)
    real(DP) :: P2(jm)
    integer :: j
    integer :: l
    real(DP) :: n

    real(DP) :: v_Lat(0:jm)
    real(DP) :: icefrac
    real(DP) :: judge_sign
    

    v_Lat(0) = - PI/2d0
    do j = 1, jm-1
       v_Lat(j) = asin( y_Lat_Weight(j) + sin(v_Lat(j-1)) )
    end do
    v_Lat(jm) = PI/2d0

    y_x(:) = sin(y_Lat)
    v_x(:) = sin(v_Lat)
    P2(:) = (3d0*y_x**2 - 1d0)*0.5d0
    
!!$    write(*,*) "iceline_lat=", iceline_lat
    
    do j=1, jm
       judge_sign =  (v_Lat(j-1) - sign(iceline_lat,y_Lat(j))) &
            &       *(v_Lat(j  ) - sign(iceline_lat,y_Lat(j)))
       
       if ( judge_sign < 0d0 ) then
          if (y_Lat(j) > 0d0) then
             icefrac = (v_x(j) - sin(iceline_lat))/y_Lat_Weight(j)
!!$             icefrac = (v_Lat(j) - iceline_lat)/(v_Lat(j) - v_Lat(j-1))
          else
             icefrac = (- sin(iceline_lat) - v_x(j-1))/y_Lat_Weight(j)
!!$             icefrac = (- iceline_lat - v_Lat(j-1))/(v_Lat(j) - v_Lat(j-1))
          end if
!!$          write(*,*) "inside"
       else if ( abs(y_Lat(j)) >= iceline_lat ) then
          icefrac = 1d0
!!$          write(*,*) "ice"
       else
          icefrac = 0d0
!!$          write(*,*) "ocean"
       end if
       y_CoAlbedo(j) =           icefrac*COALBI                        &
            &          + (1d0 - icefrac)*(COALBO_P0 + COALBO_P2*P2(j))
       
       y_S(j) = SR_P0 + SR_P2*P2(j)

       y_H(j) = y_S(j)*y_CoAlbedo(j)
!!$       write(*,*) "lat", y_Lat(j), "CoAlb", y_CoAlbedo(j), "S", y_S(j), "icefrac", icefrac
    end do

    
    l_H = l_y(y_H)
    do l=1, lm
       n = dble(l-1)
       l_H(l) = sqrt(2d0*n + 1d0)*l_H(l)
    end do
    
  end subroutine calc_Hn_AlbMod
  
  subroutine calc_Hn_exact(iceline_lat, l_H)
    real(DP), intent(in) :: iceline_lat
    real(DP), intent(out) :: l_H(lm)

    real(DP) :: y_CoAlbedo(jm)
    real(DP) :: y_S(jm)
    real(DP) :: y_H(jm)
    
    real(DP) :: x(jm)
    real(DP) :: P2(jm)
    integer :: j
    real(DP) :: G000, G002, G022, G222, G024, G004, G224
    real(DP) :: x_s

    x_s = sin(iceline_lat)
    x(:) = sin(y_Lat)
    P2(:) = (3d0*x**2 - 1d0)*0.5d0

    G000 = x_s
    G002 = 0.5d0*(x_s**3 - x_s)
    G022 = 0.25d0*(9d0/5d0*x_s**5 - 2d0*x_s**3 + x_s)
    G222 = 0.125d0*(-x_s + 3d0*x_s**3 - 0.2d0*(27d0*x_s**5) + (27d0*x_s**7)/7d0)
    G004 = 0.125d0*(7d0*x_s**5 - 10d0*x_s**3 + 3d0*x_s)
    G024 = 0.0625d0*(15d0*x_s**7 - 25d0*x_s**5 + 13d0*x_s**3 -3d0*x_s)
    G224 = x_s*(35d0/32d0*x_s**8 - 15d0/7d0*x_s**6 + 121d0/80d0*x_s**4 - 0.5d0*x_s**2 + 3d0/32d0)
!!$    write(*,*) G000, G002, G022, G222, G004, G224
!!$    stop
    
    do j=1, jm
       if ( abs(y_Lat(j)) >= iceline_lat ) then
          y_CoAlbedo(j) = COALBI
       else
          y_CoAlbedo(j) = COALBO_P0 + COALBO_P2*P2(j)
       end if
       y_S(j) = SR_P0 + SR_P2*P2(j)

       y_H(j) = y_S(j)*y_CoAlbedo(j)
!       write(*,*) "lat", y_Lat(j), "CoAlb", y_CoAlbedo(j), "S", y_S(j)
    end do

    l_H(:)= 0d0
    l_H(1) = (COALBO_P0 - COALBI)*(x_s + SR_P2*G002) &
         &  + COALBO_P2*(G002 + SR_P2*G022) + COALBI
    l_H(3) = 5d0*(COALBO_P0 - COALBI)*(G002 + SR_P2*G022) &
         &  + 5d0*COALBO_P2*(G022 + SR_P2*G222) + SR_P2*COALBI
    l_H(5) = 9d0*(COALBO_P0 - COALBI)*(G004 + SR_P2*G024) &
         &  + 9d0*COALBO_P2*(G024 + SR_P2*G224)
!!$    write(*,*) l_H(1:5)
!!$    stop
  end subroutine calc_Hn_exact
  
end program EBM_main
