!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a module to provide conservative 2nd-order remapping scheme for spherical coordinates. 
!! 
!! @author Kawai Yuta
!!
!! The scheme is based on Jones(1999). 
!!
module grid_mapping_util_jones99

  ! モジュール引用; Use statement
  !

  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  ! 宣言文; Declareration statements
  !      
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !      
  public :: gen_gridmapfile_lonlat2lonlat
  public :: set_mappingTable_interpCoef
  
contains

  subroutine gen_gridmapfile_lonlat2lonlat( &
       & filename,                          & ! (in)
       & x_LonS, y_LatS, x_LonD, y_LatD,    & ! (in)
       & x_LonIntWtS, y_LatIntWtS,          & ! (in)
       & x_LonIntWtD, y_LatIntWtD           & ! (in)
       & )

    character(*), intent(in) :: filename
    real(DP), intent(in) :: x_LonS(:)
    real(DP), intent(in) :: y_LatS(:)
    real(DP), intent(in) :: x_LonD(:)
    real(DP), intent(in) :: y_LatD(:)
    real(DP), intent(in) :: x_LonIntWtS(:)
    real(DP), intent(in) :: y_LatIntWtS(:)
    real(DP), intent(in) :: x_LonIntWtD(:)
    real(DP), intent(in) :: y_LatIntWtD(:)

    integer :: imaxS, jmaxS
    integer :: IAS, ISS, IES, JAS, JSS, JES
    integer :: imaxD, jmaxD
    integer :: IAD, ISD, IED, JAD, JSD, JED

    real(DP), allocatable :: x_LonS_(:)
    real(DP), allocatable :: u_LonS_(:)
    real(DP), allocatable :: y_LatS_(:)
    real(DP), allocatable :: v_LatS_(:)

    real(DP), allocatable :: x_LonD_(:)
    real(DP), allocatable :: u_LonD_(:)
    real(DP), allocatable :: y_LatD_(:)
    real(DP), allocatable :: v_LatD_(:)

    integer :: j
    
    imaxS = size(x_LonS); jmaxS = size(y_LatS)
    imaxD = size(x_LonD); jmaxD = size(y_LatD)    
    
    IAS = imaxS+2; ISS = 2; IES = IAS - 1
    JAS = jmaxS+2; JSS = 2; JES = JAS - 1

    IAD = imaxD+2; ISD = 2; IED = IAD - 1
    JAD = jmaxD+2; JSD = 2; JED = JAD - 1

    allocate(x_LonS_(IAS), y_LatS_(JAS), u_LonS_(IAS), v_LatS_(JAS))
    allocate(x_LonD_(IAD), y_LatD_(JAD), u_LonD_(IAD), v_LatD_(JAD))

    x_LonS_(ISS:IES) = x_LonS; y_LatS_(JSS:JES) = y_LatS
    x_LonD_(ISD:IED) = x_LonD; y_LatD_(JSD:JED) = y_LatD
    
    call calc_edge_coordinate( u_LonS_, v_LatS_,        &
         & x_LonS_, y_LatS_, x_LonIntWtS, y_LatIntWtS,  &
         & ISS, IES, JSS, JES )

    call calc_edge_coordinate( u_LonD_, v_LatD_,        &
         & x_LonD_, y_LatD_, x_LonIntWtD, y_LatIntWtD,  &
         & ISD, IED, JSD, JED )
    
    call gen_gridmapfile_lonlat2lonlatCore( &
       & filename,                              & ! (in)
       & x_LonS_, y_LatS_, x_LonD_, y_LatD_,    & ! (in)
       & u_LonS_, v_LatS_, u_LonD_, v_LatD_,    & ! (in)
       & IAS, ISS, IES, JAS, JSS, JES,          & ! (in)
       & IAD, ISD, IED, JAD, JSD, JED           & ! (in)
       & )

  contains
    subroutine calc_edge_coordinate( x_FI, y_FJ,         &
         & x_CI, y_CJ, x_IAXIS_Weight, y_JAXIS_Weight,   &
         & IS, IE, JS, JE )
      
      integer, intent(in) :: IS, IE, JS, JE
      real(DP), intent(inout) :: x_FI(:)
      real(DP), intent(inout) :: y_FJ(:)
      real(DP), intent(inout) :: x_CI(:)
      real(DP), intent(inout) :: y_CJ(:)
      real(DP), intent(in) :: x_IAXIS_Weight(IS:IE)
      real(DP), intent(in) :: y_JAXIS_Weight(JS:JE)

      integer :: i
      integer :: j
      real(DP), parameter :: PI = acos(-1d0)

      !
      do i = IS-1, IE
         x_FI(i) = 0.5d0 * (x_CI(i) + x_CI(i+1))
      end do
       
      !
      y_FJ(JS-1) = - PI/2d0
      do j = JS, JE-1
         y_FJ(j) = asin( y_JAXIS_Weight(j) + sin(y_FJ(j-1)) )
      end do
      y_FJ(JE) = PI/2d0
    end subroutine calc_edge_coordinate
    
  end subroutine gen_gridmapfile_lonlat2lonlat
  
  subroutine gen_gridmapfile_lonlat2lonlatCore( &
       & filename,                          & ! (in)
       & x_LonS, y_LatS, x_LonD, y_LatD,    & ! (in)
       & u_LonS, v_LatS, u_LonD, v_LatD,    & ! (in)
       & IAS, ISS, IES, JAS, JSS, JES,      & ! (in)
       & IAD, ISD, IED, JAD, JSD, JED       & ! (in)
       & )

    ! モジュール引用; Use statement
    !    
    use dc_iounit, only: &
         & FileOpen
    
    ! 宣言文; Declareration statements
    !                
    integer, intent(in) :: IAS, ISS, IES, JAS, JSS, JES
    integer, intent(in) :: IAD, ISD, IED, JAD, JSD, JED

    character(*), intent(in) :: filename
    real(DP), intent(in) :: x_LonS(IAS)
    real(DP), intent(in) :: y_LatS(JAS)
    real(DP), intent(in) :: x_LonD(IAD)
    real(DP), intent(in) :: y_LatD(JAD)
    real(DP), intent(in) :: u_LonS(IAS)
    real(DP), intent(in) :: v_LatS(JAS)
    real(DP), intent(in) :: u_LonD(IAD)
    real(DP), intent(in) :: v_LatD(JAD)
    
    ! 局所変数
    ! Local variables
    !            
    integer :: iD
    integer :: jD
    integer :: iS
    integer :: jS

    integer :: m
    integer :: n
    integer :: jDLat1, jDLat2
    real(DP) :: DLat
    
    integer :: NXD
    integer :: NXS
    integer :: NYD
    integer :: NYS
    
    integer :: deviceID

    real(DP) :: dlon_r!, dlat_r
    real(DP) :: dlon_s!, dlat_s
    real(DP) :: coef(4)
    logical :: extp_flag
    real(DP) :: PI = acos(-1d0)


    real(DP) :: nk
    integer :: OverWrapSrcRangeX(2)
    integer :: OverWrapSrcRangeY(2)
    
    real(DP), allocatable :: aa_rmapWt(:,:)
    real(DP), allocatable :: aa_rmapWtDLat(:,:)
    
    ! 実行文; Executable statement
    !

    NXD = IED - ISD + 1
    NXS = IES - ISS + 1
    NYD = JED - JSD + 1
    NYS = JES - JSS + 1
    
    call FileOpen(deviceID, file=trim(filename), mode='rw')

!    write(*,*) "*********************"
    do jD=JSD, JED
       do iD=ISD, IED
          write(*,*) "i,j=", iD, jD
          call search_OverwrapRange( OverWrapSrcRangeX, OverWrapSrcRangeY,  & ! (out)
               & u_LonD(iD-1), u_LonD(iD), v_LatD(jD-1), v_LatD(jD) )         ! (in)
          
          call calc_RemappingWeight( aa_rmapWt, aa_rmapWtDLat,              & ! (out)
               & OverWrapSrcRangeX, OverWrapSrcRangeY, iD, jD               & ! (in)
               & )


          do m = 1,size(aa_rmapWt,1)
             do n = 1,size(aa_rmapWt,2)
                iS = OverWrapSrcRangeX(1)+m-1
                jS = OverWrapSrcRangeY(1)+n-1

                write(deviceID,*) &
                     & iD -iSD+1, jD -jSD+1, iS -iSS+1, jS -jSS+1, aa_rmapWt(m,n)

                if (jS==jSS) then
                   jDLat1 = jS;   jDLat2 = jS+1
                else if (jS==jES) then
                   jDLat1 = jS-1; jDLat2 = jS
                else
                   jDLat1 = jS-1; jDLat2 = jS+1
                end if
                DLat = y_LatS(jDLat2) - y_LatS(jDLat1)
                
                write(deviceID,*) &
                     & iD -iSD+1, jD -jSD+1, iS -iSS+1, jDLat1 -jSS+1, -aa_rmapWtDLat(m,n)/DLat
                write(deviceID,*) &
                     & iD -iSD+1, jD -jSD+1, iS -iSS+1, jDLat2 -jSS+1, +aa_rmapWtDLat(m,n)/DLat
                
             end do
          end do          
          
       end do
    end do
    
    close(deviceID)

  contains
    subroutine calc_RemappingWeight( aa_rmapWt, aa_rmapWtDLat,  & ! (out)
         & SrcRangeX, SrcRangeY, iD, jD                         & ! (in)
         & )
      
      real(DP), intent(inout), allocatable :: aa_rmapWt(:,:)
      real(DP), intent(inout), allocatable :: aa_rmapWtDLat(:,:)
      integer, intent(in) :: SrcRangeX(2)
      integer, intent(in) :: SrcRangeY(2)
      integer ,intent(in) :: iD
      integer ,intent(in) :: jD

      integer :: i
      integer :: j
      
      integer :: NXrmap
      integer :: NYrmap

      real(DP) :: Ak
      real(DP) :: An
      real(DP) :: lat1_k, lat2_k
      real(DP) :: lat1_n, lat2_n
      real(DP) :: lat1_nk, lat2_nk
      real(DP) :: DLon_k, DLon_n, DLon_nk
      
      real(DP) :: a_AkSegLat(0:NYS)
      
      if(allocated(aa_rmapWt)) deallocate(aa_rmapWt)
      if(allocated(aa_rmapWtDLat)) deallocate(aa_rmapWtDLat)

      NXrmap = SrcRangeX(2)-SrcRangeX(1)+1
      NYrmap = SrcRangeY(2)-SrcRangeY(1)+1
      allocate( aa_rmapWt(NXrmap,NYrmap)        )
      allocate( aa_rmapWtDLat(NXrmap,NYrmap)    )


      !

      DLon_k  = 2d0*PI/dble(NXD)
      if (NXD == 1) then
         DLon_nk = 2d0*PI/dble(NXS)
         DLon_n  = DLon_nk
      else
         DLon_nk = 2d0*PI/dble(NXD)
         DLon_n  = 2d0*PI
      end if
      
      a_AkSegLat(0) = v_LatD(jD-1)
      do j=1, NYrmap-1
         a_AkSegLat(j) = v_LatS(SrcRangeY(1)+j-1)
      end do
      a_AkSegLat(NYrmap) = v_LatD(jD)

      
      lat1_k = v_LatD(jD-1); lat2_k = v_LatD(jD)
      Ak = (sin(lat2_k) - sin(lat1_k))*DLon_k

      
      ! Calculate a weight, w1nk,  of (7) in Jones(1999). 
      !
      do j = 1, NYrmap
         lat1_nk = a_AkSegLat(j-1); lat2_nk = a_AkSegLat(j)
         aa_rmapWt(:,j) = DLon_nk*(sin(lat2_nk) - sin(lat1_nk))/Ak
      end do
      
      ! Calculate a weight, w2nk,  of (7) in Jones(1999). 
      !
      do j=1, NYrmap
         lat1_nk = a_AkSegLat(j-1); lat2_nk = a_AkSegLat(j)
         lat1_n = v_LatS(SrcRangeY(1)+j-2); lat2_n = v_LatS(SrcRangeY(1)+j-1); 

         An = (sin(lat2_n) - sin(lat1_n))*DLon_n
         aa_rmapWtDLat(:,j) = &
           &   (   (cos(lat2_nk) + lat2_nk*sin(lat2_nk))                            &
           &     - (cos(lat1_nk) + lat1_nk*sin(lat1_nk)) )*DLon_nk/Ak               &
           & - (   (cos(lat2_n ) + lat2_n *sin(lat2_n ))                            &
           &     - (cos(lat1_n ) + lat1_n *sin(lat1_n )) )*DLon_n*aa_rmapWt(:,j)/An

      end do

      
      ! Check
      !
!!$      if(NXS==1)then
!!$         write(*,*) "----"
!!$         write(*,*) "SrcRangeY=", SrcRangeY, ", NYrmap=", NYrmap
!!$         write(*,*) "D: v_Lat=",  v_LatD(jD-1:jD)
!!$         write(*,*) "S: v_Lat=", v_LatS(SrcRangeY(1)-1:SrcRangeY(2))
!!$         write(*,*) "check:", sum(aa_rmapWt)
!!$         write(*,*) "segment:", a_AkSegLat(0:NYrmap)
!!$         write(*,*) "w1_nk:", aa_rmapWt(1,:)
!!$         write(*,*) "w2_nk:", aa_rmapWtDLat(1,:)         
!!$         write(*,*) "----"
!!$      end if

    end subroutine calc_RemappingWeight
    
    subroutine search_OverwrapRange( SrcRangeX, SrcRangeY, & ! (out)
         & DistX1, DistX2, DistY1, DistY2 )

      integer, intent(out) :: SrcRangeX(2)
      integer, intent(out) :: SrcRangeY(2)
      real(DP), intent(in) :: DistX1, DistX2, DistY1, DistY2

      integer :: i
      integer :: j

      SrcRangeX(1) = -1
      SrcRangeY(2) = -1
      
      do j=JSS, JES
         if( v_LatS(j-1) <= DistY1 .and. DistY1 <= v_LatS(j) ) then
            SrcRangeY(1) = j
         end if
         if( v_LatS(j-1) <= DistY2 .and. DistY2 <= v_LatS(j) ) then
            SrcRangeY(2) = j
            exit
         end if
      end do

      if ( NXD == 1 ) then
         SrcRangeX(:) = (/ ISS, IES /)
      else
         SrcRangeX(:) = (/ ISS, ISS /)
      end if

      write(*,*) "DistY", DistY1, DistY2, "SrcRange=", SrcRangeY
      
      if (SrcRangeY(1)<0.or.SrcRangeY(2)<0) then
         write(*,*) "Exception.."
         write(*,*) "v_LatS:", v_LatS(JSS-1:JES)
         stop
      end if
      
    end subroutine search_OverwrapRange
    
  end subroutine gen_gridmapfile_lonlat2lonlatCore

  !---------------------------------------------------------------------------
  
  subroutine set_mappingTable_interpCoef( gridmapfile, GNXS, GNXR, &  ! (in)
       & send_index, recv_index, coef_s                            &  ! (inout)
       & )

    ! モジュール引用; Use statement
    !    
    use dc_iounit, only: &
         & FileOpen
    
    ! 宣言文; Declareration statements
    !    
    character(*), intent(in) :: gridmapfile
    integer, intent(in) :: GNXS
    integer, intent(in) :: GNXR
    integer, intent(inout), allocatable :: send_index(:)
    integer, intent(inout), allocatable :: recv_index(:)
    real(DP), intent(inout), allocatable :: coef_s(:)

    ! 局所変数
    ! Local variables
    !              
    integer :: deviceID
    integer :: nInterpOp
    integer :: ir
    integer :: is
    integer :: jr    
    integer :: js
    real(DP) :: coef
    integer :: n

    ! 実行文; Executable statement
    !

    call FileOpen(deviceID, file=trim(gridmapfile), mode='r')

    !
    !
    nInterpOp = 0
    do
       read(deviceID, *, end=100) ir, jr, is, js, coef
       nInterpOp = nInterpOp + 1
    end do
100 continue

    if(allocated(coef_s)) deallocate(coef_s)
    allocate(send_index(nInterpOp), recv_index(nInterpOp), coef_s(nInterpOp))

    !
    !
    rewind(deviceID)
    do n=1, nInterpOp
       read(deviceID, *, end=200) ir, jr, is, js, coef
       recv_index(n) = ir + GNXR*(jr - 1)
       send_index(n) = is + GNXS*(js - 1)
       coef_s(n) = coef
    end do
200 continue

    close(deviceID)
    
  end subroutine set_mappingTable_interpCoef

  !---------------------------------------------------------------------------
  
  real(DP) function rad2deg(rad)
    real(DP), intent(in) :: rad
    rad2deg = rad*180d0/acos(-1d0)
  end function rad2deg
  
end module grid_mapping_util_jones99
