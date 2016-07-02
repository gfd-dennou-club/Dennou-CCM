module grid_mapping_util

  use dc_types
  
  implicit none
  private

  public :: gen_gridmapfile_lonlat2lonlat
  public :: set_mappingTable_interpCoef
  
contains

  subroutine gen_gridmapfile_lonlat2lonlat( &
       & filename, &
       & x_LonS, y_LatS, x_LonR, y_LatR &
       & )
    character(*), intent(in) :: filename
    real(DP), dimension(:), intent(in) :: x_LonS, y_LatS
    real(DP), dimension(:), intent(in) :: x_LonR, y_LatR

    integer :: nxr, nyr, nxs, nys
    integer :: ir, jr, is, js
    integer :: i, j
    integer :: deviceID

    real(DP) :: dlon_r!, dlat_r
    real(DP) :: dlon_s!, dlat_s
    real(DP) :: coef(4)
    logical :: extp_flag
    
    nxr = size(x_LonR); nyr = size(y_LatR)
    nxs = size(x_LonS); nys = size(y_LatS)
    
    dlon_r = 360d0/dble(nxr); !dlat_r = 180d0/(nyr - 1)
    dlon_s = 360d0/dble(nxs); !dlat_s = 180d0/(nys - 1)
    
    deviceID = 10
    open(unit=deviceID, file=trim(filename), status='replace')
    do jr=1, nyr
       !       js = int(dlat_r*(jr-1)/dlat_s) + 1
       js = get_correspondID_latS(y_LatR(jr), y_LatS, extp_flag)

       if(nxr == 1) then
          if(.not. extp_flag) then
             coef(1:2) = cal_coef_axisym(y_LatR(jr), y_LatS(js), y_LatS(js+1))
             coef(1:2) = coef(1:2)/nxs
             do is=1, nxs
                write(deviceID, *) 1, jr, is,                      js,  coef(1)
                write(deviceID, *) 1, jr, is,            mod(js,nys)+1, coef(2)
             end do
          else
             do is=1, nxs
                write(deviceID, *), 1, jr, is, js, 1d0/nxs
             end do
          end if
       else if(nxs == 1) then
          if(.not. extp_flag) then
             coef(1:2) = cal_coef_axisym(y_LatR(jr), y_LatS(js), y_LatS(js+1))
             do ir=1, nxr
                write(deviceID, *) ir, jr, 1,                      js,  coef(1)
                write(deviceID, *) ir, jr, 1,            mod(js,nys)+1, coef(2)
             end do
          else
             do ir=1, nxr
                write(deviceID, *) ir, jr, 1, js, 1d0
             end do
          end if
       else
          
          do ir=1, nxr
             is = int(dlon_r*(ir-1)/dlon_s) + 1
             coef(:) = cal_coef(x_LonR(ir), y_LatR(jr), &
                  & x_LonS(is), y_LatS(js), x_LonS(mod(is,nxs)+1), y_LatS(js+1) &
                  & )
             write(deviceID, *) ir, jr, is,                      js,  coef(1)     
             write(deviceID, *) ir, jr, mod(is,nxs)+1,           js,  coef(2)
             write(deviceID, *) ir, jr, mod(is,nxs)+1, mod(js,nys)+1, coef(3)
             write(deviceID, *) ir, jr, is,            mod(js,nys)+1, coef(4)
          end do
             
       end if
    end do

    close(deviceID)
  contains
    function get_correspondID_latS(latR, y_LatS, extp_flag) result(latSID)
      real(DP), intent(in) :: latR
      real(DP), intent(in) :: y_LatS(:)
      logical, intent(out) :: extp_flag
      integer :: latSID
      
      integer :: j, nys

      nys = size(y_LatS)
      do j=1, nys-1
         if(y_LatS(j) < latR .and. latR <= y_LatS(j+1)) then
            latSID = j; extp_flag = .false.
            return
         end if
      end do

      !
      extp_flag = .true.
      if(latR <= y_LatS(1)) latSID = 1
      if(latR >  y_LatS(nys)) latSID = nys
      
    end function get_correspondID_latS

    function cal_coef(xc, yc, x1, y1, x3, y3) result(coef)
      real(DP), intent(in) :: xc, yc, x1, y1, x3, y3
      real(DP) :: coef(4)

      real(DP) :: a1, a2, b1, b2

      a1 = (xc - x1)/(x3 - x1); a2 = 1d0 - a1
      b1 = (yc - y1)/(y3 - y1); b2 = 1d0 - b1

      coef(:) = (/ a2*b2, a1*b2, a1*b1, a2*b1 /)
      
    end function cal_coef

    function cal_coef_axisym(yc, y1, y2) result(coef)
      real(DP), intent(in) :: yc, y1, y2
      real(DP) :: coef(2)
      
      real(DP) :: b1

      b1 = (yc - y1)/(y2 - y1)
      coef(:) = (/ 1d0 - b1, b1 /)
    end function cal_coef_axisym
    
  end subroutine gen_gridmapfile_lonlat2lonlat
  
  subroutine set_mappingTable_interpCoef(gridmapfile, GNXS, GNXR, &
       & send_index, recv_index, coef_s )
    character(*), intent(in) :: gridmapfile
    integer, intent(in) :: GNXS, GNXR
    integer, intent(inout), allocatable :: &
         & send_index(:), recv_index(:)
    real(DP), intent(inout), allocatable :: coef_s(:)

    integer :: deviceID
    integer :: nInterpOp
    integer :: ir, jr, is, js
    real(DP) :: coef
    integer :: n
    
    deviceID = 10
    open(deviceID, file=trim(gridmapfile)) 

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
  
!!!!!!!!!
  
  real(DP) function rad2deg(rad)
    real(DP), intent(in) :: rad
    rad2deg = rad*180d0/acos(-1d0)
  end function rad2deg
  
end module grid_mapping_util
