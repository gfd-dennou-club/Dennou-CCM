module DSFCM_Admin_Grid_mod

  ! モジュール引用; Use statements
  !
 
  !* gtool
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  
  implicit none
  private

  integer, public :: IA
  integer, public :: IS
  integer, public :: IE
  integer, public :: IM
  integer, public :: IHALO
  
  integer, public :: JA
  integer, public :: JS
  integer, public :: JE
  integer, public :: JM
  integer, public :: JHALO

  public :: DSFCM_Admin_Grid_Init
  public :: DSFCM_Admin_Grid_Final
  
contains

  subroutine DSFCM_Admin_Grid_Init(GNX, GNY)
    integer, intent(in) :: GNX
    integer, intent(in) :: GNY

    IM = GNX
    JM = GNY
    IHALO = 1
    JHALO = 1
    
    IS = IHALO + 1
    IE = IS + IM - 1
    IA = IM + 2*IHALO
    
    JS = JHALO + 1
    JE = JS + JM - 1
    JA = JM + 2*JHALO
    
  end subroutine DSFCM_Admin_Grid_Init

  subroutine DSFCM_Admin_Grid_Final()
  end subroutine DSFCM_Admin_Grid_Final
  
end module DSFCM_Admin_Grid_mod
