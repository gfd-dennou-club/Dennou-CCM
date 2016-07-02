subroutine interpolate_data(recv_model, send_model, mapping_tag, sn1, sn2, send_data, & 
                            rn1, rn2, recv_data, num_of_data, tn, exchange_tag)
  use interpolation_data_latlon_mod, only : interpolate_data_latlon
  implicit none
  character(len=*), intent(IN) :: recv_model, send_model
  integer, intent(IN) :: mapping_tag
  integer, intent(IN) :: sn1, sn2
  real(kind=8), intent(IN) :: send_data(sn1,sn2)
  integer, intent(IN) :: rn1, rn2
  real(kind=8), intent(INOUT) :: recv_data(rn1,rn2)
  integer, intent(IN) :: num_of_data
  integer, intent(IN) :: tn
  integer, intent(IN) :: exchange_tag(tn)

  call interpolate_data_latlon(recv_model, send_model, send_data, recv_data, num_of_data, mapping_tag, exchange_tag)

end subroutine interpolate_data
