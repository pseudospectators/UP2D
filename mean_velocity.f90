subroutine mean_flow (uk)
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent(inout) :: uk
  
  select case (iMeanFlow)
  case ('none')
    ! don't do anything
  case ('constant')
    uk(0:1,0:1,1) = ux_mean
    uk(0:1,0:1,2) = uy_mean
  end select
  
end subroutine mean_flow