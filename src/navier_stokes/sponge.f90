subroutine sponge_mask (time)
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  real (kind=pr) :: R
  integer :: ix, iy

  if (use_sponge == 0) return

  mask_sponge = 0.d0

  select case (iSpongeType)
  case ('everywhere')
    mask_sponge = 1.d0
  case('none')
    mask_sponge = 0.d0
  case default
    write (*,*) "mask not defnd", iSpongeType
    stop
  end select

end subroutine sponge_mask
