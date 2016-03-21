subroutine sponge_mask (time)
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  real (kind=pr) :: R
  integer :: ix, iy

  if (use_sponge == 0) return

  mask_sponge = 0.d0

  select case (iSpongeType)
  case('frame')
      ! do ix=0,nx-1
      !   do iy==,
    case('none')
      mask = 0.d0
    case default
      write (*,*) "mask not defnd", iMask
      stop
  end select

  !$omp parallel do private(iy)
  do iy=0,ny-1
    mask(:,iy) = mask(:,iy) / eps
  enddo
  !$omp end parallel do

end subroutine sponge_mask
