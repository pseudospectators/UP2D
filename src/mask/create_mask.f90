!===============================================================================
subroutine create_mask (time, mask, us)
  use vars
  use hdf5_wrapper
  implicit none
  real(kind=pr), intent (in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us
  real(kind=pr) :: R
  integer :: ix, iy
  logical, save :: first_call = .true.

  ! if thi sis the first call to the routine, we surely have to draw the mask
  ! but if the obstacle is fixed (iMoving==0) we have to do this ONLY in the first
  ! step.
  if ((first_call .eqv. .false.).and.(iMoving==0)) return

  mask = 0.d0
  us   = 0.d0

  select case (iMask)
  case('cylinder')
    call cylinder(mask, us)
  case('ellipse')
    call ellipse(mask, us)
  case('moving_cylinder')
    call moving_cylinder(time, mask, us)
  case('from_file')
    call read_flusi_hdf5_2d_openmp( infile_mask, mask )
    ! if the maximum value of the mask is not 1.d0, then we have saved mask/eps
    ! in the previous simulation.
    if (maxval(mask) > 1.d0) then
      write(*,*) "The mask we read does not appear to be normalized.."
      write(*,*) "Previous eps=", 1.d0/maxval(mask)
      ! re-normalize to unity
      mask = mask / maxval(mask)
    endif
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


  ! we now are sure that this is no longer the first call
  first_call = .false.

end subroutine create_mask

!===============================================================================

subroutine cylinder(mask, us)
  use vars
  implicit none
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us
  integer :: ix,iy
  real(kind=pr)::R,R0,smooth

  x0 = xl/2.d0
  y0 = yl/2.d0
  R0 = 1.0d0
  smooth = 2.d0*max(dx,dy)

  !$omp parallel do private(ix,iy,R)
  do ix=0,nx-1
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 +(dble(iy)*dy-y0)**2 )
      call SmoothStep(mask(ix,iy), R, R0, smooth)
    enddo
  enddo
  !$omp end parallel do
end subroutine cylinder

!===============================================================================

subroutine ellipse(mask, us)
  use vars
  implicit none
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us
  integer :: ix,iy
  real(kind=pr)::R,R0,x,y

  x0 = xl/2.d0
  y0 = yl/2.d0

  !$omp parallel do private(ix,iy,R,x,y)
  do ix=0,nx-1
    do iy=0,ny-1
      x = dble(ix)*dx-x0
      y = dble(iy)*dy-y0

      R = (x/0.5d0)**2  +  (y/0.1d0)**2
      if (R<= 1.d0) then
        mask(ix,iy) = 1.d0
        us(ix,iy,2) = -1.d0
      endif

    enddo
  enddo
  !$omp end parallel do
end subroutine ellipse

!===============================================================================

subroutine moving_cylinder(time,mask, us)
  use vars
  implicit none
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us

  integer :: ix,iy
  real(kind=pr) :: R, R0, smooth

  R0 = 0.5d0
  smooth = 2.d0*max(dx,dy)
  x0 = xl/2.d0 + sin(2.d0*pi*time)
  y0 = yl/2.d0

  !$omp parallel do private(ix,iy,R)
  do ix=0,nx-1
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 +(dble(iy)*dy-y0)**2 )
      call SmoothStep(mask(ix,iy), R, R0, smooth)
      if (mask(ix,iy) > 0.d0) then
        us(ix,iy,1) = 2.d0*pi*cos(2.d0*pi*time)
        us(ix,iy,2) = 0.d0
      endif
    enddo
  enddo
  !$omp end parallel do
end subroutine moving_cylinder


!===============================================================================

subroutine SmoothStep (f,x,t,h)
  use vars
  implicit none
  real(kind=pr), intent (out) :: f
  real(kind=pr), intent (in)  :: x,t,h

  if (x<=t-h) then
    f = 1.d0
  elseif (((t-h)<x).and.(x<(t+h))) then
    f = 0.5d0*(1.0d0+cos((x-t+h)*pi/(2.d0*h)) )
  else
    f = 0.0d0
  endif
end subroutine SmoothStep
