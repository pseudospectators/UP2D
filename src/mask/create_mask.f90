module masks
  implicit none
contains

  !===============================================================================
  subroutine create_mask (time, mask, us)
    use share_vars
    implicit none
    real(kind=pr), intent (in) :: time
    real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask
    real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us
    real(kind=pr) :: R
    integer :: ix, iy

    mask = 0.d0
    us   = 0.d0

    select case (iMask)
    case('cylinder')
      call cylinder(mask, us)
    case('ellipse')
      call ellipse(mask, us)
    case('moving_cylinder')
      call moving_cylinder(time, mask, us)
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

  end subroutine create_mask

  !===============================================================================

  subroutine cylinder(mask, us)
    use share_vars
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
    use share_vars
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
    use share_vars
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
    use share_vars
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

end module masks
