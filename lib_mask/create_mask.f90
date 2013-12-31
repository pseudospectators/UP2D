module masks

  implicit none
  contains
!==============================================================================================================================================
subroutine create_mask (time)
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in) :: time
  real (kind=pr) :: R
  integer :: ix, iy

  mask = 0.d0 !initialize masks as 0.0
  us = 0.d0
  
  select case (iMask)
  case('cylinder')
    x0 = xl/2.d0
    y0 = yl/2.d0
    
    !$omp parallel do private(ix,iy,R)
    do ix=0,nx-1  
      do iy=0,ny-1
        R = sqrt( (real(ix)*dx-x0)**2 +(real(iy)*dy-y0)**2 )
        if (R <= 1.d0) then
          mask(ix,iy) = 1.d0
        endif
      enddo      
    enddo
    !$omp end parallel do
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




subroutine SmoothStep (f,x,t,h)
  use share_vars
  implicit none
  real (kind=pr), intent (out) :: f
  real (kind=pr), intent (in)  :: x,t,h

  if (x<=t-h) then
    f = 1.d0
  elseif (((t-h)<x).and.(x<(t+h))) then
    f = 0.5d0*(1.0d0+cos((x-t+h)*pi/(2.d0*h)) )
  else
    f = 0.0d0
  endif
end subroutine SmoothStep

end module masks