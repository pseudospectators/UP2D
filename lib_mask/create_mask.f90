module masks
  implicit none
  contains
  
!===============================================================================
subroutine create_mask (time)
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  real (kind=pr) :: R
  integer :: ix, iy

  mask = 0.d0
  us   = 0.d0
  
  select case (iMask)
  case('cylinder')
    call cylinder()
  case('lamballais')
    call lamballais_mask()
  case('couette')
    call couette_mask()
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

subroutine lamballais_mask
  use share_vars
  use fieldexport
  implicit none
  real(kind=pr) :: R,R1,R2,R3,normgrad, blend
  real(kind=pr) :: phis(1:3)
  integer :: ix,iy
  integer :: i
  ! set the inner cylinder and the outer circle
  R1=0.50d0
  R2=1.75d0
  R3=2.25d0
  x0=xl/2.d0
  y0=yl/2.d0
  
  !-- compute signed distance function   
  do ix=0,nx-1
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 + (dble(iy)*dy-y0)**2 )      
      phis(1) = (R-R1)
      phis(2) =-(R-R2)
      phis(3) = (R-R3)
      !-- closest interface, but with correct sign
      i = minloc( abs(phis), 1 )
      phi(ix,iy) = phis(i)
      
      !-- outside the cylinder, force exact solution in the ring
      if ( R > R1+dx) then
        u_BC(ix,iy,1) = uex(ix,iy,1)
        u_BC(ix,iy,2) = uex(ix,iy,2)
      endif
    enddo
  enddo  
  
  !-- compute chi function (traditional mask function)
  call phi2chi()
  
  !-- compute normal vectors
  do ix=0,nx-1
    do iy=0,ny-1
      normals(ix,iy,1) = (phi(getindex(ix+1,nx),iy) - phi(getindex(ix-1,nx),iy) )/(2.0d0*dx)
      normals(ix,iy,2) = (phi(ix,getindex(iy+1,ny)) - phi(ix,getindex(iy-1,ny)) )/(2.0d0*dy)
      !-- normalize
      normgrad = dsqrt ( normals(ix,iy,1)**2 + normals(ix,iy,2)**2 )
      if ( normgrad > 1e-10) then
        normals(ix,iy,1) = normals(ix,iy,1) / normgrad
        normals(ix,iy,2) = normals(ix,iy,2) / normgrad
      else
        normals(ix,iy,1) = 0.d0
        normals(ix,iy,2) = 0.d0      
      endif
      !-- blending: reduce normal vectors to a BL close to interface
      call SmoothStep(blend,dabs(phi(ix,iy)),4.d0*dx,3.d0*dx)
      normals(ix,iy,1) = normals(ix,iy,1) * blend
      normals(ix,iy,2) = normals(ix,iy,2) * blend
    enddo
  enddo
end subroutine lamballais_mask


!===============================================================================

subroutine phi2chi
  use share_vars
  use fieldexport
  implicit none
  integer :: ix,iy
  
  do ix=0,nx-1
    do iy=0,ny-1
      if (phi(ix,iy)<=0.0) then
        mask(ix,iy) = 1.d0
      endif
    enddo
  enddo    
  
end subroutine phi2chi

!===============================================================================

subroutine couette_mask
  use share_vars
  use fieldexport
  implicit none
  
end subroutine couette_mask

!===============================================================================

subroutine cylinder
  use share_vars
  implicit none
  integer :: ix,iy
  real(kind=pr)::R
  x0 = xl/2.d0
  y0 = yl/2.d0
  
  !$omp parallel do private(ix,iy,R)
  do ix=0,nx-1  
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 +(dble(iy)*dy-y0)**2 )
      if (R <= 1.d0) then
        mask(ix,iy) = 1.d0
      endif
    enddo      
  enddo
  !$omp end parallel do
end subroutine cylinder

!===============================================================================

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

!===============================================================================

end module masks