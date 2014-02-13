subroutine create_us ( time, u, uk )
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u, uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2) :: u_smooth
  real (kind=pr) :: R

  u_smooth = 0.d0
  us = 0.d0
  
  if ( iActive == "chantalat" ) then
    !-- use chantalat's method to compute the smooth extension
    call active_prolongation_chantalat ( u, u_smooth )
    us = u_smooth + u_BC
    
  elseif (iActive == "dave" ) then
    !-- use dave's method to compute the smooth extension
    call active_prolongation_dave ( u, u_smooth )
    us = u_smooth + u_BC
    
  elseif (iActive == "passive" ) then
    !-- us does not depend on u: classic penalization
    us = u_BC
    
  else
    write (*,*) "error: iactive is wrong."
    stop
  endif
  
end subroutine create_us



!===============================================================================



subroutine compute_beta_field (u, beta)
  ! computes the BETA field, which is the normal derivatives
  ! beta = (n \cdot \nabla) u
  use share_vars
  use FieldExport
  implicit none 
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: beta
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: ux_x,ux_y,uy_x,uy_y
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work, work2
  integer :: ix,iy
  
  !----------------------------------------------------------------------------- 
  !-- compute beta field (finite differences)
  !-----------------------------------------------------------------------------
  !$omp parallel do private(ix)
  do ix=0,nx-1
    ux_x(ix,:) = (u(getindex(ix+1,nx),:,1)-u(getindex(ix-1,nx),:,1))/(2.d0*dx)
    uy_x(ix,:) = (u(getindex(ix+1,nx),:,2)-u(getindex(ix-1,nx),:,2))/(2.d0*dx)
  enddo  
  !$omp end parallel do
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    ux_y(:,iy) = (u(:,getindex(iy+1,ny),1)-u(:,getindex(iy-1,ny),1))/(2.d0*dy)
    uy_y(:,iy) = (u(:,getindex(iy+1,ny),2)-u(:,getindex(iy-1,ny),2))/(2.d0*dy)
  enddo  
  !$omp end parallel do
  
  !----------------------------------------------------------------------------- 
  !-- compute beta field (spectral)
  !-----------------------------------------------------------------------------  
!   call coftxy(u(:,:,1),work)
!   call cofdx(work, work2)
!   call cofitxy(work2,ux_x)
!   
!   call coftxy(u(:,:,1),work)
!   call cofdy(work, work2)
!   call cofitxy(work2,ux_y)
!     
!   call coftxy(u(:,:,2),work)
!   call cofdy(work, work2)
!   call cofitxy(work2,uy_y)
!   
!   call coftxy(u(:,:,2),work)
!   call cofdx(work, work2)
!   call cofitxy(work2,uy_x)
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    beta(:,iy,1) = normals(:,iy,1)*ux_x(:,iy) + normals(:,iy,2)*ux_y(:,iy)
    beta(:,iy,2) = normals(:,iy,1)*uy_x(:,iy) + normals(:,iy,2)*uy_y(:,iy)
!     beta(:,iy,1) = beta(:,iy,1)*(1.d0-mask(:,iy)*eps)
!     beta(:,iy,2) = beta(:,iy,2)*(1.d0-mask(:,iy)*eps)
  enddo  
  !$omp end parallel do  
end subroutine compute_beta_field



!===============================================================================



subroutine active_prolongation_chantalat ( u, u_smooth )
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: u_smooth  
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2) :: beta
  real (kind=pr) :: CFL_act, umax=1.d0, Tend, dt,R,R1
  integer :: ix,iy,nt2,it
  
  u_smooth = 0.d0
  
  !-- compute the field of normal derivatives
  call compute_beta_field ( u, beta )
  
  !-----------------------------------------------------------------------------
  !-- prolongate beta field using advection/diffusion
  !-----------------------------------------------------------------------------
  CFL_act = 0.98d0
  umax = 1.d0
  dt   = CFL_act*dx/umax
  Tend = delta / umax
  nt2  = nint(Tend/dt)
  
  do it=1, nt2
    call RK4 ( beta(:,:,1), dt )
    call RK4 ( beta(:,:,2), dt )
  enddo
  
  !-----------------------------------------------------------------------------
  !-- construct u_smooth
  !-----------------------------------------------------------------------------  
  x0 = xl/2.d0
  y0 = yl/2.d0
  R1 = 0.50d0
  
  if (imask == 'lamballais') then ! this case distinguished between outer and inner cylinder
      !$omp parallel do private(iy,ix,R)
      do ix=0,nx-1
      do iy=0,ny-1
        R = dsqrt( (dble(ix)*dx-x0)**2 + (dble(iy)*dy-y0)**2 ) 
        if (R<=1.25*R1) then
          u_smooth(ix,iy,1) = (mask(ix,iy)*eps)*phi(ix,iy)*beta(ix,iy,1)
          u_smooth(ix,iy,2) = (mask(ix,iy)*eps)*phi(ix,iy)*beta(ix,iy,2)
        else
          u_smooth(ix,iy,:) = 0.d0
        endif
      enddo
      enddo
      !$omp end parallel do
  else
      !$omp parallel do private(iy,ix)
      do ix=0,nx-1
      do iy=0,ny-1
          u_smooth(ix,iy,1) = (mask(ix,iy)*eps)*phi(ix,iy)*beta(ix,iy,1)
          u_smooth(ix,iy,2) = (mask(ix,iy)*eps)*phi(ix,iy)*beta(ix,iy,2)
      enddo
      enddo
      !$omp end parallel do  
  endif
end subroutine active_prolongation_chantalat


!===============================================================================



subroutine RHS_central ( field, rhs,  dt )
  use share_vars
  use FieldExport
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent(in) :: field
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent(out) :: rhs
  real(kind=pr), intent(in) :: dt
  real(kind=pr) :: lambda, grad_x, grad_y, laplace, field_xx, field_yy
  integer :: ix,iy
  
  ! diffusion constant (experimental value)
  lambda = 0.5d0*0.5d0*dt
  
  !$omp parallel do private(iy,ix,grad_x,grad_y,field_xx,field_yy,laplace)
  do iy=0,ny-1
    do ix=0,nx-1
      if (mask(ix,iy)>0.d0) then
        !-- transport term
        grad_x = (field(getindex(ix+1,nx),iy)-field(getindex(ix-1,nx),iy))/(2.d0*dx)
        grad_y = (field(ix,getindex(iy+1,ny))-field(ix,getindex(iy-1,ny)))/(2.d0*dy)
        rhs(ix,iy) = normals(ix,iy,1)*grad_x + normals(ix,iy,2)*grad_y
        
        !-- diffusion 
        field_xx = field(getindex(ix-1,nx),iy)-2.d0*field(ix,iy)+field(getindex(ix+1,nx),iy) 
        field_yy = field(ix,getindex(iy-1,ny))-2.d0*field(ix,iy)+field(ix,getindex(iy+1,ny))
        laplace = field_xx/dx**2 + field_yy/dy**2
        
        rhs(ix,iy) = rhs(ix,iy) + lambda*laplace
      else
        rhs(ix,iy) = 0.d0
      endif      
    enddo
  enddo  
  !$omp end parallel do    
end subroutine RHS_central


!===============================================================================

subroutine RK4 ( field, dt )
  use share_vars
  use FieldExport
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent(inout) :: field
  real(kind=pr), dimension(0:nx-1,0:ny-1):: k1,k2,k3,k4
  real(kind=pr), intent(in) :: dt
  real(kind=pr) :: lambda
  integer :: ix,iy
  
  call RHS_central( field, k1, dt)
  call RHS_central( field+0.5d0*dt*k1, k2, dt)
  call RHS_central( field+0.5d0*dt*k2, k3, dt)
  call RHS_central( field+dt*k3, k4, dt)

  !$omp parallel do private(iy)
  do iy=0,ny-1
    field(:,iy)=field(:,iy)+(dt/6.d0) *( k1(:,iy)+2.d0*k2(:,iy)+2.d0*k3(:,iy)+k4(:,iy) ) 
  enddo  
  !$omp end parallel do    
end subroutine RK4


!===============================================================================

subroutine active_prolongation_dave ( u, u_smooth )
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: u_smooth
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2) :: beta
  real (kind=pr) :: s, b0, b1,ux_BC_interp,uy_BC_interp,R
  real (kind=pr) :: beta_x_interp,beta_y_interp, xi_x, xi_y,LinearInterpolation
  integer :: ix,iy
  
  u_smooth = 0.d0
  
  !-- compute the field of normal derivatives
  call compute_beta_field ( u, beta )
     
  !-----------------------------------------------------------------------------
  ! loop over points in the boundary layer where we intent to construct us
  !-----------------------------------------------------------------------------  
  !$omp parallel do private(iy,ix,xi_x,xi_y,ux_BC_interp,uy_BC_interp,beta_x_interp,beta_y_interp,s,b0,b1)
  do ix=0, nx-1
    do iy=0, ny-1
      !-- this point lies inside the obstacle and in the neighborhood of the interface
      if ((phi(ix,iy) >= -delta).and.(phi(ix,iy) <=0.d0) ) then                
        !-- coordinates of closest point on the interface
        xi_x = dble(ix)*dx - normals(ix,iy,1)*phi(ix,iy)
        xi_y = dble(iy)*dy - normals(ix,iy,2)*phi(ix,iy)
        
        !---------------------------------------
        ! interpolate values of beta
        !---------------------------------------
        ! inhomogeneous dirichlet:
        ux_BC_interp = 0.d0!LinearInterpolation ( xi_x,xi_y,u_BC(:,:,1),0.d0,0.d0,xl-dx,yl-dy )
        uy_BC_interp = 0.d0!LinearInterpolation ( xi_x,xi_y,u_BC(:,:,2),0.d0,0.d0,xl-dx,yl-dy )
        ! normal derivative
        beta_x_interp = LinearInterpolation ( xi_x,xi_y,beta(:,:,1),0.d0,0.d0,xl-dx,yl-dy )
        beta_y_interp = LinearInterpolation ( xi_x,xi_y,beta(:,:,2),0.d0,0.d0,xl-dx,yl-dy )
        
        !-- s is the dimensionless boundary layer coordinate
        s = -phi(ix,iy) / delta
        
        !-- basis functions: (they differ from dave's work)
        b0 = 3.d0*s**4 -4.d0*s**3+1.d0
        b1 = s**3 -2.d0*s**2 + s
           
        !-- construct us using both basis functions
        u_smooth(ix,iy,1) = ux_BC_interp*b0 - delta*b1*beta_x_interp
        u_smooth(ix,iy,2) = uy_BC_interp*b0 - delta*b1*beta_y_interp
      endif
    enddo
  enddo
  !$omp end parallel do 
  
  if ( iMask == "lamballais" ) then
    !$omp parallel do private(ix,iy,R)
    do ix=0,nx-1
      do iy=0,ny-1
        R = dsqrt( (dble(ix)*dx-x0)**2 + (dble(iy)*dy-y0)**2 )
        if (R > 1.0d0) then
          u_smooth(ix,iy,:) = 0.d0
        endif
      enddo
    enddo
    !$omp end parallel do
  endif  
  
  
end subroutine active_prolongation_dave



real (kind=pr) function LinearInterpolation (x_target, y_target, field2, x1_box, y1_box, x2_box, y2_box )
!  LINEAR Interpolation in a field. The field is of automatic size, indices starting with 0 both. The domain is 
!  defined by x1_box,y1_box and x2_box,y2_box. The target coordinates should lie within that box.
!  NOTE: attention on the upper point of the box. In the rest of the code, which is periodic, the grid is 0:nx-1
!        but the lattice spacing is yl/nx. This means that the point (nx-1) has NOT the coordinate yl but yl-dx
!        (otherwise this point would exist two times!)
!  NOTE3: Coordinates in the box are a constant source for errors. be careful and note that x1_box is NOT ZERO
  use share_vars
  implicit none
  integer :: i,j
  real (kind=pr) :: x,y,x_1,y_1,x_2,y_2,dx1, dy1, R1,R2
  real (kind=pr), intent (in) :: field2(0:nx-1,0:ny-1), x_target, y_target, x1_box, y1_box, x2_box, y2_box

  dx1 = (x2_box-x1_box)/real(size(field2,1)-1 )
  dy1 = (y2_box-y1_box)/real(size(field2,2)-1 )


  if ( (x_target > x2_box).or.(x_target < x1_box).or.(y_target > y2_box).or.(y_target < y1_box) ) then
    write(*,'("target: (",es11.4,"|",es11.4,") but box: (",es11.4,"|",es11.4,") x (",es11.4,"|",es11.4,")")') &
    x_target, y_target, x1_box, y1_box, x2_box, y2_box
    write (*,*) "!!! LinearInterpolation: target coordinates not in the field."
    stop
  endif

  i=int((x_target-x1_box)/dx1)  ! attention on index shift because of automatic array
  j=int((y_target-y1_box)/dy1)

  x_1= real(i)*dx1 + x1_box
  y_1= real(j)*dy1 + y1_box
  x_2= dx1*real(i+1) + x1_box
  y_2= dy1*real(j+1) + y1_box
  R1 = (x_2-x_target)*field2(i,j)/dx1   + (x_target-x_1)*field2(i+1,j)/dx1
  R2 = (x_2-x_target)*field2(i,j+1)/dx1 + (x_target-x_1)*field2(i+1,j+1)/dx1

  LinearInterpolation = (y_2-y_target)*R1/dy1 + (y_target-y_1)*R2/dy1

  return

end function LinearInterpolation