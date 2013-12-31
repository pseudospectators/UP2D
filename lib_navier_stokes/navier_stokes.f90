! ==========================================================================================================================
function timestep(time, u)
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work2
  real (kind=pr) :: dt1, u_max,timestep
  real (kind=pr), intent (in) :: time
  integer :: iy
  
  !========================================================================================
  !=              Calculate time step                                                     =
  !========================================================================================
  !-- CFL condition for the fluid
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
  work2(:,iy) = u(:,iy,1)**2 + u(:,iy,2)**2
  enddo
  !$omp end parallel do
  
  u_max = sqrt(maxval(work2))  
  dt1 = cfl * min (xl/real(nx), yl/real(ny)) / u_max

  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'u_max', status = 'unknown', access = 'append') ! Append output data file
  write (14,'(4(es11.4,1x))') time, u_max, sum(u(:,:,1))*dx*dy, sum(u(:,:,2))*dx*dy
  close (14)
  !-----------------------------------------------------------------------
  !-- u_max is very very small
  if(u_max.lt.1.e-12) then 
      dt1 = 1.e-3 
  endif

  !-----------------------------------------------------------------------
  !-- Time stepping control for volume penalization
  if ( (dt1 >= 0.9*eps) ) then
    dt1 = min(0.9*eps,dt1) ! time step is smaller than eps (Dmitry, 26 feb 08)
  endif

!   if (.not.((dt1>1e-8).and.(dt1<1.0))) then
!     write (*,*) "the time step is out of range, I give up."
!     stop
!   endif
  
  if ((Tmax - time) < dt1) then
    dt1 = Tmax - time
    write (*,*) "last time step:", dt1
  endif  

  timestep = dt1
end function  
  
  
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------  
  


subroutine cal_nlk (time, u, uk, p, vor, nlk)
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in) :: 					time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: 		vor, p
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: 	nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2
  integer :: iy  

  !-- velocity in phys. space
  call cofitxy (uk(:,:,1), u(:,:,1))
  call cofitxy (uk(:,:,2), u(:,:,2))
  
  !-- compute vorticity
  call curl (uk, work1)  
  call cofitxy(work1, vor)
    
  !-- NL term + penal
  !$omp parallel do private(iy)
  do iy=0,ny-1
    work1(:,iy) = +vor(:,iy)*u(:,iy,2) -mask(:,iy)*(u(:,iy,1)-us(:,iy,1))
    work2(:,iy) = -vor(:,iy)*u(:,iy,1) -mask(:,iy)*(u(:,iy,2)-us(:,iy,2))
  enddo
  !$omp end parallel do
  
  !-- NLK to fourier space
  call coftxy ( work1, nlk(:,:,1) )
  call coftxy ( work2, nlk(:,:,2) )  
end subroutine cal_nlk



subroutine add_pressure (nlk)
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: nlk
  real(kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2, work3
  integer :: iy
  
  ! divergence
  call cofdx ( nlk(:,:,1), work1 )
  call cofdy ( nlk(:,:,2), work2 )
  ! solve poisson eqn
  call poisson ( work1+work2, work3)
  ! gradient
  call cofdx (work3, work1)
  call cofdy (work3, work2)
  ! add gradient
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1) - work1(:,iy)    
    nlk(:,iy,2) = nlk(:,iy,2) - work2(:,iy)
  enddo
  !$omp end parallel do     
end subroutine add_pressure




subroutine add_diffusion (nlk, uk)
  ! adds the explicit diffusion term nu*laplace(u) to nlk
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (in) :: uk
  real(kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2
  integer :: iy

  ! first component
  call cofdxdx( uk(:,:,1), work1 )
  call cofdydy( uk(:,:,1), work2 )  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1) + nu*( work1(:,iy) + work2(:,iy) )    
  enddo
  !$omp end parallel do 
  
  ! second component
  call cofdxdx( uk(:,:,2), work1 )
  call cofdydy( uk(:,:,2), work2 )  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,2) = nlk(:,iy,2) + nu*( work1(:,iy) + work2(:,iy) )    
  enddo
  !$omp end parallel do   
end subroutine add_diffusion 
  


subroutine Runge (time, dt,it, u, uk, p, vort, nlk)
  use share_vars
  use FieldExport
  use PerformanceMeasurement
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (inout) :: vort, p
  real(kind=pr), intent (out) :: dt
  real(kind=pr), intent (in) :: time
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, workvis
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2) :: nlk2, uk_tmp
  integer :: iy
  integer, intent(in) :: it
  real(kind=pr) :: drag,lift,timestep

  dt = timestep(time, u)
  
  ! compute integrating factor
  call cal_vis(dt, workvis)
  
  call create_mask(time)
  call cal_nlk(time, u, uk, p, vort, nlk)
  call add_pressure(nlk)
   
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk_tmp(:,iy,1) = (uk(:,iy,1) + dt*nlk(:,iy,1) )*dealiase(:,iy)*workvis(:,iy)
    uk_tmp(:,iy,2) = (uk(:,iy,2) + dt*nlk(:,iy,2) )*dealiase(:,iy)*workvis(:,iy)
  enddo
  !$omp end parallel do 
  
  ! mean flow forcing
  call mean_flow(uk_tmp)
  
  !---------------------------------------------------------------------------------
  ! do second RK2 step (RHS evaluation with the argument defined above)
  !---------------------------------------------------------------------------------
  call create_mask(time+dt)
  call cal_nlk(time+dt, u, uk_tmp, p, vort, nlk2)  
  call add_pressure(nlk2)

  ! sum up all the terms.
  !$omp parallel do private (iy)
  do iy=0,ny-1
      uk(:,iy,1) = ( uk(:,iy,1)*workvis(:,iy) + 0.5*dt*( nlk(:,iy,1)*workvis(:,iy) + nlk2(:,iy,1) ) )*dealiase(:,iy)      
      uk(:,iy,2) = ( uk(:,iy,2)*workvis(:,iy) + 0.5*dt*( nlk(:,iy,2)*workvis(:,iy) + nlk2(:,iy,2) ) )*dealiase(:,iy)
  enddo
  !$omp end parallel do 
  
  ! mean flow forcing
  call mean_flow(uk)  
end subroutine Runge